"""
NOAA APT decoder - Python port of satilite.m (MATLAB).

Pipeline (same steps as the MATLAB version):
  FIR bandpass (64, Hamming) -> Hilbert envelope -> IIR Butterworth lowpass (zero-phase)
  -> line-rate search (max column variance) -> sync-based refinement
  -> Sync-A correlation -> Doppler correction -> optional 180 deg rotation
  -> telemetry wedge calibration -> optional temperature (deg C).

Temperature constants: NOAA KLM User's Guide, via aptdec (GPL-2.0).
"""
import numpy as np
from scipy import signal
from scipy.fft import next_fast_len
from scipy.io import wavfile

W = 2080                      # pixels per APT line
WEDGE_LEVELS = np.array([31.07, 63.02, 94.96, 126.9, 158.86, 191.1, 228.62, 255.0, 0.0])

SAT_CONSTANTS = {
    "NOAA-15": dict(
        prt=np.array([[1.36328e-06, 0.051045, 276.60157], [1.47266e-06, 0.050909, 276.62531],
                      [1.47656e-06, 0.050907, 276.67413], [1.47656e-06, 0.050966, 276.59258]]),
        rad=np.array([[925.4075, 0.337810, 0.998719], [839.8979, 0.304558, 0.999024]]),
        cor=np.array([[-4.50, 0.0004524, -0.0932, 4.76], [-3.61, 0.0002811, -0.0659, 3.83]])),
    "NOAA-18": dict(
        prt=np.array([[1.657e-06, 0.05090, 276.601], [1.482e-06, 0.05101, 276.683],
                      [1.313e-06, 0.05117, 276.565], [1.484e-06, 0.05103, 276.615]]),
        rad=np.array([[928.1460, 0.436645, 0.998607], [833.2532, 0.253179, 0.999057]]),
        cor=np.array([[-5.53, 0.00052337, -0.11069, 5.82], [-2.22, 0.00017715, -0.04360, 2.67]])),
    "NOAA-19": dict(
        prt=np.array([[1.405783e-06, 0.051111, 276.6067], [1.496037e-06, 0.051090, 276.6119],
                      [1.496990e-06, 0.051033, 276.6311], [1.493110e-06, 0.051058, 276.6268]]),
        rad=np.array([[928.9, 0.53959, 0.998534], [831.9, 0.36064, 0.998913]]),
        cor=np.array([[-5.49, 0.00054668, -0.11187, 5.70], [-3.39, 0.00024985, -0.05991, 3.58]])),
}


# ----------------------------------------------------------------- helpers
def read_wav(file):
    fs, y = wavfile.read(file)
    if y.ndim > 1:
        y = y[:, 0]
    if y.dtype == np.int16:
        y = y.astype(np.float64) / 32768.0
    elif y.dtype == np.int32:
        y = y.astype(np.float64) / 2147483648.0
    elif y.dtype == np.uint8:
        y = (y.astype(np.float64) - 128.0) / 128.0
    else:
        y = y.astype(np.float64)
    return fs, y - y.mean()


def resample_rate(x, rho, chunk=1_000_000):
    """Linear interpolation at positions 0, rho, 2*rho, ... (MATLAB: interp1(x, 1:rho:end)).
    Done in chunks to keep memory low (the signal has ~10 million samples)."""
    n = len(x)
    count = int(np.floor((n - 1) / rho + 1e-9)) + 1
    out = np.empty(count, dtype=np.float64)
    for s0 in range(0, count, chunk):
        idx = rho * np.arange(s0, min(s0 + chunk, count))
        i0 = idx.astype(np.int64)
        f = idx - i0
        i1 = np.minimum(i0 + 1, n - 1)
        out[s0:s0 + len(idx)] = x[i0] * (1.0 - f) + x[i1] * f
    return out


def to_lines(x, n_line):
    nl = len(x) // n_line
    return x[: nl * n_line].reshape(nl, n_line)


def to_width(M, width=W):
    """Linear resampling of each row to `width` pixels (MATLAB: interp1(1:N, M', linspace(1,N,W)))."""
    n = M.shape[1]
    pos = np.linspace(0, n - 1, width)
    i0 = np.minimum(pos.astype(np.int64), n - 2)
    f = pos - i0
    return M[:, i0] * (1.0 - f) + M[:, i0 + 1] * f


def line_score(env_test, T, n_line):
    M = to_lines(resample_rate(env_test, T / n_line), n_line)
    return np.max(np.var(M, axis=0, ddof=1))


def movmedian(x, k):
    h = k // 2
    n = len(x)
    out = np.empty(n)
    for i in range(n):
        out[i] = np.median(x[max(0, i - h): min(n, i + h + 1)])
    return out


def smooth_gauss(x, w=101):
    """MATLAB smoothdata(x,'gaussian',w): sigma = w/5, endpoints shrink."""
    h = w // 2
    sigma = w / 5.0
    d = np.arange(-h, h + 1)
    g = np.exp(-0.5 * (d / sigma) ** 2)
    n = len(x)
    out = np.empty(n)
    for i in range(n):
        lo, hi = max(0, i - h), min(n, i + h + 1)
        gg = g[lo - i + h: hi - i + h]
        out[i] = np.sum(gg * x[lo:hi]) / np.sum(gg)
    return out


def round_away(x):
    return np.sign(x) * np.floor(np.abs(x) + 0.5)


def prctile(x, p):
    """MATLAB-style percentile."""
    x = np.sort(np.asarray(x).ravel())
    n = len(x)
    pos = np.clip(p / 100.0 * n + 0.5, 1, n)
    lo = int(np.floor(pos))
    hi = min(lo + 1, n)
    return x[lo - 1] + (pos - lo) * (x[hi - 1] - x[lo - 1])


def sync_positions(M, tpl):
    """1-based column of the Sync-A correlation peak for every line."""
    tp = np.zeros(M.shape[1])
    tp[: len(tpl)] = tpl
    X = M - M.mean(axis=1, keepdims=True)
    xc = np.real(np.fft.ifft(np.fft.fft(X, axis=1) * np.conj(np.fft.fft(tp)), axis=1))
    return np.argmax(xc, axis=1) + 1.0


def unwrap_cols(pos, width):
    m = np.median(pos)
    d = np.mod(pos - m + width / 2.0, width) - width / 2.0
    return m + d


# ----------------------------------------------------------------- telemetry
def wedge_calibrate(img, cols, flipped):
    t = np.median(img[:, cols], axis=1)
    if flipped:
        t = t[::-1]
    lv = WEDGE_LEVELS
    tpl = np.repeat(lv, 8)
    n_l = len(t)
    nf = n_l // 128 - 1
    best, off = -np.inf, 0
    tc = tpl - tpl.mean()
    for o in range(128):
        sc = 0.0
        for f in range(nf):
            seg = t[o + f * 128: o + f * 128 + 72]
            seg = seg - seg.mean()
            sc += float(seg @ tc) / (np.linalg.norm(seg) * np.linalg.norm(tc) + 1e-16)
        sc /= nf
        if sc > best:
            best, off = sc, o
    nfr = (n_l - off) // 128
    Wf = np.zeros((nfr, 16))
    for f in range(nfr):
        for k in range(16):
            idx = off + f * 128 + k * 8 + np.arange(1, 7)      # skip transition lines
            Wf[f, k] = t[idx].mean()
    Wm = np.median(Wf, axis=0)
    a, b = np.polyfit(lv, Wm[:9], 1)
    rms = float(np.sqrt(np.mean((Wm[:9] - (a * lv + b)) ** 2)))
    mod_id = (Wm[15] - b) / a
    chid = int(np.argmin(np.abs(mod_id - np.array([31, 63, 95, 127, 159, 191])))) + 1
    return dict(t=t, offset=off, W=Wm, a=a, b=b, rms=rms, channel_id=chid, score=best)


def apply_wedge_lut(x, cal):
    lv = np.concatenate(([0.0], WEDGE_LEVELS[:8]))
    w = np.concatenate(([cal["W"][8]], cal["W"][:8]))
    return np.clip(np.interp(x, w, lv) / 255.0, 0, 1)


# ----------------------------------------------------------------- temperature
def thermal_calibrate(Ce255, cal, satellite, ch, Cs255):
    k = SAT_CONSTANTS[satellite]
    ci = ch - 4
    Vc, A, B = k["rad"][ci]
    C1, C2 = 1.1910427e-5, 1.4387752
    lv = np.concatenate(([0.0], WEDGE_LEVELS[:8]))
    w = np.concatenate(([cal["W"][8]], cal["W"][:8]))
    tele = np.clip(np.interp(cal["W"], w, lv), 0, 255)
    Tprt = np.array([k["prt"][n, 0] * (tele[9 + n] * 4) ** 2 + k["prt"][n, 1] * (tele[9 + n] * 4) + k["prt"][n, 2]
                     for n in range(4)])
    Tbb = Tprt.mean()
    Tbbs = A + Tbb * B
    Nbb = C1 * Vc ** 3 / (np.exp(C2 * Vc / Tbbs) - 1)
    Cs = Cs255 * 4
    Cb = tele[14] * 4
    Ns = k["cor"][ci, 0]
    q = k["cor"][ci, 1:4]
    Nl = Ns + (Nbb - Ns) * (Cs - 4 * Ce255) / (Cs - Cb)
    Nc = q[0] * Nl ** 2 + q[1] * Nl + q[2]
    Ne = np.maximum(Nl + Nc, 1e-3)
    Tes = C2 * Vc / np.log(C1 * Vc ** 3 / Ne + 1)
    TC = np.clip((Tes - A) / B - 273.15, -100, 60)
    par = dict(Vc=Vc, A=A, B=B, C1=C1, C2=C2, Ns=Ns, q=q, Nbb=Nbb, Cs=Cs, Cb=Cb, Tbbs=Tbbs)
    return dict(TC=TC, Tbb=Tbb, Tprt=Tprt, Nbb=Nbb, Cs=Cs / 4, Cb=Cb / 4, par=par)


def pixel_steps(par, Ce):
    """Step-by-step temperature for one pixel (used by the 'click a point' view)."""
    Nl = par["Ns"] + (par["Nbb"] - par["Ns"]) * (par["Cs"] - 4 * Ce) / (par["Cs"] - par["Cb"])
    q = par["q"]
    Nc = q[0] * Nl ** 2 + q[1] * Nl + q[2]
    Ne = max(Nl + Nc, 1e-3)
    Tes = par["C2"] * par["Vc"] / np.log(par["C1"] * par["Vc"] ** 3 / Ne + 1)
    Tk = (Tes - par["A"]) / par["B"]
    Tc = Tk - 273.15
    return dict(Ce=Ce, Nl=Nl, Nc=Nc, Ne=Ne, Tes=Tes, Tk=Tk, Tc=Tc, clipped=(Tc < -100 or Tc > 60),
                Tshow=float(np.clip(Tc, -100, 60)))


# ----------------------------------------------------------------- main decode
def decode(y, fs, rotate180=True, satellite="", progress=None):
    """Decode an APT recording. Returns a dict with images, calibration and diagnostics."""
    def prog(msg, frac):
        if progress:
            progress(msg, frac)

    out = {}
    prog("Filters", 0.05)
    n_fir = 64
    b_fir = signal.firwin(n_fir + 1, [1200, 3600], pass_zero=False, window="hamming", fs=fs)
    b_iir, a_iir = signal.butter(3, 2000 / (fs / 2), "low")
    poles = np.roots(a_iir)
    out.update(b_fir=b_fir, b_iir=b_iir, a_iir=a_iir, max_pole=float(np.max(np.abs(poles))), fs=fs)

    y_bp = signal.lfilter(b_fir, 1, np.concatenate((y, np.zeros(n_fir // 2))))[n_fir // 2:]
    prog("Hilbert demodulation", 0.15)
    # zero-pad to a fast FFT length (exact odd lengths would need Bluestein: ~2.5 GB RAM)
    n_sig = len(y_bp)
    env = np.abs(signal.hilbert(y_bp.astype(np.float32), next_fast_len(n_sig)))[:n_sig].astype(np.float64)
    del y_bp
    env = signal.filtfilt(b_iir, a_iir, env)

    # ---- line rate: max column variance
    prog("Line-rate search", 0.30)
    T_nom = fs / 2.0
    n_line = int(np.floor(T_nom))
    env_test = env[: min(len(env), int(fs * 120))]
    coarse = np.arange(-2.0, 2.0 + 1e-9, 0.05)
    sc = np.array([line_score(env_test, T_nom + d, n_line) for d in coarse])
    ic = int(np.argmax(sc))
    fine = np.arange(coarse[ic] - 0.05, coarse[ic] + 0.05 + 1e-9, 0.002)
    sf = np.array([line_score(env_test, T_nom + d, n_line) for d in fine])
    T_var = T_nom + fine[int(np.argmax(sf))]
    out.update(T_nom=T_nom, coarse=T_nom + coarse, sc=sc, fine=T_nom + fine, sf=sf, T_var=T_var)

    # ---- sync-based refinement
    prog("Sync refinement", 0.50)
    tpl = np.tile([1.0, 1.0, -1.0, -1.0], 7)
    T_ref = T_var
    pos = ok = cf = k = None
    for _ in range(3):
        rho = T_ref / n_line
        M = to_width(to_lines(resample_rate(env, rho), n_line))
        pos = unwrap_cols(sync_positions(M, tpl), W)
        k = np.arange(1, len(pos) + 1, dtype=float)
        ok = np.abs(pos - movmedian(pos, 101)) < 4
        cf = np.polyfit(k[ok], pos[ok], 1)
        T_ref = T_ref + cf[0] * n_line / W * rho
    T_true = T_ref
    sync_jitter = float(np.std(pos[ok] - np.polyval(cf, k[ok]), ddof=1))
    out.update(T_sync=T_true, sync_jitter=sync_jitter, sync_valid=float(np.mean(ok)),
               sync_pos=pos, sync_ok=ok, sync_fit=np.polyval(cf, k))

    prog("Building image", 0.65)
    rho = T_true / n_line
    img = to_width(to_lines(resample_rate(env, rho), n_line))
    del env

    # ---- Sync-A + Doppler
    prog("Sync-A and Doppler", 0.75)
    avg = img.mean(axis=0)
    avg = avg - avg.mean()
    tp = np.zeros(W)
    tp[: len(tpl)] = tpl
    xc = np.real(np.fft.ifft(np.fft.fft(avg) * np.conj(np.fft.fft(tp))))
    sync_col = int(np.argmax(xc)) + 1
    sync_quality = float(xc.max() / np.mean(np.abs(xc)))
    img = np.roll(img, -(sync_col - 1), axis=1)

    pos2 = unwrap_cols(sync_positions(img, tpl), W)
    bad = np.abs(pos2 - movmedian(pos2, 101)) >= 4
    good = np.where(~bad)[0]
    pos2 = np.interp(np.arange(len(pos2)), good, pos2[good])
    pos_s = smooth_gauss(pos2, 101)
    doppler_shift = round_away(pos_s - 1).astype(int)
    for r in range(img.shape[0]):
        img[r] = np.roll(img[r], -doppler_shift[r])
    out.update(sync_col=sync_col, sync_quality=sync_quality, doppler_max=int(np.max(np.abs(doppler_shift))))

    if rotate180:
        img = img[::-1, ::-1]
        img = np.roll(img, 1040, axis=1)

    # ---- telemetry
    prog("Telemetry calibration", 0.85)
    if rotate180:
        telA, telB = np.arange(5, 40), np.arange(1045, 1080)
    else:
        telA, telB = np.arange(1000, 1035), np.arange(2040, 2075)
    calA = wedge_calibrate(img, telA, rotate180)
    calB = wedge_calibrate(img, telB, rotate180)

    # ---- display images
    lo, hi = prctile(img, 0.5), prctile(img, 99.5)
    img_n = np.clip((img - lo) / (hi - lo), 0, 1) ** 0.75
    chA, chB = img_n[:, :1040], img_n[:, 1040:]
    c0 = 48 if rotate180 else 86                      # 0-based (MATLAB: 49 / 87)
    rawA, rawB = img[:, :1040], img[:, 1040:]
    calimgA = apply_wedge_lut(rawA[:, c0:c0 + 909], calA)
    calimgB = apply_wedge_lut(rawB[:, c0:c0 + 909], calB)

    temp = None
    if satellite and calB["channel_id"] in (4, 5):
        spc = np.arange(961, 1000) if rotate180 else np.arange(43, 82)
        Cs_meas = float(np.median(apply_wedge_lut(rawB[:, spc], calB) * 255))
        temp = thermal_calibrate(calimgB * 255, calB, satellite, calB["channel_id"], Cs_meas)

    prog("Done", 1.0)
    out.update(img_n=img_n, chA=chA, chB=chB, cropA=chA[:, c0:c0 + 909], cropB=chB[:, c0:c0 + 909],
               calimgA=calimgA, calimgB=calimgB, CeB=calimgB * 255, calA=calA, calB=calB,
               temp=temp, T_line=T_true, n_lines=img.shape[0], rotate180=rotate180, satellite=satellite)
    return out
