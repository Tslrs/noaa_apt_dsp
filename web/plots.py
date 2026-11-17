"""Plotly figures for the APT web app (mobile friendly: light payloads, no hover-only info)."""
import numpy as np
import plotly.graph_objects as go
from scipy import signal

import apt_core as ac

LIGHT = dict(template="plotly_white", margin=dict(l=40, r=10, t=40, b=40), dragmode=False)


def _single(title, height=340):
    fig = go.Figure()
    fig.update_layout(title=dict(text=title, font=dict(size=15), x=0.02, xanchor="left"), height=height, **LIGHT)
    return fig


def _zplane(title, zeros, poles, lim):
    """Z-plane with a true circle (equal axis scale)."""
    fig = _single(title, height=380)
    th = np.linspace(0, 2 * np.pi, 200)
    fig.add_trace(go.Scatter(x=np.cos(th), y=np.sin(th), mode="lines", line=dict(color="#999", dash="dot"),
                             showlegend=False, hoverinfo="skip"))
    fig.add_trace(go.Scatter(x=zeros.real, y=zeros.imag, mode="markers", marker=dict(symbol="circle-open", size=9),
                             showlegend=False, name="zeros"))
    if len(poles):
        fig.add_trace(go.Scatter(x=poles.real, y=poles.imag, mode="markers",
                                 marker=dict(symbol="x", size=10, color="red"), showlegend=False, name="poles"))
    fig.update_xaxes(range=[-lim, lim], zeroline=True)
    fig.update_yaxes(range=[-lim, lim], scaleanchor="x", scaleratio=1, zeroline=True)
    return fig


def fig_filters(r, tr):
    """Four separate figures (one per row; they do not fit side by side on a phone)."""
    fs = r["fs"]
    b_fir, b_iir, a_iir = r["b_fir"], r["b_iir"], r["a_iir"]
    zi, pi_, _ = signal.tf2zpk(b_iir, a_iir)
    zf, _, _ = signal.tf2zpk(b_fir, [1.0])
    f1 = _zplane(tr["f_z_iir"], zi, pi_, 1.5)
    out_of_view = int(np.sum(np.abs(zf) > 1.6))
    f2 = _zplane(tr["f_z_fir"], zf[np.abs(zf) <= 1.6], np.array([]), 1.6)
    if out_of_view:
        f2.add_annotation(text=tr["zeros_out"].format(n=out_of_view), xref="paper", yref="paper",
                          x=0.5, y=-0.12, showarrow=False, font=dict(size=11, color="#666"))
    w1, h1 = signal.freqz(b_iir, a_iir, 2048, fs=fs)
    w2, h2 = signal.freqz(b_fir, 1, 2048, fs=fs)
    f3 = _single(tr["f_h_iir"])
    f3.add_trace(go.Scatter(x=w1, y=20 * np.log10(np.abs(h1) + 1e-12), mode="lines", showlegend=False))
    f4 = _single(tr["f_h_fir"])
    f4.add_trace(go.Scatter(x=w2, y=20 * np.log10(np.abs(h2) + 1e-12), mode="lines", showlegend=False))
    for f in (f3, f4):
        f.update_yaxes(range=[-80, 5], title_text="|H| (dB)")
        f.update_xaxes(title_text="Hz")
    return [f1, f2, f3, f4]


def fig_stft(y, fs, tr):
    seg = y[: min(len(y), int(45 * fs))]
    f, t, S = signal.spectrogram(seg, fs, window=signal.get_window("hamming", 512), noverlap=384, nfft=1024,
                                 mode="psd")
    keep = f <= 5000
    S_db = 10 * np.log10(S[keep] + 1e-20)[:, ::4]
    fig = go.Figure(go.Heatmap(x=t[::4], y=f[keep] / 1000.0, z=np.round(S_db, 1).astype(np.float32),
                               colorscale="Viridis", colorbar=dict(title="dB/Hz")))
    fig.update_layout(title=tr["f_stft"], xaxis_title=tr["time_s"], yaxis_title="kHz", height=420, **LIGHT)
    return fig


def fig_drift(r, tr):
    a = _single(tr["f_coarse"])
    a.add_trace(go.Scatter(x=r["coarse"], y=r["sc"], mode="lines", showlegend=False))
    a.add_vline(x=r["T_var"], line=dict(color="red", dash="dash"))
    a.update_xaxes(title_text=tr["samples_line"])
    a.update_yaxes(title_text="max Var")
    b = _single(tr["f_fine"])
    b.add_trace(go.Scatter(x=r["fine"], y=r["sf"], mode="lines", showlegend=False))
    b.add_vline(x=r["T_var"], line=dict(color="red", dash="dash"))
    b.add_vline(x=r["T_sync"], line=dict(color="green"))
    b.update_xaxes(title_text=tr["samples_line"], range=[r["fine"][0], max(r["T_sync"], r["fine"][-1]) + 0.01],
                   tickformat=".2f", nticks=5)
    return [a, b]


def fig_sync(r, tr):
    ok = r["sync_ok"]
    k = np.arange(1, len(r["sync_pos"]) + 1)
    fig = _single(tr["f_sync"])
    fig.add_trace(go.Scatter(x=k[ok], y=r["sync_pos"][ok], mode="markers", marker=dict(size=3), name="Sync-A"))
    fig.add_trace(go.Scatter(x=k, y=r["sync_fit"], mode="lines", line=dict(color="red"), name="fit"))
    fig.update_layout(showlegend=False)
    fig.update_xaxes(title_text=tr["line"])
    fig.update_yaxes(title_text=tr["sync_pos"])
    return fig


def fig_wedges(r, tr):
    cB, cA = r["calB"], r["calA"]
    f1 = _single(tr["f_tele"])
    f1.add_trace(go.Scatter(y=cB["t"][:384], mode="lines", showlegend=False))
    f1.update_xaxes(title_text=tr["line"])
    f1.update_yaxes(title_text="counts")
    f2 = _single(tr["f_wedge"], height=380)
    lv = ac.WEDGE_LEVELS
    f2.add_trace(go.Scatter(x=lv, y=cB["W"][:9], mode="markers", name="B"))
    f2.add_trace(go.Scatter(x=lv, y=cA["W"][:9], mode="markers", marker=dict(symbol="square-open"), name="A"))
    xx = np.array([0, 255.0])
    f2.add_trace(go.Scatter(x=xx, y=cB["a"] * xx + cB["b"], mode="lines", line=dict(dash="dash", color="red"),
                            name=tr["lin_fit"]))
    f2.update_xaxes(title_text=tr["wedge_mod"])
    f2.update_yaxes(title_text="counts")
    f2.update_layout(legend=dict(orientation="h", y=-0.3))
    return [f1, f2]


def fig_temperature(TC, sat, tr, step=3, pick_step=6):
    """Heatmap (hover shows deg C) + an invisible scatter grid so a tap/click can pick a point."""
    sub = np.round(TC[::step, ::step], 1).astype(np.float32)
    ys = np.arange(0, TC.shape[0], step) + 1          # 1-based line/column numbers
    xs = np.arange(0, TC.shape[1], step) + 1
    fig = go.Figure()
    fig.add_trace(go.Heatmap(z=sub, x=xs, y=ys, colorscale="Turbo", zmin=-70, zmax=60,
                             colorbar=dict(title="°C"),
                             hovertemplate=tr["hover"] + "<extra></extra>"))
    gy, gx = np.meshgrid(np.arange(0, TC.shape[0], pick_step) + 1, np.arange(0, TC.shape[1], pick_step) + 1, indexing="ij")
    fig.add_trace(go.Scattergl(x=gx.ravel(), y=gy.ravel(), mode="markers",
                               marker=dict(size=7, opacity=0.0), hoverinfo="skip", name="pick", showlegend=False))
    fig.update_yaxes(autorange="reversed", scaleanchor="x", title_text=tr["line"])
    fig.update_xaxes(title_text=tr["col"])
    fig.update_layout(title=f"{sat} - {tr['f_temp']}", height=720, clickmode="event+select",
                      **LIGHT)
    return fig


def add_pins(fig, pins):
    if not pins:
        return fig
    fig.add_trace(go.Scatter(x=[p[1] + 1 for p in pins], y=[p[0] + 1 for p in pins], mode="markers+text",
                             text=[str(i + 1) for i in range(len(pins))], textposition="middle right",
                             textfont=dict(color="white", size=15),
                             marker=dict(size=14, color="rgba(0,0,0,0)", line=dict(color="white", width=3)),
                             hoverinfo="skip", showlegend=False))
    return fig
