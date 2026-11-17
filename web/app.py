"""NOAA APT decoder - web app (Streamlit). Run:  streamlit run app.py"""
import io
import zipfile
from pathlib import Path

import numpy as np
import streamlit as st
import streamlit.components.v1 as components
from matplotlib import colormaps
from PIL import Image

import apt_core as ac
import plots

TR = {
    "el": dict(
        title="NOAA APT Decoder", subtitle="Αποκωδικοποίηση δορυφορικών εικόνων NOAA APT από ηχητική λήψη (wav)",
        menu_btn="☰",
        howto_title="Πώς χρησιμοποιείται (3 βήματα)",
        howto="""1. **Ανέβασε** μια ηχογράφηση NOAA APT (wav) ή τσέκαρε «Χρήση δείγματος» αν δεν έχεις.
2. **Διάλεξε δορυφόρο** (NOAA-15/18/19) αν τον ξέρεις, για να βγουν θερμοκρασίες σε °C.
3. Πάτα **Αποκωδικοποίηση** και άνοιξε τις προβολές: *Εικόνα*, *Κανάλια*, *Θερμοκρασία* (πέρασε το δάχτυλο/ποντίκι πάνω στον χάρτη) κ.λπ.

*Το APT είναι το αναλογικό σήμα εικόνας που εκπέμπουν οι μετεωρολογικοί δορυφόροι NOAA στα 137 MHz. Η εφαρμογή το μετατρέπει σε εικόνα με ψηφιακή επεξεργασία σήματος.*""",
        lang="Γλώσσα", upload="1. Ανέβασε αρχείο ήχου (wav, 11025 Hz προτείνεται)", sat="2. Δορυφόρος (για °C)",
        sat_unknown="Άγνωστος (χωρίς °C)", rot="3. Πέρασμα προς Βορρά (περιστροφή 180°)",
        rot_help="Αν η εικόνα βγει ανάποδα, ξετσέκαρέ το.", run="Αποκωδικοποίηση",
        need_file="Ανέβασε πρώτα ένα αρχείο wav.", running="Επεξεργασία σήματος...",
        use_sample="Χρήση δείγματος (noaa3.wav, Ευρώπη, δορυφόρος άγνωστος)", err="Σφάλμα αποκωδικοποίησης", weak="Αδύναμο sync: πιθανώς δεν υπάρχει σήμα APT στο αρχείο.",
        resample="Το αρχείο έχει συχνότητα {fs} Hz. Ο αλγόριθμος υποθέτει 11025 Hz· τα αποτελέσματα μπορεί να είναι λάθος.",
        m_lines="Γραμμές", m_rate="Ρυθμός γραμμής", m_sync="Ποιότητα sync", m_chA="Κανάλι A (ID)", m_chB="Κανάλι B (ID)",
        tabs=["Εικόνα", "Κανάλια", "Θερμοκρασία", "Φίλτρα", "STFT", "Drift / Sync", "Wedges"],
        image="Τελική εικόνα (κανάλι A αριστερά, B δεξιά)", chA="Κανάλι A (ορατό/εγγύς IR) βαθμονομημένο",
        chB="Κανάλι B (θερμικό IR) - ψευδοχρώματα", download="Λήψη όλων (ZIP)",
        no_temp="Επίλεξε δορυφόρο (NOAA-15/18/19) και ξανατρέξε για να δεις θερμοκρασίες.",
        temp_help="Πέρασε το ποντίκι (ή άγγιξε) πάνω στον χάρτη για να δεις τη θερμοκρασία. Κάνε κλικ/άγγιγμα για να βάλεις σημείο, ή γράψε γραμμή/στήλη παρακάτω. Για να αφαιρέσεις σημείο, διάλεξέ το από τη λίστα παρακάτω.",
        hover="Γραμμή %{y}, στήλη %{x}: %{z:.1f} °C", f_temp="κανάλι B (θερμικό IR)", line="Γραμμή", col="Στήλη",
        manual="Προσθήκη σημείου με συντεταγμένες", remove="Αφαίρεση σημείου", remove_which="Ποιο σημείο να αφαιρεθεί;", add="Προσθήκη", clear="Καθαρισμός όλων",
        pins="Σημεία", pin_line="γραμμή {r}, στήλη {c}: {t:.1f} °C", steps="Υπολογισμός για το σημείο #{n}",
        how="Πώς υπολογίζεται η θερμοκρασία", clipped="(τιμή εκτός −100..60 °C αποκόπηκε: πιθανώς θόρυβος)",
        zeros_out="{n} μηδενικά εκτός εικόνας (−3.0 και 3.5)", f_z_iir="Z-plane IIR Butterworth 3ης τάξης", f_z_fir="Z-plane FIR bandpass (64 μηδενικά)",
        f_h_iir="IIR lowpass, fc = 2000 Hz", f_h_fir="FIR bandpass 2.4 kHz ± 1.2 kHz",
        f_stft="STFT - φέρουσα 2400 Hz", time_s="Χρόνος (s)", f_coarse="Χονδρική σάρωση", f_fine="Λεπτομερής σάρωση",
        samples_line="δείγματα/γραμμή", f_sync="Σταθερότητα συγχρονισμού (παραβολή = Doppler)",
        sync_pos="Θέση Sync-A (pixel)", f_tele="Τηλεμετρία καναλιού B", f_wedge="Wedges vs επίπεδα αναφοράς",
        lin_fit="Γραμμική προσαρμογή", wedge_mod="Διαμόρφωση wedge",
        note="Η θερμοκρασία υπολογίζεται με σταθερές NOAA KLM (μέσω aptdec). Η απόλυτη ακρίβεια δεν έχει επαληθευτεί.",
        how_text="""
1. **Ce (0-255):** φωτεινότητα pixel, βαθμονομημένη με τα 16 wedges τηλεμετρίας (αφαιρεί κέρδος/offset δέκτη).
2. **Μέλαν σώμα:** Tbb = μέσος όρος 4 αισθητήρων PRT (wedges 10-13), καθένας T = a·x² + b·x + c με x = 4·count wedge (10 bit).
   T* = A + B·Tbb = **{Tbbs:.2f} K**, Nbb = C1·Vc³/(exp(C2·Vc/T*) − 1) = **{Nbb:.3f}**
3. **Space view / back scan:** Cs = **{Cs:.1f}**, Cb = **{Cb:.1f}** (counts 10 bit)
4. **Γραμμική ακτινοβολία:** Nl = Ns + (Nbb − Ns)·(Cs − 4·Ce)/(Cs − Cb)
5. **Μη γραμμική διόρθωση:** Nc = a·Nl² + b·Nl + c, Ne = Nl + Nc
6. **Θερμοκρασία:** Te* = C2·Vc / ln(C1·Vc³/Ne + 1), T = (Te* − A)/B − 273.15 °C

Vc = {Vc:.2f} cm⁻¹, A = {A:.4f}, B = {B:.6f}, Ns = {Ns:.2f}
""",
        pt_lines=("Ce = {Ce:.1f}   (4·Ce = {Ce4:.1f} στα 10 bit)", "Nl = {Nl:.3f}   (γραμμική ακτινοβολία)",
                  "Nc = {Nc:+.3f}   (μη γραμμική διόρθωση)", "Ne = {Ne:.3f}   (διορθωμένη ακτινοβολία)",
                  "Te* = {Tes:.2f} K  (ισοδύναμη θερμοκρασία μέλανος σώματος)", "T = (Te* − A)/B = {Tk:.2f} K",
                  "ΘΕΡΜΟΚΡΑΣΙΑ = {Tshow:.1f} °C"),
    ),
    "en": dict(
        title="NOAA APT Decoder", subtitle="Decode NOAA APT weather-satellite images from an audio recording (wav)",
        menu_btn="☰",
        howto_title="How to use (3 steps)",
        howto="""1. **Upload** a NOAA APT recording (wav), or tick "Use the sample" if you have none.
2. **Choose the satellite** (NOAA-15/18/19) if you know it, to get temperatures in °C.
3. Press **Decode** and open the views: *Image*, *Channels*, *Temperature* (move your finger/mouse over the map), etc.

*APT is the analog image signal broadcast by NOAA weather satellites at 137 MHz. This app turns it into an image using digital signal processing.*""",
        lang="Language", upload="1. Upload an audio file (wav, 11025 Hz recommended)", sat="2. Satellite (for °C)",
        sat_unknown="Unknown (no °C)", rot="3. Northbound pass (rotate 180°)",
        rot_help="Untick if the image comes out upside down.", run="Decode",
        need_file="Upload a wav file first.", running="Processing signal...",
        use_sample="Use the sample recording (noaa3.wav, Europe, satellite unknown)", err="Decoding error", weak="Weak sync: the file probably contains no APT signal.",
        resample="The file sample rate is {fs} Hz. The algorithm assumes 11025 Hz; results may be wrong.",
        m_lines="Lines", m_rate="Line rate", m_sync="Sync quality", m_chA="Channel A (ID)", m_chB="Channel B (ID)",
        tabs=["Image", "Channels", "Temperature", "Filters", "STFT", "Drift / Sync", "Wedges"],
        image="Final image (channel A left, B right)", chA="Channel A (visible/near IR), calibrated",
        chB="Channel B (thermal IR) - false colour", download="Download all (ZIP)",
        no_temp="Select a satellite (NOAA-15/18/19) and decode again to get temperatures.",
        temp_help="Hover (or touch) the map to read the temperature. Click/tap to place a point, or type line/column below. To remove a point, pick it from the list below.",
        hover="Line %{y}, column %{x}: %{z:.1f} °C", f_temp="channel B (thermal IR)", line="Line", col="Column",
        manual="Add a point by coordinates", remove="Remove point", remove_which="Which point to remove?", add="Add", clear="Clear all",
        pins="Points", pin_line="line {r}, column {c}: {t:.1f} °C", steps="Calculation for point #{n}",
        how="How temperature is computed", clipped="(value outside −100..60 °C was clamped: probably noise)",
        zeros_out="{n} zeros outside the view (−3.0 and 3.5)", f_z_iir="Z-plane IIR Butterworth, order 3", f_z_fir="Z-plane FIR bandpass (64 zeros)",
        f_h_iir="IIR lowpass, fc = 2000 Hz", f_h_fir="FIR bandpass 2.4 kHz ± 1.2 kHz",
        f_stft="STFT - 2400 Hz subcarrier", time_s="Time (s)", f_coarse="Coarse scan", f_fine="Fine scan",
        samples_line="samples/line", f_sync="Sync stability (parabola = Doppler)",
        sync_pos="Sync-A position (pixel)", f_tele="Channel B telemetry", f_wedge="Wedges vs reference levels",
        lin_fit="Linear fit", wedge_mod="Wedge modulation",
        note="Temperature uses NOAA KLM constants (via aptdec). Absolute accuracy has not been verified.",
        how_text="""
1. **Ce (0-255):** pixel brightness, calibrated with the 16 telemetry wedges (removes receiver gain/offset).
2. **Blackbody:** Tbb = mean of 4 PRT sensors (wedges 10-13), each T = a·x² + b·x + c with x = 4·wedge count (10 bit).
   T* = A + B·Tbb = **{Tbbs:.2f} K**, Nbb = C1·Vc³/(exp(C2·Vc/T*) − 1) = **{Nbb:.3f}**
3. **Space view / back scan:** Cs = **{Cs:.1f}**, Cb = **{Cb:.1f}** (10-bit counts)
4. **Linear radiance:** Nl = Ns + (Nbb − Ns)·(Cs − 4·Ce)/(Cs − Cb)
5. **Non-linear correction:** Nc = a·Nl² + b·Nl + c, Ne = Nl + Nc
6. **Temperature:** Te* = C2·Vc / ln(C1·Vc³/Ne + 1), T = (Te* − A)/B − 273.15 °C

Vc = {Vc:.2f} cm⁻¹, A = {A:.4f}, B = {B:.6f}, Ns = {Ns:.2f}
""",
        pt_lines=("Ce = {Ce:.1f}   (4·Ce = {Ce4:.1f} of 10 bit)", "Nl = {Nl:.3f}   (linear radiance)",
                  "Nc = {Nc:+.3f}   (non-linear correction)", "Ne = {Ne:.3f}   (corrected radiance)",
                  "Te* = {Tes:.2f} K  (equivalent blackbody temperature)", "T = (Te* − A)/B = {Tk:.2f} K",
                  "TEMPERATURE = {Tshow:.1f} °C"),
    ),
}

SAMPLE = Path(__file__).parent / "sample" / "noaa3.wav"
st.set_page_config(page_title="NOAA APT Decoder", page_icon="🛰️", layout="centered")
ss = st.session_state

# Bigger touch targets for phones (Apple recommends >= 44 px)
st.markdown("""
<style>
/* Page width: centred column (about a phone-tablet width) on a computer, full width on a phone */
.block-container { max-width: 860px !important; }
.stButton > button, .stDownloadButton > button { min-height: 52px; font-size: 1.05rem; }
div[role="radiogroup"] { gap: 6px; flex-wrap: wrap; }
div[role="radiogroup"] > label { min-height: 46px; padding: 6px 14px; border-radius: 12px;
    border: 1px solid rgba(128,128,128,0.45); align-items: center; }
div[data-baseweb="select"] > div { min-height: 50px; }
.stCheckbox label { min-height: 44px; align-items: center; }
@media (max-width: 640px) {
  .block-container { padding-left: 0.8rem; padding-right: 0.8rem; padding-bottom: 6rem; }
}
</style>
""", unsafe_allow_html=True)

ss.setdefault("pins", [])
ss.setdefault("res", None)
ss.setdefault("last_sel", None)

# ---------------------------------------------------------------- language
lang_label = st.radio("Language / Γλώσσα", ["Ελληνικά", "English"], horizontal=True, label_visibility="collapsed")
tr = TR["el" if lang_label == "Ελληνικά" else "en"]
st.title("🛰️ " + tr["title"])
st.caption(tr["subtitle"])
with st.expander(tr["howto_title"], expanded=(ss.res is None)):
    st.markdown(tr["howto"])

# ---------------------------------------------------------------- inputs
up = st.file_uploader(tr["upload"], type=["wav"])
use_sample = st.checkbox(tr["use_sample"], value=False)
c1, c2 = st.columns(2)
sat_items = {tr["sat_unknown"]: "", "NOAA-15": "NOAA-15", "NOAA-18": "NOAA-18", "NOAA-19": "NOAA-19"}
sat_label = c1.selectbox(tr["sat"], list(sat_items))
sat = sat_items[sat_label]
rot = c2.checkbox(tr["rot"], value=True, help=tr["rot_help"])

if st.button(tr["run"], type="primary", width="stretch"):
    if up is None and not use_sample:
        st.warning(tr["need_file"])
    else:
        prog = st.progress(0.0, text=tr["running"])
        try:
            fs, y = ac.read_wav(SAMPLE if use_sample else up)
            res = ac.decode(y, fs, rotate180=rot, satellite=sat,
                            progress=lambda msg, f: prog.progress(min(f, 1.0), text=f"{tr['running']} {msg}"))
            res["y45"] = y[: int(45 * fs)].copy()
            res["file_fs"] = fs
            ss.res, ss.pins, ss.last_sel = res, [], None
        except Exception as e:  # noqa: BLE001
            st.error(f"{tr['err']}: {e}")
        finally:
            prog.empty()

res = ss.res
if res is None:
    st.stop()

# ---------------------------------------------------------------- summary
if res["file_fs"] != 11025:
    st.warning(tr["resample"].format(fs=res["file_fs"]))
if res["sync_quality"] < 8:
    st.warning(tr["weak"])
st.markdown(
    f"**{tr['m_lines']}** {res['n_lines']} · **{tr['m_rate']}** {res['T_line']:.3f} · "
    f"**{tr['m_sync']}** {res['sync_quality']:.1f} · **A** ID {res['calA']['channel_id']} · "
    f"**B** ID {res['calB']['channel_id']}"
)


def to_u8(a):
    return (np.clip(a, 0, 1) * 255).round().astype(np.uint8)


def png_bytes(arr):
    buf = io.BytesIO()
    Image.fromarray(arr).save(buf, format="PNG")
    return buf.getvalue()


turbo = colormaps["turbo"]
T = res["temp"]

# ---- view bar: [<] [menu] [>]  + swipe left/right on touch screens
st.markdown("""
<style>
.st-key-viewbar [data-testid="stHorizontalBlock"] { flex-wrap: nowrap !important; align-items: center; gap: 0.4rem; }
.st-key-viewbar [data-testid="stColumn"] { min-width: 0 !important; }
.st-key-viewbar .stButton > button, .st-key-viewbar [data-testid="stPopover"] > button { min-height: 52px; padding: 0 6px; }
.st-key-viewbar .stButton > button p { font-size: 1.7rem; line-height: 1; margin: 0; }
.st-key-viewbar [data-testid="stPopover"] > button p { font-size: 1.7rem; line-height: 1; margin: 0; }
.st-key-viewbar [data-testid="stPopoverButton"] > div:nth-child(2) { display: none !important; }   /* hide the chevron */
.st-key-viewbar [data-testid="stPopover"] > button { justify-content: center; }
/* Wide screens (computer): keep the view bar compact and centred instead of stretching across the page */
@media (min-width: 900px) {
  .st-key-viewbar { max-width: 560px; margin-left: auto; margin-right: auto; }
}
.viewlabel { background: #33709f; color: #fff; border-radius: 18px; padding: 22px 12px; text-align: center;
    font-size: 1.35rem; font-weight: 700; margin-bottom: 0.6rem; }
</style>
""", unsafe_allow_html=True)

N_VIEWS = len(tr["tabs"])


def _step(d):
    ss["view"] = min(max(int(ss.get("view", 0)) + d, 0), N_VIEWS - 1)


ss.setdefault("view", 0)


def _set(i):
    ss["view"] = i


with st.container(key="viewbar"):
    st.markdown(f"<div class='viewlabel'>{tr['tabs'][ss['view']]}</div>", unsafe_allow_html=True)   # where we are
    c_prev, c_menu, c_next = st.columns([1, 1, 1], gap="small")
    c_prev.button("◀", key="v_prev", on_click=_step, args=(-1,), width="stretch")
    with c_menu.popover(tr["menu_btn"], width="stretch"):                                    # open the list
        for _i, _name in enumerate(tr["tabs"]):
            st.button(_name, key=f"pick{_i}", on_click=_set, args=(_i,), width="stretch",
                      type="primary" if _i == ss["view"] else "secondary")
    c_next.button("▶", key="v_next", on_click=_step, args=(1,), width="stretch")
view = ss["view"]

# Swipe left/right (touch): clicks the [<] / [>] buttons. Ignores gestures that start at the screen edge
# (iOS uses those for 'back'), slow drags and mostly-vertical drags (page scrolling).
components.html("""
<script>
(function () {
  const P = window.parent;
  if (P.__aptSwipeInstalled) return;
  P.__aptSwipeInstalled = true;
  const doc = P.document;
  let x0 = null, y0 = 0, t0 = 0;
  doc.addEventListener('touchstart', function (e) {
    if (e.touches.length !== 1) { x0 = null; return; }
    const x = e.touches[0].clientX;
    if (x < 24 || x > P.innerWidth - 24) { x0 = null; return; }
    x0 = x; y0 = e.touches[0].clientY; t0 = Date.now();
  }, {passive: true});
  doc.addEventListener('touchend', function (e) {
    if (x0 === null) return;
    const t = e.changedTouches[0];
    const dx = t.clientX - x0, dy = t.clientY - y0;
    x0 = null;
    if (Date.now() - t0 > 700 || Math.abs(dx) < 90 || Math.abs(dy) > Math.abs(dx) * 0.5) return;
    const want = dx < 0 ? '\u25B6' : '\u25C0';
    const btns = doc.querySelectorAll('.st-key-viewbar button');
    for (const b of btns) { if (b.innerText.trim() === want) { b.click(); break; } }
  }, {passive: true});
})();
</script>
""", height=0)

PLOT_CFG = {"displayModeBar": False, "scrollZoom": False, "doubleClick": False}   # no toolbar, page still scrolls

# ---- image
if view == 0:
    st.image(to_u8(res["img_n"]), caption=tr["image"], use_container_width=True)

# ---- channels + downloads
if view == 1:
    B_color = (turbo(res["calimgB"])[..., :3] * 255).astype(np.uint8)
    st.image(to_u8(res["calimgA"]), caption=tr["chA"], use_container_width=True)
    st.image(B_color, caption=tr["chB"], use_container_width=True)
    if "_zip" not in res:                       # build the ZIP only once per decoded file
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
            z.writestr("apt_full.png", png_bytes(to_u8(res["img_n"])))
            z.writestr("apt_A_cal.png", png_bytes(to_u8(res["calimgA"])))
            z.writestr("apt_B_cal.png", png_bytes(to_u8(res["calimgB"])))
            z.writestr("apt_B_cal_color.png", png_bytes(B_color))
            if T is not None:
                z.writestr("apt_B_temp_gray.png", png_bytes(to_u8((T["TC"] + 100) / 160)))
                z.writestr("apt_B_temp_color.png",
                           png_bytes((turbo((T["TC"] + 70) / 130)[..., :3] * 255).astype(np.uint8)))
        res["_zip"] = buf.getvalue()
    st.download_button(tr["download"], res["_zip"], file_name="noaa_apt_results.zip", mime="application/zip")

# ---- temperature (hover + tap/click pins)
if view == 2:
    if T is None:
        st.info(tr["no_temp"])
    else:
        st.caption(tr["temp_help"])
        fig = plots.add_pins(plots.fig_temperature(T["TC"], res["satellite"], tr), ss.pins)
        event = st.plotly_chart(fig, use_container_width=True, config=PLOT_CFG, on_select="rerun", selection_mode="points",
                                key="tempmap")
        pts = (event.selection.points if event and event.selection else []) or []
        if pts and pts[0].get("curve_number") == 1:
            row, col = int(round(pts[0]["y"])) - 1, int(round(pts[0]["x"])) - 1
            row = min(max(row, 0), T["TC"].shape[0] - 1)
            col = min(max(col, 0), T["TC"].shape[1] - 1)
            if ss.last_sel != (row, col):
                ss.last_sel = (row, col)
                ss.pins.append((row, col))
                st.rerun()

        n_r, n_c = T["TC"].shape
        with st.expander(tr["manual"]):
            a, b, c = st.columns([2, 2, 1])
            r_in = a.number_input(tr["line"], 1, n_r, 900)
            c_in = b.number_input(tr["col"], 1, n_c, 450)
            if c.button(tr["add"]):
                ss.pins.append((int(r_in) - 1, int(c_in) - 1))
                st.rerun()

        if ss.pins:
            st.subheader(tr["pins"])
            for i, (r_, c_) in enumerate(ss.pins):
                st.markdown(f"**#{i + 1}**  " + tr["pin_line"].format(r=r_ + 1, c=c_ + 1, t=T["TC"][r_, c_]))
            rm = st.selectbox(tr["remove_which"], list(range(1, len(ss.pins) + 1)), key="rm_sel")
            b1, b2 = st.columns(2)
            if b1.button(tr["remove"], width="stretch"):
                ss.pins.pop(int(rm) - 1)
                st.rerun()
            if b2.button(tr["clear"], width="stretch"):
                ss.pins = []
                ss.last_sel = None
                st.rerun()
            r_, c_ = ss.pins[-1]
            s = ac.pixel_steps(T["par"], float(res["CeB"][r_, c_]))
            st.markdown(f"**{tr['steps'].format(n=len(ss.pins))}**")
            lines = [ln.format(Ce=s["Ce"], Ce4=4 * s["Ce"], Nl=s["Nl"], Nc=s["Nc"], Ne=s["Ne"], Tes=s["Tes"],
                               Tk=s["Tk"], Tshow=s["Tshow"]) for ln in tr["pt_lines"]]
            if s["clipped"]:
                lines.append(tr["clipped"])
            st.code("\n".join(lines), language=None)

        with st.expander(tr["how"], expanded=False):
            st.markdown(tr["how_text"].format(**T["par"]))
            st.caption(tr["note"])

# ---- figures
if view == 3:
    for f in plots.fig_filters(res, tr):
        st.plotly_chart(f, use_container_width=True, config=PLOT_CFG)
if view == 4:
    st.plotly_chart(plots.fig_stft(res["y45"], res["fs"], tr), use_container_width=True, config=PLOT_CFG)
if view == 5:
    for f in plots.fig_drift(res, tr):
        st.plotly_chart(f, use_container_width=True, config=PLOT_CFG)
    st.plotly_chart(plots.fig_sync(res, tr), use_container_width=True, config=PLOT_CFG)
if view == 6:
    for f in plots.fig_wedges(res, tr):
        st.plotly_chart(f, use_container_width=True, config=PLOT_CFG)
