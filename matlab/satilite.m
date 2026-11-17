function result = satilite(wav_file, rotate180, out_dir, satellite, show_figs, lang)
% =========================================================================
%  NOAA APT DECODER (MATLAB) - Φάση 1
%  Αποδιαμόρφωση AM (Hilbert) + FIR bandpass + IIR lowpass (zero-phase)
%  + STFT + βέλτιστη αντιστάθμιση clock drift + συγχρονισμός Sync-A (xcorr)
%
%  Χρήση:   satilite                       (data/noaa3.wav, με περιστροφή 180)
%           satilite('path/to/file.wav')
%           satilite('file.wav', false)    (χωρίς περιστροφή)
%           satilite('file.wav', true, 'φάκελος_εξόδου', 'NOAA-18')   (+ θερμοκρασία σε °C)
%           satilite(..., show_figs, 'en')   (γλώσσα μηνυμάτων/σχημάτων: 'el' ή 'en')
%  Έξοδος:  struct με εικόνες, ρυθμό γραμμής, παραμέτρους ποιότητας.
%  Απαιτεί: Signal Processing Toolbox. Αποθηκεύει PNG στο ../output/
% =========================================================================

here = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(wav_file)
    wav_file = fullfile(here, '..', 'web', 'sample', 'noaa3.wav');
end
if nargin < 2 || isempty(rotate180), rotate180 = true; end
if nargin < 4, satellite = ''; end
if nargin < 6 || isempty(lang), lang = 'el'; end
L = @(el, en) tr_(lang, el, en);
if nargin < 3 || isempty(out_dir), out_dir = fullfile(here, '..', 'output'); end
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
vis = 'on'; if ~usejava('desktop'), vis = 'off'; end
if nargin >= 5 && ~isempty(show_figs) && ~show_figs, vis = 'off'; end

disp('===================================================');
disp(L(' NOAA APT decoder - Φάση 1', ' NOAA APT decoder'));
disp('===================================================');

%% 1. Διάβασμα
[y, fs] = audioread(wav_file);
y = y(:, 1);
y = y - mean(y);
fprintf(L('1. Αρχείο: %s | fs = %d Hz | διάρκεια = %.1f s\n', '1. File: %s | fs = %d Hz | duration = %.1f s\n'), wav_file, fs, numel(y)/fs);

%% 2. Φίλτρα
% 2a. FIR bandpass N=64, Hamming, 2.4kHz +- 1.2kHz, γραμμική φάση (τg = 32 δείγματα)
N_fir = 64;
b_fir = fir1(N_fir, [1200 3600]/(fs/2), 'bandpass', hamming(N_fir+1));
% 2b. IIR Butterworth 3ης τάξης, fc = 2000 Hz
[b_iir, a_iir] = butter(3, 2000/(fs/2), 'low');
poles = roots(a_iir);
fprintf(L('2. Πόλοι IIR: max|p| = %.4f (<1 => BIBO ευσταθές)\n', '2. IIR poles: max|p| = %.4f (<1 => BIBO stable)\n'), max(abs(poles)));

% FIR: φιλτράρισμα + αντιστάθμιση καθυστέρησης ομάδας N/2
y_bp = filter(b_fir, 1, [y; zeros(N_fir/2, 1)]);
y_bp = y_bp(N_fir/2 + 1 : end);

%% 3. Αποδιαμόρφωση Hilbert + zero-phase lowpass
env = abs(hilbert(y_bp));
env = filtfilt(b_iir, a_iir, env);
disp(L('3. Αποδιαμόρφωση Hilbert + filtfilt ολοκληρώθηκε', '3. Hilbert demodulation + filtfilt done'));

%% 4. Clock drift: βελτιστοποίηση μήκους γραμμής (max κατακόρυφης διακύμανσης)
%  Ονομαστικά T_nom = fs/2 = 5512.5 δείγματα/γραμμή. Επαναδειγματοληψία ώστε
%  κάθε γραμμή να έχει ακριβώς N_line = 5512 δείγματα (ρ = T/5512).
T_nom  = fs/2;
N_line = floor(T_nom);
test_len = min(numel(env), fs*120);
env_test = env(1:test_len);

score = @(T) line_score(env_test, T, N_line);
coarse = -2.0 : 0.05 : 2.0;
sc = arrayfun(@(d) score(T_nom + d), coarse);
[~, i] = max(sc);
fine = coarse(i) - 0.05 : 0.002 : coarse(i) + 0.05;
sf = arrayfun(@(d) score(T_nom + d), fine);
[~, j] = max(sf);
T_true = T_nom + fine(j);
fprintf(L('4. Βέλτιστος ρυθμός γραμμής: %.3f δείγματα/γραμμή (ονομαστικά %.1f)\n', '4. Optimal line rate: %.3f samples/line (nominal %.1f)\n'), T_true, T_nom);

rho = T_true / N_line;
env_fixed = interp1(1:numel(env), env, 1:rho:numel(env), 'linear');
env_fixed(isnan(env_fixed)) = [];
num_lines = floor(numel(env_fixed) / N_line);
img = reshape(env_fixed(1:num_lines*N_line), N_line, num_lines)';

%% 4b. Εκλέπτυνση drift με τον Sync-A (γραμμική παλινδρόμηση θέσης sync ανά γραμμή)
%  Η μετρική max-variance έχει πλατό/θόρυβο. Εδώ μετράμε τη θέση του sync σε
%  κάθε γραμμή: αν ο ρυθμός είναι σωστός, η θέση είναι σταθερή (κλίση = 0).
W = 2080;
tpl = repmat([1 1 -1 -1], 1, 7);
T_var = T_true;                               % ρυθμός από max-variance (πριν την εκλέπτυνση)
T_ref = T_true;
for it = 1:3
    rho = T_ref / N_line;
    ef = interp1(1:numel(env), env, 1:rho:numel(env), 'linear');
    ef(isnan(ef)) = [];
    nl = floor(numel(ef) / N_line);
    M = reshape(ef(1:nl*N_line), N_line, nl)';
    M = interp1(1:N_line, M', linspace(1, N_line, W), 'linear')';
    pos = sync_positions(M, tpl);
    pos = unwrap_cols(pos, W);
    k = (1:numel(pos))';
    ok = abs(pos - movmedian(pos, 101)) < 4;      % απόρριψη γραμμών με θόρυβο
    cf = polyfit(k(ok), pos(ok), 1);              % pixels/γραμμή (ανθεκτική)
    T_ref = T_ref + cf(1) * N_line / W * rho;
    fprintf(L('   sync-refine %d: κλίση = %+.4f px/γραμμή -> T = %.4f\n', '   sync-refine %d: slope = %+.4f px/line -> T = %.4f\n'), it, cf(1), T_ref);
end
T_sync = T_ref;
T_true = T_sync;
sync_jitter = std(pos(ok) - polyval(cf, k(ok)));
fprintf(L('   έγκυρες γραμμές sync: %.1f%%\n', '   valid sync lines: %.1f%%\n'), 100*mean(ok));
fprintf(L('4b. Τελικός ρυθμός (sync): %.4f δείγματα/γραμμή | jitter sync = %.2f px\n', '4b. Final rate (sync): %.4f samples/line | sync jitter = %.2f px\n'), T_true, sync_jitter);
rho = T_true / N_line;
env_fixed = interp1(1:numel(env), env, 1:rho:numel(env), 'linear');
env_fixed(isnan(env_fixed)) = [];
num_lines = floor(numel(env_fixed) / N_line);
img = reshape(env_fixed(1:num_lines*N_line), N_line, num_lines)';

%% 5. Κανονικοποίηση πλάτους σε 2080 pixels/γραμμή
img = interp1(1:N_line, img', linspace(1, N_line, W), 'linear')';
fprintf(L('5. Εικόνα: %d γραμμές x %d pixels\n', '5. Image: %d lines x %d pixels\n'), size(img,1), W);

%% 6. Συγχρονισμός με συσχέτιση με το Sync-A (7 παλμοί 1040 Hz)
%  Sync-A: 7 περίοδοι (2 pixels υψηλά, 2 χαμηλά) = 28 pixels
avg_line = mean(img, 1);
avg_line = avg_line - mean(avg_line);
tpl_pad = zeros(1, W); tpl_pad(1:numel(tpl)) = tpl;
xc = real(ifft(fft(avg_line) .* conj(fft(tpl_pad))));   % κυκλική συσχέτιση
[pk, sync_col] = max(xc);
sync_quality = pk / mean(abs(xc));
fprintf(L('6. Sync-A στη στήλη %d (ποιότητα κορυφής = %.1f)\n', '6. Sync-A at column %d (peak quality = %.1f)\n'), sync_col, sync_quality);
img = circshift(img, -(sync_col - 1), 2);   % το Sync-A ξεκινά στη στήλη 1

% Διόρθωση Doppler: η θέση του sync ακολουθεί παραβολή (το περιθώριο "καμπυλώνει").
% Μετράμε τη θέση ανά γραμμή, την εξομαλύνουμε και μετακινούμε κάθε γραμμή.
pos2 = unwrap_cols(sync_positions(img, tpl), W);
bad  = abs(pos2 - movmedian(pos2, 101)) >= 4;
pos2(bad) = NaN;
pos2 = fillmissing(pos2, 'linear', 'EndValues', 'nearest');
pos_s = smoothdata(pos2, 'gaussian', 101);
doppler_shift = round(pos_s - 1);
for r = 1:size(img, 1)
    img(r, :) = circshift(img(r, :), -doppler_shift(r));
end
fprintf(L('   Doppler: μέγιστη μετατόπιση γραμμής = %d px\n', '   Doppler: max line shift = %d px\n'), max(abs(doppler_shift)));
% Τώρα: [SyncA | SpaceA | Εικόνα A | TelemA | SyncB | SpaceB | Εικόνα B | TelemB]

%% 7. Προσανατολισμός (North-Up) - για πέρασμα προς Βορρά
if rotate180
    img = rot90(img, 2);
    img = circshift(img, 1040, 2);   % κανάλι A ξανά αριστερά
    disp(L('7. Περιστροφή 180 (North-Up) εφαρμόστηκε', '7. 180 deg rotation (North-Up) applied'));
end

%% 7b. Βαθμονόμηση από την τηλεμετρία (16 wedges x 8 γραμμές = πλαίσιο 128 γραμμών)
%  Wedges 1-8: 31,63,95,...,255 (% διαμόρφωσης), wedge 9: 0, wedge 16: ID καναλιού.
%  Επίπεδα αναφοράς: 10.6%..87% του 255 => 31.07,63.02,94.96,126.9,158.86,191.1,228.62,255.
%  (πηγή: aptdec / NOAA KLM). Με αντιστοίχιση counts -> modulation αφαιρούμε κέρδος/offset.
if rotate180
    telA = 6:40;     telB = 1046:1080;
else
    telA = 1001:1035; telB = 2041:2075;
end
calA = wedge_calibrate(img, telA, rotate180);
calB = wedge_calibrate(img, telB, rotate180);
fprintf(L('7b. Τηλεμετρία A: offset γραμμής %d, ID καναλιού = %d, RMS προσαρμογής = %.2f\n', '7b. Telemetry A: line offset %d, channel ID = %d, fit RMS = %.2f\n'), ...
    calA.offset, calA.channel_id, calA.rms);
fprintf(L('    Τηλεμετρία B: offset γραμμής %d, ID καναλιού = %d, RMS προσαρμογής = %.2f\n', '    Telemetry B: line offset %d, channel ID = %d, fit RMS = %.2f\n'), ...
    calB.offset, calB.channel_id, calB.rms);

%% 8. Κανονικοποίηση & εμφάνιση
lo = prctile(img(:), 0.5); hi = prctile(img(:), 99.5);
img_n = min(max((img - lo) / (hi - lo), 0), 1) .^ 0.75;
chA = img_n(:, 1:1040);
chB = img_n(:, 1041:2080);
% Καθαρές εικόνες (χωρίς sync/space/τηλεμετρία): 909 pixels ανά κανάλι.
% Διάταξη APT: Sync(39) Space(47) Εικόνα(909) Τηλεμετρία(45). Μετά από 180 ° αντιστρέφεται.
if rotate180, c0 = 49; else, c0 = 87; end
cropA = chA(:, c0 : c0+908);
cropB = chB(:, c0 : c0+908);
imwrite(cropA, fullfile(out_dir, 'apt_A_crop.png'));
imwrite(cropB, fullfile(out_dir, 'apt_B_crop.png'));
% Ψευδοχρώματα για το θερμικό κανάλι (ψυχρά νέφη = λευκά/κόκκινα)
imwrite(ind2rgb(uint8(round(255*cropB)), turbo(256)), fullfile(out_dir, 'apt_B_color.png'));
% Βαθμονομημένες εικόνες: αντιστοίχιση counts -> % διαμόρφωσης με πίνακα (LUT) από τα wedges,
% ώστε να διορθώνεται και η μη γραμμικότητα (συμπίεση) του δέκτη στις υψηλές τιμές.
rawA = img(:, 1:1040); rawB = img(:, 1041:2080);
calimgA = apply_wedge_lut(rawA(:, c0:c0+908), calA);
calimgB = apply_wedge_lut(rawB(:, c0:c0+908), calB);
% ---- Θερμοκρασία (°C) για το θερμικό κανάλι, αν δηλωθεί δορυφόρος ----
temp = [];
if ~isempty(satellite) && ismember(calB.channel_id, [4 5])
    if rotate180, spc = 962:1000; else, spc = 44:82; end      % στήλες Space-data μέσα στο μπλοκ του B
    Cs_meas = median(median(apply_wedge_lut(rawB(:, spc), calB) * 255));
    temp = thermal_calibrate(calimgB * 255, calB, satellite, calB.channel_id, Cs_meas);
    fprintf(L('7c. Θερμοκρασία (%s, κανάλι %d): Tbb = %.2f K | Cs(μετρημένο) = %.1f | Cb = %.1f | Nbb = %.3f\n', '7c. Temperature (%s, channel %d): Tbb = %.2f K | Cs(measured) = %.1f | Cb = %.1f | Nbb = %.3f\n'), ...
        satellite, calB.channel_id, temp.Tbb, Cs_meas, temp.Cb, temp.Nbb);
    fprintf(L('    Εύρος θερμοκρασίας εικόνας (1%%-99%%): %.1f ... %.1f °C\n', '    Image temperature range (1%%-99%%): %.1f ... %.1f °C\n'), ...
        prctile(temp.TC(:), 1), prctile(temp.TC(:), 99));
    imwrite(uint8(min(max((temp.TC + 100) / 160, 0), 1) * 255), fullfile(out_dir, 'apt_B_temp_gray.png'));
    f7 = lightfig('Name','Temperature','Visible',vis,'Color','w','Position',[100 100 600 1000]);
    imagesc(temp.TC, [-70 60]); axis image off; colormap(turbo(256)); cb = colorbar; cb.Label.String = L('Θερμοκρασία (°C)', 'Temperature (°C)');
    title(sprintf(L('%s - κανάλι %d (θερμικό IR)', '%s - channel %d (thermal IR)'), satellite, calB.channel_id));
    exportgraphics(f7, fullfile(out_dir, 'fig7_temperature.png'), 'Resolution', 150);
end
imwrite(calimgA, fullfile(out_dir, 'apt_A_cal.png'));
imwrite(calimgB, fullfile(out_dir, 'apt_B_cal.png'));
imwrite(ind2rgb(uint8(round(255*calimgB)), turbo(256)), fullfile(out_dir, 'apt_B_cal_color.png'));

imwrite(img_n, fullfile(out_dir, 'apt_full.png'));
imwrite(chA,   fullfile(out_dir, 'apt_channel_A.png'));
imwrite(chB,   fullfile(out_dir, 'apt_channel_B.png'));

%% 9. Διαγράμματα για την αναφορά
% (α) Z-plane + Bode
f1 = lightfig('Name','Z-plane & Bode','Visible',vis,'Color','w','Position',[100 100 1100 800]);
subplot(2,2,1); zplane(b_iir, a_iir); title(L('(α) Z-plane IIR Butterworth 3ης τάξης', '(a) Z-plane IIR Butterworth, order 3'));
subplot(2,2,2); zplane(b_fir, 1);     title(L('(β) Z-plane FIR bandpass (64 μηδενικά)', '(b) Z-plane FIR bandpass (64 zeros)'));
[H1, w1] = freqz(b_iir, a_iir, 2048, fs);
subplot(2,2,3); plot(w1, 20*log10(abs(H1)), 'LineWidth', 1.4); grid on; ylim([-80 5]);
xlabel(L('Συχνότητα (Hz)', 'Frequency (Hz)')); ylabel('|H| (dB)'); title(L('(γ) IIR Lowpass, f_c = 2000 Hz', '(c) IIR lowpass, f_c = 2000 Hz'));
[H2, w2] = freqz(b_fir, 1, 2048, fs);
subplot(2,2,4); plot(w2, 20*log10(abs(H2)), 'LineWidth', 1.4); grid on; ylim([-80 5]);
xlabel(L('Συχνότητα (Hz)', 'Frequency (Hz)')); ylabel('|H| (dB)'); title(L('(δ) FIR bandpass 2.4kHz \pm 1.2kHz', '(d) FIR bandpass 2.4kHz \pm 1.2kHz'));
exportgraphics(f1, fullfile(out_dir, 'fig1_filters.png'), 'Resolution', 150);

% (β) STFT: 45 s, Hamming 512, επικάλυψη 75% (384)
f2 = lightfig('Name','STFT','Visible',vis,'Color','w','Position',[100 100 900 500]);
seg = y(1 : min(numel(y), 45*fs));
spectrogram(seg, hamming(512), 384, 1024, fs, 'yaxis');
ylim([0 5]); title(L('STFT - φέρουσα 2400 Hz', 'STFT - 2400 Hz subcarrier'));
exportgraphics(f2, fullfile(out_dir, 'fig2_stft.png'), 'Resolution', 150);

% (γ) Βελτιστοποίηση drift
f3 = lightfig('Name','Drift','Visible',vis,'Color','w','Position',[100 100 900 400]);
subplot(1,2,1); plot(T_nom + coarse, sc, 'LineWidth', 1.2); grid on; xline(T_var, 'r--', 'LineWidth', 1.2);
xlabel(L('T (δείγματα/γραμμή)', 'T (samples/line)')); ylabel('max Var'); title(L('Χονδρική σάρωση', 'Coarse scan'));
subplot(1,2,2); plot(T_nom + fine, sf, 'LineWidth', 1.2); grid on; hold on;
h1 = xline(T_var, 'r--', 'LineWidth', 1.2); h2 = xline(T_sync, '-', 'Color', [0 0.6 0], 'LineWidth', 1.5);
xlim([T_nom + fine(1), max(T_sync, T_nom + fine(end)) + 0.01]);
legend([h1 h2], {'max-variance', L('τελικό (sync)', 'final (sync)')}, 'Location', 'southeast');
xticks(5512.36:0.04:5512.48); xtickformat('%.2f');
xlabel(L('T (δείγματα/γραμμή)', 'T (samples/line)')); title(L('Λεπτομερής σάρωση', 'Fine scan'));
exportgraphics(f3, fullfile(out_dir, 'fig3_drift.png'), 'Resolution', 150);

% (δ) Τελική εικόνα
f4 = lightfig('Name','NOAA APT','Visible',vis,'Color','w','Position',[100 100 1000 900]);
imshow(img_n); title(sprintf(L('NOAA APT - %d γραμμές | T = %.3f | sync q = %.1f', 'NOAA APT - %d lines | T = %.3f | sync q = %.1f'), ...
    size(img_n,1), T_true, sync_quality));
exportgraphics(f4, fullfile(out_dir, 'fig4_image.png'), 'Resolution', 150);

f6 = lightfig('Name','Wedges','Visible',vis,'Color','w','Position',[100 100 1000 380]);
subplot(1,2,1); plot(calB.t(1:min(384,end)), 'LineWidth', 1); grid on;
xlabel(L('Γραμμή', 'Line')); ylabel(L('Τηλεμετρία (counts)', 'Telemetry (counts)')); title(L('Τηλεμετρία καναλιού B (3 πλαίσια)', 'Channel B telemetry (3 frames)'));
subplot(1,2,2); lv = wedge_levels();
plot(lv, calB.W(1:9), 'o', 'LineWidth', 1.2); hold on; xx = [0 255];
plot(xx, calB.a*xx + calB.b, 'r--'); plot(lv, calA.W(1:9), 's'); grid on;
xlabel(L('Διαμόρφωση wedge (%)', 'Wedge modulation (%)')); ylabel('Counts'); legend('B (wedges)', L('Γραμμική προσαρμογή', 'Linear fit'), 'A (wedges)', 'Location', 'northwest');
title(sprintf(L('Wedges vs επίπεδα αναφοράς (RMS = %.4f)', 'Wedges vs reference levels (RMS = %.4f)'), calB.rms));
exportgraphics(f6, fullfile(out_dir, 'fig6_wedges.png'), 'Resolution', 150);

f5 = lightfig('Name','Sync','Visible',vis,'Color','w','Position',[100 100 900 350]);
plot(k(ok), pos(ok), '.', 'MarkerSize', 4); hold on; plot(polyval(cf, k), 'r', 'LineWidth', 1.2); grid on;
xlabel(L('Γραμμή', 'Line')); ylabel(L('Θέση Sync-A (pixel)', 'Sync-A position (pixel)')); title(L('Σταθερότητα συγχρονισμού μετά τη διόρθωση drift', 'Sync stability after drift correction'));
exportgraphics(f5, fullfile(out_dir, 'fig5_sync.png'), 'Resolution', 150);

result = struct('CeB', calimgB * 255, 'temp', temp, 'calA', calA, 'calB', calB, 'sync_jitter', sync_jitter, 'img', img_n, 'chA', chA, 'chB', chB, 'T_line', T_true, ...
    'sync_col', sync_col, 'sync_quality', sync_quality, 'max_pole', max(abs(poles)));
disp(L('-> ΟΛΟΚΛΗΡΩΘΗΚΕ. Αποτελέσματα στον φάκελο output/', '-> DONE. Results saved in the output folder'));
end

function s = line_score(env_test, T, N_line)
    rho = T / N_line;
    x = interp1(1:numel(env_test), env_test, 1:rho:numel(env_test), 'linear');
    x(isnan(x)) = [];
    nl = floor(numel(x) / N_line);
    M = reshape(x(1:nl*N_line), N_line, nl)';
    s = max(var(M, 0, 1));
end

function pos = sync_positions(M, tpl)
    W = size(M, 2);
    tp = zeros(1, W); tp(1:numel(tpl)) = tpl;
    Tf = conj(fft(tp));
    X = M - mean(M, 2);
    xc = real(ifft(fft(X, [], 2) .* Tf, [], 2));
    [~, pos] = max(xc, [], 2);
    pos = pos(:);
end

function p = unwrap_cols(pos, W)
    % κεντράρισμα γύρω από τη διάμεσο ώστε να μη σπάει στο όριο της στήλης
    m = median(pos);
    d = mod(pos - m + W/2, W) - W/2;
    p = m + d;
end

function f = lightfig(varargin)
    f = figure(varargin{:});
    try, theme(f, 'light'); catch, end
end

function c = wedge_calibrate(img, cols, flipped)
    t = median(img(:, cols), 2);                 % τηλεμετρία ανά γραμμή
    if flipped, t = flipud(t); end               % η περιστροφή 180 αντιστρέφει τη σειρά των wedges
    lv = wedge_levels();
    tpl = repelem(lv, 8)';                       % wedges 1-9 (72 γραμμές)
    nL = numel(t); nf = floor(nL / 128) - 1;
    best = -inf; off = 0;
    for o = 0:127
        sc = 0;
        for f = 0:nf-1
            seg = t(o + f*128 + (1:72));
            sc = sc + corrcoef_(seg, tpl);
        end
        sc = sc / nf;
        if sc > best, best = sc; off = o; end
    end
    nfr = floor((nL - off) / 128);
    Wf = zeros(nfr, 16);
    for f = 0:nfr-1
        for k = 1:16
            idx = off + f*128 + (k-1)*8 + (2:7);  % παράλειψη γραμμών μετάβασης
            Wf(f+1, k) = mean(t(idx));
        end
    end
    W = median(Wf, 1);
    pf = polyfit(lv, W(1:9), 1);                 % counts = a*mod + b
    a = pf(1); b = pf(2);
    rms = sqrt(mean((W(1:9) - polyval(pf, lv)).^2));
    modID = (W(16) - b) / a;                     % wedge 16 -> % διαμόρφωσης
    [~, chid] = min(abs(modID - [31 63 95 127 159 191]));
    c = struct('t', t, 'offset', off, 'W', W, 'a', a, 'b', b, 'rms', rms, ...
               'channel_id', chid, 'score', best);
end

function r = corrcoef_(x, y)
    x = x(:) - mean(x); y = y(:) - mean(y);
    r = (x' * y) / (norm(x) * norm(y) + eps);
end

function out = apply_wedge_lut(x, c)
    lv = [0 wedge_levels(1:8)];
    w  = [c.W(9) c.W(1:8)];                      % counts των wedges 9,1..8 (αύξουσα σειρά)
    out = interp1(w, lv, x, 'linear', 'extrap') / 255;
    out = min(max(out, 0), 1);
end

function lv = wedge_levels(idx)
    % Επίπεδα διαμόρφωσης wedges 1-8 και wedge 9 (μηδενική διαμόρφωση), κλίμακα 0-255
    lv = [31.07 63.02 94.96 126.9 158.86 191.1 228.62 255.0 0.0];
    if nargin > 0, lv = lv(idx); end
end

function t = thermal_calibrate(Ce255, cal, satellite, ch, Cs255)
    % Θερμοκρασία από APT τηλεμετρία (αλγόριθμος NOAA KLM User's Guide, όπως στο aptdec).
    % Ce255: pixels σε κλίμακα 0-255 μετά τη βαθμονόμηση wedges. Τελεστές 10-bit = 4*counts.
    k = sat_constants(satellite);
    ci = ch - 3;                                   % 4 -> 1, 5 -> 2
    Vc = k.rad(ci,1); A = k.rad(ci,2); B = k.rad(ci,3);
    C1 = 1.1910427e-5; C2 = 1.4387752;
    tele = [];                                     % τηλεμετρία σε κλίμακα 0-255 (wedges 10-15)
    lv = [0 wedge_levels(1:8)]; w = [cal.W(9) cal.W(1:8)];
    tele = min(max(interp1(w, lv, cal.W, 'linear', 'extrap'), 0), 255);
    Tprt = zeros(1, 4);
    for n = 1:4
        x = tele(9 + n) * 4;                       % wedges 10-13 = PRT 1-4
        Tprt(n) = k.prt(n,1)*x^2 + k.prt(n,2)*x + k.prt(n,3);
    end
    Tbb  = mean(Tprt);
    Tbbs = A + Tbb * B;
    Nbb  = C1 * Vc^3 / (exp(C2 * Vc / Tbbs) - 1);
    Cs   = Cs255 * 4;
    Cb   = tele(15) * 4;                           % back scan (wedge 15)
    Ns   = k.cor(ci,1); q = k.cor(ci,2:4);
    Nl   = Ns + (Nbb - Ns) * (Cs - 4*Ce255) / (Cs - Cb);
    Nc   = q(1)*Nl.^2 + q(2)*Nl + q(3);
    Ne   = max(Nl + Nc, 1e-3);                     % αποφυγή αρνητικής ακτινοβολίας (θόρυβος)
    Tes  = C2 * Vc ./ log(C1 * Vc^3 ./ Ne + 1);
    par = struct('Vc', Vc, 'A', A, 'B', B, 'C1', C1, 'C2', C2, 'Ns', Ns, 'q', q, 'Nbb', Nbb, 'Cs', Cs, 'Cb', Cb, 'Tbbs', Tbbs, 'prt', k.prt);
    t = struct('TC', min(max((Tes - A) / B - 273.15, -100), 60), 'Tbb', Tbb, 'Tprt', Tprt, 'Nbb', Nbb, 'Cs', Cs/4, 'Cb', Cb/4, 'par', par);
end

function k = sat_constants(name)
    % Σταθερές βαθμονόμησης (NOAA KLM User's Guide, μέσω aptdec/src/calibration.c, GPL-2.0).
    % prt: [a b c] τετραγωνικό για PRT1-4. rad: [vc A B] (κανάλια 4,5). cor: [Ns a b c].
    switch upper(strrep(name, '-', ''))
        case 'NOAA15'
            k.prt = [1.36328e-06 0.051045 276.60157; 1.47266e-06 0.050909 276.62531; ...
                     1.47656e-06 0.050907 276.67413; 1.47656e-06 0.050966 276.59258];
            k.rad = [925.4075 0.337810 0.998719; 839.8979 0.304558 0.999024];
            k.cor = [-4.50 0.0004524 -0.0932 4.76; -3.61 0.0002811 -0.0659 3.83];
        case 'NOAA18'
            k.prt = [1.657e-06 0.05090 276.601; 1.482e-06 0.05101 276.683; ...
                     1.313e-06 0.05117 276.565; 1.484e-06 0.05103 276.615];
            k.rad = [928.1460 0.436645 0.998607; 833.2532 0.253179 0.999057];
            k.cor = [-5.53 0.00052337 -0.11069 5.82; -2.22 0.00017715 -0.04360 2.67];
        case 'NOAA19'
            k.prt = [1.405783e-06 0.051111 276.6067; 1.496037e-06 0.051090 276.6119; ...
                     1.496990e-06 0.051033 276.6311; 1.493110e-06 0.051058 276.6268];
            k.rad = [928.9 0.53959 0.998534; 831.9 0.36064 0.998913];
            k.cor = [-5.49 0.00054668 -0.11187 5.70; -3.39 0.00024985 -0.05991 3.58];
        otherwise
            error('Άγνωστος δορυφόρος: %s (NOAA-15, NOAA-18, NOAA-19)', name);
    end
end

function out = tr_(lang, el, en)
    if strcmpi(lang, 'en'), out = en; else, out = el; end
end
