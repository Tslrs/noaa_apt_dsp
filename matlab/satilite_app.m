function ui = satilite_app(varargin)
% =========================================================================
%  NOAA APT DECODER - Γραφική εφαρμογή / Graphical app (MATLAB)
%  Τρέξε / Run:  satilite_app            (ή satilite_app('Lang','en'))
%  1) Διάλεξε wav (11025 Hz mono, λήψη NOAA APT) / Choose a wav recording
%  2) Διάλεξε δορυφόρο για θερμοκρασία °C / Choose a satellite for °C
%  3) Αποκωδικοποίηση / Decode (~1-2 min for 15 min of signal)
%  Στην καρτέλα Θερμοκρασία: κλικ πάνω στην εικόνα => θερμοκρασία στο σημείο
%  + ανάλυση του υπολογισμού / Temperature tab: click the map => temperature
%  at that point + step-by-step calculation.
%  Χρησιμοποιεί το satilite.m (ίδιος φάκελος). Απαιτεί Signal Processing Toolbox.
% =========================================================================
%#function satilite
p = inputParser;
addParameter(p, 'Visible', 'on');
addParameter(p, 'Lang', 'el');
parse(p, varargin{:});

here = fileparts(mfilename('fullpath'));
addpath(here);

lang = lower(p.Results.Lang);
S = strs(lang);
state = struct('wav', '', 'outdir', '', 'result', [], 'sat', '', 'marker', gobjects(0), 'pins', zeros(0, 2));

fig = uifigure('Name', 'NOAA APT Decoder', 'Position', [60 60 1320 800], ...
    'Visible', p.Results.Visible);
root = uigridlayout(fig, [1 2]);
root.ColumnWidth = {290, '1x'};
root.Padding = [8 8 8 8];

%% ---- Αριστερά / Left ----
left = uigridlayout(root, [13 1]);
left.RowHeight = {32, 28, 22, 40, 28, 22, 28, 22, 36, 36, 28, 150, '1x'};
left.Padding = [0 0 0 0];
titleLbl  = uilabel(left, 'Text', S.title, 'FontSize', 18, 'FontWeight', 'bold');
langDD    = uidropdown(left, 'Items', {'Ελληνικά', 'English'}, 'ItemsData', {'el', 'en'}, ...
    'Value', lang, 'ValueChangedFcn', @onLang);
fileHdr   = uilabel(left, 'Text', S.file_hdr);
fileLbl   = uilabel(left, 'Text', S.file_none, 'WordWrap', 'on', 'FontColor', [0.5 0.5 0.5]);
pickBtn   = uibutton(left, 'Text', S.pick, 'ButtonPushedFcn', @onPick);
satHdr    = uilabel(left, 'Text', S.sat_hdr);
satDD     = uidropdown(left, 'Items', S.sat_items, 'ItemsData', {'', 'NOAA-15', 'NOAA-18', 'NOAA-19'}, 'Value', '');
orientHdr = uilabel(left, 'Text', S.orient_hdr);
rotCB     = uicheckbox(left, 'Text', S.rot, 'Value', true, 'WordWrap', 'on');
runBtn    = uibutton(left, 'Text', S.run, 'FontWeight', 'bold', ...
    'BackgroundColor', [0.20 0.45 0.80], 'FontColor', [1 1 1], 'ButtonPushedFcn', @onRun);
saveBtn   = uibutton(left, 'Text', S.save, 'Enable', 'off', 'ButtonPushedFcn', @onSave);
summary   = uitextarea(left, 'Editable', 'off', 'Value', {S.summary_init});
hintLbl   = uilabel(left, 'Text', S.hint, 'WordWrap', 'on', 'FontColor', [0.5 0.5 0.5], ...
    'VerticalAlignment', 'top');

%% ---- Δεξιά: καρτέλες / Right: tabs ----
tabs = uitabgroup(root);
keys  = {'image', 'chA', 'chB', 'temp', 'filters', 'stft', 'drift', 'sync', 'wedges', 'log'};
files = {'apt_full.png', 'apt_A_cal.png', 'apt_B_cal_color.png', '', 'fig1_filters.png', ...
         'fig2_stft.png', 'fig3_drift.png', 'fig5_sync.png', 'fig6_wedges.png', ''};
tabObj = gobjects(1, numel(keys));
imgs   = gobjects(1, numel(keys));
for i = 1:numel(keys)
    tabObj(i) = uitab(tabs, 'Title', S.tabs{i});
    switch keys{i}
        case 'temp'
            tg = uigridlayout(tabObj(i), [3 2]);
            tg.RowHeight = {30, '1x', 250}; tg.ColumnWidth = {'1x', '1x'}; tg.Padding = [4 4 4 4];
            hoverLbl = uilabel(tg, 'Text', S.hover_idle, 'FontSize', 15, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center');
            hoverLbl.Layout.Row = 1; hoverLbl.Layout.Column = [1 2];
            ax = uiaxes(tg); ax.Layout.Row = 2; ax.Layout.Column = [1 2];
            ax.Interactions = []; ax.Toolbar.Visible = 'off';
            ptArea = uitextarea(tg, 'Editable', 'off', 'FontName', 'Menlo', 'FontSize', 11);
            ptArea.Layout.Row = 3; ptArea.Layout.Column = 1;
            howArea = uitextarea(tg, 'Editable', 'off', 'FontName', 'Menlo', 'FontSize', 11);
            howArea.Layout.Row = 3; howArea.Layout.Column = 2;
        case 'log'
            gl = uigridlayout(tabObj(i), [1 1]);
            logArea = uitextarea(gl, 'Editable', 'off', 'FontName', 'Menlo');
        otherwise
            g = uigridlayout(tabObj(i), [1 1]); g.Padding = [4 4 4 4];
            imgs(i) = uiimage(g, 'ScaleMethod', 'fit');
    end
end
showTempPlaceholder();
ptArea.Value = {S.click_hint};
howArea.Value = splitlines(string(howText(lang, [])));
fig.WindowButtonMotionFcn = @onMotion;

ui = struct('fig', fig, 'run', @runDecode, 'pick', @setFile, 'click', @pointInfo, 'hover', @hoverInfo, 'remove', @removePin, 'clear', @clearPins, 'getPins', @getPins, 'hoverLbl', hoverLbl, ...
            'setLang', @setLang, 'summary', summary, 'images', imgs, 'sat', satDD, ...
            'rot', rotCB, 'ptArea', ptArea, 'howArea', howArea, 'tabs', tabs, ...
            'tabObj', tabObj, 'langDD', langDD, 'ax', ax);

%% ---- callbacks ----
    function onLang(~, ~)
        setLang(langDD.Value);
    end

    function setLang(newLang)
        lang = newLang; S = strs(lang);
        if ~strcmp(langDD.Value, lang), langDD.Value = lang; end
        titleLbl.Text = S.title; fileHdr.Text = S.file_hdr; pickBtn.Text = S.pick;
        satHdr.Text = S.sat_hdr; orientHdr.Text = S.orient_hdr; rotCB.Text = S.rot;
        runBtn.Text = S.run; saveBtn.Text = S.save; hintLbl.Text = S.hint;
        if isempty(state.wav), fileLbl.Text = S.file_none; end
        cur = satDD.Value;
        satDD.Items = S.sat_items; satDD.ItemsData = {'', 'NOAA-15', 'NOAA-18', 'NOAA-19'}; satDD.Value = cur;
        for k = 1:numel(tabObj), tabObj(k).Title = S.tabs{k}; end
        hoverLbl.Text = S.hover_idle;
        if isempty(state.result)
            summary.Value = {S.summary_init};
            ptArea.Value = {S.click_hint};
            howArea.Value = splitlines(string(howText(lang, [])));
            showTempPlaceholder();
        else
            summary.Value = buildSummary(state.result, state.sat);
            howArea.Value = splitlines(string(howText(lang, parOf(state.result))));
            ptArea.Value = {S.click_hint};
            showTemp();
        end
    end

    function onPick(~, ~)
        [f, d] = uigetfile({'*.wav', S.pick_filter}, S.pick_title);
        figure(fig);
        if isequal(f, 0), return; end
        setFile(fullfile(d, f));
    end

    function setFile(path)
        state.wav = path;
        fileLbl.Text = path;
        fileLbl.FontColor = [0 0 0];
    end

    function onRun(~, ~)
        runDecode(state.wav, rotCB.Value, satDD.Value);
    end

    function runDecode(wav, rot, sat)
        if isempty(wav) || ~isfile(wav)
            uialert(fig, S.err_nofile, S.err_nofile_title);
            return;
        end
        outdir = fullfile(tempdir, ['satilite_app_' datestr(now, 'yyyymmdd_HHMMSS')]);
        runBtn.Enable = 'off'; saveBtn.Enable = 'off';
        dlg = uiprogressdlg(fig, 'Title', S.progress_title, 'Indeterminate', 'on', 'Message', S.progress_msg);
        cleanup = onCleanup(@() cleanupRun(dlg)); %#ok<NASGU>
        r = []; logtxt = ''; curLang = lang;
        try
            logtxt = evalc('r = satilite(wav, rot, outdir, sat, false, curLang);');
        catch err
            uialert(fig, err.message, S.err_title);
            logArea.Value = splitlines(string(err.message));
            return;
        end
        state.result = r; state.outdir = outdir; state.sat = sat;
        logArea.Value = splitlines(string(logtxt));
        for k = 1:numel(files)
            if isempty(files{k}), continue; end
            f = fullfile(outdir, files{k});
            if isfile(f), imgs(k).ImageSource = f; else, imgs(k).ImageSource = ''; end
        end
        summary.Value = buildSummary(r, sat);
        howArea.Value = splitlines(string(howText(lang, parOf(r))));
        ptArea.Value = {S.click_hint};
        showTemp();
        tabs.SelectedTab = tabObj(1);
        saveBtn.Enable = 'on';
    end

    function cleanupRun(dlg)
        if isvalid(dlg), close(dlg); end
        if isvalid(runBtn), runBtn.Enable = 'on'; end
    end

    function par = parOf(r)
        par = [];
        if ~isempty(r.temp), par = r.temp.par; end
    end

    function showTempPlaceholder()
        cla(ax); state.marker = gobjects(0); state.pins = zeros(0, 2);
        ax.XTick = []; ax.YTick = []; ax.Box = 'on';
        text(ax, 0.5, 0.5, S.no_temp, 'Units', 'normalized', 'HorizontalAlignment', 'center', ...
            'FontSize', 14, 'Color', [0.5 0.5 0.5]);
        title(ax, '');
    end

    function showTemp()
        r = state.result;
        if isempty(r.temp), showTempPlaceholder(); return; end
        cla(ax); state.marker = gobjects(0); state.pins = zeros(0, 2);
        h = imagesc(ax, r.temp.TC, [-70 60]);
        h.ButtonDownFcn = @(~, ev) onClickImage(ev);
        ax.ButtonDownFcn = @(~, ev) onClickImage(ev);
        axis(ax, 'image'); colormap(ax, turbo(256));
        ax.XTickMode = 'auto'; ax.YTickMode = 'auto';
        cb = colorbar(ax); cb.Label.String = S.temp_label;
        title(ax, sprintf(S.temp_title, state.sat));
        xlabel(ax, S.col); ylabel(ax, S.row);
    end

    function onMotion(~, ~)
        % Ζωντανή θερμοκρασία κάτω από τον δείκτη (χωρίς κλικ)
        if isempty(state.result) || isempty(state.result.temp), return; end
        if tabs.SelectedTab ~= tabObj(4), return; end
        cp = ax.CurrentPoint;
        x = cp(1, 1); y = cp(1, 2);
        if x < ax.XLim(1) || x > ax.XLim(2) || y < ax.YLim(1) || y > ax.YLim(2)
            hoverLbl.Text = S.hover_idle; return;
        end
        hoverInfo(round(y), round(x));
    end

    function hoverInfo(row, col)
        r = state.result;
        if isempty(r) || isempty(r.temp), return; end
        [nr, nc] = size(r.temp.TC);
        row = min(max(row, 1), nr); col = min(max(col, 1), nc);
        hoverLbl.Text = sprintf(S.hover_fmt, row, col, r.temp.TC(row, col));
    end

    function onClickImage(ev)
        pt = ev.IntersectionPoint;
        if strcmp(fig.SelectionType, 'alt')
            removePin(round(pt(2)), round(pt(1)));     % δεξί κλικ (ή ctrl+κλικ)
        else
            pointInfo(round(pt(2)), round(pt(1)));     % αριστερό κλικ
        end
    end

    function pointInfo(row, col)
        % Προσθήκη καρφίτσας στο (row, col)
        r = state.result;
        if isempty(r) || isempty(r.temp), return; end
        [nr, nc] = size(r.temp.TC);
        row = min(max(row, 1), nr); col = min(max(col, 1), nc);
        state.pins(end+1, :) = [row col];
        redrawPins(r.temp.par, r);
    end

    function removePin(row, col)
        % Αφαίρεση της πλησιέστερης καρφίτσας (εντός 40 pixel)
        r = state.result;
        if isempty(r) || isempty(r.temp) || isempty(state.pins), return; end
        d = hypot(state.pins(:,1) - row, state.pins(:,2) - col);
        [dm, k] = min(d);
        if dm > 40, return; end
        state.pins(k, :) = [];
        redrawPins(r.temp.par, r);
    end

    function p = getPins()
        p = state.pins;
    end

    function clearPins()
        state.pins = zeros(0, 2);
        r = state.result;
        if ~isempty(r) && ~isempty(r.temp), redrawPins(r.temp.par, r); end
    end

    function redrawPins(par, r)
        if ~isempty(state.marker), delete(state.marker(isvalid(state.marker))); end
        state.marker = gobjects(0);
        hold(ax, 'on');
        for k = 1:size(state.pins, 1)
            rr = state.pins(k, 1); cc = state.pins(k, 2);
            state.marker(end+1) = plot(ax, cc, rr, 'wo', 'MarkerSize', 13, 'LineWidth', 2.5, 'HitTest', 'off');
            state.marker(end+1) = plot(ax, cc, rr, 'k+', 'MarkerSize', 13, 'LineWidth', 1.5, 'HitTest', 'off');
            state.marker(end+1) = text(ax, cc + 14, rr, num2str(k), 'Color', 'w', 'FontWeight', 'bold', ...
                'FontSize', 13, 'BackgroundColor', [0 0 0 0.55], 'HitTest', 'off', 'Margin', 1);
        end
        hold(ax, 'off');
        if isempty(state.pins)
            ptArea.Value = {S.click_hint};
        else
            ptArea.Value = splitlines(string(pinsText(lang, par, r.CeB, r.temp.TC, state.pins)));
        end
    end

    function lines = buildSummary(r, sat)
        lines = {sprintf(S.s_lines, size(r.img, 1)), ...
                 sprintf(S.s_rate, r.T_line), ...
                 sprintf(S.s_sync, r.sync_quality), ...
                 sprintf(S.s_chA, r.calA.channel_id), ...
                 sprintf(S.s_chB, r.calB.channel_id)};
        if r.sync_quality < 8
            lines{end+1} = S.s_weak;
        end
        if ~isempty(r.temp)
            lines{end+1} = sprintf(S.s_temp, sat, prctile(r.temp.TC(:), 1), prctile(r.temp.TC(:), 99));
        elseif isempty(sat)
            lines{end+1} = S.s_notemp;
        end
    end

    function onSave(~, ~)
        d = uigetdir(pwd, S.savedir_title);
        figure(fig);
        if isequal(d, 0), return; end
        copyfile(state.outdir, d);
        uialert(fig, [S.saved_msg d], S.saved_title, 'Icon', 'success');
    end
end

%% ================== Κείμενα / Texts ==================
function S = strs(lang)
if strcmpi(lang, 'en')
    S.title = 'NOAA APT Decoder';
    S.file_hdr = '1. Audio file (wav)';
    S.file_none = '(no file selected)';
    S.pick = 'Choose file...';
    S.sat_hdr = '2. Satellite (for °C)';
    S.sat_items = {'Unknown (no °C)', 'NOAA-15', 'NOAA-18', 'NOAA-19'};
    S.orient_hdr = '3. Orientation';
    S.rot = 'Northbound pass (rotate 180°)';
    S.run = 'Decode';
    S.save = 'Save results...';
    S.summary_init = 'Results will appear here.';
    S.hint = 'If the image is upside down, untick the orientation box. Figures are drawn in the language used at decode time.';
    S.tabs = {'Image', 'Channel A', 'Channel B (colour)', 'Temperature', 'Filters', 'STFT', ...
              'Drift', 'Sync/Doppler', 'Wedges', 'Log'};
    S.progress_title = 'Decoding'; S.progress_msg = 'Processing signal... (1-2 minutes)';
    S.err_title = 'Decoding error'; S.err_nofile_title = 'No file';
    S.err_nofile = 'Choose a valid wav file first.';
    S.saved_title = 'Done'; S.saved_msg = 'Saved to: ';
    S.pick_title = 'Choose a NOAA APT recording'; S.pick_filter = 'WAV audio files';
    S.savedir_title = 'Save folder';
    S.s_lines = 'Lines: %d'; S.s_rate = 'Line rate: %.3f samples'; S.s_sync = 'Sync quality: %.1f';
    S.s_chA = 'Channel A: ID %d'; S.s_chB = 'Channel B: ID %d';
    S.s_weak = 'WARNING: weak sync - there may be no APT signal.';
    S.s_temp = 'Temperature (%s): %.0f ... %.0f °C (1%%-99%%)';
    S.s_notemp = 'No °C: no satellite selected.';
    S.click_hint = 'Left-click the map to place a point (see the step-by-step calculation). Right-click near a point to remove it.';
    S.hover_idle = 'Move the mouse over the map to read the temperature';
    S.hover_fmt = 'Line %d, column %d:  %.1f °C';
    S.no_temp = 'Select a satellite (NOAA-15/18/19) and decode to get temperatures.';
    S.temp_label = 'Temperature (°C)'; S.temp_title = '%s - channel B (thermal IR)';
    S.col = 'Pixel column'; S.row = 'Line';
else
    S.title = 'NOAA APT Decoder';
    S.file_hdr = '1. Αρχείο ήχου (wav)';
    S.file_none = '(δεν έχει επιλεγεί αρχείο)';
    S.pick = 'Επιλογή αρχείου...';
    S.sat_hdr = '2. Δορυφόρος (για °C)';
    S.sat_items = {'Άγνωστος (χωρίς °C)', 'NOAA-15', 'NOAA-18', 'NOAA-19'};
    S.orient_hdr = '3. Προσανατολισμός';
    S.rot = 'Πέρασμα προς Βορρά (περιστροφή 180°)';
    S.run = 'Αποκωδικοποίηση';
    S.save = 'Αποθήκευση αποτελεσμάτων...';
    S.summary_init = 'Τα αποτελέσματα θα εμφανιστούν εδώ.';
    S.hint = 'Αν η εικόνα είναι ανάποδα, ξετσέκαρε τον προσανατολισμό. Τα σχήματα σχεδιάζονται στη γλώσσα που ήταν επιλεγμένη κατά την αποκωδικοποίηση.';
    S.tabs = {'Εικόνα', 'Κανάλι A', 'Κανάλι B (χρώμα)', 'Θερμοκρασία', 'Φίλτρα', 'STFT', ...
              'Drift', 'Sync/Doppler', 'Wedges', 'Καταγραφή'};
    S.progress_title = 'Αποκωδικοποίηση'; S.progress_msg = 'Επεξεργασία σήματος... (1-2 λεπτά)';
    S.err_title = 'Σφάλμα αποκωδικοποίησης'; S.err_nofile_title = 'Χωρίς αρχείο';
    S.err_nofile = 'Διάλεξε πρώτα ένα έγκυρο αρχείο wav.';
    S.saved_title = 'Έτοιμο'; S.saved_msg = 'Αποθηκεύτηκαν στο: ';
    S.pick_title = 'Διάλεξε λήψη NOAA APT'; S.pick_filter = 'Αρχεία ήχου WAV';
    S.savedir_title = 'Φάκελος αποθήκευσης';
    S.s_lines = 'Γραμμές: %d'; S.s_rate = 'Ρυθμός γραμμής: %.3f δείγματα'; S.s_sync = 'Ποιότητα sync: %.1f';
    S.s_chA = 'Κανάλι A: ID %d'; S.s_chB = 'Κανάλι B: ID %d';
    S.s_weak = 'ΠΡΟΣΟΧΗ: αδύναμο sync - πιθανώς δεν υπάρχει σήμα APT.';
    S.s_temp = 'Θερμοκρασία (%s): %.0f ... %.0f °C (1%%-99%%)';
    S.s_notemp = 'Χωρίς °C: δεν επιλέχθηκε δορυφόρος.';
    S.click_hint = 'Αριστερό κλικ στον χάρτη: βάζεις σημείο (βλέπεις τον υπολογισμό βήμα βήμα). Δεξί κλικ κοντά σε σημείο: το βγάζεις.';
    S.hover_idle = 'Πέρασε το ποντίκι πάνω από τον χάρτη για να δεις τη θερμοκρασία';
    S.hover_fmt = 'Γραμμή %d, στήλη %d:  %.1f °C';
    S.no_temp = 'Επίλεξε δορυφόρο (NOAA-15/18/19) και αποκωδικοποίησε για να δεις θερμοκρασίες.';
    S.temp_label = 'Θερμοκρασία (°C)'; S.temp_title = '%s - κανάλι B (θερμικό IR)';
    S.col = 'Στήλη pixel'; S.row = 'Γραμμή';
end
end

function t = howText(lang, par)
% Επεξήγηση του υπολογισμού (με τις πραγματικές τιμές του τρεξίματος όταν υπάρχουν)
en = strcmpi(lang, 'en');
if isempty(par)
    v = struct('Tbbs', NaN, 'Nbb', NaN, 'Cs', NaN, 'Cb', NaN, 'Vc', NaN, 'A', NaN, 'B', NaN, 'Ns', NaN);
else
    v = par;
end
if en
    L = {'HOW TEMPERATURE IS COMPUTED (channel B, thermal IR)', ...
         '1. Count Ce (0-255): pixel brightness, calibrated with the', ...
         '   16 telemetry wedges (removes receiver gain/offset).', ...
         '2. Blackbody: Tbb = mean of 4 PRT sensors (wedges 10-13),', ...
         '   each T = a*x^2 + b*x + c, x = 4*wedge count (10 bit).', ...
         sprintf('   T* = A + B*Tbb = %.2f K ; Nbb = C1*Vc^3/(exp(C2*Vc/T*)-1) = %.3f', v.Tbbs, v.Nbb), ...
         '3. Space view Cs and back scan Cb:', ...
         sprintf('   Cs = %.1f , Cb = %.1f   (10-bit counts)', v.Cs, v.Cb), ...
         '4. Linear radiance: Nl = Ns + (Nbb-Ns)*(Cs-4*Ce)/(Cs-Cb)', ...
         '5. Non-linear correction: Nc = a*Nl^2 + b*Nl + c ; Ne = Nl + Nc', ...
         '6. Te* = C2*Vc / ln(C1*Vc^3/Ne + 1) ; T = (Te*-A)/B - 273.15 [°C]', ...
         sprintf('   Vc = %.2f cm^-1 , A = %.4f , B = %.6f , Ns = %.2f', v.Vc, v.A, v.B, v.Ns), ...
         'Constants: NOAA KLM User''s Guide (via aptdec). Absolute accuracy not verified.'};
else
    L = {'ΠΩΣ ΥΠΟΛΟΓΙΖΕΤΑΙ Η ΘΕΡΜΟΚΡΑΣΙΑ (κανάλι B, θερμικό IR)', ...
         '1. Count Ce (0-255): φωτεινότητα pixel, βαθμονομημένη με τα', ...
         '   16 wedges τηλεμετρίας (αφαιρεί κέρδος/offset δέκτη).', ...
         '2. Μέλαν σώμα: Tbb = μέσος όρος 4 αισθητήρων PRT (wedges 10-13),', ...
         '   καθένας T = a*x^2 + b*x + c, x = 4*count wedge (10 bit).', ...
         sprintf('   T* = A + B*Tbb = %.2f K ; Nbb = C1*Vc^3/(exp(C2*Vc/T*)-1) = %.3f', v.Tbbs, v.Nbb), ...
         '3. Space view Cs και back scan Cb:', ...
         sprintf('   Cs = %.1f , Cb = %.1f   (counts 10 bit)', v.Cs, v.Cb), ...
         '4. Γραμμική ακτινοβολία: Nl = Ns + (Nbb-Ns)*(Cs-4*Ce)/(Cs-Cb)', ...
         '5. Μη γραμμική διόρθωση: Nc = a*Nl^2 + b*Nl + c ; Ne = Nl + Nc', ...
         '6. Te* = C2*Vc / ln(C1*Vc^3/Ne + 1) ; T = (Te*-A)/B - 273.15 [°C]', ...
         sprintf('   Vc = %.2f cm^-1 , A = %.4f , B = %.6f , Ns = %.2f', v.Vc, v.A, v.B, v.Ns), ...
         'Σταθερές: NOAA KLM User''s Guide (μέσω aptdec). Η απόλυτη ακρίβεια δεν έχει επαληθευτεί.'};
end
t = strjoin(L, newline);
end

function t = pointText(lang, par, Ce, row, col)
% Αναλυτικός υπολογισμός για το επιλεγμένο pixel
Nl = par.Ns + (par.Nbb - par.Ns) * (par.Cs - 4*Ce) / (par.Cs - par.Cb);
Nc = par.q(1)*Nl^2 + par.q(2)*Nl + par.q(3);
Ne = max(Nl + Nc, 1e-3);
Tes = par.C2 * par.Vc / log(par.C1 * par.Vc^3 / Ne + 1);
Tk = (Tes - par.A) / par.B;
Tc = Tk - 273.15;
clipped = Tc < -100 || Tc > 60;
Tshow = min(max(Tc, -100), 60);
if strcmpi(lang, 'en')
    L = {sprintf('SELECTED POINT: line %d, column %d', row, col), ...
         sprintf('Ce  = %.1f   (4*Ce = %.1f of 10 bit)', Ce, 4*Ce), ...
         sprintf('Nl  = %.3f   (linear radiance)', Nl), ...
         sprintf('Nc  = %+.3f  (non-linear correction)', Nc), ...
         sprintf('Ne  = %.3f   (corrected radiance)', Ne), ...
         sprintf('Te* = %.2f K (equivalent blackbody temp.)', Tes), ...
         sprintf('T   = (Te*-A)/B = %.2f K', Tk), ...
         sprintf('TEMPERATURE = %.1f °C', Tshow)};
    if clipped, L{end+1} = '(value outside -100..60 °C was clamped: probably noise)'; end
else
    L = {sprintf('ΕΠΙΛΕΓΜΕΝΟ ΣΗΜΕΙΟ: γραμμή %d, στήλη %d', row, col), ...
         sprintf('Ce  = %.1f   (4*Ce = %.1f στα 10 bit)', Ce, 4*Ce), ...
         sprintf('Nl  = %.3f   (γραμμική ακτινοβολία)', Nl), ...
         sprintf('Nc  = %+.3f  (μη γραμμική διόρθωση)', Nc), ...
         sprintf('Ne  = %.3f   (διορθωμένη ακτινοβολία)', Ne), ...
         sprintf('Te* = %.2f K (ισοδύναμη θερμ. μέλανος σώματος)', Tes), ...
         sprintf('T   = (Te*-A)/B = %.2f K', Tk), ...
         sprintf('ΘΕΡΜΟΚΡΑΣΙΑ = %.1f °C', Tshow)};
    if clipped, L{end+1} = '(τιμή εκτός -100..60 °C αποκόπηκε: πιθανώς θόρυβος)'; end
end
t = strjoin(L, newline);
end

function t = pinsText(lang, par, CeB, TC, pins)
% Λίστα όλων των σημείων + αναλυτικός υπολογισμός για το τελευταίο
en = strcmpi(lang, 'en');
n = size(pins, 1);
if en, hdr = 'POINTS (right-click near a point removes it):'; else, hdr = 'ΣΗΜΕΙΑ (δεξί κλικ κοντά σε σημείο το αφαιρεί):'; end
L = {hdr};
for k = 1:n
    if en
        L{end+1} = sprintf('#%d  line %d, col %d :  %.1f °C', k, pins(k,1), pins(k,2), TC(pins(k,1), pins(k,2)));
    else
        L{end+1} = sprintf('#%d  γραμμή %d, στήλη %d :  %.1f °C', k, pins(k,1), pins(k,2), TC(pins(k,1), pins(k,2)));
    end
end
L{end+1} = '';
last = pins(end, :);
lastTxt = pointText(lang, par, CeB(last(1), last(2)), last(1), last(2));
if en, L{end+1} = sprintf('--- calculation for point #%d ---', n); else, L{end+1} = sprintf('--- υπολογισμός για το σημείο #%d ---', n); end
L{end+1} = lastTxt;
t = strjoin(L, newline);
end
