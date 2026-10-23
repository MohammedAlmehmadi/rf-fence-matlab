function varargout = RF_Drone_Fence_App(mode)
%RF_DRONE_FENCE_APP Live 3-D demonstration of the RF Drone Fence forward-scatter simulation.
%
% USAGE
%   RF_Drone_Fence_App                  % open the app
%   RF_Drone_Fence_App('selftest')      % headless checks of the app logic (no window)
%   T = RF_Drone_Fence_App('selftest'); % same, returns the result table
%
% WHAT IS SHOWN (nothing in the plots is drawn "for effect")
%   * Received signal r(t): fsr_detector_comparison's simulateTrial, i.e. the paraxial
%     Babinet/Fresnel field of fsr_scattered_field through fsr_channel_model
%     (complex white noise at the chosen CNR, 3 Hz residual CFO). Multipath is off.
%   * Detector scores: the SAME scoreA / scoreB functions used by the tested protocol,
%     reached through fsr_detector_comparison('api'). Both are causal, so drawing the
%     score up to the cursor shows exactly what a real-time detector would have known.
%   * Threshold: calibrated on target-free recordings only, with the protocol's own seeds
%     and calibrateThreshold, to cfg.target_fa_per_h. Cached on disk after the first run.
%   * "TARGET DETECTED" appears only after an event onset from eventOnsets(score, thr),
%     i.e. when the score has really crossed the calibrated threshold.
%   * Ground-truth evaluation (hit / miss / false alarm) uses t0 ONLY for scoring, with
%     the protocol's scoreCrossing hit window; the decision never sees t0.
%   * Pd and FA/h come from a saved protocol run (fsr_comparison_results.mat), never
%     from the single crossing on screen.
%   The whole recording and both scores are computed once per setting (~0.2 s); the
%   animation only reveals precomputed samples, so motion and data stay synchronised:
%   drone position y(t) = v (t - t0) is the same trajectory the physics used.
%   ALL SIGNALS ARE SYNTHETIC SIMULATION, NOT FIELD MEASUREMENTS. The drone is drawn as
%   an enlarged icon; the scatterer actually modelled is the red 0.30 x 0.20 m flat screen.
%
% REQUIREMENTS: MATLAB R2020b or newer (uifigure, uigridlayout, uilamp, uiprogressdlg,
%   timer). Base MATLAB only, no toolboxes.
% FILES (same folder): fsr_detector_comparison.m (v3, with 'api'), fsr_channel_model.m,
%   fsr_scattered_field.m, fsr_babinet_paraxial.m, fsr_fresnel_cs.m, fsr_wavelength.m.
% OPTIONAL: fsr_comparison_results.mat containing R from
%   R = fsr_detector_comparison(); save fsr_comparison_results R
%   (enables the Pd / FA/h panel and makes the app use the same calibration length).

if nargin >= 1
    switch lower(char(mode))
        case 'selftest'
            T = appSelfTest();
            if nargout > 0, varargout{1} = T; end
            return
        otherwise
            error('rfapp:arg', 'Unknown option ''%s''. Use no argument or ''selftest''.', char(mode));
    end
end
checkDependencies(true);
api = fsr_detector_comparison('api');

%% ------------------------------------------------------------------ state
C = theme();
appDir = fileparts(mfilename('fullpath'));
S = struct();
S.api = api;
S.cfg0 = api.defaultConfig();
S.scen = scenarioList();
S.cache = loadThresholdCache(appDir);
S.R = loadResults(fullfile(appDir, 'fsr_comparison_results.mat'));
S.D = [];
S.tau = 0;
S.playing = false;
S.wall0 = [];
S.tau0 = 0;
S.seedOffset = 0;
S.dlg = [];
S.lastText = struct('decision', '', 'eval', '', 'score', '', 'time', '');
S.az0 = -38; S.el0 = 20;
S.Ysc = 22;                         % half-width of the scene along track [m]
S.iconScale = 2.5;                  % drone icon enlargement (display only)

%% ------------------------------------------------------------------ figure
fig = uifigure('Name', 'RF Drone Fence - forward-scatter simulation (synthetic)', ...
    'Color', C.bg, 'Position', [30 40 1600 920]);
fig.CloseRequestFcn = @(~,~) onClose();

main = uigridlayout(fig, [3 2]);
main.RowHeight = {58, '1x', 210};
main.ColumnWidth = {'1.35x', '1x'};
main.BackgroundColor = C.bg;
main.Padding = [10 10 10 10]; main.RowSpacing = 8; main.ColumnSpacing = 10;

% header
hdr = uigridlayout(main, [1 2]);
hdr.Layout.Row = 1; hdr.Layout.Column = [1 2];
hdr.ColumnWidth = {'1x', 520}; hdr.Padding = [16 6 16 6]; hdr.BackgroundColor = C.panel;
uilabel(hdr, 'Text', 'RF DRONE FENCE  ·  Forward-Scatter Detection Simulator', ...
    'FontSize', 24, 'FontWeight', 'bold', 'FontColor', C.fg);
uilabel(hdr, 'Text', 'SYNTHETIC SIMULATION  ·  NOT FIELD MEASUREMENTS', ...
    'FontSize', 15, 'FontWeight', 'bold', 'FontColor', C.amber, 'HorizontalAlignment', 'right');

% 3-D scene
ax3 = uiaxes(main);
ax3.Layout.Row = 2; ax3.Layout.Column = 1;
styleAxes(ax3, C);
title(ax3, ' ', 'Color', C.fg, 'FontSize', 14);

% signal + score
right = uigridlayout(main, [2 1]);
right.Layout.Row = 2; right.Layout.Column = 2;
right.RowHeight = {'1x', '1x'}; right.Padding = [0 0 0 0]; right.RowSpacing = 8;
right.BackgroundColor = C.bg;
axAmp = uiaxes(right);   styleAxes(axAmp, C);
axScore = uiaxes(right); styleAxes(axScore, C);

% controls
ctrlPanel = uipanel(main, 'Title', 'SIMULATION CONTROLS', 'BackgroundColor', C.panel, ...
    'ForegroundColor', C.dim, 'FontWeight', 'bold', 'FontSize', 12);
ctrlPanel.Layout.Row = 3; ctrlPanel.Layout.Column = 1;
g = uigridlayout(ctrlPanel, [3 7]);
g.RowHeight = {20, 46, 40}; g.ColumnWidth = repmat({'1x'}, 1, 7);
g.Padding = [12 6 12 8]; g.ColumnSpacing = 14; g.RowSpacing = 6; g.BackgroundColor = C.panel;

% ONE creation per control (the old file created every control twice)
W.scen  = addCtl(g, 1, 'Scenario', uidropdown(g, 'Items', {S.scen.name}, ...
    'ItemsData', 1:numel(S.scen), 'Value', 1));
W.det   = addCtl(g, 2, 'Detector', uidropdown(g, ...
    'Items', {'A: amplitude + noise est.', 'B: matched-filter bank'}, 'ItemsData', [1 2], 'Value', 2));
W.freq  = addCtl(g, 3, 'Carrier', uidropdown(g, 'Items', {'2.4 GHz', '5.8 GHz'}, ...
    'ItemsData', [2.4e9 5.8e9], 'Value', 2.4e9));
W.cnr   = addCtl(g, 4, 'CNR (dB)', uispinner(g, 'Limits', [5 45], 'Step', 1, 'Value', 30, ...
    'ValueDisplayFormat', '%.0f dB'));
W.speed = addCtl(g, 5, 'Drone speed (m/s)', uislider(g, 'Limits', [2 12], 'Value', 5, ...
    'MajorTicks', 2:2:12, 'MinorTicks', []));
W.height = addCtl(g, 6, 'Drone height (m)  [LOS 1.5 m]', uislider(g, 'Limits', [0.5 4], 'Value', 1.5, ...
    'MajorTicks', [0.5 1 1.5 2 2.5 3 3.5 4], 'MinorTicks', []));
W.rate  = addCtl(g, 7, 'Playback speed', uidropdown(g, 'Items', {'0.25x', '0.5x', '1x', '2x', '4x'}, ...
    'ItemsData', [0.25 0.5 1 2 4], 'Value', 1));

W.play  = addBtn(g, 1, 'Play',  C.green, @(~,~) onPlay());
W.pause = addBtn(g, 2, 'Pause', C.amber, @(~,~) onPause());
W.reset = addBtn(g, 3, 'Reset', C.dim,   @(~,~) onReset());
W.noise = addBtn(g, 4, 'New noise realisation', C.blue, @(~,~) onNewNoise());
W.cam   = uidropdown(g, 'Items', {'3-D perspective', 'From receiver', 'Side view', 'Top view'}, ...
    'Value', '3-D perspective');
W.cam.Layout.Row = 3; W.cam.Layout.Column = 5;
W.orbit = uicheckbox(g, 'Text', 'Orbit camera', 'Value', false, 'FontColor', C.fg);
W.orbit.Layout.Row = 3; W.orbit.Layout.Column = 6;
W.load  = addBtn(g, 7, 'Load results...', C.dim, @(~,~) onLoadResults());
styleControls(W, C);

% status
stPanel = uipanel(main, 'Title', 'DETECTOR STATUS', 'BackgroundColor', C.panel, ...
    'ForegroundColor', C.dim, 'FontWeight', 'bold', 'FontSize', 12);
stPanel.Layout.Row = 3; stPanel.Layout.Column = 2;
sg = uigridlayout(stPanel, [5 2]);
sg.ColumnWidth = {54, '1x'}; sg.RowHeight = {44, 22, 22, 22, '1x'};
sg.Padding = [12 6 12 6]; sg.RowSpacing = 3; sg.BackgroundColor = C.panel;
W.lamp = uilamp(sg, 'Color', C.green); W.lamp.Layout.Row = 1; W.lamp.Layout.Column = 1;
W.decision = uilabel(sg, 'Text', 'MONITORING', 'FontSize', 26, 'FontWeight', 'bold', 'FontColor', C.green);
W.decision.Layout.Row = 1; W.decision.Layout.Column = 2;
W.scoreTxt = statusLabel(sg, 2, C);
W.timeTxt  = statusLabel(sg, 3, C);
W.evalTxt  = statusLabel(sg, 4, C);
W.protoTxt = uilabel(sg, 'Text', '', 'FontSize', 12, 'FontColor', C.dim, 'WordWrap', 'on', ...
    'VerticalAlignment', 'top');
W.protoTxt.Layout.Row = 5; W.protoTxt.Layout.Column = [1 2];

%% ------------------------------------------------------------------ graphics (created once)
H = buildScene(ax3, S, C);
H = buildSignalPlots(axAmp, axScore, H, C);

% callbacks (set after all controls exist)
W.scen.ValueChangedFcn   = @(~,~) onScenario();
W.det.ValueChangedFcn    = @(~,~) onDetector();
W.freq.ValueChangedFcn   = @(~,~) onParam();
W.cnr.ValueChangedFcn    = @(~,~) onParam();
W.speed.ValueChangedFcn  = @(~,~) onParam();
W.height.ValueChangedFcn = @(~,~) onParam();
W.rate.ValueChangedFcn   = @(~,~) onRate();
W.cam.ValueChangedFcn    = @(~,~) applyCamera(0);

tmr = timer('ExecutionMode', 'fixedRate', 'Period', 0.04, 'BusyMode', 'drop', ...
    'TimerFcn', @(~,~) onTick(), 'Name', 'RFDroneFenceAnimation');

applyScenario(1);
prepare();

%% ================================================================== nested: actions
    function onScenario()
        applyScenario(W.scen.Value);
        S.seedOffset = 0;
        prepare();
    end

    function applyScenario(k)
        sc = S.scen(k);
        W.cnr.Value = sc.cnr;
        W.speed.Value = sc.speed;
        W.height.Value = sc.height;
    end

    function onParam()
        prepare();
    end

    function onDetector()
        S.lastText = struct('decision', '', 'eval', '', 'score', '', 'time', '');
        updateScoreAxes();
        render(S.tau);
    end

    function onRate()
        if S.playing
            S.tau0 = S.tau; S.wall0 = tic;
        end
    end

    function onNewNoise()
        S.seedOffset = S.seedOffset + 1;
        prepare();
    end

    function onPlay()
        if isempty(S.D), return; end
        if S.tau >= S.D.t(end) - 1e-9, S.tau = S.D.t(1); end
        S.tau0 = S.tau; S.wall0 = tic; S.playing = true;
        if strcmp(tmr.Running, 'off'), start(tmr); end
    end

    function onPause()
        S.playing = false;
        if strcmp(tmr.Running, 'on'), stop(tmr); end
    end

    function onReset()
        onPause();
        if isempty(S.D), return; end
        S.tau = S.D.t(1);
        S.lastText = struct('decision', '', 'eval', '', 'score', '', 'time', '');
        render(S.tau);
    end

    function onLoadResults()
        [f, p] = uigetfile('*.mat', 'Select a .mat file containing R from fsr_detector_comparison');
        figure(fig);
        if isequal(f, 0), return; end
        R = loadResults(fullfile(p, f));
        if isempty(R)
            uialert(fig, 'That file has no variable R from fsr_detector_comparison.', 'Results');
            return
        end
        S.R = R;
        prepare();
    end

    function onClose()
        try
            if isvalid(tmr), stop(tmr); delete(tmr); end
        catch
        end
        delete(fig);
    end

%% ================================================================== nested: compute
    function prepare()
        onPause();
        p = currentParams();
        if ~hasThreshold(S.cache, p, S.cfg0, S.R)
            dlg = uiprogressdlg(fig, 'Title', 'Calibrating thresholds (target-free data only)', ...
                'Message', sprintf(['First use of %.1f GHz, CNR %g dB: scoring target-free calibration ' ...
                'recordings with both detectors. Cached afterwards.'], p.fc/1e9, p.cnr), 'Value', 0);
            S.dlg = dlg;
            progress = @(x) setProgress(x);
        else
            dlg = []; progress = [];
        end
        try
            [S.D, S.cache, isNew] = computeScenario(p, S.api, S.cache, S.R, progress);
            if isNew, saveThresholdCache(appDir, S.cache); end
        catch ME
            if ~isempty(dlg), close(dlg); end
            uialert(fig, ME.message, 'Simulation error');
            return
        end
        if ~isempty(dlg), close(dlg); end
        S.dlg = [];
        S.tau = S.D.t(1);
        S.lastText = struct('decision', '', 'eval', '', 'score', '', 'time', '');
        updateStaticGraphics();
        render(S.tau);
    end

    function setProgress(x)
        if ~isempty(S.dlg) && isvalid(S.dlg), S.dlg.Value = min(1, x); end
    end

    function p = currentParams()
        sc = S.scen(W.scen.Value);
        p = struct('fc', W.freq.Value, 'cnr', W.cnr.Value, 'speed', W.speed.Value, ...
            'height', W.height.Value, 'target', sc.target, 'seed', sc.seed + S.seedOffset, ...
            'yHalf', S.Ysc, 'name', sc.name);
    end

%% ================================================================== nested: graphics
    function updateStaticGraphics()
        D = S.D;
        % Fresnel zone and flight path depend on fc / height
        [FX, FY, FZ] = fresnelSurface(D.cfg, D.lambda);
        set(H.fres, 'XData', FX, 'YData', FY, 'ZData', FZ);
        zf = D.cfg.hAnt + D.z0;
        set(H.path, 'XData', [D.cfg.xc D.cfg.xc], 'YData', [-S.Ysc S.Ysc], 'ZData', [zf zf]);
        set(H.trail, 'XData', nan, 'YData', nan, 'ZData', nan);
        vis = onOff(D.target);
        set([H.path H.trail H.body H.rotors H.arms H.screen H.shadow H.drop H.droneLbl], 'Visible', vis);
        set(H.noTarget, 'Visible', onOff(~D.target));
        % amplitude axes
        set(H.cleanLine, 'XData', nan, 'YData', nan);
        set(H.ampLine, 'XData', nan, 'YData', nan);
        yl = [min([D.amp_db; D.clean_db]) max([D.amp_db; D.clean_db])];
        yl = yl + [-0.15 0.15]*max(0.5, diff(yl));
        axAmp.XLim = [D.t(1) D.t(end)]; axAmp.YLim = yl;
        set(H.cursorA, 'YData', yl);
        if D.target
            set(H.t0A, 'XData', [D.t0 D.t0], 'YData', yl, 'Visible', 'on');
        else
            set(H.t0A, 'Visible', 'off');
        end
        title(axAmp, sprintf('Received amplitude |r(t)|  ·  %.1f GHz  ·  CNR %g dB  ·  seed %d', ...
            D.cfg.fc/1e9, D.cfg.cnr_db, D.seed), 'Color', C.fg);
        updateScoreAxes();
        W.protoTxt.Text = protocolText(S.R, D);
    end

    function updateScoreAxes()
        D = S.D; k = W.det.Value;
        s = D.score{k}; thr = D.thr(k);
        yl = [0, max(1.6*thr, 1.08*max(s))];
        axScore.XLim = [D.t(1) D.t(end)]; axScore.YLim = yl;
        set(H.thrLine, 'XData', [D.t(1) D.t(end)], 'YData', [thr thr]);
        set(H.cursorS, 'YData', yl);
        set(H.scoreLine, 'Color', C.det{k});
        if D.target
            set(H.t0S, 'XData', [D.t0 D.t0], 'YData', yl, 'Visible', 'on');
        else
            set(H.t0S, 'Visible', 'off');
        end
        names = {'A (amplitude + noise estimate)', 'B (matched-filter bank)'};
        title(axScore, sprintf('Detector %s score vs calibrated threshold %.2f', names{k}, thr), 'Color', C.fg);
        k1 = frameIndex(D, S.tau);
        set(H.scoreLine, 'XData', D.t(1:k1), 'YData', s(1:k1));
        on = D.onsets{k}; on = on(on <= S.tau);
        set(H.onsetMk, 'XData', on, 'YData', thr*ones(size(on)));
    end

    function onTick()
        try
            if ~isvalid(fig), stop(tmr); return; end
            if ~S.playing, return; end
            tau = S.tau0 + W.rate.Value * toc(S.wall0);
            if tau >= S.D.t(end)
                tau = S.D.t(end);
                S.playing = false;
                stop(tmr);
            end
            S.tau = tau;
            render(tau);
        catch ME
            S.playing = false;
            try, stop(tmr); catch, end
            if isvalid(fig), uialert(fig, ME.message, 'Animation stopped'); end
        end
    end

    function render(tau)
        D = S.D; if isempty(D), return; end
        k = frameIndex(D, tau);
        det = W.det.Value;
        % ---- 3-D scene
        if D.target
            pos = dronePosition(D, tau);
            inScene = abs(pos(2)) <= S.Ysc;
            set([H.body H.rotors H.arms H.screen H.shadow H.drop H.droneLbl], 'Visible', onOff(inScene));
            if inScene
                set(H.body,   'Vertices', H.bodyV0  + pos);
                set(H.rotors, 'Vertices', H.rotorV0 + pos);
                set(H.arms,   'XData', H.armX0 + pos(1), 'YData', H.armY0 + pos(2), 'ZData', H.armZ0 + pos(3));
                set(H.screen, 'Vertices', H.screenV0 + pos);
                set(H.shadow, 'Vertices', H.shadowV0 + [pos(1) pos(2) 0.01]);
                set(H.drop,   'XData', [pos(1) pos(1)], 'YData', [pos(2) pos(2)], 'ZData', [0.02 pos(3)]);
                set(H.droneLbl, 'Position', pos + [0 0 0.9], 'String', ...
                    sprintf('y = %+.1f m   h = %.2f m', pos(2), pos(3)));
            end
            yTrail = linspace(-S.Ysc, max(-S.Ysc, min(S.Ysc, pos(2))), 60);
            set(H.trail, 'XData', D.cfg.xc*ones(1,60), 'YData', yTrail, 'ZData', pos(3)*ones(1,60));
            sceneTitle = sprintf('t = %6.2f s   ·   along-track y = %+6.2f m   ·   %.1f m above LOS', ...
                tau, pos(2), D.z0);
        else
            sceneTitle = sprintf('t = %6.2f s   ·   target-free recording', tau);
        end
        st = decisionState(D, det, k);
        if st.active, set(H.los, 'Color', C.red, 'LineWidth', 4);
        else, set(H.los, 'Color', C.cyan, 'LineWidth', 2.5); end
        ax3.Title.String = sceneTitle;
        if W.orbit.Value, applyCamera(tau - D.t(1)); end
        % ---- signal plots (precomputed samples revealed up to the cursor)
        set(H.ampLine,   'XData', D.t(1:k), 'YData', D.amp_db(1:k));
        set(H.cleanLine, 'XData', D.t(1:k), 'YData', D.clean_db(1:k));
        set(H.scoreLine, 'XData', D.t(1:k), 'YData', D.score{det}(1:k));
        set(H.cursorA, 'XData', [tau tau]);
        set(H.cursorS, 'XData', [tau tau]);
        on = D.onsets{det}; on = on(on <= D.t(k));
        set(H.onsetMk, 'XData', on, 'YData', D.thr(det)*ones(size(on)));
        % ---- status
        if st.detected
            dTxt = 'TARGET DETECTED'; col = C.red;
        elseif st.active
            dTxt = 'ALARM'; col = C.red;
        else
            dTxt = 'MONITORING'; col = C.green;
        end
        if setText(W.decision, 'decision', dTxt)
            W.decision.FontColor = col; W.lamp.Color = col;
        end
        setText(W.scoreTxt, 'score', sprintf('Detector %s   score %6.2f   threshold %.2f   %s', ...
            char('A' + det - 1), D.score{det}(k), D.thr(det), ternary(st.active, '> THRESHOLD', '')));
        if st.detected
            tt = sprintf('First alarm at t = %.2f s (event onset from eventOnsets)', st.tFirst);
        else
            tt = sprintf('No alarm yet  ·  calibrated on %.2f h target-free data to %.1f FA/h', ...
                D.calHours, D.cfg.target_fa_per_h);
        end
        setText(W.timeTxt, 'time', tt);
        setText(W.evalTxt, 'eval', evaluationText(S.api, D, det, D.t(k)));
        drawnow limitrate
    end

    function applyCamera(elapsed)
        switch W.cam.Value
            case '3-D perspective', az = S.az0; el = S.el0;
            case 'From receiver',   az = 90;    el = 10;
            case 'Side view',       az = 0;     el = 6;
            otherwise,              az = 0;     el = 90;
        end
        if W.orbit.Value, az = az + 6*elapsed; end
        view(ax3, az, el);
    end

    function changed = setText(h, key, txt)
        changed = ~strcmp(S.lastText.(key), txt);
        if changed
            h.Text = txt;
            S.lastText.(key) = txt;
        end
    end
end

%% =====================================================================
%% Engine (no UI): used by the app and by the self-test
%% =====================================================================
function sc = scenarioList()
% Seeds 900000+ are disjoint from the protocol's calibration (1e5..), empty (2e5..) and
% target (3e5..) seed ranges, so no demo recording is reused for calibration.
sc = struct( ...
    'name',   {'1  Clear crossing', '2  Difficult crossing (low CNR, offset)', '3  No target'}, ...
    'target', {true, true, false}, ...
    'cnr',    {30, 15, 20}, ...
    'height', {1.5, 2.5, 1.5}, ...
    'speed',  {5, 7, 5}, ...
    'seed',   {900001, 900101, 900201});
end

function [D, cache, isNew] = computeScenario(p, api, cache, R, progress)
% One full recording + both detector scores, computed once. p: fc, cnr, speed, height,
% target, seed, yHalf.
cfg = api.defaultConfig();
cfg.plot = false; cfg.verbose = false;
cfg.fc = p.fc; cfg.cnr_db = p.cnr;
cfg.testSpeedRange = [p.speed p.speed];     % speed fixed by the user (rand is still drawn)
cfg.calTrials = calTrialsFor(R, cfg);
z0 = p.height - cfg.hAnt;
[thr, bank, cache, isNew] = getThresholds(cfg, api, cache, progress);

dur = 120;
tr = api.simulateTrial(cfg, dur, p.target, z0, p.seed);
sA = api.scoreA(tr.r, cfg);
sB = api.scoreB(tr.r, bank, cfg);
onA = api.eventOnsets(sA, thr(1), cfg);
onB = api.eventOnsets(sB, thr(2), cfg);

if p.target, tc = tr.t0; else, tc = 75; end
half = max(p.yHalf / p.speed, 5);           % >= 5 s so the hit window [-1.5, 3] s fits
ts = max(tc - half, cfg.warmup);
te = min(tc + half, tr.t(end));
idx = find(tr.t >= ts & tr.t <= te);

D = struct();
D.cfg = cfg; D.seed = p.seed; D.target = p.target; D.lambda = fsr_wavelength(cfg.fc);
D.t0 = tr.t0; D.v = p.speed; D.z0 = z0; D.fs = cfg.fs;
D.t = tr.t(idx);
D.amp_db   = 20*log10(max(abs(tr.r(idx)), 1e-6));
D.clean_db = 20*log10(max(abs(tr.clean(idx)), 1e-6));
D.score = {sA(idx), sB(idx)};
D.thr = thr;
D.onsets = {onA(onA >= D.t(1) & onA <= D.t(end)), onB(onB >= D.t(1) & onB <= D.t(end))};
D.onsetsFull = {onA, onB};
D.calHours = cfg.calTrials*(cfg.calDur - cfg.warmup)/3600;
D.scenario = p;
end

function n = calTrialsFor(R, cfg)
% Same calibration length as the loaded protocol results, so thresholds match them.
n = cfg.calTrials;
if ~isempty(R) && isfield(R, 'calHours')
    n = round(R.calHours*3600/(cfg.calDur - cfg.warmup));
end
end

function key = thresholdKey(cfg)
key = sprintf('v3|fc=%.6g|cnr=%.6g|cal=%dx%g|fa=%g|drift=%g|interf=%g|cfo=%g|mp=%d', ...
    cfg.fc, cfg.cnr_db, cfg.calTrials, cfg.calDur, cfg.target_fa_per_h, ...
    cfg.drift_db, cfg.interf_db, cfg.cfo_hz, numel(cfg.paths));
end

function tf = hasThreshold(cache, p, cfg0, R)
cfg = cfg0; cfg.fc = p.fc; cfg.cnr_db = p.cnr; cfg.calTrials = calTrialsFor(R, cfg0);
tf = isKey(cache, thresholdKey(cfg));
end

function [thr, bank, cache, isNew] = getThresholds(cfg, api, cache, progress)
% Calibration exactly as in runComparison: target-free recordings, protocol seeds.
bank = api.buildBank(cfg);
key = thresholdKey(cfg);
isNew = ~isKey(cache, key);
if ~isNew
    thr = cache(key);
    return
end
iC = find(abs(cfg.cnrList - cfg.cnr_db) < 1e-9, 1);
if isempty(iC), iC = 9; end          % CNR not in the protocol grid: own seed block
sA = cell(cfg.calTrials, 1); sB = sA;
for k = 1:cfg.calTrials
    tr = api.simulateTrial(cfg, cfg.calDur, false, 0, api.seedOf(cfg, 'cal', iC, 0, k));
    sA{k} = api.scoreA(tr.r, cfg);
    sB{k} = api.scoreB(tr.r, bank, cfg);
    if ~isempty(progress), progress(k/cfg.calTrials); end
end
thr = [api.calibrateThreshold(sA, cfg), api.calibrateThreshold(sB, cfg)];
cache(key) = thr;
end

function k = frameIndex(D, tau)
k = floor((tau - D.t(1))*D.fs + 1e-9) + 1;
k = max(1, min(numel(D.t), k));
end

function pos = dronePosition(D, tau)
% Same trajectory as the physics: y(t) = v (t - t0), x = xc, z = hAnt + z0.
pos = [D.cfg.xc, D.v*(tau - D.t0), D.cfg.hAnt + D.z0];
end

function st = decisionState(D, det, k)
% Decision uses ONLY the causal score and the calibrated threshold (via eventOnsets).
s = D.score{det};
st.active = s(k) > D.thr(det);
on = D.onsets{det};
past = on(on <= D.t(k) + 1e-9);
st.detected = ~isempty(past);
if st.detected, st.tFirst = past(1); else, st.tFirst = NaN; end
end

function txt = evaluationText(api, D, det, tNow)
% Scoring against ground truth (display of the evaluation only; never fed back).
on = D.onsets{det}; on = on(on <= tNow + 1e-9);
if D.target
    hw = D.cfg.hitWin;
    [hit, lat] = api.scoreCrossing(on, D.t0, D.cfg);
    if hit
        txt = sprintf('Ground truth: HIT  ·  latency %+.2f s from LOS crossing', lat);
    elseif tNow > D.t0 + hw(2)
        txt = 'Ground truth: MISS  ·  no alarm inside the hit window';
    elseif ~isempty(on)
        txt = 'Ground truth: alarm before the hit window (would count as false alarm)';
    else
        txt = sprintf('Ground truth: drone reaches the LOS at t = %.2f s (scoring only)', D.t0);
    end
else
    if isempty(on)
        txt = 'Ground truth: no target  ·  no alarm so far (correct)';
    else
        txt = sprintf('Ground truth: no target  ·  FALSE ALARM (%d event%s)', numel(on), ternary(numel(on) > 1, 's', ''));
    end
end
end

function txt = protocolText(R, D)
if isempty(R)
    txt = sprintf(['Pd / FA-per-hour: no protocol results loaded. Run once:  R = fsr_detector_comparison(); ' ...
        'save fsr_comparison_results R   then restart the app or press "Load results...". ' ...
        'Threshold here: %.2f h target-free calibration.'], D.calHours);
    return
end
fcR = 2.4e9; if isfield(R, 'cfg'), fcR = R.cfg.fc; end
iC = find(abs(R.cnr - D.cfg.cnr_db) < 1e-9, 1);
iZ = find(abs(R.z0 - D.z0) < 1e-6, 1);
mpOK = ~isfield(R, 'multipath') || ~isempty(strfind(R.multipath, 'off'));
if abs(fcR - D.cfg.fc) > 1 || isempty(iC) || isempty(iZ) || ~mpOK
    txt = sprintf(['Protocol Pd / FA-h: no saved cell for %.1f GHz, CNR %g dB, z0 %.1f m. ' ...
        'Saved grid: %.1f GHz, CNR %s dB, z0 %s m.'], D.cfg.fc/1e9, D.cfg.cnr_db, D.z0, ...
        fcR/1e9, mat2str(R.cnr), mat2str(R.z0));
    return
end
spd = '';
if isfield(R, 'cfg'), spd = sprintf(', speeds U(%g, %g) m/s', R.cfg.testSpeedRange); end
txt = sprintf(['PROTOCOL (synthetic, %d crossings%s; %.1f h independent target-free):  ' ...
    'Pd A %.2f [%.2f-%.2f]  ·  Pd B %.2f [%.2f-%.2f]  ·  FA/h A %.2f [%.2f-%.2f]  ·  ' ...
    'FA/h B %.2f [%.2f-%.2f]  ·  median latency A %.2f s, B %.2f s'], ...
    R.crossings, spd, R.emptyHours, ...
    R.pdA(iC,iZ), R.pdA_ci(iC,iZ,1), R.pdA_ci(iC,iZ,2), ...
    R.pdB(iC,iZ), R.pdB_ci(iC,iZ,1), R.pdB_ci(iC,iZ,2), ...
    R.faA(iC), R.faA_ci(iC,1), R.faA_ci(iC,2), R.faB(iC), R.faB_ci(iC,1), R.faB_ci(iC,2), ...
    R.latA(iC,iZ), R.latB(iC,iZ));
if isfield(R, 'thrA') && abs(R.thrA(iC) - D.thr(1)) < 1e-12 && abs(R.thrB(iC) - D.thr(2)) < 1e-12
    txt = [txt '  ·  thresholds identical to the protocol run'];
end
end

%% =====================================================================
%% Persistence
%% =====================================================================
function cache = loadThresholdCache(appDir)
cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
f = fullfile(appDir, 'fsr_app_threshold_cache.mat');
if exist(f, 'file') == 2
    try
        L = load(f, 'keys', 'vals');
        for i = 1:numel(L.keys), cache(L.keys{i}) = L.vals{i}; end
    catch
    end
end
end

function saveThresholdCache(appDir, cache)
keys = cache.keys(); vals = cache.values(); %#ok<NASGU>
try
    save(fullfile(appDir, 'fsr_app_threshold_cache.mat'), 'keys', 'vals');
catch
end
end

function R = loadResults(f)
R = [];
if exist(f, 'file') ~= 2, return; end
try
    L = load(f);
    if isfield(L, 'R') && isstruct(L.R) && isfield(L.R, 'pdA'), R = L.R; end
catch
end
end

function ok = checkDependencies(throwErr)
need = {'fsr_detector_comparison', 'fsr_channel_model', 'fsr_scattered_field', ...
        'fsr_babinet_paraxial', 'fsr_fresnel_cs', 'fsr_wavelength'};
miss = need(cellfun(@(n) exist(n, 'file') ~= 2, need));
rel = version('-release');
oldMatlab = str2double(rel(1:4)) < 2020 || strcmp(rel, '2020a');
ok = isempty(miss) && ~oldMatlab;
if ~ok && throwErr
    if ~isempty(miss)
        error('rfapp:deps', 'Missing on the MATLAB path: %s.m', strjoin(miss, '.m, '));
    end
    error('rfapp:version', 'MATLAB R2020b or newer is required (this is R%s).', rel);
end
try
    api = fsr_detector_comparison('api'); %#ok<NASGU>
catch
    if throwErr
        error('rfapp:deps', 'fsr_detector_comparison.m has no ''api'' option: use the v3 file.');
    end
    ok = false;
end
end

%% =====================================================================
%% Scene construction (called once)
%% =====================================================================
function C = theme()
C.bg    = [0.055 0.063 0.082];
C.panel = [0.098 0.110 0.137];
C.axes  = [0.075 0.086 0.110];
C.fg    = [0.88 0.91 0.95];
C.dim   = [0.58 0.63 0.71];
C.grid  = [0.35 0.42 0.52];
C.cyan  = [0.20 0.82 1.00];
C.blue  = [0.35 0.60 1.00];
C.amber = [1.00 0.74 0.22];
C.red   = [1.00 0.30 0.30];
C.green = [0.30 0.86 0.48];
C.det   = {[0.35 0.65 1.00], [1.00 0.60 0.20]};
end

function styleAxes(ax, C)
ax.Color = C.axes; ax.XColor = C.dim; ax.YColor = C.dim; ax.ZColor = C.dim;
ax.GridColor = C.grid; ax.GridAlpha = 0.25; ax.FontSize = 12;
ax.Title.Color = C.fg; ax.XLabel.Color = C.dim; ax.YLabel.Color = C.dim;
grid(ax, 'on'); box(ax, 'on');
end

function w = addCtl(g, col, labelText, w)
C = theme();
lb = uilabel(g, 'Text', labelText, 'FontSize', 12, 'FontColor', C.dim);
lb.Layout.Row = 1; lb.Layout.Column = col;
w.Layout.Row = 2; w.Layout.Column = col;
end

function b = addBtn(g, col, txt, col3, fcn)
b = uibutton(g, 'Text', txt, 'ButtonPushedFcn', fcn, 'FontSize', 14, 'FontWeight', 'bold', ...
    'BackgroundColor', col3*0.55, 'FontColor', [1 1 1]);
b.Layout.Row = 3; b.Layout.Column = col;
end

function h = statusLabel(parent, row, C)
h = uilabel(parent, 'Text', '', 'FontSize', 14, 'FontColor', C.fg);
h.Layout.Row = row; h.Layout.Column = [1 2];
end

function styleControls(W, C)
for h = [W.scen W.det W.freq W.rate W.cam]
    h.BackgroundColor = C.axes; h.FontColor = C.fg; h.FontSize = 13;
end
W.cnr.BackgroundColor = C.axes; W.cnr.FontColor = C.fg; W.cnr.FontSize = 13;
W.speed.FontColor = C.dim; W.height.FontColor = C.dim;
end

function H = buildScene(ax, S, C)
cfg = S.cfg0; L = cfg.L; hA = cfg.hAnt; Y = S.Ysc;
hold(ax, 'on');
% ground
[gx, gy] = meshgrid(-4:2:L+4, -Y:2:Y);
surf(ax, gx, gy, zeros(size(gx)), 'FaceColor', [0.085 0.10 0.12], 'EdgeColor', [0.16 0.22 0.28], ...
    'FaceAlpha', 1, 'HitTest', 'off');
% masts and antenna panels
[cx, cy, cz] = cylinder(0.05, 14);
for xm = [0 L]
    surf(ax, cx + xm, cy, cz*hA, 'FaceColor', [0.62 0.65 0.70], 'EdgeColor', 'none');
end
drawBox(ax, [0 0 hA], [0.10 0.45 0.55], C.blue);
drawBox(ax, [L 0 hA], [0.10 0.45 0.55], C.red*0.9);
text(ax, 0, 0, hA + 0.75, 'Tx', 'Color', C.fg, 'FontSize', 15, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
text(ax, L, 0, hA + 0.75, 'Rx', 'Color', C.fg, 'FontSize', 15, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
text(ax, L/2, 0, hA - 0.55, 'line of sight', 'Color', C.cyan, 'FontSize', 11, 'HorizontalAlignment', 'center');
% line of sight and first Fresnel zone
H.los = plot3(ax, [0 L], [0 0], [hA hA], '-', 'Color', C.cyan, 'LineWidth', 2.5);
H.fres = surf(ax, nan(2), nan(2), nan(2), 'FaceColor', C.cyan, 'FaceAlpha', 0.08, ...
    'EdgeColor', 'none', 'HitTest', 'off');
% flight path and trail
H.path  = plot3(ax, nan, nan, nan, '--', 'Color', [0.55 0.58 0.68], 'LineWidth', 1);
H.trail = plot3(ax, nan, nan, nan, '-', 'Color', C.amber, 'LineWidth', 2.2);
% drone icon (enlarged, pitched 10 deg forward along +y), true-size modelled screen
s = S.iconScale;
[H.bodyV0, bodyF] = boxGeometry([0 0 0], [0.24 0.24 0.07]*s);
armEnds = [ 0.19  0.19; -0.19  0.19; -0.19 -0.19;  0.19 -0.19]*s;
ax0 = []; ay0 = []; az0 = [];
for i = 1:4
    ax0 = [ax0 0 armEnds(i,1) NaN]; ay0 = [ay0 0 armEnds(i,2) NaN]; az0 = [az0 0 0.03*s NaN]; %#ok<AGROW>
end
th = linspace(0, 2*pi, 17); th(end) = [];
rotV = zeros(64, 3); rotF = zeros(4, 16);
for i = 1:4
    rows = (i-1)*16 + (1:16);
    rotV(rows, :) = [armEnds(i,1) + 0.10*s*cos(th'), armEnds(i,2) + 0.10*s*sin(th'), 0.05*s*ones(16,1)];
    rotF(i, :) = rows;
end
Rp = pitchMatrix(-10);
H.bodyV0 = H.bodyV0 * Rp';
rotV = rotV * Rp';
armP = [ax0' ay0' az0'] * Rp';
H.armX0 = armP(:,1)'; H.armY0 = armP(:,2)'; H.armZ0 = armP(:,3)';
H.rotorV0 = rotV;
H.body   = patch(ax, 'Vertices', H.bodyV0, 'Faces', bodyF, 'FaceColor', [0.86 0.88 0.92], ...
    'EdgeColor', [0.25 0.27 0.32]);
H.rotors = patch(ax, 'Vertices', H.rotorV0, 'Faces', rotF, 'FaceColor', C.cyan, 'FaceAlpha', 0.35, ...
    'EdgeColor', C.cyan);
H.arms   = plot3(ax, H.armX0, H.armY0, H.armZ0, '-', 'Color', [0.80 0.82 0.86], 'LineWidth', 3);
w = cfg.w; h = cfg.h;
H.screenV0 = [0 -w/2 -h/2; 0 w/2 -h/2; 0 w/2 h/2; 0 -w/2 h/2];
H.screen = patch(ax, 'Vertices', H.screenV0, 'Faces', 1:4, 'FaceColor', C.red, 'FaceAlpha', 0.55, ...
    'EdgeColor', C.red, 'LineWidth', 1.5);
H.shadowV0 = [0.30*s*cos(th') 0.30*s*sin(th') zeros(16,1)];
H.shadow = patch(ax, 'Vertices', H.shadowV0, 'Faces', 1:16, 'FaceColor', [0 0 0], 'FaceAlpha', 0.45, ...
    'EdgeColor', 'none');
H.drop = plot3(ax, nan, nan, nan, ':', 'Color', [0.7 0.7 0.75], 'LineWidth', 1);
H.droneLbl = text(ax, 0, 0, 0, '', 'Color', C.amber, 'FontSize', 12, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center');
H.noTarget = text(ax, L/2, 0, hA + 2.2, 'NO TARGET IN SCENE  (target-free recording)', ...
    'Color', C.dim, 'FontSize', 15, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Visible', 'off');
text(ax, -3.5, -Y + 1, 0.05, sprintf('Icon enlarged x%.1f; red rectangle = modelled 0.30 x 0.20 m scatterer', s), ...
    'Color', C.dim, 'FontSize', 10);
hold(ax, 'off');
xlim(ax, [-4 L+4]); ylim(ax, [-Y Y]); zlim(ax, [0 6]);
daspect(ax, [1 1 1]);
ax.Projection = 'perspective';
xlabel(ax, 'x along baseline (m)'); ylabel(ax, 'y along track (m)'); zlabel(ax, 'height (m)');
ax.ZLabel.Color = C.dim;
view(ax, S.az0, S.el0);
end

function H = buildSignalPlots(axAmp, axScore, H, C)
hold(axAmp, 'on');
H.cleanLine = plot(axAmp, nan, nan, '--', 'Color', [0.75 0.78 0.85], 'LineWidth', 1.2);
H.ampLine   = plot(axAmp, nan, nan, '-', 'Color', C.cyan, 'LineWidth', 1.3);
H.t0A       = plot(axAmp, nan, nan, ':', 'Color', C.dim, 'LineWidth', 1);
H.cursorA   = plot(axAmp, [0 0], [0 1], '-', 'Color', C.amber, 'LineWidth', 2);
hold(axAmp, 'off');
xlabel(axAmp, 'recording time (s)'); ylabel(axAmp, '20 log_{10}|r| (dB re free space)');
lg = legend(axAmp, [H.ampLine H.cleanLine H.t0A], {'received (with noise)', 'noise-free physics', 'LOS crossing t_0 (truth)'}, ...
    'Location', 'southwest');
lg.TextColor = C.fg; lg.Color = C.panel; lg.EdgeColor = C.grid;

hold(axScore, 'on');
H.thrLine   = plot(axScore, nan, nan, '--', 'Color', C.red, 'LineWidth', 1.6);
H.scoreLine = plot(axScore, nan, nan, '-', 'Color', C.det{2}, 'LineWidth', 1.4);
H.t0S       = plot(axScore, nan, nan, ':', 'Color', C.dim, 'LineWidth', 1);
H.onsetMk   = plot(axScore, nan, nan, 'v', 'Color', C.red, 'MarkerFaceColor', C.red, 'MarkerSize', 9);
H.cursorS   = plot(axScore, [0 0], [0 1], '-', 'Color', C.amber, 'LineWidth', 2);
hold(axScore, 'off');
xlabel(axScore, 'recording time (s)'); ylabel(axScore, 'detector score');
lg = legend(axScore, [H.scoreLine H.thrLine H.onsetMk], {'score (causal)', 'calibrated threshold', 'alarm onset'}, ...
    'Location', 'northwest');
lg.TextColor = C.fg; lg.Color = C.panel; lg.EdgeColor = C.grid;
end

function [FX, FY, FZ] = fresnelSurface(cfg, lambda)
% First Fresnel zone around the LOS: r1(x) = sqrt(lambda x (L - x) / L).
xs = linspace(0, cfg.L, 60)'; th = linspace(0, 2*pi, 36);
r1 = sqrt(lambda .* xs .* (cfg.L - xs) ./ cfg.L);
FX = repmat(xs, 1, numel(th));
FY = r1 * cos(th);
FZ = cfg.hAnt + r1 * sin(th);
end

function drawBox(ax, c, sz, col)
[V, F] = boxGeometry(c, sz);
patch(ax, 'Vertices', V, 'Faces', F, 'FaceColor', col, 'EdgeColor', col*0.6, 'FaceAlpha', 0.95);
end

function [V, F] = boxGeometry(c, sz)
d = sz/2;
V = [-1 -1 -1; 1 -1 -1; 1 1 -1; -1 1 -1; -1 -1 1; 1 -1 1; 1 1 1; -1 1 1] .* d + c;
F = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
end

function Rm = pitchMatrix(deg)
a = deg*pi/180;
Rm = [1 0 0; 0 cos(a) -sin(a); 0 sin(a) cos(a)];
end

function s = onOff(tf)
if tf, s = 'on'; else, s = 'off'; end
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end

%% =====================================================================
%% Self-test (headless): checks the app logic, not the physics or detectors
%% =====================================================================
function T = appSelfTest()
fprintf('\nRF_Drone_Fence_App - self-test (no window; synthetic data)\n');
T = struct('name', {}, 'pass', {}, 'detail', {});
ok = checkDependencies(false);
T = rec(T, 'D1 dependencies, v3 api and MATLAB release (R2020b+)', ok, version('-release'));
if ~ok
    fprintf('Stopping: fix dependencies first.\n'); return
end
api = fsr_detector_comparison('api');
cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
sc = scenarioList();
calSmall = 4;                                       % short calibration for test speed
Rfake = struct('calHours', calSmall*(600 - 35)/3600);  % forces calTrials = 4

% D2 thresholds equal the protocol's own calibration (same seeds, same functions)
p = struct('fc', 2.4e9, 'cnr', 30, 'speed', 5, 'height', 1.5, 'target', true, ...
    'seed', sc(1).seed, 'yHalf', 22);
[D1, cache] = computeScenario(p, api, cache, Rfake, []);
cP = struct('cnrList', 30, 'z0List', 0, 'calTrials', calSmall, 'emptyTrials', 1, ...
    'targetTrials', 1, 'plot', false, 'verbose', false);
Rp = fsr_detector_comparison(cP);
T = rec(T, 'D2 app thresholds = fsr_detector_comparison thresholds (CNR 30 dB)', ...
    isequal(D1.thr, [Rp.thrA(1) Rp.thrB(1)]), ...
    sprintf('app [%.4f %.4f], protocol [%.4f %.4f]', D1.thr, Rp.thrA(1), Rp.thrB(1)));

% D3 no decision without a real threshold crossing, for every frame and both detectors
okD3 = true; nFrames = 0;
for det = 1:2
    s = D1.score{det}; thr = D1.thr(det);
    for k = 1:numel(D1.t)
        st = decisionState(D1, det, k);
        if st.detected && ~any(s(1:k) > thr), okD3 = false; end
        if st.active ~= (s(k) > thr), okD3 = false; end
        nFrames = nFrames + 1;
    end
    on = D1.onsets{det};
    for i = 1:numel(on)
        j = round((on(i) - D1.t(1))*D1.fs) + 1;
        okD3 = okD3 && s(j) > thr && (j == 1 || s(j-1) <= thr);
    end
end
T = rec(T, 'D3 "detected" only after score > threshold; onsets are real crossings', okD3, ...
    sprintf('%d frames checked', nFrames));

% D4 motion and data are synchronised with the physics trajectory
okD4 = true;
for k = 1:25:numel(D1.t)
    tau = D1.t(k) + 0.3/D1.fs;
    kk = frameIndex(D1, tau);
    pos = dronePosition(D1, tau);
    okD4 = okD4 && kk == k && abs(pos(2) - D1.v*(tau - D1.t0)) < 1e-12 && ...
        abs(pos(3) - (D1.cfg.hAnt + D1.z0)) < 1e-12;
end
% the drone trajectory shown, fed independently to the physics, reproduces the plotted signal
yShown = D1.v*(D1.t - D1.t0);
Eind = 1 - fsr_scattered_field('paraxial', api.makeGeom(D1.cfg), yShown, D1.z0, D1.lambda);
eSync = max(abs(20*log10(max(abs(Eind), 1e-6)) - D1.clean_db));
T = rec(T, 'D4 drone y = v(t - t0) and frame index are synchronised with the signal', ...
    okD4 && eSync < 1e-9, sprintf('physics(trajectory shown) vs plotted noise-free signal: %.1e dB', eSync));

% D5 recomputation is deterministic; new noise realisation changes the noise only
D1b = computeScenario(p, api, cache, Rfake, []);
p2 = p; p2.seed = p.seed + 1;
D1c = computeScenario(p2, api, cache, Rfake, []);
T = rec(T, 'D5 same settings -> identical data; new seed -> different noise', ...
    isequal(D1.amp_db, D1b.amp_db) && isequal(D1.score, D1b.score) && ~isequal(D1.amp_db, D1c.amp_db), '');

% D6 all three scenarios compute; the no-target one carries no truth crossing
okD6 = true; det6 = '';
for i = 1:3
    pi_ = struct('fc', 2.4e9, 'cnr', sc(i).cnr, 'speed', sc(i).speed, 'height', sc(i).height, ...
        'target', sc(i).target, 'seed', sc(i).seed, 'yHalf', 22);
    [Di, cache] = computeScenario(pi_, api, cache, Rfake, []);
    okD6 = okD6 && all(isfinite(Di.amp_db)) && all(isfinite(Di.score{1})) && all(isfinite(Di.score{2})) ...
        && (Di.target == ~isnan(Di.t0)) && Di.t(1) >= Di.cfg.warmup;
    stA = decisionState(Di, 1, numel(Di.t)); stB = decisionState(Di, 2, numel(Di.t));
    det6 = [det6 sprintf('%s: A %s, B %s; ', strtrim(sc(i).name(3:end)), ...
        ternary(stA.detected, 'alarm', 'no alarm'), ternary(stB.detected, 'alarm', 'no alarm'))]; %#ok<AGROW>
end
T = rec(T, 'D6 three scenarios compute (clear / difficult / no target)', okD6, det6);

fprintf('\n================ APP SELF-TEST SUMMARY ================\n');
for i = 1:numel(T)
    if T(i).pass, tag = 'PASS'; else, tag = 'FAIL'; end
    fprintf('[%s] %s\n       %s\n', tag, T(i).name, T(i).detail);
end
nF = sum(~[T.pass]);
fprintf('\n%d tests, %d failed. Alarm outcomes above are single synthetic recordings, not Pd.\n', numel(T), nF);
if nF > 0, error('rfapp:selftest', '%d app self-tests failed.', nF); end
end

function T = rec(T, name, pass, detail)
T(end+1) = struct('name', name, 'pass', logical(pass), 'detail', detail);
if pass, tag = 'PASS'; else, tag = 'FAIL'; end
fprintf('  [%s] %s\n', tag, name);
end
