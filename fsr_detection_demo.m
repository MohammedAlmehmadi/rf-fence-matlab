function varargout = fsr_detection_demo(arg)
%FSR_DETECTION_DEMO Synthetic detection layer on top of the verified FSR physics.
%
% USAGE
%   fsr_detection_demo                    % run the default experiment + figure
%   fsr_detection_demo('tests')           % run the basic tests (no target / weak / strong noise)
%   cfg.cnr_db = 20; fsr_detection_demo(cfg)   % override any field of the default config
%   r = fsr_detection_demo(cfg);          % also return all numbers in a struct
%
% WHAT IT DOES
%   1. Signal: E/E0 = 1 - s(t) from fsr_scattered_field('paraxial', ...) with
%      the target centre moving along y(t) = v (t - t0) through the line of
%      sight (LOS) at x = xc. Physics files are only CALLED, never changed.
%   2. Channel (all SIMULATED, documented in defaultConfig):
%        r(t) = g(t) * [E(t) + c(t)] * exp(j 2 pi f_cfo t) + n(t)
%      n(t) complex white Gaussian noise, E|n|^2 = 10^(-CNR/10) (|E0| = 1);
%      g(t) optional slow gain drift (dB, Ornstein-Uhlenbeck);
%      c(t) optional slow complex interference (simple multipath-like wander);
%      f_cfo  residual carrier offset between independent Tx/Rx clocks.
%   3. Detector (amplitude only, so independent of carrier offset/phase):
%        a = |r|, b = slow running mean of a (tau_base),  e = a/b - 1,
%        e_hp = e - causal mean of e over hp_win      (removes slow drift),
%        e_lp = causal mean of e_hp over lp_win       (smooths noise),
%        sigma = causal RMS of e_lp over cfar_win, delayed by guard (CFAR),
%        score = |e_lp| / sigma.
%   4. Threshold: calibrated ONLY on target-free calibration recordings (own
%      seeds) to cfg.target_fa_per_h; evaluated on separate empty and target
%      recordings (different seeds).
%   5. Alarm events: threshold exceedances closer than mergeGap seconds are one
%      event; the event time is its first sample above threshold (onset).
%   6. Metrics: Pd per crossing (Wilson 95% CI), false alarms per hour on the
%      independent empty recordings (exact Poisson 95% CI), detection latency.
%
% DEPENDENCIES (must be in the same folder, unchanged):
%   fsr_scattered_field.m, fsr_babinet_paraxial.m, fsr_fresnel_cs.m, fsr_wavelength.m
% Base MATLAB only: no toolbox functions (no awgn, medfilt1, butter, prctile,
% chi2inv, norminv). Uses randn, filter, cumsum, gammaincinv, erfinv.
%
% SCOPE: every number produced here is SYNTHETIC. It shows how this detector
% behaves on this simulated channel. It is NOT a field result, NOT a measured
% Pd or false-alarm rate, and says nothing about birds, people or real drones.

if nargin < 1, arg = struct(); end
if ischar(arg) || isstring(arg)
    switch lower(char(arg))
        case 'tests'
            results = runTests();
            if nargout > 0, varargout{1} = results; end
            return
        case 'demo'
            arg = struct();
        otherwise
            error('fsr:arg','Unknown option ''%s''. Use ''tests'' or a config struct.', char(arg));
    end
end
cfg = mergeConfig(defaultConfig(), arg);
results = runExperiment(cfg);
if nargout > 0, varargout{1} = results; end
end

%% =====================================================================
%% Configuration
%% =====================================================================
function cfg = defaultConfig()
% Geometry (same as the verified physics tests)
cfg.fc        = 2.4e9;      % carrier [Hz] (2.4e9 or 5.8e9)
cfg.L         = 20.0;       % baseline [m]
cfg.hAnt      = 1.5;        % Tx/Rx antenna height [m]
cfg.xc        = 10.0;       % crossing point along the baseline [m]
cfg.w         = 0.30;       % screen width (along track) [m]
cfg.h         = 0.20;       % screen height [m]
cfg.opacity   = 1.0;        % 1 = fully opaque screen
cfg.z0        = 0.0;        % vertical offset of the target centre from the LOS [m]
cfg.speedRange = [5 5];     % crossing speed [m/s]; [vmin vmax] draws uniformly per crossing
% Sampling
cfg.fs        = 250;        % baseband rate seen by the detector [Hz] (signature <= v/w ~ 17 Hz)
% Simulated channel
cfg.cnr_db    = 30;         % direct-carrier-to-noise ratio per complex sample [dB]
cfg.drift_db  = 0;          % slow gain drift, std in dB (0 = off), e.g. 0.3
cfg.drift_tau = 60;         % drift correlation time [s]
cfg.interf_db = -Inf;       % slow complex interference power vs direct path [dB] (-Inf = off), e.g. -30
cfg.interf_tau = 5;         % interference correlation time [s]
cfg.cfo_hz    = 3;          % residual carrier offset [Hz] (amplitude detector ignores it)
% Detector
cfg.tau_base  = 3.0;        % running-mean time constant for normalisation [s]
cfg.hp_win    = 0.25;       % high-pass: subtract causal mean over this window [s]
cfg.lp_win    = 0.04;       % low-pass: causal mean over this window [s]
cfg.cfar_win  = 30;         % noise-estimate window [s]
cfg.guard     = 2;          % gap between noise window and current sample [s]
cfg.warmup    = 35;         % first seconds ignored (detector not settled) [s]
% Calibration / evaluation protocol
cfg.target_fa_per_h = 0.5;  % calibration target (half of 1/h, as margin)
cfg.calTrials   = 24;  cfg.calDur    = 600;   % 24 x 10 min = 4 h target-free calibration
cfg.emptyTrials = 24;  cfg.emptyDur  = 600;   % 4 h independent target-free test
cfg.targetTrials = 150; cfg.targetDur = 120;  % 150 crossings, one per 2-min recording
cfg.seedCal   = 10000;      % seeds: calibration 10000+, empty test 20000+, target test 30000+
cfg.seedEmpty = 20000;
cfg.seedTarget = 30000;
cfg.mergeGap  = 1.0;        % exceedances closer than this are one event [s]
cfg.hitWin    = [-1.5 3.0]; % event onset in [t0-1.5, t0+3] s counts as a detection
cfg.assocWin  = 6.0;        % events within +-6 s of t0 are attributed to the target
cfg.plot      = true;
cfg.verbose   = true;
end

function cfg = mergeConfig(cfg, user)
if ~isstruct(user), error('fsr:arg','Config must be a struct.'); end
f = fieldnames(user);
for i = 1:numel(f)
    if ~isfield(cfg, f{i})
        error('fsr:arg','Unknown config field ''%s''.', f{i});
    end
    cfg.(f{i}) = user.(f{i});
end
if numel(cfg.speedRange) == 1, cfg.speedRange = [cfg.speedRange cfg.speedRange]; end
if any(cfg.speedRange <= 0), error('fsr:arg','speedRange must be positive.'); end
if cfg.targetDur < cfg.warmup + 30
    error('fsr:arg','targetDur must be at least warmup + 30 s.');
end
end

%% =====================================================================
%% Experiment: calibrate on empty data, test on independent data
%% =====================================================================
function R = runExperiment(cfg)
tic;
lambda = fsr_wavelength(cfg.fc);
deff = cfg.xc*(cfg.L-cfg.xc)/cfg.L;
if cfg.verbose
    fprintf('\nRF Drone Fence - synthetic detection experiment (NOT field data)\n');
    fprintf('  f = %.1f GHz, L = %.0f m, xc = %.0f m, screen %.2f x %.2f m, z0 = %.2f m, v = %g-%g m/s\n', ...
        cfg.fc/1e9, cfg.L, cfg.xc, cfg.w, cfg.h, cfg.z0, cfg.speedRange(1), cfg.speedRange(2));
    fprintf('  CNR = %g dB, drift = %g dB, interference = %g dB, fs = %g Hz\n', ...
        cfg.cnr_db, cfg.drift_db, cfg.interf_db, cfg.fs);
end

% ---- 1. calibration on target-free recordings only
calScores = cell(cfg.calTrials,1);
for k = 1:cfg.calTrials
    tr = simulateTrial(cfg, cfg.calDur, false, cfg.seedCal + k);
    calScores{k} = detectorScore(tr.r, cfg);
end
[thr, calFA, calHours] = calibrateThreshold(calScores, cfg);
if cfg.verbose
    fprintf('  Calibration: %d target-free recordings (%.1f h) -> threshold %.3f (%d events, %.2f /h)\n', ...
        cfg.calTrials, calHours, thr, calFA, calFA/calHours);
end

% ---- 2. independent target-free test
nFA = 0; hoursEmpty = 0; exampleEmpty = [];
for k = 1:cfg.emptyTrials
    tr = simulateTrial(cfg, cfg.emptyDur, false, cfg.seedEmpty + k);
    sc = detectorScore(tr.r, cfg);
    on = eventOnsets(sc, thr, cfg);
    nFA = nFA + numel(on);
    hoursEmpty = hoursEmpty + (cfg.emptyDur - cfg.warmup)/3600;
    if k == 1, exampleEmpty = struct('t',tr.t,'score',sc,'onsets',on); end
end

% ---- 3. independent target test (one crossing per recording)
hit = false(cfg.targetTrials,1); lat = nan(cfg.targetTrials,1);
vAll = zeros(cfg.targetTrials,1); nFAt = 0; hoursT = 0; exampleT = [];
for k = 1:cfg.targetTrials
    tr = simulateTrial(cfg, cfg.targetDur, true, cfg.seedTarget + k);
    sc = detectorScore(tr.r, cfg);
    on = eventOnsets(sc, thr, cfg);
    t0 = tr.t0; vAll(k) = tr.v;
    inHit = on >= t0 + cfg.hitWin(1) & on <= t0 + cfg.hitWin(2);
    if any(inHit), hit(k) = true; lat(k) = on(find(inHit,1)) - t0; end
    nFAt = nFAt + sum(abs(on - t0) > cfg.assocWin);
    hoursT = hoursT + (cfg.targetDur - cfg.warmup - 2*cfg.assocWin)/3600;
    if k == 1
        exampleT = tr; exampleT.score = sc; exampleT.onsets = on;
    end
end

% ---- 4. metrics
R.config = cfg;
R.threshold = thr;
R.calibration = struct('events',calFA,'hours',calHours);
R.Pd = mean(hit);
[ciLo, ciHi] = wilsonCI(sum(hit), numel(hit)); R.Pd_ci = [ciLo ciHi];
R.nCrossings = numel(hit);
R.FA_per_h = nFA/hoursEmpty;
[ciLo, ciHi] = poissonCI(nFA, hoursEmpty); R.FA_ci = [ciLo ciHi];
R.nFA = nFA; R.hoursEmpty = hoursEmpty;
R.FA_per_h_inTargetRecordings = nFAt/hoursT;
R.latency_median_s = medianOrNaN(lat(hit));
R.latency_p90_s = percentileOrNaN(lat(hit), 90);
% signature size vs noise (for interpretation)
lobe = 2*lambda*deff/cfg.w/mean(cfg.speedRange);
geom = makeGeom(cfg);
yy = linspace(-6,6,2001)';
sP = fsr_scattered_field('paraxial', geom, yy, cfg.z0, lambda);
R.peakRelAmpChange = max(abs(abs(1 - sP) - 1));
R.ampNoiseStd = sqrt(10^(-cfg.cnr_db/10)/2);
R.signatureSNR_db = 20*log10(R.peakRelAmpChange / R.ampNoiseStd);
R.mainLobeDuration_s = lobe;
R.runtime_s = toc;
R.label = 'SYNTHETIC simulation result - not a field measurement';

if cfg.verbose
    fprintf('\n  ---- SYNTHETIC RESULTS (simulated channel, not measured) ----\n');
    fprintf('  Peak relative amplitude change of the clean signal: %.4f (%.2f dB)\n', ...
        R.peakRelAmpChange, 20*log10(1+R.peakRelAmpChange));
    fprintf('  Amplitude noise std per sample: %.4f -> signature SNR %.1f dB per sample\n', ...
        R.ampNoiseStd, R.signatureSNR_db);
    fprintf('  Pd = %.3f  [95%% CI %.3f-%.3f]  (%d crossings)\n', R.Pd, R.Pd_ci(1), R.Pd_ci(2), R.nCrossings);
    fprintf('  False alarms = %.2f /h  [95%% CI %.2f-%.2f]  (%d events in %.1f h of independent empty data)\n', ...
        R.FA_per_h, R.FA_ci(1), R.FA_ci(2), R.nFA, R.hoursEmpty);
    fprintf('  Unassociated alarms in target recordings: %.2f /h\n', R.FA_per_h_inTargetRecordings);
    fprintf('  Detection latency (onset - t0): median %.2f s, 90th percentile %.2f s\n', ...
        R.latency_median_s, R.latency_p90_s);
    fprintf('  (negative latency = alarm before the target centre reaches the LOS; Fresnel side lobes)\n');
    fprintf('  Runtime %.1f s\n', R.runtime_s);
end
if cfg.plot
    plotResults(cfg, R, exampleT, exampleEmpty, lambda);
end
end

%% =====================================================================
%% Signal and channel simulation
%% =====================================================================
function geom = makeGeom(cfg)
geom = struct('src',[0 0 cfg.hAnt],'rcv',[cfg.L 0 cfg.hAnt],'xc',cfg.xc, ...
              'w',cfg.w,'h',cfg.h,'opacity',cfg.opacity);
end

function tr = simulateTrial(cfg, dur, withTarget, seed)
% One recording. Order of random draws is fixed, so a seed reproduces it exactly.
rng(seed, 'twister');
fs = cfg.fs;
n = round(dur*fs);
t = (0:n-1)'/fs;
lambda = fsr_wavelength(cfg.fc);
tr.t = t; tr.t0 = NaN; tr.v = NaN; tr.y = nan(n,1);
if withTarget
    v  = cfg.speedRange(1) + rand*(cfg.speedRange(2) - cfg.speedRange(1));
    t0 = cfg.warmup + 10 + rand*(dur - cfg.warmup - 10 - 15);   % crossing well after warm-up
    y  = v*(t - t0);                                             % along-track offset from the LOS
    s  = fsr_scattered_field('paraxial', makeGeom(cfg), y, cfg.z0, lambda);
    tr.t0 = t0; tr.v = v; tr.y = y;
else
    rand; rand;                                                  % keep the draw order identical
    s = zeros(n,1);
end
E = 1 - s;                                                       % E/E0 = 1 - s (verified physics)
% slow complex interference (off when interf_db = -Inf)
if isfinite(cfg.interf_db)
    c = (ou(n, fs, cfg.interf_tau) + 1i*ou(n, fs, cfg.interf_tau))/sqrt(2) * 10^(cfg.interf_db/20);
else
    c = zeros(n,1);
end
% slow gain drift (off when drift_db = 0)
if cfg.drift_db > 0
    g = 10.^(cfg.drift_db*ou(n, fs, cfg.drift_tau)/20);
else
    g = ones(n,1);
end
rot = exp(1i*2*pi*cfg.cfo_hz*t);
if isfinite(cfg.cnr_db)
    sig = sqrt(10^(-cfg.cnr_db/10)/2);
    noise = sig*(randn(n,1) + 1i*randn(n,1));
else
    noise = zeros(n,1);
end
tr.clean = E;                    % clean normalised field (no channel effects)
tr.r = g.*(E + c).*rot + noise;  % simulated received baseband
end

function x = ou(n, fs, tau)
% Unit-variance Ornstein-Uhlenbeck (AR(1)) sequence, correlation time tau [s].
rho = exp(-1/(fs*tau));
e = randn(n,1)*sqrt(1 - rho^2);
x = filter(1, [1 -rho], e, rho*randn);
end

%% =====================================================================
%% Detector
%% =====================================================================
function [score, eLP, sigma] = detectorScore(r, cfg)
fs = cfg.fs;
a = abs(r);                                   % amplitude: immune to CFO and common phase
al = exp(-1/(fs*cfg.tau_base));
b = filter(1-al, [1 -al], a, al*a(1));        % slow running mean (tracks drift)
e = a./b - 1;                                 % relative amplitude deviation
eHP = e - causalMean(e, max(1,round(cfg.hp_win*fs)));
eLP = causalMean(eHP, max(1,round(cfg.lp_win*fs)));
ms = causalMean(eLP.^2, max(1,round(cfg.cfar_win*fs)));
gN = round(cfg.guard*fs);
if gN > 0
    msd = [repmat(ms(min(gN+1,numel(ms))), gN, 1); ms(1:end-gN)];
else
    msd = ms;
end
sigma = sqrt(max(msd, 1e-12));
score = abs(eLP)./sigma;
score(1:min(numel(score), round(cfg.warmup*fs))) = 0;   % detector not settled yet
end

function m = causalMean(x, n)
% Mean of the last n samples (fewer at the start). Base MATLAB, no toolbox.
c = cumsum([0; x(:)]);
idx = (1:numel(x))';
lo = max(0, idx - n);
m = (c(idx+1) - c(lo+1)) ./ (idx - lo);
m = reshape(m, size(x));
end

function on = eventOnsets(score, thr, cfg)
% Group exceedances: samples above thr closer than mergeGap form one event.
idx = find(score(:) > thr);
if isempty(idx), on = zeros(0,1); return; end
keep = [true; diff(idx) > round(cfg.mergeGap*cfg.fs)];
on = (idx(keep) - 1)/cfg.fs;      % event onset times [s]
end

function [thr, nEv, hours] = calibrateThreshold(scores, cfg)
% Smallest threshold whose event rate on the CALIBRATION data <= target.
hours = numel(scores)*(cfg.calDur - cfg.warmup)/3600;
lo = 0;
hi = 0;
for k = 1:numel(scores), hi = max(hi, max(scores{k})); end
hi = hi*1.001 + 1e-9;
for it = 1:50
    mid = 0.5*(lo + hi);
    nEv = 0;
    for k = 1:numel(scores), nEv = nEv + numel(eventOnsets(scores{k}, mid, cfg)); end
    if nEv/hours > cfg.target_fa_per_h, lo = mid; else, hi = mid; end
end
thr = hi;
nEv = 0;
for k = 1:numel(scores), nEv = nEv + numel(eventOnsets(scores{k}, thr, cfg)); end
end

%% =====================================================================
%% Statistics (base MATLAB only)
%% =====================================================================
function [lo, hi] = wilsonCI(k, n)
if n == 0, lo = NaN; hi = NaN; return; end
z = sqrt(2)*erfinv(0.95);         % 1.95996...
p = k/n; den = 1 + z^2/n;
c = (p + z^2/(2*n))/den;
hw = z*sqrt(p*(1-p)/n + z^2/(4*n^2))/den;
lo = max(0, c - hw); hi = min(1, c + hw);
end

function [lo, hi] = poissonCI(k, hours)
% Exact (Garwood) 95% interval for a Poisson rate, events per hour.
% chi2inv(p, 2k)/2 = gammaincinv(p, k)  (base MATLAB)
if k == 0, lo = 0; else, lo = gammaincinv(0.025, k); end
hi = gammaincinv(0.975, k + 1);
lo = lo/hours; hi = hi/hours;
end

function m = medianOrNaN(x)
if isempty(x), m = NaN; else, m = median(x); end
end

function p = percentileOrNaN(x, q)
% Linear-interpolation percentile without the Statistics Toolbox.
x = sort(x(:));
if isempty(x), p = NaN; return; end
pos = 1 + (numel(x)-1)*q/100;
i0 = floor(pos); f = pos - i0;
if i0 >= numel(x), p = x(end); else, p = x(i0) + f*(x(i0+1) - x(i0)); end
end

%% =====================================================================
%% Figures
%% =====================================================================
function plotResults(cfg, R, ex, exE, lambda)
fig = figure('Name','RF Drone Fence - synthetic detection (simulation)','Color','w', ...
    'Position',[60 60 1250 820]);
tl = tiledlayout(fig, 3, 2, 'TileSpacing','compact','Padding','compact');
title(tl, sprintf('SIMULATED detection on a synthetic channel - %.1f GHz, CNR %g dB (not field data)', ...
    cfg.fc/1e9, cfg.cnr_db), 'FontWeight','bold');
tWin = ex.t0 + [-4 4];
m = ex.t >= tWin(1) & ex.t <= tWin(2);

% (1) geometry
ax = nexttile(tl, 1); hold(ax,'on'); grid(ax,'on'); view(ax, 35, 22);
plot3(ax, 0, 0, cfg.hAnt, 'r^', 'MarkerSize', 9, 'MarkerFaceColor','r');
plot3(ax, cfg.L, 0, cfg.hAnt, 'bs', 'MarkerSize', 9, 'MarkerFaceColor','b');
plot3(ax, [0 cfg.L], [0 0], [cfg.hAnt cfg.hAnt], 'k-', 'LineWidth', 1.2);
yTr = [-8 8];
plot3(ax, [cfg.xc cfg.xc], yTr, [1 1]*(cfg.hAnt + cfg.z0), 'm--', 'LineWidth', 1.2);
yc = 0; zc = cfg.hAnt + cfg.z0;
fill3(ax, cfg.xc*[1 1 1 1], yc + cfg.w/2*[-1 1 1 -1], zc + cfg.h/2*[-1 -1 1 1], ...
    [0.85 0.4 0.1], 'FaceAlpha', 0.8, 'EdgeColor','k');
xlabel(ax,'x along baseline [m]'); ylabel(ax,'y along track [m]'); zlabel(ax,'z [m]');
legend(ax, {'Tx','Rx','LOS baseline','target path','screen at LOS crossing'}, 'Location','northeast');
title(ax, sprintf('Geometry: L = %.0f m, h_{ant} = %.1f m, crossing at x = %.0f m', cfg.L, cfg.hAnt, cfg.xc));
zlim(ax, [0 3]); axis(ax,'tight');

% (2) target path vs time
ax = nexttile(tl, 2); hold(ax,'on'); grid(ax,'on');
plot(ax, ex.t, ex.y, 'm-', 'LineWidth', 1.2);
xline(ax, ex.t0, 'k--', 't_0 (on LOS)');
yline(ax, 0, 'k:');
xlim(ax, ex.t0 + [-10 10]); ylim(ax, [-40 40]);
xlabel(ax,'time [s]'); ylabel(ax,'y(t) = v (t - t_0) [m]');
title(ax, sprintf('Example crossing: v = %.2f m/s, t_0 = %.2f s', ex.v, ex.t0));

% (3) clean vs noisy signal
ax = nexttile(tl, 3); hold(ax,'on'); grid(ax,'on');
rn = abs(ex.r(m)); rn = rn/mean(rn);
plot(ax, ex.t(m), 20*log10(rn), 'Color',[0.6 0.6 0.6], 'LineWidth', 0.7);
plot(ax, ex.t(m), 20*log10(abs(ex.clean(m))), 'b-', 'LineWidth', 1.6);
xline(ax, ex.t0, 'k--');
xlabel(ax,'time [s]'); ylabel(ax,'amplitude [dB re mean]');
legend(ax, {'simulated noisy |r| (normalised)','clean |E/E_0| (paraxial model)'}, 'Location','southwest');
title(ax, 'Received amplitude: clean model vs simulated noisy signal');

% (4) detection score, threshold and events (example crossing)
ax = nexttile(tl, 4); hold(ax,'on'); grid(ax,'on');
plot(ax, ex.t(m), ex.score(m), 'Color',[0.1 0.35 0.6]);
yline(ax, R.threshold, 'r-', sprintf('threshold %.2f', R.threshold), 'LineWidth', 1.2);
xline(ax, ex.t0, 'k--');
onIn = ex.onsets(ex.onsets >= tWin(1) & ex.onsets <= tWin(2));
for k = 1:numel(onIn), xline(ax, onIn(k), 'g-', 'LineWidth', 2); end
patch(ax, ex.t0 + [cfg.hitWin(1) cfg.hitWin(2) cfg.hitWin(2) cfg.hitWin(1)], ...
    [0 0 1 1]*max([ex.score(m); R.threshold])*1.05, [0.2 0.8 0.2], 'FaceAlpha', 0.08, 'EdgeColor','none');
xlabel(ax,'time [s]'); ylabel(ax,'score = |e_{lp}| / \sigma');
title(ax, 'Detection score near the crossing (green = event onset, band = hit window)');

% (5) score on an independent target-free recording
ax = nexttile(tl, 5); hold(ax,'on'); grid(ax,'on');
plot(ax, exE.t, exE.score, 'Color',[0.45 0.45 0.45]);
yline(ax, R.threshold, 'r-', 'LineWidth', 1.2);
for k = 1:numel(exE.onsets), xline(ax, exE.onsets(k), 'm-', 'LineWidth', 1.5); end
xlabel(ax,'time [s]'); ylabel(ax,'score');
title(ax, sprintf('Independent empty recording (%d false-alarm events shown in magenta)', numel(exE.onsets)));

% (6) metrics summary
ax = nexttile(tl, 6); axis(ax,'off');
txt = {
 'SYNTHETIC RESULTS - simulated channel, not field-validated'
 ''
 sprintf('Threshold (from %.1f h target-free calibration): %.3f', R.calibration.hours, R.threshold)
 sprintf('Pd = %.3f   [95%% CI %.3f - %.3f], %d crossings', R.Pd, R.Pd_ci(1), R.Pd_ci(2), R.nCrossings)
 sprintf('False alarms = %.2f /h   [95%% CI %.2f - %.2f], %.1f h empty test', R.FA_per_h, R.FA_ci(1), R.FA_ci(2), R.hoursEmpty)
 sprintf('Latency (onset - t_0): median %.2f s, p90 %.2f s', R.latency_median_s, R.latency_p90_s)
 sprintf('Peak clean change %.4f; signature SNR %.1f dB/sample', R.peakRelAmpChange, R.signatureSNR_db)
 sprintf('Main-lobe duration ~ %.2f s at %.1f m/s', R.mainLobeDuration_s, mean(cfg.speedRange))
 ''
 'Assumptions: flat opaque screen, scalar paraxial diffraction, straight'
 'perpendicular crossing, complex white Gaussian noise, simple slow drift.'
 'No birds, people, wind clutter, multipath or hardware effects.'
 };
text(ax, 0, 1, txt, 'VerticalAlignment','top', 'FontName','FixedWidth', 'FontSize', 9);
end

%% =====================================================================
%% Basic tests
%% =====================================================================
function T = runTests()
fprintf('\nfsr_detection_demo - basic tests (synthetic data only)\n');
T = struct('name',{},'pass',{},'detail',{},'limits',{});
base = defaultConfig();
base.plot = false; base.verbose = false;
base.calTrials = 12; base.emptyTrials = 12; base.targetTrials = 60;   % ~2 h + 2 h + 60 crossings

% T0: noise-free, target-free -> received field equals the direct field exactly
c0 = base; c0.cnr_db = Inf; c0.cfo_hz = 0;
tr = simulateTrial(c0, 120, false, 1);
dev = max(abs(tr.r - 1));
sc = detectorScore(tr.r, c0);
T = rec(T, 'T0 no target, no noise: E/E0 = 1 and score = 0', dev < 1e-12 && max(sc) < 1e-6, ...
    sprintf('max|r-1| = %.1e, max score = %.1e', dev, max(sc)), ...
    'Checks bookkeeping only (E/E0 = 1 - s with s = 0); says nothing about noise.');

% T0b: clean target signal reproduces the verified physics value
tr = simulateTrial(c0, 120, true, 2);
pk = max(abs(abs(tr.clean) - 1));
T = rec(T, 'T0b clean crossing: peak |dA|/A matches the verified 0.0856', abs(pk - 0.0856) < 1.5e-3, ...
    sprintf('peak = %.4f (sampled at %g Hz)', pk, c0.fs), ...
    'Consistency with fsr_validate_physics Stage 7 only; sampling can lower the peak slightly.');

% T1: no target, noise present -> false-alarm rate controlled on independent data
c1 = base; c1.cnr_db = 30;
R1 = runExperiment(c1);
T = rec(T, 'T1 no target (CNR 30 dB): independent FA rate <= 2 /h', R1.FA_per_h <= 2.0, ...
    sprintf('%.2f /h [%.2f-%.2f], %d events in %.1f h; threshold %.2f', R1.FA_per_h, R1.FA_ci(1), R1.FA_ci(2), R1.nFA, R1.hoursEmpty, R1.threshold), ...
    ['Statistical: with ~2 h of test data a correct detector can still exceed 2 /h by chance ' ...
     '(rare). Gaussian noise only; real clutter is not Gaussian or stationary.']);

% T2: target with weak noise -> almost always detected
c2 = base; c2.cnr_db = 40;
R2 = runExperiment(c2);
T = rec(T, 'T2 target + weak noise (CNR 40 dB): Pd >= 0.95 and FA <= 2 /h', R2.Pd >= 0.95 && R2.FA_per_h <= 2.0, ...
    sprintf('Pd = %.3f [%.3f-%.3f], FA = %.2f /h, median latency %.2f s, signature SNR %.1f dB', ...
    R2.Pd, R2.Pd_ci(1), R2.Pd_ci(2), R2.FA_per_h, R2.latency_median_s, R2.signatureSNR_db), ...
    'Nominal opaque 0.30 x 0.20 m screen on the LOS at mid-baseline only; best case.');

% T3: target with strong noise -> FA still controlled, Pd degrades
c3 = base; c3.cnr_db = 10;
R3 = runExperiment(c3);
T = rec(T, 'T3 target + strong noise (CNR 10 dB): FA <= 2 /h and Pd below the weak-noise Pd', ...
    R3.FA_per_h <= 2.0 && R3.Pd < R2.Pd, ...
    sprintf('Pd = %.3f [%.3f-%.3f], FA = %.2f /h, signature SNR %.1f dB', ...
    R3.Pd, R3.Pd_ci(1), R3.Pd_ci(2), R3.FA_per_h, R3.signatureSNR_db), ...
    ['Shows that calibration keeps false alarms in check while Pd falls; ' ...
     'the exact Pd value is not a requirement and is not field-relevant.']);

% T4 (info): slow drift + slow interference, CNR 30 dB
c4 = base; c4.cnr_db = 30; c4.drift_db = 0.3; c4.interf_db = -30;
R4 = runExperiment(c4);
T = rec(T, 'T4 (info) CNR 30 dB + 0.3 dB drift + -30 dB slow interference', true, ...
    sprintf('Pd = %.3f, FA = %.2f /h', R4.Pd, R4.FA_per_h), ...
    'Reported only: a simple synthetic disturbance, not a model of a real site.');

fprintf('\n================ TEST SUMMARY ================\n');
for i = 1:numel(T)
    if T(i).pass, tag = 'PASS'; else, tag = 'FAIL'; end
    fprintf('[%s] %s\n       %s\n       limits: %s\n', tag, T(i).name, T(i).detail, T(i).limits);
end
nF = sum(~[T.pass]);
fprintf('\n%d tests, %d failed. All results are synthetic, not field validation.\n', numel(T), nF);
if nF > 0, error('fsr:tests','%d detection tests failed.', nF); end
end

function T = rec(T, name, pass, detail, limits)
T(end+1) = struct('name',name,'pass',logical(pass),'detail',detail,'limits',limits);
if pass, tag = 'PASS'; else, tag = 'FAIL'; end
fprintf('  [%s] %s  (%s)\n', tag, name, detail);
end
