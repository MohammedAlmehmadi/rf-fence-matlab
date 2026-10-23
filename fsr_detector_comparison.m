function varargout = fsr_detector_comparison(arg)
%FSR_DETECTOR_COMPARISON Compare two crossing detectors on identical synthetic data.
%
% USAGE
%   fsr_detector_comparison                 % full grid: CNR {30,20,15} dB x z0 {0,0.5,1,2} m
%   fsr_detector_comparison('quick')        % smaller run (fewer hours / crossings) to try it out
%   fsr_detector_comparison('tests')        % leakage / causality / sanity / channel-integration tests
%   cfg.targetTrials = 200; fsr_detector_comparison(cfg)   % override config fields
%   R = fsr_detector_comparison(...);       % return all numbers
%   api = fsr_detector_comparison('api');   % handles to the protocol functions (used by the app)
%
% DETECTORS
%   A  "amplitude deviation + noise estimate": identical to fsr_detection_demo.m
%      score_A = |e_lp| / sigma, where e = |r|/slow mean - 1, e_hp = e - mean(e, 0.25 s),
%      e_lp = mean(e_hp, 0.04 s), sigma = causal RMS of e_lp over 30 s (2 s guard).
%   B  "matched-filter bank": sliding normalised correlation of e_hp with templates
%      built ONLY from fsr_scattered_field('paraxial', ...).
%        * Unknown crossing time  -> the correlation slides over every sample (causal).
%        * Unknown speed          -> a bank of templates for cfg.bankSpeeds; max over bank.
%        * Unknown vertical offset -> in the paraxial model s(t;z0) = G_z(z0)/G_z(0) * s(t;0)
%          (a constant complex factor), and to first order |1-s|-1 ~ -Re(s). So every z0
%          lies in the 2-D subspace {Re s0(t), Im s0(t)}; the statistic is the energy of
%          the projection onto that subspace (non-coherent, phase unknown).
%      Each template column passes through the SAME high-pass as the data, is Hann-tapered
%      and orthonormalised. Output CFAR: T_k / causal mean of T_k over 30 s (2 s guard).
%      score_B = max_k of that ratio.
%   Neither detector receives t0, the true speed or z0. Both use only the received r(t).
%
% FAIR PROTOCOL (same for A and B)
%   * Same recordings: every simulated r(t) is scored by both detectors.
%   * Thresholds: per CNR, calibrated on target-free calibration recordings only, each
%     detector to cfg.target_fa_per_h. Then frozen.
%   * Test: independent target-free recordings (FA/h) and target recordings (Pd), with
%     seeds disjoint from calibration. Same event merging, hit window and association.
%   * Reported: Pd with Wilson 95% CI, FA/h with exact Poisson 95% CI, hours of
%     target-free data, number of crossings, latency.
%
% ASSUMPTIONS / LIMITS (also printed)
%   Geometry fixed: L = 20 m, antennas 1.5 m, crossing at xc = 10 m (the bank assumes this
%   xc; an unknown xc would need templates over xc too). Flat opaque 0.30 x 0.20 m screen,
%   paraxial scalar diffraction, straight perpendicular crossings, speeds uniform in
%   cfg.testSpeedRange, complex white Gaussian noise, 3 Hz residual CFO. No clutter,
%   birds, people or hardware effects.
%   ALL METRICS ARE FROM SYNTHETIC SIMULATION, NOT FIELD MEASUREMENTS.
%
% CHANNEL (v3)
%   The received signal is produced by fsr_channel_model.m (same equations, same random
%   draw order as the v2 channel). Multipath is OFF by default (cfg.paths empty); with it
%   off the recordings are bit-identical to v2 (tests C1, C2). To enable reflected paths:
%     cfg = struct(); cfg.paths = struct('gain',-0.7,'excess_m',0.2237,'modulation','static');
%     fsr_detector_comparison(cfg)
%   (assign cfg.paths with dot syntax). Path parameters are USER ASSUMPTIONS, not measured.
%   Calibration, empty-test and target recordings all use the same cfg, hence the same
%   channel. The detector-B templates still come from the direct-path physics only.
%   cfg.channelImpl = 'legacy' runs the frozen v2 channel code; it exists for regression
%   tests only.
%
% DEPENDENCIES (same folder): fsr_channel_model.m, fsr_scattered_field.m,
%   fsr_babinet_paraxial.m, fsr_fresnel_cs.m, fsr_wavelength.m. Base MATLAB only.
%   fsr_detection_demo.m is NOT called or modified; detector A is copied here verbatim
%   so this file is self-contained.

if nargin < 1, arg = struct(); end
if ischar(arg) || isstring(arg)
    switch lower(char(arg))
        case 'tests'
            out = runTests();
            if nargout > 0, varargout{1} = out; end
            return
        case 'api'
            % Read-only access to the SAME functions used by the protocol (for the app).
            varargout{1} = apiHandles();
            return
        case 'quick'
            arg = struct('calTrials',12,'emptyTrials',12,'targetTrials',60);
        case 'full'
            arg = struct();
        otherwise
            error('fsr:arg','Unknown option ''%s''. Use ''quick'', ''tests'' or a config struct.', char(arg));
    end
end
cfg = mergeConfig(defaultConfig(), arg);
R = runComparison(cfg);
if nargout > 0, varargout{1} = R; end
end

%% =====================================================================
%% Configuration
%% =====================================================================
function cfg = defaultConfig()
% geometry (verified physics setup)
cfg.fc = 2.4e9; cfg.L = 20.0; cfg.hAnt = 1.5; cfg.xc = 10.0;
cfg.w = 0.30; cfg.h = 0.20; cfg.opacity = 1.0;
cfg.fs = 250;
% comparison grid
cfg.cnrList = [30 20 15];          % dB
cfg.z0List  = [0 0.5 1 2];         % m, vertical offset of target centre from the LOS
cfg.testSpeedRange = [3 8];        % m/s, true speed drawn uniformly per crossing
cfg.bankSpeeds = [2.5 3.5 5 7 10]; % m/s, templates of detector B (fixed BEFORE any test data)
% channel (same as fsr_detection_demo defaults)
cfg.drift_db = 0; cfg.drift_tau = 60; cfg.interf_db = -Inf; cfg.interf_tau = 5; cfg.cfo_hz = 3;
cfg.paths = struct([]);            % reflected paths for fsr_channel_model (empty = off = v2)
cfg.channelImpl = 'model';         % 'model' (fsr_channel_model) | 'legacy' (frozen v2, tests only)
% detector A parameters (same as fsr_detection_demo)
cfg.tau_base = 3.0; cfg.hp_win = 0.25; cfg.lp_win = 0.04; cfg.cfar_win = 30; cfg.guard = 2;
cfg.warmup = 35;
% detector B parameters
cfg.tplSpan = 1.5;                 % template covers |y| <= tplSpan * lambda*d_eff/w
% protocol
cfg.target_fa_per_h = 0.5;
cfg.calTrials = 24;   cfg.calDur = 600;     % per CNR: 24 x 10 min = 3.8 h usable
cfg.emptyTrials = 24; cfg.emptyDur = 600;   % per CNR: independent 3.8 h usable
cfg.targetTrials = 120; cfg.targetDur = 120;% per (CNR, z0): 120 crossings
cfg.seedBaseCal = 100000; cfg.seedBaseEmpty = 200000; cfg.seedBaseTarget = 300000;
cfg.mergeGap = 1.0; cfg.hitWin = [-1.5 3.0]; cfg.assocWin = 6.0;
cfg.plot = true; cfg.verbose = true;
end

function cfg = mergeConfig(cfg, user)
if ~isstruct(user), error('fsr:arg','Config must be a struct.'); end
f = fieldnames(user);
for i = 1:numel(f)
    if ~isfield(cfg, f{i}), error('fsr:arg','Unknown config field ''%s''.', f{i}); end
    cfg.(f{i}) = user.(f{i});
end
if max([cfg.calTrials cfg.emptyTrials cfg.targetTrials]) >= 1000
    error('fsr:arg','Trial counts must be < 1000 so seed ranges stay disjoint.');
end
if ~any(strcmp(char(cfg.channelImpl), {'model','legacy'}))
    error('fsr:arg','cfg.channelImpl must be ''model'' or ''legacy''.');
end
if ~isempty(cfg.paths) && ~isstruct(cfg.paths)
    error('fsr:arg','cfg.paths must be empty or a struct array (see fsr_channel_model).');
end
end

function s = seedOf(cfg, kind, iCnr, iZ, k)
% Disjoint seed ranges: calibration 1e5.., empty test 2e5.., target test 3e5..
switch kind
    case 'cal',    s = cfg.seedBaseCal    + 1000*iCnr + k;
    case 'empty',  s = cfg.seedBaseEmpty  + 1000*iCnr + k;
    case 'target', s = cfg.seedBaseTarget + 10000*iCnr + 1000*iZ + k;
end
end

%% =====================================================================
%% Main comparison
%% =====================================================================
function R = runComparison(cfg)
tic;
bank = buildBank(cfg);           % built from physics + cfg only, before any data exists
nC = numel(cfg.cnrList); nZ = numel(cfg.z0List);
calH = cfg.calTrials*(cfg.calDur - cfg.warmup)/3600;
empH = cfg.emptyTrials*(cfg.emptyDur - cfg.warmup)/3600;
if cfg.verbose
    fprintf('\nRF Drone Fence - detector comparison on SYNTHETIC data (not field measurements)\n');
    fprintf('  Geometry: L = %.0f m, h_ant = %.1f m, xc = %.0f m, screen %.2f x %.2f m, %.1f GHz\n', ...
        cfg.L, cfg.hAnt, cfg.xc, cfg.w, cfg.h, cfg.fc/1e9);
    fprintf('  True speed ~ U(%g, %g) m/s; bank speeds %s m/s; t0, v, z0 unknown to both detectors\n', ...
        cfg.testSpeedRange(1), cfg.testSpeedRange(2), mat2str(cfg.bankSpeeds));
    fprintf('  Per CNR: calibration %.2f h target-free, independent test %.2f h target-free\n', calH, empH);
    fprintf('  Per (CNR, z0): %d crossings\n', cfg.targetTrials);
    fprintf('  Channel (%s): CFO %g Hz, drift %g dB, interference %g dB, %s\n\n', ...
        char(cfg.channelImpl), cfg.cfo_hz, cfg.drift_db, cfg.interf_db, multipathText(cfg));
end
R = struct();
R.cnr = cfg.cnrList; R.z0 = cfg.z0List;
R.thrA = zeros(1,nC); R.thrB = zeros(1,nC);
R.faA = zeros(1,nC); R.faB = zeros(1,nC); R.faA_ci = zeros(nC,2); R.faB_ci = zeros(nC,2);
R.nfaA = zeros(1,nC); R.nfaB = zeros(1,nC);
R.pdA = zeros(nC,nZ); R.pdB = zeros(nC,nZ);
R.pdA_ci = zeros(nC,nZ,2); R.pdB_ci = zeros(nC,nZ,2);
R.latA = nan(nC,nZ); R.latB = nan(nC,nZ);
R.calHours = calH; R.emptyHours = empH; R.crossings = cfg.targetTrials;
R.channelImpl = char(cfg.channelImpl); R.multipath = multipathText(cfg);
R.paths = cfg.paths;
R.cfg = cfg;                     % full configuration, for reproducibility and the app
example = [];

for iC = 1:nC
    c = cfg; c.cnr_db = cfg.cnrList(iC);
    % ---- calibration: target-free only
    sA = cell(c.calTrials,1); sB = cell(c.calTrials,1);
    for k = 1:c.calTrials
        tr = simulateTrial(c, c.calDur, false, 0, seedOf(cfg,'cal',iC,0,k));
        sA{k} = scoreA(tr.r, c); sB{k} = scoreB(tr.r, bank, c);
    end
    R.thrA(iC) = calibrateThreshold(sA, c);
    R.thrB(iC) = calibrateThreshold(sB, c);
    % ---- independent target-free test
    nA = 0; nB = 0;
    for k = 1:c.emptyTrials
        tr = simulateTrial(c, c.emptyDur, false, 0, seedOf(cfg,'empty',iC,0,k));
        nA = nA + numel(eventOnsets(scoreA(tr.r, c), R.thrA(iC), c));
        nB = nB + numel(eventOnsets(scoreB(tr.r, bank, c), R.thrB(iC), c));
    end
    R.nfaA(iC) = nA; R.nfaB(iC) = nB;
    R.faA(iC) = nA/empH; R.faB(iC) = nB/empH;
    [lo, hi] = poissonCI(nA, empH); R.faA_ci(iC,:) = [lo hi];
    [lo, hi] = poissonCI(nB, empH); R.faB_ci(iC,:) = [lo hi];
    % ---- target test, per z0
    for iZ = 1:nZ
        c.z0 = cfg.z0List(iZ);
        hA = false(c.targetTrials,1); hB = hA; lA = nan(size(hA)); lB = lA;
        for k = 1:c.targetTrials
            tr = simulateTrial(c, c.targetDur, true, c.z0, seedOf(cfg,'target',iC,iZ,k));
            onA = eventOnsets(scoreA(tr.r, c), R.thrA(iC), c);   % decisions use r only ...
            onB = eventOnsets(scoreB(tr.r, bank, c), R.thrB(iC), c);
            [hA(k), lA(k)] = scoreCrossing(onA, tr.t0, c);        % ... truth used only here
            [hB(k), lB(k)] = scoreCrossing(onB, tr.t0, c);
            if iC == 1 && iZ == 1 && k == 1
                example = tr; example.sA = scoreA(tr.r, c); example.sB = scoreB(tr.r, bank, c);
            end
        end
        R.pdA(iC,iZ) = mean(hA); R.pdB(iC,iZ) = mean(hB);
        [lo, hi] = wilsonCI(sum(hA), numel(hA)); R.pdA_ci(iC,iZ,:) = [lo hi];
        [lo, hi] = wilsonCI(sum(hB), numel(hB)); R.pdB_ci(iC,iZ,:) = [lo hi];
        R.latA(iC,iZ) = medianOrNaN(lA(hA)); R.latB(iC,iZ) = medianOrNaN(lB(hB));
    end
    if cfg.verbose
        fprintf('  CNR %2d dB done: thr A = %.2f, thr B = %.2f (%.0f s elapsed)\n', ...
            c.cnr_db, R.thrA(iC), R.thrB(iC), toc);
    end
end
R.runtime_s = toc;
R.label = 'SYNTHETIC simulation results - not field measurements';
if cfg.verbose, printTable(cfg, R); end
if cfg.plot, plotComparison(cfg, R, bank, example); end
end

function [hit, lat] = scoreCrossing(on, t0, cfg)
% Scoring step: the ONLY place where the true crossing time is used.
inHit = on >= t0 + cfg.hitWin(1) & on <= t0 + cfg.hitWin(2);
hit = any(inHit);
if hit, lat = on(find(inHit,1)) - t0; else, lat = NaN; end
end

%% =====================================================================
%% Channel: target signature from the physics, channel from fsr_channel_model
%% =====================================================================
function geom = makeGeom(cfg)
geom = struct('src',[0 0 cfg.hAnt],'rcv',[cfg.L 0 cfg.hAnt],'xc',cfg.xc, ...
              'w',cfg.w,'h',cfg.h,'opacity',cfg.opacity);
end

function api = apiHandles()
% Handles to local functions of this file. Nothing is copied: the app calls exactly the
% code that the protocol and its tests use.
api = struct('defaultConfig', @defaultConfig, 'seedOf', @seedOf, ...
    'simulateTrial', @simulateTrial, 'scoreA', @scoreA, 'scoreB', @scoreB, ...
    'buildBank', @buildBank, 'calibrateThreshold', @calibrateThreshold, ...
    'eventOnsets', @eventOnsets, 'scoreCrossing', @scoreCrossing, 'makeGeom', @makeGeom);
end

function txt = multipathText(cfg)
if isempty(cfg.paths)
    txt = 'multipath off (direct path only)';
else
    txt = sprintf('multipath ON: %d reflected path(s), ASSUMED parameters', numel(cfg.paths));
end
end

function tr = simulateTrial(cfg, dur, withTarget, z0, seed)
% Draw order: target (2 x rand), then fsr_channel_model draws interference, drift,
% noise from the same stream (its seed field is left empty on purpose). Noise and
% interference are added ONLY inside fsr_channel_model.
if strcmp(char(cfg.channelImpl), 'legacy')
    tr = legacySimulateTrial(cfg, dur, withTarget, z0, seed);
    return
end
rng(seed, 'twister');
fs = cfg.fs; n = round(dur*fs); t = (0:n-1)'/fs;
lambda = fsr_wavelength(cfg.fc);
tr.t = t; tr.t0 = NaN; tr.v = NaN; tr.z0 = NaN; tr.seed = seed;
if withTarget
    v  = cfg.testSpeedRange(1) + rand*(cfg.testSpeedRange(2) - cfg.testSpeedRange(1));
    t0 = cfg.warmup + 10 + rand*(dur - cfg.warmup - 10 - 15);
    s  = fsr_scattered_field('paraxial', makeGeom(cfg), v*(t - t0), z0, lambda);
    tr.t0 = t0; tr.v = v; tr.z0 = z0;
else
    rand; rand;
    s = zeros(n,1);
end
E = 1 - s;
ch = struct('fs', fs, 'fc', cfg.fc, 'drift_db', cfg.drift_db, 'drift_tau', cfg.drift_tau, ...
    'interf_db', cfg.interf_db, 'interf_tau', cfg.interf_tau, 'cfo_hz', cfg.cfo_hz, ...
    'cnr_db', cfg.cnr_db);
ch.paths = cfg.paths;                       % dot assignment: keeps ch a scalar struct
tr.clean = E;                               % direct-path field before the channel
tr.r = fsr_channel_model(E, t, ch);
end

function tr = legacySimulateTrial(cfg, dur, withTarget, z0, seed)
% FROZEN copy of the v2 simulateTrial (copied from fsr_detection_demo.m). Do not edit:
% it is the reference that tests C1 and C2 compare the new channel against.
rng(seed, 'twister');
fs = cfg.fs; n = round(dur*fs); t = (0:n-1)'/fs;
lambda = fsr_wavelength(cfg.fc);
tr.t = t; tr.t0 = NaN; tr.v = NaN; tr.z0 = NaN; tr.seed = seed;
if withTarget
    v  = cfg.testSpeedRange(1) + rand*(cfg.testSpeedRange(2) - cfg.testSpeedRange(1));
    t0 = cfg.warmup + 10 + rand*(dur - cfg.warmup - 10 - 15);
    s  = fsr_scattered_field('paraxial', makeGeom(cfg), v*(t - t0), z0, lambda);
    tr.t0 = t0; tr.v = v; tr.z0 = z0;
else
    rand; rand;
    s = zeros(n,1);
end
E = 1 - s;
if isfinite(cfg.interf_db)
    c = (ou(n, fs, cfg.interf_tau) + 1i*ou(n, fs, cfg.interf_tau))/sqrt(2)*10^(cfg.interf_db/20);
else
    c = zeros(n,1);
end
if cfg.drift_db > 0, g = 10.^(cfg.drift_db*ou(n, fs, cfg.drift_tau)/20); else, g = ones(n,1); end
rot = exp(1i*2*pi*cfg.cfo_hz*t);
if isfinite(cfg.cnr_db)
    sig = sqrt(10^(-cfg.cnr_db/10)/2); noise = sig*(randn(n,1) + 1i*randn(n,1));
else
    noise = zeros(n,1);
end
tr.clean = E;
tr.r = g.*(E + c).*rot + noise;
end

function x = ou(n, fs, tau)
rho = exp(-1/(fs*tau));
e = randn(n,1)*sqrt(1 - rho^2);
x = filter(1, [1 -rho], e, rho*randn);
end

%% =====================================================================
%% Shared front end
%% =====================================================================
function eHP = frontEnd(r, cfg)
fs = cfg.fs;
a = abs(r);
al = exp(-1/(fs*cfg.tau_base));
b = filter(1-al, [1 -al], a, al*a(1));
e = a./b - 1;
eHP = highPass(e, cfg);
end

function y = highPass(x, cfg)
y = x - causalMean(x, max(1, round(cfg.hp_win*cfg.fs)));
end

function m = cfarMean(x, cfg)
% causal mean over cfar_win, delayed by guard
ms = causalMean(x, max(1, round(cfg.cfar_win*cfg.fs)));
gN = round(cfg.guard*cfg.fs);
if gN > 0, m = [repmat(ms(min(gN+1,numel(ms))), gN, 1); ms(1:end-gN)]; else, m = ms; end
end

%% =====================================================================
%% Detector A (identical to fsr_detection_demo.m)
%% =====================================================================
function score = scoreA(r, cfg)
eHP = frontEnd(r, cfg);
eLP = causalMean(eHP, max(1, round(cfg.lp_win*cfg.fs)));
sigma = sqrt(max(cfarMean(eLP.^2, cfg), 1e-12));
score = abs(eLP)./sigma;
score(1:min(numel(score), round(cfg.warmup*cfg.fs))) = 0;
end

%% =====================================================================
%% Detector B: matched-filter bank from the physics model
%% =====================================================================
function bank = buildBank(cfg)
% Uses only cfg (geometry, bank speeds) and fsr_scattered_field. No trial data.
lambda = fsr_wavelength(cfg.fc);
deff = cfg.xc*(cfg.L - cfg.xc)/cfg.L;
yNull = lambda*deff/cfg.w;
geom = makeGeom(cfg);
bank = struct('v',{},'Q',{},'len',{});
for k = 1:numel(cfg.bankSpeeds)
    v = cfg.bankSpeeds(k);
    tt = (-cfg.tplSpan*yNull/v : 1/cfg.fs : cfg.tplSpan*yNull/v)';
    s0 = fsr_scattered_field('paraxial', geom, v*tt, 0, lambda);   % z0 = 0 reference shape
    B = [real(s0) imag(s0)];
    pad = zeros(round(cfg.hp_win*cfg.fs), 2);
    Bp = [pad; B];                           % let the high-pass settle, as in the data
    Bp = [highPass(Bp(:,1), cfg) highPass(Bp(:,2), cfg)];
    B = Bp(size(pad,1)+1:end, :);
    nT = size(B,1);
    win = 0.5 - 0.5*cos(2*pi*(0:nT-1)'/(nT-1));   % Hann taper
    B = B .* win;
    B = B - mean(B,1);
    [Q, ~] = qr(B, 0);                        % orthonormal 2-D signal subspace
    bank(k).v = v; bank(k).Q = Q; bank(k).len = nT;
end
end

function score = scoreB(r, bank, cfg)
eHP = frontEnd(r, cfg);
score = zeros(size(eHP));
for k = 1:numel(bank)
    Q = bank(k).Q;
    c1 = causalCorr(eHP, Q(:,1));
    c2 = causalCorr(eHP, Q(:,2));
    T = c1.^2 + c2.^2;                        % energy in the signal subspace
    Tn = T ./ max(cfarMean(T, cfg), 1e-20);   % output CFAR normalisation
    score = max(score, Tn);
end
score(1:min(numel(score), round(cfg.warmup*cfg.fs))) = 0;
end

function y = causalCorr(x, q)
% y(n) = sum_m q(m) x(n - M + m): correlation with the template ending at sample n.
% Linear (zero-padded) FFT convolution, so no future samples are used.
h = flipud(q(:));
nOut = numel(x) + numel(h) - 1;
N = 2^nextpow2(nOut);
y = real(ifft(fft(x(:), N) .* fft(h, N)));
y = y(1:numel(x));
end

%% =====================================================================
%% Events, calibration, statistics
%% =====================================================================
function m = causalMean(x, n)
c = cumsum([0; x(:)]);
idx = (1:numel(x))';
lo = max(0, idx - n);
m = (c(idx+1) - c(lo+1)) ./ (idx - lo);
m = reshape(m, size(x));
end

function on = eventOnsets(score, thr, cfg)
idx = find(score(:) > thr);
if isempty(idx), on = zeros(0,1); return; end
keep = [true; diff(idx) > round(cfg.mergeGap*cfg.fs)];
on = (idx(keep) - 1)/cfg.fs;
end

function thr = calibrateThreshold(scores, cfg)
% Target-free calibration recordings ONLY: smallest threshold with rate <= target.
hours = numel(scores)*(cfg.calDur - cfg.warmup)/3600;
lo = 0; hi = 0;
for k = 1:numel(scores), hi = max(hi, max(scores{k})); end
hi = hi*1.001 + 1e-9;
for it = 1:50
    mid = 0.5*(lo + hi); nEv = 0;
    for k = 1:numel(scores), nEv = nEv + numel(eventOnsets(scores{k}, mid, cfg)); end
    if nEv/hours > cfg.target_fa_per_h, lo = mid; else, hi = mid; end
end
thr = hi;
end

function [lo, hi] = wilsonCI(k, n)
if n == 0, lo = NaN; hi = NaN; return; end
z = sqrt(2)*erfinv(0.95);
p = k/n; den = 1 + z^2/n; c = (p + z^2/(2*n))/den;
hw = z*sqrt(p*(1-p)/n + z^2/(4*n^2))/den;
lo = max(0, c - hw); hi = min(1, c + hw);
end

function [lo, hi] = poissonCI(k, hours)
if k == 0, lo = 0; else, lo = gammaincinv(0.025, k); end
hi = gammaincinv(0.975, k + 1);
lo = lo/hours; hi = hi/hours;
end

function m = medianOrNaN(x)
if isempty(x), m = NaN; else, m = median(x); end
end

%% =====================================================================
%% Output
%% =====================================================================
function printTable(cfg, R)
fprintf('\n=========== SYNTHETIC RESULTS (simulation only, not field data) ===========\n');
fprintf('Target-free data per CNR: calibration %.2f h, independent test %.2f h. Crossings per cell: %d.\n', ...
    R.calHours, R.emptyHours, R.crossings);
fprintf('Thresholds calibrated to %.2f FA/h on target-free calibration data only.\n\n', cfg.target_fa_per_h);
fprintf('%-5s %-5s | %-22s | %-22s | %-15s\n', 'CNR', 'z0', 'Pd A [95% CI]', 'Pd B [95% CI]', 'latency A / B');
fprintf('%s\n', repmat('-', 1, 80));
for iC = 1:numel(R.cnr)
    for iZ = 1:numel(R.z0)
        fprintf('%-5d %-5.1f | %.3f [%.3f, %.3f]    | %.3f [%.3f, %.3f]    | %5.2f / %5.2f s\n', ...
            R.cnr(iC), R.z0(iZ), R.pdA(iC,iZ), R.pdA_ci(iC,iZ,1), R.pdA_ci(iC,iZ,2), ...
            R.pdB(iC,iZ), R.pdB_ci(iC,iZ,1), R.pdB_ci(iC,iZ,2), R.latA(iC,iZ), R.latB(iC,iZ));
    end
end
fprintf('\n%-5s | %-28s | %-28s\n', 'CNR', 'FA/h A [95% CI] (events)', 'FA/h B [95% CI] (events)');
fprintf('%s\n', repmat('-', 1, 68));
for iC = 1:numel(R.cnr)
    fprintf('%-5d | %.2f [%.2f, %.2f] (%d)        | %.2f [%.2f, %.2f] (%d)\n', R.cnr(iC), ...
        R.faA(iC), R.faA_ci(iC,1), R.faA_ci(iC,2), R.nfaA(iC), ...
        R.faB(iC), R.faB_ci(iC,1), R.faB_ci(iC,2), R.nfaB(iC));
end
fprintf(['\nReading: where the two Pd confidence intervals overlap, this run does not show a\n' ...
         'difference between the detectors. Latency = first event onset - t0 (negative = before\n' ...
         'the target centre reaches the LOS). Runtime %.0f s.\n'], R.runtime_s);
fprintf('All numbers are synthetic; none is a field measurement.\n');
end

function plotComparison(cfg, R, bank, ex)
fig = figure('Name','RF Drone Fence - detector comparison (synthetic)','Color','w','Position',[50 50 1300 820]);
tl = tiledlayout(fig, 2, 3, 'TileSpacing','compact','Padding','compact');
title(tl, 'Detector A (amplitude + noise estimate) vs B (matched-filter bank) - SYNTHETIC data, not field measurements', ...
    'FontWeight','bold');
colA = [0.10 0.35 0.65]; colB = [0.80 0.40 0.05];
for iC = 1:numel(R.cnr)
    ax = nexttile(tl, iC); hold(ax,'on'); grid(ax,'on');
    dz = 0.04;
    pa = squeeze(R.pdA_ci(iC,:,:)); pb = squeeze(R.pdB_ci(iC,:,:));
    if isvector(pa), pa = reshape(pa, [], 2); pb = reshape(pb, [], 2); end
    errorbar(ax, R.z0 - dz, R.pdA(iC,:), R.pdA(iC,:) - pa(:,1)', pa(:,2)' - R.pdA(iC,:), 'o-', ...
        'Color', colA, 'LineWidth', 1.4, 'MarkerFaceColor', colA);
    errorbar(ax, R.z0 + dz, R.pdB(iC,:), R.pdB(iC,:) - pb(:,1)', pb(:,2)' - R.pdB(iC,:), 's-', ...
        'Color', colB, 'LineWidth', 1.4, 'MarkerFaceColor', colB);
    ylim(ax, [-0.02 1.05]); xlim(ax, [min(R.z0)-0.2 max(R.z0)+0.2]);
    xlabel(ax, 'vertical offset z_0 [m]'); ylabel(ax, 'Pd per crossing (95% CI)');
    title(ax, sprintf('CNR %d dB: FA/h A %.2f, B %.2f', R.cnr(iC), R.faA(iC), R.faB(iC)));
    legend(ax, {'A amplitude + noise est.','B matched-filter bank'}, 'Location','southwest');
end
% templates
ax = nexttile(tl, 4); hold(ax,'on'); grid(ax,'on');
cm = lines(numel(bank));
for k = 1:numel(bank)
    tt = ((0:bank(k).len-1) - (bank(k).len-1)/2)/cfg.fs;
    plot(ax, tt, bank(k).Q(:,1), '-', 'Color', cm(k,:), 'LineWidth', 1.1);
end
xlabel(ax, 'time from template centre [s]'); ylabel(ax, 'template (1st basis vector)');
title(ax, 'Detector B templates from fsr\_scattered\_field');
legend(ax, arrayfun(@(b) sprintf('v = %.1f m/s', b.v), bank, 'UniformOutput', false), 'Location','northeast');
% example scores
if ~isempty(ex)
    m = ex.t >= ex.t0 - 6 & ex.t <= ex.t0 + 6;
    ax = nexttile(tl, 5); hold(ax,'on'); grid(ax,'on');
    plot(ax, ex.t(m) - ex.t0, ex.sA(m), 'Color', colA);
    yline(ax, R.thrA(1), '--', 'Color', colA, 'LineWidth', 1.2);
    xline(ax, 0, 'k:');
    xlabel(ax, 'time from true crossing t_0 [s] (for display only)'); ylabel(ax, 'score A');
    title(ax, sprintf('Example crossing, CNR %d dB, z_0 = %.1f m, v = %.2f m/s', R.cnr(1), R.z0(1), ex.v));
    ax = nexttile(tl, 6); hold(ax,'on'); grid(ax,'on');
    plot(ax, ex.t(m) - ex.t0, ex.sB(m), 'Color', colB);
    yline(ax, R.thrB(1), '--', 'Color', colB, 'LineWidth', 1.2);
    xline(ax, 0, 'k:');
    xlabel(ax, 'time from true crossing t_0 [s] (for display only)'); ylabel(ax, 'score B');
    title(ax, 'Same recording, detector B (dashed: calibrated threshold)');
end
end

%% =====================================================================
%% Tests: no target leakage, causality, sanity
%% =====================================================================
function T = runTests()
fprintf('\nfsr_detector_comparison - leakage, causality, sanity and channel-integration tests (synthetic)\n');
T = struct('name',{},'pass',{},'detail',{});
cfg = defaultConfig(); cfg.plot = false; cfg.verbose = false; cfg.cnr_db = 20;

% L1 seed ranges are disjoint between calibration, empty test and target test
sc = []; se = []; st = [];
for iC = 1:numel(cfg.cnrList)
    sc = [sc, arrayfun(@(k) seedOf(cfg,'cal',iC,0,k), 1:cfg.calTrials)]; %#ok<AGROW>
    se = [se, arrayfun(@(k) seedOf(cfg,'empty',iC,0,k), 1:cfg.emptyTrials)]; %#ok<AGROW>
    for iZ = 1:numel(cfg.z0List)
        st = [st, arrayfun(@(k) seedOf(cfg,'target',iC,iZ,k), 1:cfg.targetTrials)]; %#ok<AGROW>
    end
end
ok = isempty(intersect(sc,se)) && isempty(intersect(sc,st)) && isempty(intersect(se,st)) && ...
     numel(unique(st)) == numel(st) && numel(unique(sc)) == numel(sc);
T = rec(T, 'L1 calibration / empty-test / target-test seeds are disjoint', ok, ...
    sprintf('%d cal, %d empty, %d target seeds', numel(sc), numel(se), numel(st)));

% L2 calibration recordings contain no target: t0 is NaN and the clean field is exactly 1
okL2 = true;
for k = 1:3
    tr = simulateTrial(cfg, 120, false, 0, seedOf(cfg,'cal',1,0,k));
    okL2 = okL2 && isnan(tr.t0) && max(abs(tr.clean - 1)) == 0;
end
T = rec(T, 'L2 calibration recordings are target-free (s = 0, no t0)', okL2, 'checked 3 recordings');

% L3 thresholds do not depend on anything about targets (speed range, z0, target count)
c1 = cfg; c1.calTrials = 4;
c2 = c1;  c2.testSpeedRange = [1 2]; c2.targetTrials = 7; c2.z0 = 1.7;
thr1 = calibrateBoth(c1); thr2 = calibrateBoth(c2);
T = rec(T, 'L3 thresholds identical when target settings change', isequal(thr1, thr2), ...
    sprintf('A %.4f vs %.4f, B %.4f vs %.4f', thr1(1), thr2(1), thr1(2), thr2(2)));

% L4 the template bank does not depend on test targets
b1 = buildBank(cfg); c3 = cfg; c3.testSpeedRange = [9 9.5]; c3.z0List = 3; b2 = buildBank(c3);
same = numel(b1) == numel(b2);
for k = 1:min(numel(b1),numel(b2)), same = same && isequal(b1(k).Q, b2(k).Q); end
T = rec(T, 'L4 matched-filter bank is fixed by cfg.bankSpeeds and geometry only', same, ...
    sprintf('%d templates, speeds %s m/s', numel(b1), mat2str([b1.v])));

% L5 causality: changing samples after n0 never changes scores at or before n0
tr = simulateTrial(cfg, 120, true, 0, 999);
n0 = round(80*cfg.fs);
r2 = tr.r; r2(n0+1:end) = r2(n0+1:end) + 0.3*exp(1i*0.7);   % large change in the future
bank = buildBank(cfg);
sA1 = scoreA(tr.r, cfg); sA2 = scoreA(r2, cfg);
dA = max(abs(sA1(1:n0) - sA2(1:n0)));
sB1 = scoreB(tr.r, bank, cfg); sB2 = scoreB(r2, bank, cfg);
dB = max(abs(sB1(1:n0) - sB2(1:n0)));
T = rec(T, 'L5 causality: future samples do not affect current scores (A and B)', dA < 1e-9 && dB < 1e-6*max(1,max(abs(sB1(1:n0)))), ...
    sprintf('max change before n0: A %.1e, B %.1e', dA, dB));

% L6 no hidden use of truth: on TARGET-FREE recordings, scoring against random fake t0
% must give a "Pd" close to the chance level (alarms do not follow t0)
c6 = cfg; c6.calTrials = 8;
[thr6, bank6] = calibrateBoth(c6);
hitsA = 0; hitsB = 0; nF = 40; rng(4242,'twister');
for k = 1:nF
    tr = simulateTrial(c6, 120, false, 0, 777000 + k);
    fake = c6.warmup + 10 + rand*(120 - c6.warmup - 25);
    hitsA = hitsA + scoreCrossing(eventOnsets(scoreA(tr.r,c6), thr6(1), c6), fake, c6);
    hitsB = hitsB + scoreCrossing(eventOnsets(scoreB(tr.r,bank6,c6), thr6(2), c6), fake, c6);
end
T = rec(T, 'L6 fake-t0 control on target-free data: "Pd" stays at chance (<= 0.10)', ...
    hitsA/nF <= 0.10 && hitsB/nF <= 0.10, sprintf('A %.3f, B %.3f over %d recordings', hitsA/nF, hitsB/nF, nF));

% S1 template sanity: the bank reproduces the physics signature at its own speed
v = cfg.bankSpeeds(3); lam = fsr_wavelength(cfg.fc);
tt = (-6:1/cfg.fs:6)'; s0 = fsr_scattered_field('paraxial', makeGeom(cfg), v*tt, 0.5, lam);
dAmp = abs(1 - s0) - 1;
B = [real(s0) imag(s0)];
resid = norm(dAmp - B*(B\dAmp))/norm(dAmp);
T = rec(T, 'S1 amplitude signature at z0 = 0.5 m lies in {Re s0, Im s0} (first order)', resid < 0.05, ...
    sprintf('relative residual %.3f', resid));

% S2 detector B fires on a clean strong crossing and stays quiet before it
c7 = cfg; c7.cnr_db = 40;
tr = simulateTrial(c7, 120, true, 0, 4321);
sB = scoreB(tr.r, bank, c7);
i0 = round(tr.t0*c7.fs);
pre = max(sB(round(40*c7.fs):i0 - round(5*c7.fs)));
peak = max(sB(i0 - round(1.5*c7.fs):i0 + round(3*c7.fs)));
T = rec(T, 'S2 detector B: peak near the crossing >> pre-crossing maximum', peak > 3*pre, ...
    sprintf('peak %.1f, pre-crossing max %.1f', peak, pre));

% A1 the 'api' handles run the same functions as the protocol
api = fsr_detector_comparison('api');
trA1 = simulateTrial(cfg, 120, true, 0, 1234);
bA1 = api.buildBank(cfg);
okA1 = isequal(api.simulateTrial(cfg, 120, true, 0, 1234), trA1) && ...
       isequal(api.scoreA(trA1.r, cfg), scoreA(trA1.r, cfg)) && ...
       isequal(api.scoreB(trA1.r, bA1, cfg), scoreB(trA1.r, bank, cfg));
T = rec(T, 'A1 api handles reproduce the protocol functions exactly (for the app)', okA1, ...
    'simulateTrial, scoreA, scoreB compared with isequal');

% ---------------- channel integration (v3) ----------------
% C1 default channel (multipath off) gives bit-identical recordings to the frozen v2 code
okC1 = true; nC1 = 0;
variants = {struct(), struct('drift_db',1.5), struct('interf_db',-20), ...
            struct('drift_db',1,'interf_db',-25,'cnr_db',15), struct('cnr_db',Inf,'cfo_hz',0)};
for iv = 1:numel(variants)
    c = cfg; fv = fieldnames(variants{iv});
    for i = 1:numel(fv), c.(fv{i}) = variants{iv}.(fv{i}); end
    cL = c; cL.channelImpl = 'legacy';
    for withTarget = [false true]
        for sd = [100001 300123 777]
            a = simulateTrial(cL, 120, withTarget, 0.5, sd);
            b = simulateTrial(c,  120, withTarget, 0.5, sd);
            okC1 = okC1 && isequal(a.r, b.r) && isequal(a.clean, b.clean) && ...
                   isequaln(a.t0, b.t0) && isequaln(a.v, b.v);
            nC1 = nC1 + 1;
        end
    end
end
T = rec(T, 'C1 multipath off: recordings bit-identical to the v2 channel', okC1, ...
    sprintf('%d recordings (5 channel settings x target/none x 3 seeds), isequal', nC1));

% C2 the whole protocol (calibration, thresholds, FA/h, Pd, CI, latency) is unchanged
cs = cfg; cs.cnrList = [20 15]; cs.z0List = [0 1];
cs.calTrials = 3; cs.calDur = 200; cs.emptyTrials = 3; cs.emptyDur = 200;
cs.targetTrials = 8; cs.targetDur = 120;
csL = cs; csL.channelImpl = 'legacy';
Rm = runComparison(cs); Rl = runComparison(csL);
flds = {'thrA','thrB','faA','faB','nfaA','nfaB','faA_ci','faB_ci','pdA','pdB','pdA_ci','pdB_ci','latA','latB'};
sameAll = true;
for i = 1:numel(flds), sameAll = sameAll && isequaln(Rm.(flds{i}), Rl.(flds{i})); end
T = rec(T, 'C2 multipath off: thresholds, FA/h, Pd, CIs, latency identical to v2', sameAll, ...
    sprintf('reduced grid 2 CNR x 2 z0, %d fields compared with isequaln; thr A %s', ...
    numel(flds), mat2str(Rm.thrA, 4)));

% C3 noise and interference are added exactly once
c3 = cfg; c3.cnr_db = 20;                                  % drift off, interference off
tr = simulateTrial(c3, 600, false, 0, 4242);
rot = exp(1i*2*pi*c3.cfo_hz*tr.t);
Pn = mean(abs(tr.r - tr.clean.*rot).^2); P0 = 10^(-c3.cnr_db/10); nS = numel(tr.r);
c3b = c3; c3b.cnr_db = Inf;                                % no noise at all
trb = simulateTrial(c3b, 120, true, 0, 4243);
eNo = max(abs(trb.r - trb.clean.*exp(1i*2*pi*c3b.cfo_hz*trb.t)));
c3c = c3b; c3c.interf_db = -10; c3c.interf_tau = 0.2;      % interference only
trc = simulateTrial(c3c, 600, false, 0, 4244);
Pc = mean(abs(trc.r./exp(1i*2*pi*c3c.cfo_hz*trc.t) - trc.clean).^2);
rho = exp(-1/(c3c.fs*c3c.interf_tau));
okC3 = abs(Pn/P0 - 1) < 6/sqrt(nS) && eNo < 1e-12 && ...
       abs(Pc/0.1 - 1) < 6*sqrt((1 + rho^2)/(1 - rho^2)/numel(trc.r));
T = rec(T, 'C3 noise and interference added once (power = configured, not doubled)', okC3, ...
    sprintf('noise P/P0 = %.4f (tol %.4f); no-noise residual %.1e; interf P = %.4f (0.1)', ...
    Pn/P0, 6/sqrt(nS), eNo, Pc));

% C4 multipath ON reaches every recording type, including target-free calibration
c4 = cfg; c4.cnr_db = Inf; c4.cfo_hz = 0;
c4.paths = struct('gain', -0.7, 'excess_m', 0.2237, 'modulation', 'static', ...
    'label', 'test value only, not a site model');
lam4 = fsr_wavelength(c4.fc);
cyc = 0.2237/lam4; pf = -0.7*exp(-1i*2*pi*(cyc - floor(cyc)));
trCal = simulateTrial(c4, 120, false, 0, seedOf(c4,'cal',1,0,1));
trTgt = simulateTrial(c4, 120, true, 0.5, 5151);
e4a = max(abs(trCal.r - (1 + pf)));
e4b = max(abs(trTgt.r - (trTgt.clean + pf)));
okC4 = e4a < 1e-12 && e4b < 1e-12;
T = rec(T, 'C4 multipath on: calibration and target recordings both get r = E + path', okC4, ...
    sprintf('calibration err %.1e, target err %.1e, |1+path| = %.3f', e4a, e4b, abs(1 + pf)));

% C5 the full protocol runs end to end with multipath on (smoke test, not a result)
c5 = cs; c5.paths = c4.paths;
R5 = runComparison(c5);
okC5 = all(isfinite([R5.thrA R5.thrB])) && all(R5.pdA(:) >= 0 & R5.pdA(:) <= 1) && ...
       all(R5.pdB(:) >= 0 & R5.pdB(:) <= 1) && ~isempty(strfind(R5.multipath, 'ON'));
T = rec(T, 'C5 protocol runs with multipath on (thresholds re-calibrated with it)', okC5, ...
    sprintf('thr A %s vs off %s; thr B %s vs off %s', mat2str(R5.thrA,4), mat2str(Rm.thrA,4), ...
    mat2str(R5.thrB,4), mat2str(Rm.thrB,4)));

fprintf('\n================ TEST SUMMARY ================\n');
for i = 1:numel(T)
    if T(i).pass, tag = 'PASS'; else, tag = 'FAIL'; end
    fprintf('[%s] %s  (%s)\n', tag, T(i).name, T(i).detail);
end
nF = sum(~[T.pass]);
fprintf('\n%d tests, %d failed. Tests check the protocol, not real-world performance.\n', numel(T), nF);
if nF > 0, error('fsr:tests','%d comparison tests failed.', nF); end
end

function [thr, bank] = calibrateBoth(cfg)
bank = buildBank(cfg);
sA = cell(cfg.calTrials,1); sB = sA;
for k = 1:cfg.calTrials
    tr = simulateTrial(cfg, cfg.calDur, false, 0, seedOf(cfg,'cal',1,0,k));
    sA{k} = scoreA(tr.r, cfg); sB{k} = scoreB(tr.r, bank, cfg);
end
thr = [calibrateThreshold(sA, cfg), calibrateThreshold(sB, cfg)];
end

function T = rec(T, name, pass, detail)
T(end+1) = struct('name',name,'pass',logical(pass),'detail',detail);
if pass, tag = 'PASS'; else, tag = 'FAIL'; end
fprintf('  [%s] %s  (%s)\n', tag, name, detail);
end