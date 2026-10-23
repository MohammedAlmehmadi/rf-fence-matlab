function T = fsr_channel_model_tests()
%FSR_CHANNEL_MODEL_TESTS Deterministic and statistical tests of fsr_channel_model.m
%
%   >> fsr_channel_model_tests
%   >> T = fsr_channel_model_tests;      % also returns the result table
%
% Needs in the same folder: fsr_channel_model.m, fsr_wavelength.m, and (for test 9)
% fsr_scattered_field.m, fsr_babinet_paraxial.m, fsr_fresnel_cs.m. Base MATLAB only.
%
% What these tests establish: the code implements equations (1)-(2) in the header
% of fsr_channel_model.m, is causal, is reproducible, rejects bad input, and
% reproduces the legacy channel of fsr_detector_comparison.m bit for bit.
% What they do NOT establish: that any channel parameter (reflection gain, drift,
% interference, CNR) matches a real site. Nothing here is a measurement.
%
% Statistical tests use fixed seeds AND tolerances of 6 standard deviations of the
% estimator, so a pass does not depend on a lucky draw and a fail is meaningful.

fprintf('\nfsr_channel_model - tests (synthetic, mathematical consistency only)\n');
need = {'fsr_channel_model','fsr_wavelength','fsr_scattered_field','fsr_babinet_paraxial','fsr_fresnel_cs'};
for i = 1:numel(need)
    if exist(need{i}, 'file') ~= 2
        error('fsr_channel_tests:deps', 'Missing %s.m on the MATLAB path.', need{i});
    end
end
T = struct('name', {}, 'pass', {}, 'detail', {});
c0 = 299792458.0;
fs = 250; fc = 2.4e9; lam = fsr_wavelength(fc);

%% 1 Identity
n = 3000; t = (0:n-1)'/fs;
rng(1, 'twister');
E = 1 - 0.1*(randn(n,1) + 1i*randn(n,1));
r = fsr_channel_model(E, t, struct('fs', fs));
ch = fsr_channel_model('defaults'); ch.fs = fs;
r2 = fsr_channel_model(E, t, ch);
T = rec(T, '1 Identity: all effects off, no noise -> r == E exactly', isequal(r, E) && isequal(r2, E), ...
    sprintf('max|r-E| = %g', max(abs(r - E))));

%% 2 Single reflected path, static: r = E + G*exp(-1i*k*dL) exactly, two-ray formula for |r|
G = 0.6*exp(1i*2.2); dL = 0.2237;               % arbitrary test values, not a site model
cfg = struct('fs', fs, 'fc', fc);
cfg.paths = struct('gain', G, 'excess_m', dL, 'modulation', 'static');
E1 = ones(n,1);
[r, info] = fsr_channel_model(E1, t, cfg);
expect = 1 + G*exp(-1i*2*pi*dL/lam);
psi = angle(G) - 2*pi*dL/lam;
twoRay = sqrt(1 + abs(G)^2 + 2*abs(G)*cos(psi));
e2 = max(abs(r - expect)); e2b = max(abs(abs(r) - twoRay));
T = rec(T, '2 Single static path = analytic 1 + G e^{-jk dL}; |r| = two-ray formula', e2 < 1e-12 && e2b < 1e-12, ...
    sprintf('field err %.1e, |r| err %.1e; delay %.2e samples (unresolvable at %g Hz)', e2, e2b, info.delay_samples, fs));

%% 3 Delay: artificial large delays to exercise the interpolator
% Gaussian pulse, sigma = 20 samples; linear interpolation error bound h^2/8*max|x''| = 1/(8*sigma^2)
sig = 20; nn = 1000; tt = (0:nn-1)'/fs; k0 = 400;
x = exp(-((0:nn-1)' - k0).^2/(2*sig^2));
ok3 = true; worst = 0; worstC = 0;
for D = [1 2.25 3.4 7.9]
    cfg = struct('fs', fs, 'fc', fc);
    cfg.paths = struct('gain', 1, 'excess_m', D*c0/fs, 'modulation', x);
    [~, info] = fsr_channel_model(zeros(nn,1), tt, cfg);
    y = (info.Emp) / info.path_factor;           % remove the complex path factor
    exact = exp(-((0:nn-1)' - k0 - D).^2/(2*sig^2));
    errW = max(abs(y - exact));
    nIdx = (0:nn-1)';
    Dhat = sum(nIdx.*real(y))/sum(real(y)) - sum(nIdx.*x)/sum(x);   % centroid shift
    worst = max(worst, errW / (1/(8*sig^2))); worstC = max(worstC, abs(Dhat - D));
    ok3 = ok3 && errW <= 1/(8*sig^2) + 1e-12 && abs(Dhat - D) < 1e-6 && abs(info.delay_samples - D) < 1e-9;
end
T = rec(T, '3 Delay: shape within h^2/8 max|x''''| bound, centroid recovers D (1..7.9 samples)', ok3, ...
    sprintf('worst shape err = %.2f of bound; worst centroid err %.1e samples', worst, worstC));

%% 4 Phase: dL = lambda/4 -> -pi/2, lambda/2 -> pi, 1.125 lambda -> -pi/4
ok4 = true; det4 = '';
for q = [0.25 0.5 1.125]
    cfg = struct('fs', fs, 'fc', fc);
    cfg.paths = struct('gain', 0.5, 'excess_m', q*lam, 'modulation', 'static');
    [~, info] = fsr_channel_model(ones(10,1), (0:9)'/fs, cfg);
    want = 0.5*exp(-1i*2*pi*q);
    ok4 = ok4 && abs(info.path_factor - want) < 1e-12;
    det4 = [det4 sprintf('%.3g lam: arg %.4f rad; ', q, angle(info.path_factor))]; %#ok<AGROW>
end
cfg = struct('fs', fs, 'fc', fc);                 % phase of the gain itself also carried
cfg.paths = struct('gain', 0.3*exp(1i*0.7), 'excess_m', lam, 'modulation', 'direct');
Eh = exp(1i*0.1*(1:10)');
r = fsr_channel_model(Eh, (0:9)'/fs, cfg);
D = lam/c0*fs;  % tiny
ok4 = ok4 && max(abs(r - (Eh + 0.3*exp(1i*0.7)*((1-D)*Eh + D*[Eh(1); Eh(1:end-1)])))) < 1e-12;
T = rec(T, '4 Phase: -2*pi*dL/lambda and arg(gain) reproduced exactly', ok4, det4);

%% 5 Multipath: hand-computed 5-sample example in artificial units
% fs = 1 Hz and fc = 1 Hz, so 1 sample = 1 wavelength = c0 metres. Path 1 ('direct'):
% gain 0.5, dL = 1.5 c0 -> D = 1.5 samples, phase -2*pi*1.5 -> factor -0.5.
% Path 2 ('static'): gain 0.2, dL = 0.75 c0 -> phase -2*pi*0.75 -> factor +0.2i.
% Delayed path 1: y = 0.5*E(n-1) + 0.5*E(n-2), held at E(1) before the start:
%   y = [1; 1; 1.5; 1+1.5i; 2+1.5i];  r = E - 0.5*y + 0.2i
Eh = [1; 2; 3i; 4; 5];
cfg = struct('fs', 1, 'fc', 1);
cfg.paths = struct('gain', {0.5, 0.2}, 'excess_m', {1.5*c0, 0.75*c0}, 'modulation', {'direct', 'static'});
[r, info] = fsr_channel_model(Eh, (0:4)', cfg);
hand = [0.5+0.2i; 1.5+0.2i; -0.75+3.2i; 3.5-0.55i; 4-0.55i];
e5 = max(abs(r - hand)); e5f = max(abs(info.path_factor - [-0.5, 0.2i]));
T = rec(T, '5 Two-path example equals a hand calculation', e5 < 1e-12 && e5f < 1e-12, ...
    sprintf('field err %.1e, path-factor err %.1e', e5, e5f));

%% 6 Noise power, drift and interference statistics (fixed seeds, 6-sigma tolerances)
N = 400000; tN = (0:N-1)'/fs; ok6 = true; det6 = '';
for cnr = [30 20 10]
    [~, info] = fsr_channel_model(ones(N,1), tN, struct('fs', fs, 'cnr_db', cnr, 'seed', 61));
    Pm = mean(abs(info.noise).^2); P0 = 10^(-cnr/10);
    tol = 6/sqrt(N);
    balance = abs(var(real(info.noise)) - var(imag(info.noise)))/P0;
    ok6 = ok6 && abs(Pm/P0 - 1) < tol && balance < 6/sqrt(N) && abs(mean(info.noise))^2 < 36*P0/N;
    det6 = [det6 sprintf('CNR %d: P/P0-1 = %+.4f (tol %.4f); ', cnr, Pm/P0 - 1, tol)]; %#ok<AGROW>
end
tau = 0.2; rho = exp(-1/(fs*tau));
[~, info] = fsr_channel_model(ones(N,1), tN, struct('fs', fs, 'drift_db', 2, 'drift_tau', tau, 'seed', 62));
u = 20*log10(info.g);
relTolStd = 6*0.5*sqrt(2/N*(1 + rho^2)/(1 - rho^2));    % std of a sample std, AR(1)
r1 = sum(u(1:end-1).*u(2:end))/sum(u.^2);
okD = abs(std(u)/2 - 1) < relTolStd && abs(r1 - rho) < 6*sqrt((1 - rho^2)/N);
[~, info] = fsr_channel_model(ones(N,1), tN, struct('fs', fs, 'interf_db', -10, 'interf_tau', tau, 'seed', 63));
Pc = mean(abs(info.c).^2);
okI = abs(Pc/0.1 - 1) < 6*sqrt((1 + rho^2)/(1 - rho^2)/N);
ok6 = ok6 && okD && okI;
det6 = [det6 sprintf('drift std %.3f dB (2), lag-1 %.4f (%.4f); interf P %.4f (0.1)', std(u), r1, rho, Pc)];
T = rec(T, '6 Noise power -> 10^(-CNR/10); drift std/correlation; interference power', ok6, det6);

%% 7 Reproducibility
cfgR = struct('fs', fs, 'fc', fc, 'cnr_db', 15, 'drift_db', 1, 'interf_db', -20, 'cfo_hz', 3, 'seed', 7);
cfgR.paths = struct('gain', 0.4, 'excess_m', 0.3, 'modulation', 'direct');
Er = 1 - 0.05*exp(-((0:n-1)' - 1500).^2/2000);
ra = fsr_channel_model(Er, t, cfgR);
rng(12345, 'twister'); before = rng;
rb = fsr_channel_model(Er, t, cfgR);
after = rng;
cfgR.seed = 8; rc = fsr_channel_model(Er, t, cfgR);
cfgR.seed = [];
rng(99, 'twister'); rd = fsr_channel_model(Er, t, cfgR);
rng(99, 'twister'); re = fsr_channel_model(Er, t, cfgR);
ok7 = isequal(ra, rb) && ~isequal(ra, rc) && isequal(before.State, after.State) && isequal(rd, re);
T = rec(T, '7 Same seed -> identical; other seed -> different; caller RNG state restored', ok7, ...
    'seed mode and global-stream mode both checked');

%% 8 Causality: changing E after n0 never changes r(1:n0)
rng(3, 'twister');
Ec = 1 - 0.1*(randn(n,1) + 1i*randn(n,1));
cfgC = struct('fs', fs, 'fc', fc, 'cnr_db', 20, 'drift_db', 1, 'interf_db', -15, 'cfo_hz', 3, 'seed', 81);
cfgC.paths = struct('gain', {0.5, 0.3, 0.2}, 'excess_m', {2.7*c0/fs, 5*c0/fs, 0.2}, ...
    'modulation', {'direct', 'direct', 'static'});
n0 = 1700;
E2 = Ec; E2(n0+1:end) = E2(n0+1:end) + 0.5*exp(1i*0.4);
ra = fsr_channel_model(Ec, t, cfgC);
rb = fsr_channel_model(E2, t, cfgC);
T = rec(T, '8 Causality: future input changes leave past output identical', ...
    isequal(ra(1:n0), rb(1:n0)) && ~isequal(ra(n0+1:end), rb(n0+1:end)), ...
    sprintf('max change before n0 = %g; after = %.3g', max(abs(ra(1:n0) - rb(1:n0))), max(abs(ra(n0+1:end) - rb(n0+1:end)))));

%% 9 Backward compatibility with simulateTrial of fsr_detector_comparison.m
base = struct('fc', 2.4e9, 'L', 20, 'hAnt', 1.5, 'xc', 10, 'w', 0.30, 'h', 0.20, 'opacity', 1, ...
    'fs', 250, 'testSpeedRange', [3 8], 'warmup', 35, ...
    'drift_db', 0, 'drift_tau', 60, 'interf_db', -Inf, 'interf_tau', 5, 'cfo_hz', 3, 'cnr_db', 20);
variants = {struct(), struct('drift_db', 1.5), struct('interf_db', -20), ...
            struct('drift_db', 1, 'interf_db', -25, 'cnr_db', 15), struct('cnr_db', Inf, 'cfo_hz', 0)};
ok9 = true; nCmp = 0;
for v = 1:numel(variants)
    c = base; fv = fieldnames(variants{v});
    for i = 1:numel(fv), c.(fv{i}) = variants{v}.(fv{i}); end
    for withTarget = [false true]
        for seed = [100001 300123 777]
            a = legacySimulateTrial(c, 120, withTarget, 0.5, seed);
            b = newSimulateTrial(c, 120, withTarget, 0.5, seed);
            ok9 = ok9 && isequal(a.r, b.r) && isequal(a.clean, b.clean) && isequaln(a.t0, b.t0);
            nCmp = nCmp + 1;
        end
    end
end
T = rec(T, '9 Compatibility: paths off -> bit-identical to legacy simulateTrial', ok9, ...
    sprintf('%d recordings (5 channel settings x target/no target x 3 seeds), isequal', nCmp));

%% 10 Input validation
okC = struct('fs', fs, 'fc', fc);
bad = {
  'negative excess',    @() fsr_channel_model(ones(10,1), (0:9)'/fs, withPath(okC, struct('gain',0.5,'excess_m',-1,'modulation','static'))), 'fsr_channel:path'
  'zero excess',        @() fsr_channel_model(ones(10,1), (0:9)'/fs, withPath(okC, struct('gain',0.5,'excess_m',0,'modulation','direct'))), 'fsr_channel:path'
  'delay > record',     @() fsr_channel_model(ones(10,1), (0:9)'/fs, withPath(okC, struct('gain',0.5,'excess_m',20*c0/fs,'modulation','direct'))), 'fsr_channel:path'
  'gain not scalar',    @() fsr_channel_model(ones(10,1), (0:9)'/fs, withPath(okC, struct('gain',[1 2],'excess_m',1,'modulation','static'))), 'fsr_channel:path'
  'gain NaN',           @() fsr_channel_model(ones(10,1), (0:9)'/fs, withPath(okC, struct('gain',NaN,'excess_m',1,'modulation','static'))), 'fsr_channel:path'
  'typo in path field', @() fsr_channel_model(ones(10,1), (0:9)'/fs, withPath(okC, struct('gian',0.5,'excess_m',1,'modulation','static'))), 'fsr_channel:path'
  'bad modulation',     @() fsr_channel_model(ones(10,1), (0:9)'/fs, withPath(okC, struct('gain',0.5,'excess_m',1,'modulation','reflect'))), 'fsr_channel:path'
  'modulation length',  @() fsr_channel_model(ones(10,1), (0:9)'/fs, withPath(okC, struct('gain',0.5,'excess_m',1,'modulation',ones(9,1)))), 'fsr_channel:path'
  'paths without fc',   @() fsr_channel_model(ones(10,1), (0:9)'/fs, withPath(struct('fs',fs), struct('gain',0.5,'excess_m',1,'modulation','static'))), 'fsr_channel:path'
  'typo in cfg field',  @() fsr_channel_model(ones(10,1), (0:9)'/fs, struct('fs',fs,'cnr_dB',20)), 'fsr_channel:config'
  'missing fs',         @() fsr_channel_model(ones(10,1), (0:9)'/fs, struct('cnr_db',20)), 'fsr_channel:config'
  'struct-array cfg',   @() fsr_channel_model(ones(10,1), (0:9)'/fs, struct('fs',{fs,fs})), 'fsr_channel:config'
  'negative drift',     @() fsr_channel_model(ones(10,1), (0:9)'/fs, struct('fs',fs,'drift_db',-1)), 'fsr_channel:config'
  'NaN cnr',            @() fsr_channel_model(ones(10,1), (0:9)'/fs, struct('fs',fs,'cnr_db',NaN)), 'fsr_channel:config'
  'bad seed',           @() fsr_channel_model(ones(10,1), (0:9)'/fs, struct('fs',fs,'seed',1.5)), 'fsr_channel:config'
  'row E',              @() fsr_channel_model(ones(1,10), (0:9)'/fs, struct('fs',fs)), 'fsr_channel:signal'
  'E/t length',         @() fsr_channel_model(ones(10,1), (0:8)'/fs, struct('fs',fs)), 'fsr_channel:signal'
  'NaN in E',           @() fsr_channel_model([1;NaN;1], (0:2)'/fs, struct('fs',fs)), 'fsr_channel:signal'
  't wrong rate',       @() fsr_channel_model(ones(10,1), (0:9)'/100, struct('fs',fs)), 'fsr_channel:signal'
  };
ok10 = true; failed = {};
for i = 1:size(bad, 1)
    try
        bad{i,2}();
        ok10 = false; failed{end+1} = [bad{i,1} ' (no error)']; %#ok<AGROW>
    catch ME
        if ~strcmp(ME.identifier, bad{i,3})
            ok10 = false; failed{end+1} = sprintf('%s (got %s)', bad{i,1}, ME.identifier); %#ok<AGROW>
        end
    end
end
if ok10, d10 = sprintf('%d invalid inputs, each rejected with the expected error id', size(bad,1));
else, d10 = strjoin(failed, '; '); end
T = rec(T, '10 Input validation: bad delays, sizes, gains, fields -> clear errors', ok10, d10);

%% Summary
fprintf('\n================ TEST SUMMARY ================\n');
for i = 1:numel(T)
    if T(i).pass, tag = 'PASS'; else, tag = 'FAIL'; end
    fprintf('[%s] %s\n       %s\n', tag, T(i).name, T(i).detail);
end
nF = sum(~[T.pass]);
fprintf('\n%d tests, %d failed. These check the implementation of the stated equations,\n', numel(T), nF);
fprintf('not the realism of any parameter. All signals are synthetic.\n');
if nF > 0, error('fsr_channel_tests:failed', '%d channel-model tests failed.', nF); end
end

%% =====================================================================
function c = withPath(c, p)
c.paths = p;
end

function T = rec(T, name, pass, detail)
T(end+1) = struct('name', name, 'pass', logical(pass), 'detail', detail);
if pass, tag = 'PASS'; else, tag = 'FAIL'; end
fprintf('  [%s] %s\n', tag, name);
end

function geom = makeGeom(cfg)
geom = struct('src',[0 0 cfg.hAnt],'rcv',[cfg.L 0 cfg.hAnt],'xc',cfg.xc, ...
              'w',cfg.w,'h',cfg.h,'opacity',cfg.opacity);
end

function tr = legacySimulateTrial(cfg, dur, withTarget, z0, seed)
% VERBATIM copy of simulateTrial + ou from fsr_detector_comparison.m (reference only)
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

function tr = newSimulateTrial(cfg, dur, withTarget, z0, seed)
% Proposed form of simulateTrial: target part unchanged, channel delegated.
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
    'interf_db', cfg.interf_db, 'interf_tau', cfg.interf_tau, 'cfo_hz', cfg.cfo_hz, 'cnr_db', cfg.cnr_db);
tr.clean = E;
tr.r = fsr_channel_model(E, t, ch);     % seed = [] -> continues the stream seeded above
end
