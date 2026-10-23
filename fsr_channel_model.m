function [r, info] = fsr_channel_model(E, t, cfg)
%FSR_CHANNEL_MODEL Apply channel effects to a clean complex-baseband FSR field.
%
% USAGE
%   r          = fsr_channel_model(E, t, cfg)
%   [r, info]  = fsr_channel_model(E, t, cfg)
%   ch         = fsr_channel_model('defaults')     % default (identity) channel struct
%
% INPUTS
%   E    n-by-1 complex column. Clean, normalised DIRECT-path field E/E0 = 1 - s(t),
%        normally from fsr_scattered_field (E = 1 means free space, no target).
%   t    n-by-1 real column, uniformly sampled: t(k+1) - t(k) = 1/cfg.fs  [s].
%   cfg  scalar struct with ONLY the channel fields below (unknown fields are an
%        error, so a typo cannot silently disable an effect). Missing fields take
%        the IDENTITY defaults shown, i.e. an empty cfg (plus fs) leaves E unchanged.
%
%   field        default  unit   meaning
%   fs           (req.)   Hz     complex sample rate of E and t
%   fc           []       Hz     carrier; required only when paths is non-empty
%   paths        empty    -      reflected paths (struct array, see MULTIPATH)
%   drift_db     0        dB     std of slow log-normal gain drift (0 = off)
%   drift_tau    60       s      drift correlation time
%   interf_db    -Inf     dB     additive slow complex interference power rel. |E0|^2
%   interf_tau   5        s      interference correlation time
%   cfo_hz       0        Hz     residual carrier frequency offset
%   cnr_db       Inf      dB     carrier-to-noise ratio (Inf = no noise)
%   seed         []       -      [] = draw from the CURRENT global stream (legacy
%                                behaviour: the caller seeds); integer = use
%                                rng(seed,'twister') inside and restore the
%                                caller's global RNG state on exit.
%
%   NOTE: the identity defaults are NOT the demo defaults (cnr 20 dB, cfo 3 Hz in
%   fsr_detector_comparison.m). The caller must pass its own values.
%   NOTE: assign a struct array with dot syntax (cfg.paths = P). struct('paths',P)
%   with a struct array P would create a struct ARRAY cfg, which is rejected.
%
% SIGNAL MODEL (complex baseband, time convention exp(+1i*2*pi*fc*t), same as
% fsr_kirchhoff_numeric: a path that is dL metres longer is multiplied by
% exp(-1i*2*pi*dL/lambda))
%
%   Emp(t) = E(t) + sum_k  gain_k * exp(-1i*2*pi*dL_k/lambda) * x_k(t - tau_k)    (1)
%   r(t)   = g(t) .* ( Emp(t) + c(t) ) .* exp(1i*2*pi*cfo_hz*t) + n(t)            (2)
%
%   Equation (2) with paths = empty is, operation by operation and in the same
%   random-draw order, the channel of simulateTrial in fsr_detector_comparison.m
%   (and, per that file's header, fsr_detection_demo.m). Order of draws:
%   interference (real OU, imag OU), drift OU, noise (real, imag).
%
% MULTIPATH  (phenomenological tapped-delay model, engineering approximation)
%   The direct path is E itself (gain 1, zero excess length). It is NEVER listed in
%   paths, so it cannot be double-counted. Each element of cfg.paths is ONE
%   reflected path with fields
%     gain       complex scalar  amplitude/phase of the path relative to the direct
%                path, EXCLUDING the propagation phase (i.e. reflection coefficient
%                times spreading ratio L_direct/L_path). Dimensionless.
%     excess_m   real scalar > 0 path length minus direct length dL_k [m]
%                (a reflected path is always strictly longer than the direct one).
%     modulation 'static'  path is not crossed by the target: x_k(t) = 1
%                'direct'  path carries the same target signature: x_k = E
%                n-by-1    user-supplied clean field of that path (e.g. a later
%                          geometry-based model). Normalised like E.
%     label      optional char/string, for documentation only.
%   Delay  tau_k = dL_k / c0,  c0 = 299792458 m/s.  In samples D_k = tau_k*fs.
%   Carrier phase  phi_k = -2*pi*dL_k/lambda (evaluated modulo one wavelength).
%   Envelope delay x_k(t - tau_k): causal linear interpolation between samples,
%     y(n) = (1-f)*x(n-m) + f*x(n-m-1),  m = floor(D), f = D - m.
%   Exact for integer D; for a smooth x the error is <= h^2/8*max|x''| (h = 1
%   sample). Before the first sample the input is held at x(1) (the recording is
%   assumed to start in the same state it was in just before).
%   RESOLUTION LIMIT: at fs = 250 Hz one sample is c0/fs ~ 1199 km of path. A
%   physical reflection (cm .. km) therefore has D << 1 and is NOT resolvable in
%   time; its only effects are the complex factor gain*exp(-1i*k*dL) and the choice
%   of x_k. Larger D values are accepted for testing the interpolator, not because
%   they are physical at this sample rate. info.delay_samples reports D_k.
%   Not modelled here: a geometry-based reflection model (surfaces, Fresnel
%   coefficients, per-path FSR signature), per-path Doppler, diffuse scattering.
%
% SLOW GAIN DRIFT (engineering approximation of slow fading / hardware drift)
%   g(t) = 10^( drift_db * u(t) / 20 ),  u = discrete Ornstein-Uhlenbeck (AR(1))
%   with unit variance and time constant drift_tau:
%   u(n) = rho*u(n-1) + sqrt(1-rho^2)*w(n),  rho = exp(-1/(fs*drift_tau)),
%   u(0) ~ N(0,1) (stationary start). So 20*log10(g) ~ N(0, drift_db^2)
%   (log-normal gain), lag-1 correlation rho.
%
% INTERFERENCE (additive, slowly varying; engineering approximation)
%   c(t) = 10^(interf_db/20) * (u1(t) + 1i*u2(t))/sqrt(2), u1, u2 independent
%   unit OU processes with time constant interf_tau, so E|c|^2 = 10^(interf_db/10)
%   relative to the free-space carrier power |E0|^2 = 1. It passes through g and the
%   CFO like the carrier (legacy behaviour).
%
% NOISE (convention unchanged)
%   n(t) complex white Gaussian, E|n|^2 = 10^(-cnr_db/10), real and imaginary parts
%   each of variance 10^(-cnr_db/10)/2. cnr_db is the per-sample ratio of the
%   FREE-SPACE direct carrier power (|E0|^2 = 1) to noise power in the complex
%   sample bandwidth fs. It is NOT referenced to the multipath-modified carrier:
%   noise is a receiver property, so destructive multipath lowers the effective
%   CNR. info.effective_cnr_db = cnr_db + 10*log10(mean|Emp|^2) reports that.
%
% OUTPUT
%   r    n-by-1 complex received baseband samples.
%   info struct: Emp, c, g, rot, noise, delay_samples, path_factor (complex factor
%        gain*exp(-1i*k*dL) per path), carrier_power_db, effective_cnr_db, fs.
%
% REFERENCES (for the model forms, not for any parameter value)
%   T. S. Rappaport, Wireless Communications: Principles and Practice, 2nd ed.,
%     Prentice Hall, 2002 - multipath impulse response; two-ray ground reflection;
%     log-normal shadowing.
%   T. I. Laakso, V. Valimaki, M. Karjalainen, U. K. Laine, "Splitting the unit
%     delay", IEEE Signal Processing Magazine 13(1):30-60, 1996 - fractional delay
%     by (Lagrange/linear) interpolation.
%   D. T. Gillespie, "Exact numerical simulation of the Ornstein-Uhlenbeck process
%     and its integral", Phys. Rev. E 54:2084, 1996 - exact AR(1) update.
%   All parameter values are user assumptions; none is a measured value.
%
% DEPENDENCIES: fsr_wavelength.m (only when paths is non-empty). Base MATLAB only.

if nargin == 1 && (ischar(E) || isstring(E))
    if strcmpi(char(E), 'defaults'), r = defaultChannel(); info = []; return; end
    error('fsr_channel:usage', 'Single text argument must be ''defaults''.');
end
if nargin < 3
    error('fsr_channel:usage', 'Usage: r = fsr_channel_model(E, t, cfg). cfg.fs is required.');
end
ch = checkConfig(cfg);
[E, t] = checkSignal(E, t, ch.fs);
n  = numel(E);
fs = ch.fs;
P  = checkPaths(ch.paths, n, ch.fc, fs);

% ---- random-number stream
if ~isempty(ch.seed)
    prevRng = rng;
    restoreRng = onCleanup(@() rng(prevRng)); %#ok<NASGU>
    rng(ch.seed, 'twister');
end

% ---- (1) multipath: deterministic, makes NO random draws
Emp = E;
delaySamples = zeros(1, numel(P));
pathFactor = complex(zeros(1, numel(P)));
for k = 1:numel(P)
    lambda = fsr_wavelength(ch.fc);
    cycles = P(k).excess_m / lambda;                     % excess length in wavelengths
    pathFactor(k) = P(k).gain * exp(-1i*2*pi*(cycles - floor(cycles)));
    D = P(k).excess_m / 299792458.0 * fs;                % delay in samples
    delaySamples(k) = D;
    switch P(k).kind
        case 'static', x = ones(n, 1);
        case 'direct', x = fracDelay(E, D);
        case 'custom', x = fracDelay(P(k).x, D);
    end
    Emp = Emp + pathFactor(k) * x;
end

% ---- (2) legacy channel, same expressions and draw order as simulateTrial
if isfinite(ch.interf_db)
    c = (ouProcess(n, fs, ch.interf_tau) + 1i*ouProcess(n, fs, ch.interf_tau))/sqrt(2)*10^(ch.interf_db/20);
else
    c = zeros(n, 1);
end
if ch.drift_db > 0, g = 10.^(ch.drift_db*ouProcess(n, fs, ch.drift_tau)/20); else, g = ones(n, 1); end
rot = exp(1i*2*pi*ch.cfo_hz*t);
if isfinite(ch.cnr_db)
    sig = sqrt(10^(-ch.cnr_db/10)/2); noise = sig*(randn(n,1) + 1i*randn(n,1));
else
    noise = zeros(n, 1);
end
r = g.*(Emp + c).*rot + noise;

if nargout > 1
    info = struct();
    info.Emp = Emp; info.c = c; info.g = g; info.rot = rot; info.noise = noise;
    info.delay_samples = delaySamples;
    info.path_factor = pathFactor;
    info.carrier_power_db = 10*log10(mean(abs(Emp).^2));
    info.effective_cnr_db = ch.cnr_db + info.carrier_power_db;
    info.fs = fs;
end
end

%% =====================================================================
function ch = defaultChannel()
ch = struct('fs', [], 'fc', [], 'paths', struct([]), ...
    'drift_db', 0, 'drift_tau', 60, 'interf_db', -Inf, 'interf_tau', 5, ...
    'cfo_hz', 0, 'cnr_db', Inf, 'seed', []);
end

function ch = checkConfig(cfg)
if ~isstruct(cfg) || ~isscalar(cfg)
    error('fsr_channel:config', ['cfg must be a SCALAR struct (if you used struct(''paths'',P) ' ...
        'with a struct array P, assign cfg.paths = P instead).']);
end
ch = defaultChannel();
f = fieldnames(cfg);
for i = 1:numel(f)
    if ~isfield(ch, f{i})
        error('fsr_channel:config', 'Unknown channel field ''%s''. Allowed fields: %s.', ...
            f{i}, strjoin(fieldnames(ch)', ', '));
    end
    ch.(f{i}) = cfg.(f{i});
end
if isempty(ch.fs), error('fsr_channel:config', 'cfg.fs (sample rate, Hz) is required.'); end
needScalar(ch.fs, 'fs', isfinite(ch.fs) && ch.fs > 0, 'a finite positive number');
if ~isempty(ch.fc)
    needScalar(ch.fc, 'fc', isfinite(ch.fc) && ch.fc > 0, 'a finite positive number');
end
needScalar(ch.drift_db,   'drift_db',   isfinite(ch.drift_db) && ch.drift_db >= 0, 'finite and >= 0');
needScalar(ch.drift_tau,  'drift_tau',  isfinite(ch.drift_tau) && ch.drift_tau > 0, 'finite and > 0');
needScalar(ch.interf_db,  'interf_db',  ~isnan(ch.interf_db) && ch.interf_db < Inf, '-Inf (off) or finite');
needScalar(ch.interf_tau, 'interf_tau', isfinite(ch.interf_tau) && ch.interf_tau > 0, 'finite and > 0');
needScalar(ch.cfo_hz,     'cfo_hz',     isfinite(ch.cfo_hz), 'finite');
needScalar(ch.cnr_db,     'cnr_db',     ~isnan(ch.cnr_db) && ch.cnr_db > -Inf, 'finite or Inf (no noise)');
if ~isempty(ch.seed)
    needScalar(ch.seed, 'seed', isfinite(ch.seed) && ch.seed >= 0 && ch.seed == round(ch.seed) && ...
        ch.seed < 2^32, 'an integer in [0, 2^32)');
end
end

function needScalar(v, name, ok, what)
if ~(isnumeric(v) && isreal(v) && isscalar(v)) || ~ok
    error('fsr_channel:config', 'cfg.%s must be a real numeric scalar, %s.', name, what);
end
end

function [E, t] = checkSignal(E, t, fs)
if ~isnumeric(E) || ~iscolumn(E) || isempty(E)
    error('fsr_channel:signal', 'E must be a non-empty numeric n-by-1 column (got %s).', mat2str(size(E)));
end
if ~all(isfinite(E))
    error('fsr_channel:signal', 'E contains NaN or Inf.');
end
if ~isnumeric(t) || ~isreal(t) || ~iscolumn(t) || numel(t) ~= numel(E)
    error('fsr_channel:signal', 't must be a real n-by-1 column with the same length as E (E %s, t %s).', ...
        mat2str(size(E)), mat2str(size(t)));
end
if ~all(isfinite(t))
    error('fsr_channel:signal', 't contains NaN or Inf.');
end
if numel(t) > 1 && max(abs(diff(t)*fs - 1)) > 1e-6
    error('fsr_channel:signal', 't is not uniformly sampled at cfg.fs = %g Hz (max spacing error %.3g samples).', ...
        fs, max(abs(diff(t)*fs - 1)));
end
E = double(E); t = double(t);
end

function P = checkPaths(paths, n, fc, fs)
P = struct('gain', {}, 'excess_m', {}, 'kind', {}, 'x', {});
if isempty(paths), return; end
if ~isstruct(paths)
    error('fsr_channel:path', 'cfg.paths must be a struct array (or empty).');
end
allowed = {'gain', 'excess_m', 'modulation', 'label'};
f = fieldnames(paths);
bad = setdiff(f, allowed);
if ~isempty(bad)
    error('fsr_channel:path', 'Unknown path field ''%s''. Allowed: %s.', bad{1}, strjoin(allowed, ', '));
end
for req = {'gain', 'excess_m', 'modulation'}
    if ~isfield(paths, req{1})
        error('fsr_channel:path', 'Every path needs the field ''%s''.', req{1});
    end
end
if isempty(fc)
    error('fsr_channel:path', 'cfg.fc (carrier, Hz) is required when cfg.paths is non-empty.');
end
for k = 1:numel(paths)
    p = paths(k);
    if ~(isnumeric(p.gain) && isscalar(p.gain) && isfinite(p.gain))
        error('fsr_channel:path', 'paths(%d).gain must be a finite (complex) numeric scalar.', k);
    end
    if ~(isnumeric(p.excess_m) && isreal(p.excess_m) && isscalar(p.excess_m) && isfinite(p.excess_m))
        error('fsr_channel:path', 'paths(%d).excess_m must be a finite real scalar [m].', k);
    end
    if p.excess_m <= 0
        error('fsr_channel:path', ['paths(%d).excess_m = %g m: a reflected path must be strictly longer ' ...
            'than the direct path (excess_m > 0). Negative values would be non-causal; zero would ' ...
            'duplicate the direct path.'], k, p.excess_m);
    end
    D = p.excess_m / 299792458.0 * fs;
    if D > n - 1
        error('fsr_channel:path', 'paths(%d) delay %.3g samples exceeds the record length (%d samples).', k, D, n);
    end
    m = p.modulation;
    if ischar(m) || isstring(m)
        m = lower(char(m));
        if ~ismember(m, {'static', 'direct'})
            error('fsr_channel:path', 'paths(%d).modulation must be ''static'', ''direct'' or an n-by-1 vector.', k);
        end
        P(k).kind = m; P(k).x = [];
    elseif isnumeric(m) && iscolumn(m) && numel(m) == n && all(isfinite(m))
        P(k).kind = 'custom'; P(k).x = double(m);
    else
        error('fsr_channel:path', ['paths(%d).modulation must be ''static'', ''direct'' or a finite ' ...
            '%d-by-1 numeric column (got %s %s).'], k, n, class(m), mat2str(size(m)));
    end
    P(k).gain = double(p.gain);
    P(k).excess_m = double(p.excess_m);
end
end

function y = fracDelay(x, D)
% Causal linear-interpolation delay by D >= 0 samples, pre-history held at x(1).
n = numel(x);
m = floor(D); f = D - m;
idx = (1:n)' - m;
i0 = max(idx, 1);
i1 = max(idx - 1, 1);
y = (1 - f)*x(i0) + f*x(i1);
end

function x = ouProcess(n, fs, tau)
% Identical to ou() in fsr_detector_comparison.m (same draws, same order).
rho = exp(-1/(fs*tau));
e = randn(n,1)*sqrt(1 - rho^2);
x = filter(1, [1 -rho], e, rho*randn);
end
