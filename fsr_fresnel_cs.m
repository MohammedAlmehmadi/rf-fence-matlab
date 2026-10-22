function [C,S] = fsr_fresnel_cs(x)
%FSR_FRESNEL_CS Fresnel cosine/sine integrals for real x, toolbox-free.
% C(x) = integral_0^x cos(pi*t^2/2) dt
% S(x) = integral_0^x sin(pi*t^2/2) dt
%
% Uses stable rational approximations in the small-argument region and
% asymptotic auxiliary functions for larger arguments. Coefficients/form are
% based on the Cephes Math Library Fresnel routine by Stephen L. Moshier.
% Original routine notice: Copyright 1984, 1987, 1989 Stephen L. Moshier.
% See README_MATLAB.txt for source attribution.

validateattributes(x, {'numeric'}, {'real','finite'});
origSize = size(x);
x = double(x(:));
ax = abs(x);
x2 = ax.^2;
C = zeros(size(x));
S = zeros(size(x));

% Coefficients for the small-argument rational approximations.
SN = [-2.99181919401019853726E3, 7.08840045257738576863E5, ...
      -6.29741486205862506537E7, 2.54890880573376359104E9, ...
      -4.42979518059697779103E10, 3.18016297876567817986E11];
SD = [2.81376268889994315696E2, 4.55847810806532581675E4, ...
      5.17343888770096400730E6, 4.19320245898111231129E8, ...
      2.24411795645340920940E10, 6.07366389490084639049E11];
CN = [-4.98843114573573548651E-8, 9.50428062829859605134E-6, ...
      -6.45191435683965050962E-4, 1.88843319396703850064E-2, ...
      -2.05525900955013891793E-1, 9.99999999999999998822E-1];
CD = [3.99982968972495980367E-12, 9.15439215774657478799E-10, ...
      1.25001862479598821474E-7, 1.22262789024179030997E-5, ...
      8.68029542941784300606E-4, 4.12142090722199792936E-2, ...
      1.00000000000000000118E0];

% Coefficients for asymptotic auxiliary functions f and g.
FN = [4.21543555043677546506E-1, 1.43407919780758885261E-1, ...
      1.15220955073585758835E-2, 3.45017939782574027900E-4, ...
      4.63613749287867322088E-6, 3.05568983790257605827E-8, ...
      1.02304514164907233465E-10, 1.72010743268161828879E-13, ...
      1.34283276233062758925E-16, 3.76329711269987889006E-20];
FD = [7.51586398353378947175E-1, 1.16888925859191382142E-1, ...
      6.44051526508858611005E-3, 1.55934409164153020873E-4, ...
      1.84627567348930545870E-6, 1.12699224763999035261E-8, ...
      3.60140029589371370404E-11, 5.88754533621578410010E-14, ...
      4.52001434074129701496E-17, 1.25443237090011264384E-20];
GN = [5.04442073643383265887E-1, 1.97102833525523411709E-1, ...
      1.87648584092575249293E-2, 6.84079380915393090172E-4, ...
      1.15138826111884280931E-5, 9.82852443688422223854E-8, ...
      4.45344415861750144738E-10, 1.08268041139020870318E-12, ...
      1.37555460633261799868E-15, 8.36354435630677421531E-19, ...
      1.86958710162783235106E-22];
GD = [1.47495759925128324529E0, 3.37748989120019970451E-1, ...
      2.53603741420338795122E-2, 8.14679107184306179049E-4, ...
      1.27545075667729118702E-5, 1.04314589657571990585E-7, ...
      4.60680728146520428211E-10, 1.10273215066240270757E-12, ...
      1.38796531259578871258E-15, 8.39158816283118707363E-19, ...
      1.86958710162783236342E-22];

small = (x2 < 2.5625);
if any(small)
    xx = ax(small);
    xx2 = x2(small);
    t = xx2.^2;
    S(small) = xx .* xx2 .* polyval(SN,t) ./ polyval([1,SD],t);
    C(small) = xx .* polyval(CN,t) ./ polyval(CD,t);
end

large = ~small;
if any(large)
    xx = ax(large);
    cc = zeros(size(xx));
    ss = zeros(size(xx));
    huge = xx > 36974;
    cc(huge) = 0.5;
    ss(huge) = 0.5;
    use = ~huge;
    if any(use)
        xu = xx(use);
        t = pi .* xu.^2;
        u = 1 ./ (t.^2);
        f = 1 - u .* polyval(FN,u) ./ polyval([1,FD],u);
        g = (1./t) .* polyval(GN,u) ./ polyval([1,GD],u);
        phase = (pi/2) .* xu.^2;
        cc(use) = 0.5 + (f.*sin(phase) - g.*cos(phase)) ./ (pi.*xu);
        ss(use) = 0.5 - (f.*cos(phase) + g.*sin(phase)) ./ (pi.*xu);
    end
    C(large) = cc;
    S(large) = ss;
end

sgn = sign(x);
C = reshape(C .* sgn, origSize);
S = reshape(S .* sgn, origSize);
end
