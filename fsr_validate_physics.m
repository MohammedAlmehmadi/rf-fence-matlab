%% RF Drone Fence - staged physics verification (v2)
% Run from the folder that contains the fsr_*.m files:  >> fsr_validate_physics
%
% What this script CAN establish: internal consistency of the scalar
% Babinet/Fresnel model, numerical convergence of the Kirchhoff quadrature, and
% how far the paraxial closed form agrees with the exact-geometry scalar integral.
% What it CANNOT establish: that scalar Kirchhoff diffraction describes a real
% drone, or any real-world detection performance. Nothing here is measured data.
%
% Categories (kept separate in the printout):
%   [UNIT]  formula limits and code behaviour
%   [CONV]  numerical convergence of the quadrature
%   [MODEL] agreement between two mathematical models
%   [INFO]  reported, not asserted (sensitivity, regions outside the main lobe)
clear; clc; close all;

R = struct('cat',{},'name',{},'pass',{},'detail',{});

fcList = [2.4e9, 5.8e9];
L  = 20.0;  xc = 10.0;  hAnt = 1.5;          % baseline and antenna height [m]
w  = 0.30;  h = 0.20;   v = 5.0;            % nominal screen and speed
d1 = xc;    d2 = L - xc;  deff = d1*d2/(d1+d2);
geom = struct('src',[0 0 hAnt],'rcv',[L 0 hAnt],'xc',xc,'w',w,'h',h,'opacity',1);

fprintf('RF Drone Fence - staged physics verification (v2)\n');
fprintf('Scalar Babinet/Fresnel model of a flat opaque %.2f x %.2f m screen\n', w, h);
fprintf('Baseline %.0f m at height %.1f m, crossing at x = %.0f m (d_eff = %.2f m)\n\n', L, hAnt, xc, deff);

%% Stage 0 - dimensions and regime (printed, not asserted)
fprintf('Stage 0 - dimensions and validity regime\n');
for fc = fcList
    lam = fsr_wavelength(fc);
    fprintf('  %.1f GHz: lambda = %.4f m | w/lambda = %.2f, h/lambda = %.2f | ', fc/1e9, lam, w/lam, h/lam);
    fprintf('Fresnel radius sqrt(lambda*deff) = %.3f m | N_F = w^2/(lambda*deff) = %.3f\n', ...
        sqrt(lam*deff), w^2/(lam*deff));
end
fprintf(['  Note: the screen is only ~1.6-5.8 wavelengths across. Kirchhoff (physical-optics)\n' ...
         '  boundary values are approximate at this size; both models share that approximation.\n\n']);

%% Stage 1 - unit tests of limits
fprintf('Stage 1 - limits\n');
lam = fsr_wavelength(2.4e9);
yT = linspace(-3,3,41);
s0p = fsr_babinet_paraxial(yT, 0, w, h, d1, d2, lam, 0);
s0k = fsr_scattered_field('kirchhoff', setfield(geom,'opacity',0), yT, 0, lam, 8); %#ok<SFLD>
R = rec(R,'UNIT','No target (opacity 0) gives E/E0 = 1 exactly, both models', ...
    max(abs(s0p))==0 && max(abs(s0k))==0, sprintf('max|s| = %g / %g', max(abs(s0p)), max(abs(s0k))));
sBlock = fsr_babinet_paraxial(0,0,1e4,1e4,d1,d2,lam,1);
R = rec(R,'UNIT','Very large screen blocks the link (|1-s| -> 0)', abs(1-sBlock) < 1e-3, sprintf('|1-s| = %.3g', abs(1-sBlock)));
sFar = fsr_babinet_paraxial(500,0,w,h,d1,d2,lam,1);
R = rec(R,'UNIT','Screen far from the LOS is transparent (|s| -> 0)', abs(sFar) < 1e-3, sprintf('|s| = %.3g', abs(sFar)));
sSmall = fsr_babinet_paraxial(0,0,0.05,0.05,d1,d2,lam,1);
eSmall = abs(abs(sSmall) - 0.05*0.05/(lam*deff))/(0.05*0.05/(lam*deff));
R = rec(R,'UNIT','Small-target limit |s| = A/(lambda*d_eff)', eSmall < 0.02, sprintf('rel. error = %.3g', eSmall));
r24 = abs(fsr_babinet_paraxial(0,0,0.1,0.1,d1,d2,fsr_wavelength(2.4e9),1));
r58 = abs(fsr_babinet_paraxial(0,0,0.1,0.1,d1,d2,fsr_wavelength(5.8e9),1));
eRat = abs(r58/r24 - 5.8/2.4)/(5.8/2.4);
R = rec(R,'UNIT','Small-target |s| scales as 1/lambda', eRat < 0.03, sprintf('ratio %.4f vs %.4f', r58/r24, 5.8/2.4));

%% Stage 2 - symmetry, reciprocity, geometry bookkeeping
fprintf('Stage 2 - symmetry and geometry\n');
yS = linspace(-2,2,41);
for fc = fcList
    lam = fsr_wavelength(fc);
    sp  = fsr_scattered_field('paraxial',  geom,  yS, 0.3, lam);
    spm = fsr_scattered_field('paraxial',  geom, -yS, 0.3, lam);
    spz = fsr_scattered_field('paraxial',  geom,  yS,-0.3, lam);
    sk  = fsr_scattered_field('kirchhoff', geom,  yS, 0.3, lam, 8);
    skm = fsr_scattered_field('kirchhoff', geom, -yS, 0.3, lam, 8);
    skz = fsr_scattered_field('kirchhoff', geom,  yS,-0.3, lam, 8);
    e1 = max([abs(sp-spm), abs(sp-spz)])/max(abs(sp));
    e2 = max([abs(sk-skm), abs(sk-skz)])/max(abs(sk));
    R = rec(R,'UNIT',sprintf('%.1f GHz: s(y,z) = s(-y,z) = s(y,-z), both models', fc/1e9), ...
        e1 < 1e-12 && e2 < 1e-10, sprintf('paraxial %.1e, Kirchhoff %.1e', e1, e2));
    gA = geom; gA.xc = 6;  gB = geom; gB.xc = 14;
    kA = fsr_scattered_field('kirchhoff', gA, yS, 0.3, lam, 8);
    kB = fsr_scattered_field('kirchhoff', gB, yS, 0.3, lam, 8);
    eR = max(abs(kA-kB))/max(abs(kA));
    R = rec(R,'UNIT',sprintf('%.1f GHz: reciprocity, crossing at 6 m = crossing at 14 m', fc/1e9), ...
        eR < 1e-10, sprintf('rel. diff %.1e', eR));
end
% The coordinate trap the wrapper prevents: raw Kirchhoff call with an offset as if absolute
lam = fsr_wavelength(2.4e9);
yNull = lam*deff/w;  yM = linspace(-yNull, yNull, 61);
sOK  = fsr_scattered_field('kirchhoff', geom, yM, 0.4, lam, 8);
sBad = fsr_kirchhoff_numeric(geom.src, geom.rcv, xc, yM, 0.4, w, h, lam, 8, 1, 'exact');  % z=0.4 m ABSOLUTE
sPx  = fsr_scattered_field('paraxial', geom, yM, 0.4, lam);
eOK  = max(abs(sPx - sOK))/max(abs(sOK));
eBad = max(abs(sPx - sBad))/max(abs(sOK));
R = rec(R,'UNIT','Antennas at 1.5 m: wrapper keeps both models on the same target', eOK < 0.05, ...
    sprintf('wrapper %.4f; raw call with z as offset would differ by %.2f', eOK, eBad));
t0 = 3.0;  tt = linspace(0,6,601);  yt = v*(tt - t0);
[~,i0] = min(abs(yt));
R = rec(R,'UNIT','Trajectory y(t) = v(t - t0) is on the LOS plane at t0', abs(tt(i0)-t0) < 1e-12 && abs(yt(i0)) < 1e-12, ...
    sprintf('y(t0) = %.1e m; main-lobe duration 2*lambda*deff/(w*v) = %.3f s at 2.4 GHz', yt(i0), 2*yNull/v));
try
    fsr_kirchhoff_numeric([0 0 0],[20 1 0],10,0,0,w,h,lam);   % tilted baseline must be refused
    refused = false;
catch
    refused = true;
end
R = rec(R,'UNIT','Tilted baseline is refused (screen orientation assumes x-normal)', refused, '');

%% Stage 3 - independent sign/normalisation check: opaque disk on axis (Arago spot)
fprintf('Stage 3 - Arago spot (opaque disk on the LOS)\n');
for fc = fcList
    lam = fsr_wavelength(fc);
    for a = [0.15, 1.0]
        [sE, sP] = fsr_disk_axial(a, d1, d2, lam);
        R = rec(R,'UNIT',sprintf('%.1f GHz, disk a = %.2f m: paraxial |E/E0| = 1 exactly', fc/1e9, a), ...
            abs(abs(1-sP)-1) < 1e-12, sprintf('|1-s_parax| - 1 = %.1e', abs(1-sP)-1));
        if a < 0.5
            % small disk (a/d = 0.015): paraxial regime, the two must agree closely
            R = rec(R,'MODEL',sprintf('%.1f GHz, disk a = %.2f m: exact 1-D integral = paraxial', fc/1e9, a), ...
                abs(sE-sP) < 1e-3, sprintf('|s_exact - s_parax| = %.2e, |E/E0|_exact = %.5f', abs(sE-sP), abs(1-sE)));
        else
            % large disk (a/d = 0.1): non-paraxial terms are physically expected, so report only
            R = rec(R,'INFO',sprintf('%.1f GHz, disk a = %.2f m (non-paraxial)', fc/1e9, a), true, ...
                sprintf('|s_exact - s_parax| = %.2e, |E/E0|_exact = %.5f', abs(sE-sP), abs(1-sE)));
        end
    end
end

%% Stage 4 - bridge test: same paraxial integrand, closed form vs quadrature
fprintf('Stage 4 - closed form vs numerically integrated paraxial integrand\n');
for fc = fcList
    lam = fsr_wavelength(fc);
    yNull = lam*deff/w;  yM = linspace(-yNull, yNull, 61);
    sC = fsr_scattered_field('paraxial',           geom, yM, 0, lam);
    sQ = fsr_scattered_field('kirchhoff_paraxial', geom, yM, 0, lam, 8);
    eB = max(abs(sC-sQ))/max(abs(sQ));
    R = rec(R,'MODEL',sprintf('%.1f GHz main lobe: closed form = quadrature of the same integrand', fc/1e9), ...
        eB < 2e-3, sprintf('max rel. diff %.2e (tests Fresnel routine + quadrature, not physics)', eB));
end

%% Stage 5 - convergence of the exact Kirchhoff quadrature, near and away from the LOS
fprintf('Stage 5 - quadrature convergence (reference: 64 points per wavelength)\n');
pplList = [4 8 16 32 64];
for fc = fcList
    lam = fsr_wavelength(fc);
    yNull = lam*deff/w;  yC = [0, yNull, 2*yNull, 5];
    S = zeros(numel(pplList), numel(yC));
    for p = 1:numel(pplList)
        S(p,:) = fsr_scattered_field('kirchhoff', geom, yC, 0, lam, pplList(p));
    end
    rel = abs(S(1:end-1,:) - S(end,:)) ./ abs(S(end,:));
    fprintf('  %.1f GHz  y = %s m\n', fc/1e9, mat2str(round(yC,3)));
    for p = 1:numel(pplList)-1
        fprintf('    %2d ppl: rel. change vs 64 ppl = %s\n', pplList(p), mat2str(rel(p,:),3));
    end
    R = rec(R,'CONV',sprintf('%.1f GHz: 8 ppl within 1%% of 64 ppl for |y| <= y_null', fc/1e9), ...
        all(rel(2,1:2) < 1e-2), sprintf('8 ppl: %s at y = 0, y_null', mat2str(rel(2,1:2),3)));
    R = rec(R,'INFO',sprintf('%.1f GHz: 8 ppl at y = 5 m', fc/1e9), true, ...
        sprintf('rel. change %.3g -> use >= 32 ppl beyond the main lobe', rel(2,4)));
end

%% Stage 6 - paraxial closed form vs exact-geometry Kirchhoff
fprintf('Stage 6 - paraxial vs exact-geometry Kirchhoff (32 ppl)\n');
fig1 = figure('Name','FSR v2: paraxial vs exact Kirchhoff','Color','w');
tiledlayout(fig1,2,2,'TileSpacing','compact','Padding','compact');
fig2 = figure('Name','FSR v2: model difference vs offset','Color','w');
ax2 = axes(fig2); hold(ax2,'on'); set(ax2,'YScale','log'); grid(ax2,'on');
cols = {[0.10 0.35 0.80],[0.85 0.33 0.10]};
for jj = 1:numel(fcList)
    fc = fcList(jj);  lam = fsr_wavelength(fc);
    yNull = lam*deff/w;
    % (a) inside the main lobe: asserted
    yM = linspace(-yNull, yNull, 61);
    sP = fsr_scattered_field('paraxial',  geom, yM, 0, lam);
    sK = fsr_scattered_field('kirchhoff', geom, yM, 0, lam, 16);
    eS = max(abs(sP-sK))/max(abs(sK));
    eF = max(abs((1-sP)-(1-sK)));            % field error relative to |E0| = 1
    R = rec(R,'MODEL',sprintf('%.1f GHz main lobe |y| <= %.2f m: paraxial vs exact', fc/1e9, yNull), ...
        eS < 0.05, sprintf('scattered-term error %.4f; field error %.5f of |E0|', eS, eF));
    % (b) wider range: reported only
    yW = linspace(-5, 5, 201);
    sPw = fsr_scattered_field('paraxial',           geom, yW, 0, lam);
    sQw = fsr_scattered_field('kirchhoff_paraxial', geom, yW, 0, lam, 32);
    sKw = fsr_scattered_field('kirchhoff',          geom, yW, 0, lam, 32);
    for Y = [2 3 5]
        m = abs(yW) <= Y + 1e-12;
        R = rec(R,'INFO',sprintf('%.1f GHz |y| <= %d m', fc/1e9, Y), true, sprintf( ...
            'paraxial vs exact %.4f (of which quadrature %.1e)', ...
            max(abs(sPw(m)-sKw(m)))/max(abs(sKw(m))), max(abs(sPw(m)-sQw(m)))/max(abs(sKw(m)))));
    end
    % plots
    figure(fig1);
    nexttile((jj-1)*2+1);
    plot(yW, 20*log10(abs(1-sPw)), '-', 'Color', cols{jj}, 'LineWidth', 1.4); hold on;
    plot(yW, 20*log10(abs(1-sKw)), 'k:', 'LineWidth', 1.2);
    xline([-yNull yNull], '--', 'Color', [0.5 0.5 0.5]); hold off; grid on;
    xlabel('Along-track offset y (m)'); ylabel('|E/E_0| (dB)');
    title(sprintf('%.1f GHz amplitude (dashed: main lobe)', fc/1e9));
    legend('Paraxial closed form','Exact-geometry Kirchhoff','Location','southeast');
    nexttile((jj-1)*2+2);
    plot(yW, rad2deg(angle(1-sPw)), '-', 'Color', cols{jj}, 'LineWidth', 1.4); hold on;
    plot(yW, rad2deg(angle(1-sKw)), 'k:', 'LineWidth', 1.2);
    xline([-yNull yNull], '--', 'Color', [0.5 0.5 0.5]); hold off; grid on;
    xlabel('Along-track offset y (m)'); ylabel('arg(E/E_0) (deg)');
    title(sprintf('%.1f GHz phase', fc/1e9));
    plot(ax2, yW, abs(sPw-sKw)/max(abs(sKw)), '-', 'Color', cols{jj}, 'LineWidth', 1.4, ...
        'DisplayName', sprintf('%.1f GHz: paraxial vs exact', fc/1e9));
    plot(ax2, yW, max(abs(sPw-sQw),1e-12)/max(abs(sKw)), ':', 'Color', cols{jj}, 'LineWidth', 1.2, ...
        'DisplayName', sprintf('%.1f GHz: quadrature only', fc/1e9));
end
xlabel(ax2,'Along-track offset y (m)'); ylabel(ax2,'|s_{paraxial} - s| / max|s_{exact}|');
title(ax2,'Where the paraxial model departs from the exact-geometry integral');
legend(ax2,'Location','southeast');

%% Stage 7 - parameter sensitivity (closed form; reported, a few physical trends asserted)
fprintf('Stage 7 - sensitivity of the peak relative amplitude change (closed form)\n');
yP = linspace(-6, 6, 4001);
peakDev = @(fc,ww,hh,x1,z0) max(abs(abs(1 - fsr_babinet_paraxial(yP, z0, ww, hh, x1, L-x1, fsr_wavelength(fc), 1)) - 1));
fprintf('  %-8s %-11s %-6s %-6s %-14s\n','f [GHz]','w x h [m]','xc [m]','z0 [m]','peak |dA|/A');
cases = [2.4 0.30 0.20 10 0; 2.4 0.30 0.20 5 0; 2.4 0.30 0.20 2 0; 2.4 0.30 0.20 10 0.5; ...
         2.4 0.30 0.20 10 1.0; 2.4 0.20 0.10 10 0; 5.8 0.30 0.20 10 0; 5.8 0.20 0.10 10 0];
pd = zeros(size(cases,1),1);
for c = 1:size(cases,1)
    pd(c) = peakDev(cases(c,1)*1e9, cases(c,2), cases(c,3), cases(c,4), cases(c,5));
    fprintf('  %-8.1f %-11s %-6.0f %-6.1f %-14.4f\n', cases(c,1), sprintf('%.2fx%.2f',cases(c,2),cases(c,3)), ...
        cases(c,4), cases(c,5), pd(c));
end
R = rec(R,'UNIT','Signature grows toward an antenna (xc 10 -> 5 -> 2 m)', pd(1) < pd(2) && pd(2) < pd(3), ...
    sprintf('%.4f < %.4f < %.4f', pd(1), pd(2), pd(3)));
R = rec(R,'UNIT','Signature grows with frequency for the nominal screen', pd(7) > pd(1), sprintf('%.4f > %.4f', pd(7), pd(1)));
R = rec(R,'UNIT','Smaller screen gives a smaller signature', pd(6) < pd(1), sprintf('%.4f < %.4f', pd(6), pd(1)));
fprintf('  Speed does not change s(y); it only rescales time: main-lobe duration = %.2f s at 5 m/s, %.2f s at 10 m/s (2.4 GHz).\n\n', ...
    2*fsr_wavelength(2.4e9)*deff/w/5, 2*fsr_wavelength(2.4e9)*deff/w/10);

%% Summary
fprintf('\n================ SUMMARY ================\n');
cats = {'UNIT','CONV','MODEL','INFO'};
for c = 1:numel(cats)
    idx = find(strcmp({R.cat}, cats{c}));
    for i = idx
        if strcmp(cats{c},'INFO'), tag = ' INFO '; elseif R(i).pass, tag = ' PASS '; else, tag = '*FAIL*'; end
        fprintf('[%s][%-5s] %s', tag, R(i).cat, R(i).name);
        if ~isempty(R(i).detail), fprintf('  (%s)', R(i).detail); end
        fprintf('\n');
    end
end
nFail = sum(~[R(~strcmp({R.cat},'INFO')).pass]);
fprintf('\n%d asserted checks, %d failed.\n', sum(~strcmp({R.cat},'INFO')), nFail);
fprintf(['These are unit, convergence and model-agreement checks only. There is no\n' ...
         'validation against measurements and no experimental field validation.\n']);
if nFail > 0
    error('fsr:verification','%d physics checks failed - do not integrate into the app yet.', nFail);
end

%% Local function
function R = rec(R, cat, name, pass, detail)
R(end+1) = struct('cat',cat,'name',name,'pass',logical(pass),'detail',detail);
if strcmp(cat,'INFO'), tag = ' INFO '; elseif pass, tag = ' PASS '; else, tag = '*FAIL*'; end
fprintf('  [%s] %s', tag, name);
if ~isempty(detail), fprintf('  (%s)', detail); end
fprintf('\n');
end
