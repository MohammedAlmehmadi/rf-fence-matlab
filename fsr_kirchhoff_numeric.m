function s = fsr_kirchhoff_numeric(src,rcv,xc,y,z,w,h,lambda,ptsPerLambda,opacity,mode)
%FSR_KIRCHHOFF_NUMERIC Numerical scalar Fresnel-Kirchhoff (Babinet aperture) integral.
%
%   s = fsr_kirchhoff_numeric(src,rcv,xc,y,z,w,h,lambda)
%   s = fsr_kirchhoff_numeric(...,ptsPerLambda,opacity,mode)
%
% Returns the complex scattered fraction s such that E/E0 = 1 - s, where E0 is
% the free-space direct field exp(-1i*k*L)/L (time convention exp(+1i*w*t)).
%
% GEOMETRY (v2: now checked, not assumed)
%   src, rcv : [x y z] in metres. The baseline MUST be parallel to the x-axis
%              (same y and z for src and rcv), because the screen is placed in
%              the plane x = xc with its normal along x.
%   xc       : x-position of the screen, strictly between src(1) and rcv(1).
%   y, z     : ABSOLUTE coordinates of the screen centre (not offsets from the
%              line of sight). Use fsr_scattered_field.m if you want to pass
%              offsets from the line of sight; it converts them for you.
%
% MODE
%   'exact'    (default) exact path lengths r1+r2, spreading L/(r1*r2) and the
%              obliquity factor (cos t1 + cos t2)/2.
%   'paraxial' the SAME integral with the paraxial integrand: phase
%              rho^2/2*(1/d1+1/d2), spreading L/(d1*d2), obliquity 1. Its exact
%              value is the closed form in fsr_babinet_paraxial.m, so comparing
%              the two tests the Fresnel routine and the quadrature only.
%
% Quadrature: midpoint rule on an ny-by-nz grid, ptsPerLambda points per
% wavelength along each side (default 8). Accuracy degrades as the screen moves
% away from the line of sight (more phase variation across the screen); see
% the convergence test in fsr_validate_physics.m before trusting |y| > a few m.
%
% This is scalar physical-optics diffraction for a flat opaque screen. It is
% NOT a Maxwell/full-wave solution and is not ground truth for a real drone.

if nargin < 9  || isempty(ptsPerLambda), ptsPerLambda = 8; end
if nargin < 10 || isempty(opacity), opacity = 1.0; end
if nargin < 11 || isempty(mode), mode = 'exact'; end
validateattributes(src, {'numeric'}, {'real','finite','vector','numel',3});
validateattributes(rcv, {'numeric'}, {'real','finite','vector','numel',3});
validateattributes(xc, {'numeric'}, {'real','finite','scalar'});
validateattributes(y, {'numeric'}, {'real','finite'});
validateattributes(z, {'numeric'}, {'real','finite'});
validateattributes(w, {'numeric'}, {'real','finite','scalar','positive'});
validateattributes(h, {'numeric'}, {'real','finite','scalar','positive'});
validateattributes(lambda, {'numeric'}, {'real','finite','scalar','positive'});
validateattributes(ptsPerLambda, {'numeric'}, {'real','finite','scalar','>=',1});
validateattributes(opacity, {'numeric'}, {'real','finite','scalar','>=',0,'<=',1});
mode = validatestring(mode, {'exact','paraxial'});
if ~isscalar(z) && ~isequal(size(y),size(z))
    error('fsr:size','z must be scalar or have the same size as y.');
end

src = reshape(double(src),1,3);
rcv = reshape(double(rcv),1,3);
if abs(src(2)-rcv(2)) > 1e-9 || abs(src(3)-rcv(3)) > 1e-9
    error('fsr:geometry', ['Baseline must be parallel to the x-axis ' ...
        '(src and rcv need equal y and z). Got src=[%g %g %g], rcv=[%g %g %g].'], src, rcv);
end
if ~(xc > src(1) && xc < rcv(1))
    error('fsr:geometry','xc=%g must lie strictly between src(1)=%g and rcv(1)=%g.', xc, src(1), rcv(1));
end

L  = rcv(1) - src(1);
d1 = xc - src(1);
d2 = rcv(1) - xc;
k  = 2*pi/lambda;

ny = max(4, ceil(w/lambda*ptsPerLambda));
nz = max(4, ceil(h/lambda*ptsPerLambda));
gy = ((1:ny)-0.5)/ny*w - w/2;
gz = ((1:nz)-0.5)/nz*h - h/2;
dA = (w/ny)*(h/nz);
[GY,GZ] = ndgrid(gy,gz);
GY = GY(:); GZ = GZ(:);

yShape = size(y);
yv = double(y(:));
if isscalar(z), zv = repmat(double(z),size(yv)); else, zv = double(z(:)); end
out = complex(zeros(size(yv)));

for ii = 1:numel(yv)
    py = yv(ii) + GY;               % absolute y of each panel point
    pz = zv(ii) + GZ;               % absolute z of each panel point
    rho2 = (py-src(2)).^2 + (pz-src(3)).^2;   % squared distance from the LOS
    if strcmp(mode,'exact')
        r1 = sqrt(d1^2 + rho2);
        r2 = sqrt(d2^2 + rho2);
        obliquity = 0.5 .* (d1./r1 + d2./r2);
        integrand = L ./ (r1.*r2) .* exp(-1i*k.*(r1 + r2 - L)) .* obliquity;
    else
        integrand = (L/(d1*d2)) .* exp(-1i*k.*rho2.*(0.5*(1/d1 + 1/d2)));
    end
    out(ii) = sum(integrand) * dA;
end
s = reshape(opacity .* (1i/lambda) .* out, yShape);
end
