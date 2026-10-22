function s = fsr_scattered_field(model, geom, y0, z0, lambda, ptsPerLambda)
%FSR_SCATTERED_FIELD One entry point for both models, using LOS-relative offsets.
%
%   s = fsr_scattered_field('paraxial', geom, y0, z0, lambda)
%   s = fsr_scattered_field('kirchhoff', geom, y0, z0, lambda, ptsPerLambda)
%   s = fsr_scattered_field('kirchhoff_paraxial', geom, y0, z0, lambda, ptsPerLambda)
%
% geom is a struct with fields
%   src, rcv : [x y z] antenna positions in metres (baseline parallel to x)
%   xc       : crossing point along x, metres
%   w, h     : screen width (along y) and height (along z), metres
%   opacity  : 0..1 effective-area factor (1 = fully opaque)
% y0, z0     : target-centre OFFSETS from the line of sight, metres
%              (y0 = along-track, z0 = vertical). y0 may be a vector.
%
% Why this file exists: fsr_babinet_paraxial takes offsets from the LOS, while
% fsr_kirchhoff_numeric takes absolute coordinates. With antennas at z = 1.5 m,
% passing the same z to both silently puts the two targets 1.5 m apart. This
% wrapper does the conversion in one place so the app and tests cannot mix them.
%
% Returns s with E/E0 = 1 - s (E0 = free-space direct field).
if nargin < 6 || isempty(ptsPerLambda), ptsPerLambda = 8; end
model = char(model);
if ~ismember(model, {'paraxial','kirchhoff','kirchhoff_paraxial'})
    error('fsr:model','model must be ''paraxial'', ''kirchhoff'' or ''kirchhoff_paraxial''.');
end
src = reshape(geom.src,1,3);
rcv = reshape(geom.rcv,1,3);
if abs(src(2)-rcv(2)) > 1e-9 || abs(src(3)-rcv(3)) > 1e-9
    error('fsr:geometry','Baseline must be parallel to the x-axis.');
end
d1 = geom.xc - src(1);
d2 = rcv(1) - geom.xc;
if d1 <= 0 || d2 <= 0
    error('fsr:geometry','Crossing point must lie strictly between Tx and Rx.');
end
switch model
    case 'paraxial'
        s = fsr_babinet_paraxial(y0, z0, geom.w, geom.h, d1, d2, lambda, geom.opacity);
    case 'kirchhoff'
        s = fsr_kirchhoff_numeric(src, rcv, geom.xc, src(2)+y0, src(3)+z0, ...
            geom.w, geom.h, lambda, ptsPerLambda, geom.opacity, 'exact');
    case 'kirchhoff_paraxial'
        s = fsr_kirchhoff_numeric(src, rcv, geom.xc, src(2)+y0, src(3)+z0, ...
            geom.w, geom.h, lambda, ptsPerLambda, geom.opacity, 'paraxial');
end
end
