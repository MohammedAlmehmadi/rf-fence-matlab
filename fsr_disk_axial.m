function [sExact, sParax] = fsr_disk_axial(a, d1, d2, lambda)
%FSR_DISK_AXIAL Opaque circular disk centred on the line of sight (Arago spot test).
%
%   [sExact, sParax] = fsr_disk_axial(a, d1, d2, lambda)
%
% a      : disk radius in metres;  d1, d2 : Tx-disk and disk-Rx distances (m).
% sExact : the 'exact'-mode Kirchhoff integral reduced to one dimension by
%          rotational symmetry and evaluated with MATLAB integral():
%            s = (1i/lambda) * 2*pi * int_0^a L/(r1 r2) exp(-1i k (r1+r2-L))
%                                   * (d1/r1 + d2/r2)/2 * rho d rho
% sParax : the same integral with the paraxial integrand, which has the closed form
%            s = 1 - exp(-1i*k*a^2/(2*deff)),  deff = d1*d2/(d1+d2).
%
% Physics check: behind an opaque disk on axis E/E0 = 1 - s = exp(-1i*k*a^2/(2*deff))
% in the paraxial limit, so |E/E0| = 1 (Poisson-Arago spot). This tests the sign,
% the 1i/lambda prefactor and the normalisation with a geometry that is
% independent of the rectangular Fresnel-integral code.
validateattributes(a, {'numeric'}, {'real','finite','scalar','positive'});
validateattributes(d1, {'numeric'}, {'real','finite','scalar','positive'});
validateattributes(d2, {'numeric'}, {'real','finite','scalar','positive'});
validateattributes(lambda, {'numeric'}, {'real','finite','scalar','positive'});
k = 2*pi/lambda;
L = d1 + d2;
deff = d1*d2/L;
f = @(rho) (L ./ (sqrt(d1^2+rho.^2).*sqrt(d2^2+rho.^2))) ...
    .* exp(-1i*k.*(sqrt(d1^2+rho.^2) + sqrt(d2^2+rho.^2) - L)) ...
    .* 0.5 .* (d1./sqrt(d1^2+rho.^2) + d2./sqrt(d2^2+rho.^2)) .* rho;
sExact = (1i/lambda) * 2*pi * integral(f, 0, a, 'RelTol', 1e-10, 'AbsTol', 1e-14);
sParax = 1 - exp(-1i*k*a^2/(2*deff));
end
