function s = fsr_babinet_paraxial(y,z,w,h,d1,d2,lambda,opacity)
%FSR_BABINET_PARAXIAL Scalar Babinet/Fresnel model for a rectangular screen.
%
% Returns the complex scattered fraction s, where the normalised field is
% E/E0 = 1 - s. y and z are target-centre offsets from the LOS, in metres.
%
% This is a scalar, paraxial, opaque-rectangle approximation; it is NOT a
% full-wave electromagnetic simulation and does not model a real drone shape.
if nargin < 8 || isempty(opacity), opacity = 1.0; end
validateattributes(y, {'numeric'}, {'real','finite'});
validateattributes(z, {'numeric'}, {'real','finite'});
validateattributes(w, {'numeric'}, {'real','finite','scalar','positive'});
validateattributes(h, {'numeric'}, {'real','finite','scalar','positive'});
validateattributes(d1, {'numeric'}, {'real','finite','scalar','positive'});
validateattributes(d2, {'numeric'}, {'real','finite','scalar','positive'});
validateattributes(lambda, {'numeric'}, {'real','finite','scalar','positive'});
validateattributes(opacity, {'numeric'}, {'real','finite','scalar','>=',0,'<=',1});
if ~isscalar(z) && ~isequal(size(y),size(z))
    error('z must be scalar or have the same size as y.');
end

deff = d1*d2/(d1+d2);
scale = sqrt(2/(lambda*deff));
gy = fresnelFactor(y - w/2, y + w/2, scale);
gz = fresnelFactor(z - h/2, z + h/2, scale);
s = opacity .* (0.5i) .* gy .* gz;
end

function g = fresnelFactor(a,b,scale)
[Ca,Sa] = fsr_fresnel_cs(a .* scale);
[Cb,Sb] = fsr_fresnel_cs(b .* scale);
g = (Cb-Ca) - 1i.*(Sb-Sa);
end
