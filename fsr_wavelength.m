function lambda = fsr_wavelength(fc_hz)
%FSR_WAVELENGTH Wavelength in metres from carrier frequency in Hz.
% Uses the exact speed-of-light constant used by the Python source model.
validateattributes(fc_hz, {'numeric'}, {'real','finite','scalar','positive'});
c0 = 299792458.0;
lambda = c0 / fc_hz;
end
