RF Drone Fence MATLAB physics layer - v2 (verification stage)
==============================================================

Put all seven .m files in ONE folder (you can create each file in MATLAB Online
with  >> edit <name>.m  and paste its contents). Then run:

    >> fsr_validate_physics

Files
-----
unchanged from v1:
  fsr_wavelength.m         wavelength from frequency
  fsr_fresnel_cs.m         toolbox-free Fresnel C/S (Cephes). Checked: matches the
                           reference Fresnel integrals to 4e-15 for |x| <= 60
  fsr_babinet_paraxial.m   closed-form paraxial Babinet/Fresnel model; y, z are
                           OFFSETS from the line of sight
replaced:
  fsr_kirchhoff_numeric.m  same exact-geometry integral as v1, plus
                           - 'paraxial' mode (same integral, paraxial integrand) for
                             the bridge test
                           - refuses tilted baselines and xc outside the baseline
                           - documents that y, z are ABSOLUTE coordinates
  fsr_validate_physics.m   staged checks: UNIT / CONV / MODEL / INFO
new:
  fsr_scattered_field.m    one entry point taking LOS-relative offsets for every
                           model (prevents the z-reference mismatch when the
                           antennas are at 1.5 m)
  fsr_disk_axial.m         opaque disk on axis (Poisson-Arago spot) reference

Expected key numbers (from a line-by-line NumPy transcription of these files;
MATLAB itself has NOT been run on them - please report any difference)
---------------------------------------------------------------------------
Stage 0  2.4 GHz: lambda 0.1249 m, w/lambda 2.40, N_F 0.144
         5.8 GHz: lambda 0.0517 m, w/lambda 5.80, N_F 0.348
Stage 1  large screen |1-s| ~ 1.0e-4; far |s| ~ 6.3e-5; small-target err ~ 9e-7;
         1/lambda ratio 2.4165 (expected 2.4167)
Stage 2  symmetry / reciprocity errors ~ 1e-15 or 0; wrapper check 0.0271
         (a raw call passing an LOS offset as an absolute z differs by ~0.91)
         main-lobe duration 0.833 s at 2.4 GHz, 5 m/s
Stage 3  disk a = 0.15 m: |s_exact - s_parax| ~ 2e-5 (2.4 GHz), ~ 5e-5 (5.8 GHz)
         disk a = 1.00 m (INFO): ~ 1.8e-2 and ~ 3.2e-2 (non-paraxial, expected)
Stage 4  bridge (quadrature only): 7.6e-4 (2.4 GHz), 1.3e-4 (5.8 GHz)
Stage 5  8 ppl vs 64 ppl:  2.4 GHz  y = 0, 2.08, 4.16, 5 m -> ~2e-4, 3.5e-3, 1.4e-2, 1.9e-2
                           5.8 GHz  y = 0, 0.86, 1.72, 5 m -> ~8e-5, 7e-4, 2.6e-3, 2.0e-2
Stage 6  main lobe, paraxial vs exact: 0.0242 (2.4 GHz), 0.0028 (5.8 GHz)
         field error (fraction of |E0|): 0.0023 and 0.00065
         INFO |y|<=2, 3, 5 m: 2.4 GHz 0.024, 0.193, 0.253 ; 5.8 GHz 0.054, 0.153, 0.160
         (quadrature part of these <= 2e-4, so the rest is the paraxial approximation)
Stage 7  peak |dA|/A: 0.0856 (nominal), 0.110 (xc 5 m), 0.197 (xc 2 m),
         0.0904 (z0 0.5 m), 0.0756 (z0 1 m), 0.0304 (0.2x0.1 m), 0.1775 (5.8 GHz)
Summary: every asserted check should PASS.

Scope
-----
Scalar physical-optics diffraction of a flat opaque rectangle. Not a full-wave
solver; not validated against measurements; no experimental field validation.
