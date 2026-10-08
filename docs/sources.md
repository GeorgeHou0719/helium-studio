# Property and model sources

References below identify the data and equations used by the implementation.
Publication PDFs, unrelated research repositories, and local run histories are
not bundled. Dependency licenses and source-data terms remain with their owners.

| Source | Use in this project | Limits / implementation notes |
| --- | --- | --- |
| [Donnelly & Barenghi (1998), DOI 10.1063/1.556028](https://doi.org/10.1063/1.556028) | Saturated-helium pressure, latent heat, and liquid enthalpy splines. | Supports pot equilibrium quantities; not a general superfluid transport EOS. |
| [Hurly & Moldover (2000), NIST, Table A1](https://nvlpubs.nist.gov/nistpubs/jres/105/5/j55hur.pdf) | Dilute helium-4 viscosity and thermal conductivity, with published derivative interpolation. | Zero-density gas transport does not justify continuum no-slip flow at low pressure. |
| [NIST Cryogenic Materials: SS304](https://trc.nist.gov/cryogenics/materials/304Stainless/304Stainless_rev.htm) | Stainless-steel conductivity fit. | Equation 1–300 K; underlying data 4–300 K. |
| [NIST IR 8474 (2023)](https://nvlpubs.nist.gov/nistpubs/ir/2023/NIST.IR.8474.pdf) | Normal-fluid helium EOS scope and independent saturation check values. | He II is excluded. Tests compare Table 4 saturation quantities. |
| [CoolProp helium documentation](https://coolprop.org/fluid_properties/fluids/Helium.html) | HEOS normal-fluid thermodynamics and transport. | Enthalpy is shifted to match DB98 saturated liquid at 4.2 K; joining references is approximate. |
| [HiPace 300 M manufacturer page](https://www.pfeiffervacuum.com/ca/en/product/hipace-300-m.html) | Nominal helium preset: 240 L/s and compression ratio 10⁸. | A nominal preset is not a full pump operating map. |
| [HiPace 300 M / TM700 manual](https://www.pfeiffervacuum.com/global/assets/134568/1/Instruction%20Manual%20HiPace%20300%20M%20with%20TM%20700_EN_en.pdf) | Documents the configuration-dependent 215 L/s / >10⁸ values. | Confirm the revision; no curve is invented from conflicting nominal data. |
| [Berg (2005), DOI 10.1088/0026-1394/42/1/002](https://doi.org/10.1088/0026-1394/42/1/002) | Capillary-flow Poiseuille baseline. | Does not validate extrapolation to cryogenic flashing flow. |
| [NASA Glenn: Mass Flow Rate](https://www.grc.nasa.gov/www/k-12/BGP/mflchk.html) | Legacy ideal-gas nozzle/orifice choking relation. | Discharge coefficient remains empirical. |
| [NBS Technical Note 633 (1973)](https://nvlpubs.nist.gov/nistpubs/Legacy/TN/nbstechnicalnote633.pdf) | Homogeneous-equilibrium critical-flow framework. | Phase equilibrium and zero slip are engineering assumptions. |
| [NASA CR-202662, Arauz (1996/1997)](https://ntrs.nasa.gov/api/citations/19970005355/downloads/19970005355.pdf) | Pressure/enthalpy mixture descriptions and transport closures. | General model background, not helium microcapillary validation. |
| [COMSOL: Internal Film Resistance](https://doc.comsol.com/6.3/doc/com.comsol.help.pipe/pipe_ug_heattransfer.06.08.html) | Fully developed laminar circular-tube baseline `Nu=3.66`. | Its extension to homogeneous-mixture condensation is an explicit project assumption. |

The friction implementation uses Churchill (1977), *Friction-factor equation
spans all fluid-flow regimes*, Chemical Engineering 84(24), 91–92, with the Darcy
factor convention and the `64/Re` laminar limit. The
[two-phase document](two-phase-model.md) gives the full closure and reference list.

The CuNi capillary dimensions and measured impedance were supplied as a
literature-based example. Their calibration temperature and uncertainty are not
specified. Treat these as editable example inputs, not manufacturer-certified
performance of the user's hardware.

The optional legacy gas-conductivity interpolation preserves the earlier
engineering anchor values in `models/comsol_pipe/legacy_properties.py`. It has no
primary experimental citation in the original implementation and is not the
default NIST conductivity model.
