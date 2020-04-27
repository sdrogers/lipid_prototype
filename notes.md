## Lipid data - readme

mzMine steps

- Mass detection: exact mass, noise level = 1e5
- ADAP Chromatogram builder: had to set scan numbers (some obviously had no peaks in them, right at the end), min group = 5, group intensity thresh  1e5, min highest intensity = 1e5
- Peak detection: chromatogram deconvolution, wavelets (ADAP). S/N = 3 min feat height = 10, coeff = 10, duration range = 0.00 to 0.8, RT wavelet range = 0.00 - 0.2
- CSV export: provide full boxes as we will need to make ROIs to correlate isotopes

	