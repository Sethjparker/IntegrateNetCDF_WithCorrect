# IntegrateNetCDF_WithCorrect
How to integrate GC-MS data and correct for natural isotope abundance

Requirements:
- GC-MS data exported as .CDF format.
- IntegrateNetCDF_WithCorrect.m
- Library_template.m
- Matlab (R2016 or higher)

The "IntegrateNetCDF_WithCorrect" script accounts for the natural abundance of isotopes for correction of stable-isotope tracing experiments. The current script is configured to account for incorporation of 13C atoms. If a different tracer is used (e.g., 2H, 15N), then the script must be adjusted by navigating to "function [x,M] = mscorrect(x,y)" in the *.m file and adjusting the z-value to either  "1" (for 2H) or "7" (for 15N).

To integrate GC-MS data, call the script and enter three commands:

>>IntegrateNetCDF_WithCorrect('A',B,C)

A: Mapped directory location where the data is stored (e.g., 'R:\GCMS Analysis\Sample Data').

B: Library consisting of individual metabolite fragments, retention times, and chemical formulae. See 'Library_template.m' for an example of how to format this.

C: 0 (for no isotope correction) or 1 (for natural isotope correction).

>>IntegrateNetCDF_WithCorrect('R:\GCMS Analysis\Sample Data', Library_template, 1)



If the script is used for publication, we politely request that the following statement and citation be included in the methods:

"Mass isotopologue distributions of selected metabolite ion fragments were quantified and corrected for natural isotope abundance using algorithms adapted from (Fernandez et al. 1996). The code is available on https://github.com/Sethjparker/IntegrateNetCDF_WithCorrect (Accessed on xxx) under MIT license.".

Fernandez, C.A., Des Rosiers, C., Previs, S.F., David, F. & Brunengraber, H. Correction of 13C mass isotopomer distributions for natural stable isotope abundance. J Mass Spectrom 31, 255-62 (1996).
