# Combine bathymetry with IFS lake information

`ecwam-combine_bathy_laked-dp` adjusts the depths in a GRIB bathymetry table
produced by `ecwam-preproc-dp`, using IFS land/sea mask, lake cover, and lake depth.
The single-precision executable, when enabled, is `ecwam-combine_bathy_laked-sp`.
The program takes no command-line arguments and reads these files from its
working directory:

| GRIB File         | GRIB parameter | Contents                               |
| ----------------- | -------------- | -------------------------------------- |
| `wam_grid_tables` | 140219         | WAM bathymetry; missing points are dry |
| `lsm`             | 172            | IFS land/sea mask                      |
| `clake`           | 26             | IFS lake cover                         |
| `lakedl`          | 228007         | IFS lake depth in metres               |

All inputs must have compatible grids, scanning order, and point counts.
At existing wet points, lake cover of at least 0.65 selects replacement by
the IFS lake depth. Points with `lsm <= 0.5` and `clake <= 0.01` instead receive
a weighted average of WAM and IFS depths. IFS depths are capped at 999 m;
the averaging weights depend on whether either depth is at most 3 m.
The output is a GRIB file named `new_bathy`; original inputs are not modified.
In IFS processing, this output replaces the table installed as
`wam_grid_tables_${WAMNFRE}_${WAMNANG}`. For a fixed horizontal grid and
bathymetry, those frequency/direction counts do not enter coupling-mask
extraction; they identify the spectral configuration and table filename.

**Mask preservation:** this program must change depths only, never wet/dry
classification. Wet means non-missing bathymetry, not positive depth;
a valid zero depth remains wet. The original input mask is saved in memory
and compared point by point with the values decoded from the re-encoded output
GRIB. A change in either direction aborts with the first changed point,
its original wet/dry status, output depth, and IFS input values. The abort
happens before writing output.

**Lake-depth validation:** before changing any depths, the input lake depth `lakedl`
is validated. The program rejects missing, negative, or non-finite lake depths
selected for replacement or averaging at existing wet points.
It reports the first affected point and its input values, then aborts without writing output.
Zero remains valid; invalid lake depths at dry points or points not selected for an update are
ignored. This also prevents invalid averages that would remain non-missing
and therefore escape a mask-only check.
