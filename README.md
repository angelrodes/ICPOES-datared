# MATLAB: ICPOES datared
ICPOES data reduction from Qtegra raw data.
*Ángel Rodés*, SUERC (2021)

## Requirements

these files must be in the same folder:

* ICP data in ```.csv```
* ```ICPOESdatared_v2_3.m```
* ```linear_regression_chisq_fn.m```

## How to make it work with Qtegra

1. Export data including ```Raw.Average```, ```Raw.STD``` and ```ExtCal.StandardConcentration``` in a ```.csv``` file with semicolon (```;```) as separator:

``` csv
Raw.Average;(...);Raw.STD;(...);ExtCal.StandardConcentration
```

See example file: ```MMR 10 2 20 Al CONC_20200224.csv```

You can export other data too. They will be ignored.

2. Fill the missing data in ```ExtCal.StandardConcentration``` columns: ```0``` for blanks, and the corresponding values for the standards not used by Qtegra in the calibration (the standards in between samples or at the end of the run).

3. Run ```ICPOESdatared_v2_3.m``` in MATLAB/Octave and select your file. A lot of graphical output will appear, and a text file (```.txt```) containing a few tables will be created (concentrations, calibrations, LOD/LOQ; by analyte, by element, etc.).
