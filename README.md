# MATLAB: ICPOES datared

ICPOES data reduction from Qtegra raw data.

*Ángel Rodés*, SUERC (2021)\
[angelrodes.com](https://angelrodes.wordpress.com/)

## Requirements

These files must be in the same folder:

* ICP data in ```.csv```
* ```ICPOESdatared_v2_4.m```
* ```linear_regression_chisq_fn.m```

## How to make it work with Qtegra

1. After running the [Thermo ICP-OES](https://github.com/angelrodes/ICPOES-datared/blob/main/ICPOES_checklist.md), export data including ```Raw.Average```, ```Raw.STD``` and ```ExtCal.StandardConcentration``` in a ```.csv``` file with semicolon (```;```) as separator:

``` csv
Raw.Average;(...);Raw.STD;(...);ExtCal.StandardConcentration
```

See example file: ```MMR 10 2 20 Al CONC_20200224.csv```

You can export other data too. They will be ignored.

2. Fill the missing data in ```ExtCal.StandardConcentration``` columns: ```0``` for blanks, and the corresponding values for the standards not used by Qtegra in the calibration (the standards in between samples or at the end of the run).

3. Run ```ICPOESdatared_v2_4.m``` in MATLAB/Octave and select your file. A lot of graphical output will appear, and a text file (```.txt```) containing a few tables will be created (concentrations, calibrations, LOD/LOQ; by analyte, by element, etc.).

![Graphical output](https://user-images.githubusercontent.com/53089531/124753740-45a09a00-df21-11eb-9fcc-508e1f4e7713.jpg)

## Changelog

v2.4 (2021.09.08)\
If `Raw.STD` is not included, `Raw.RSD` is used instead (`Raw.STD=Raw.Average*Raw.RSD/100`)
