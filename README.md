# ICPOES datared

MATLAB scripts for ICPOES data reduction from Qtegra raw data.

*Ángel Rodés*, SUERC (2021)\
[angelrodes.com](https://angelrodes.wordpress.com/)

## Requirements

These files must be in the same folder:

* ICP data in ```.csv```
* ```ICPOESdatared_vX.m```
* ```linear_regression_chisq_fn.m```
* ```normrnd_BoxMuller.m```

## How to make it work with Qtegra

1. After running the [Thermo ICP-OES](https://github.com/angelrodes/ICPOES-datared/blob/main/ICPOES_checklist.md), export data including ```Raw.Average```, ```Raw.STD``` and ```ExtCal.StandardConcentration``` in a ```.csv``` file with semicolon (```;```) as separator:

``` csv
Raw.Average;(...);Raw.STD;(...);ExtCal.StandardConcentration
```

See example file: ```AR20200205Aliquots.csv```

You can export other data too. They will be ignored.

2. Modify data in ```ExtCal.StandardConcentration``` columns if needed. E.g. if you want to use some standard data not stated in Qtegra.

3. Run ```ICPOESdatared_vX.m``` in MATLAB/Octave and select your file. A lot of graphical output will appear, and a text file (```.txt```) containing a few tables will be created (concentrations, calibrations, LOD/LOQ; by analyte, by element, etc.).

![Graphical output](https://user-images.githubusercontent.com/53089531/124753740-45a09a00-df21-11eb-9fcc-508e1f4e7713.jpg)

## Changelog

v2

- Some code was added to transform Qtegra format to standard ICP format (comma separated values with one analyte per line)

v2.4 (2021.09.08)

- If `Raw.STD` is not included, `Raw.RSD` is used instead (`Raw.STD=Raw.Average*Raw.RSD/100`)
- Samples named 'blank" or 'blk' are considered blanks
- Standard's nominal values ara copied to all locations if they have the same name.

v3 (2022.01.21)

- Lots of comments added.
- Two new functions defined (normrnd_BoxMuller and normpdf_local) to replace normrnd and normpdf just in case the statistic package is not available in Octave.
- Some redundant code was deleted.
- New dialog allows the user selecting standards and blanks.
- Maximum chi-squared formula was simplified to min-chi-squared+DOF in ```linear_regression_chisq_fn.m``` to avoid using functions form the statistic package. (original code is commented).
