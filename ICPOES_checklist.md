# How to run the Thermo ICP-OES

* Check the air extraction valve (ICP-MS off, ICP-OES on)
* Start Qtegra.
* Check that "Purge Gas Flow" is set to "Trickle".
* Check the rest of indicators: all should be green excepting "Detector water flow", "Detector temperature" and "Plasma".
* Check that the pressure of the gas PQII+ at the wall indicator is >80 (5.5).
* Start chiller and wait until all the indicators except "Plasma" are on.
* Set up tubing and check they are in the right direction.
* Press button "GET READY" in Qtegra (top of the screen).
* Press OK.
* The plasma should switch on now...
* Prepare samples (check positions in Qtegra>any lab book> ASX-S60).
* Go to "Lab Books">"Create new".
* Choose a name (usually AR+date+description).
* From existing labbook: choose "MMR 26 09 2018 Standards" or other newer.
* Go to ICAPOES: check peaks.
* Go to Acquisition parameters: check that all are "Radial".
* Go to Standards: Select or create parent material and check or input the right concentrations.
* Go to Sample list: input and check rack, position, repeats, etc.
* Press the green triangle (Play) button to start.
* See Evaluation Results and wait until all the samples are measured.
* Home Page>Dashboard>ASX-S60: go back to H2O for 5 minutes.
* Prob up.
* ICAPOES>GET READY>Shutdown.
* Release tubing.
* Wait 5 minutes while removing the sample tubes.
* Switch off water chiller.
* Export data as csv for MATLAB analysis. Make sure that you have the following headers:  Raw.Average, Raw.STD & ExtCal.StandardConcentration
* Write any incident in the ICP book.

*Ángel Rodés, 2020*\
[**angelrodes.com**](https://angelrodes.wordpress.com/)
