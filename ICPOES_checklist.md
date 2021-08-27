# How to run the Thermo ICP-OES

1. Check the air extraction valve (ICP-MS off, ICP-OES on)
1. Start Qtegra.
1. Check that "Purge Gas Flow" is set to "Trickle".
1. Check the rest of indicators: all should be green excepting "Detector water flow", "Detector temperature" and "Plasma".
1. Check that the pressure of the gas PQII+ at the wall indicator is >80 (5.5).
1. Start chiller and wait until all the indicators except "Plasma" are on.
1. Set up tubing and check they are in the right direction.
1. Press button "GET READY" in Qtegra (top of the screen).
1. Press OK.
1. The plasma should switch on now...
1. Prepare samples (check positions in Qtegra>any lab book> ASX-S60).
1. Go to "Lab Books">"Create new".
1. Choose a name (usually AR+date+description).
1. From existing labbook: choose "MMR 26 09 2018 Standards" or other newer.
1. Go to ICAPOES: check peaks.
1. Go to Acquisition parameters: check that all are "Radial".
1. Go to Standards: Select or create parent material and check or input the right concentrations.
1. Go to Sample list: input and check rack, position, repeats, etc.
1. Press the green triangle (Play) button to start.
1. See Evaluation Results and wait until all the samples are measured.
1. Home Page>Dashboard>ASX-S60: go back to H2O for 5 minutes.
1. Prob up.
1. ICAPOES>GET READY>Shutdown.
1. Release tubing.
1. Wait 5 minutes while removing the sample tubes.
1. Switch off water chiller.
1. Export data as csv for MATLAB analysis. Make sure that you have the following headers:  Raw.Average, Raw.STD & ExtCal.StandardConcentration
1. Write any incident in the ICP book.

*Ángel Rodés, 2020*\
[**angelrodes.com**](https://angelrodes.wordpress.com/)
