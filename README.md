# Phasor Analysis GUI

## What is it? 

This GUI example is designed to perform single point time-correlated single photon counting (TCSPC) measurements and display the resulting point on the phasor plot in real-time.

The Phasor Analysis GUI offers a complete pipeline for single point fluorescence lifetime analysis, starting from the calibration on a reference fluorophore for then analyzing the fluorescence lifetime decay of unknown samples using the TCSPC histogram and the phasor plot.

The TCSPC histogram is obtained by dividing the laser period in 256 time bins  and the single fluorescence photon events falling in each bin, recorded by a detector and time-tagged with tens/hundreds of picoseconds precision by the data acquisition card, are counted and passed to a 2D histogram to reconstruct the profile of the fluorescence lifetime decay curve.

For instance, if the pulsed laser is set with a repetition frequency of 80 MHz, that corresponds to a laser period of 12.5 nanoseconds, then time bins of 0.048 nanoseconds (48 picoseconds) are created and the number of photons falling in each 48 picoseconds time bin will be passed to the histogram.

The phasor analysis is performed in real-time applying a cosine and sine transformation to the TCSPC data, in order to associate to a TCSPC decay histogram a point of coordinates  G (cosine component) and S (sine component) in the phasor space.

The GUI allows to perform both phasor analysis on a single batch of data, applying the phasor transformation to the entire data acquired, or on different batches, whose number can be obtained by dividing the total acquisition time of the experiment for a given number decided by the experimenters. In this way it is possible to get both a single point in the phasor plot at the end of the acquisition or a cloud of points corresponding to different batches of TCSPC data.

For immediate reference, the code used for designing this GUI is reported and commented in the folder [Single-point-spectroscopy-phasor-analysis](/Single-point-spectroscopy-phasor-analysis) 

This is how the GUI looks like and how it shows data during the experiments 

![Phasor analysis](/images/phasor_plot.png "Photons_tracing")

## How to get the requirements

You can install requested dependencies with the following *pip* command:

```
pip install flim-labs-api

```