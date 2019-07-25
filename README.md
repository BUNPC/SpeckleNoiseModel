# SpeckleNoiseModel
Code to model the noise of a fiber based NIRS system, as well as to plot the figures from associated papers.

All .m files need to be in the same directory for this to work.

A good starting point is example.m, it calculates the noise model for a common specific device (OPT101) and a liquid breast phantom, using a high coherence laser diode to display the increased noise due to speckle. This code is a good illustration of how the structures for the medium, source and detector need to be created for the functions to work.
The main code is noisecalcfunc.m. It is a function that performs all the calculations based on source detector separation and the source, medium and detector specified. It uses two different auxiliary functions, calcbeta1.m to calculate the speckle contrast, and getRr to calculate the expected power per unit are at a given source detector separation.


Most code by Antonio Ortega (Boston University) with auxiliary code from Xiaojun Cheng (BU) and Steven Jacques (https://omlc.org/news/apr08/skinspectra/index.html)
