DDAVLLock Code HowTo

***questions or bugs? -> mailto thomas.stolz@tum.de***

1. Analyzing spectra with “spectrum.py”


>> from spectrum import *


This file contains the class “SpecFitter", which analyzes transmission spectra, e.g. as recorded with digilock, 
and provides different options to fit them.


>> s = SpecFitter(data, numpeaks=6, FModel=Hg199_204, BModel=Linear)


'data' must be a two dimensional array with x and y data, 'numpeaks' gives the number of peaks to detect (if this is 
None the peak detection algorithm uses default preferences and all detected peaks are used), FModel & BModel are 
the fitmodels for the foreground and background. Fitmodels are defined in the file 'fitmodels.py' and have to inherit 
from the 'FitModel' superclass. 

When a SpecFitter instance is created, a number of internal functions are executed that detect peaks in the spectrum, 
estimate the noiselevel and instantiate the fitmodels. 

If the peak detection was not successful (e.g if less than 'numpeaks' were found), the variable 'self.alive' is set 
to 'False'. This should be checked by calling 'self.isAlive' before trying to perform a fit.
If the peak detection was successful, the height, position, and width are stored in 'self.peaks'. It is noted that,
if more than numpeaks were found, the smallest peaks are deleted to match the specification.


>> s.fit()


The spectrum is fitted with ROOT, the fit result is processed by the FitModel, and the residuals are calculated. The
fit parameters and their errors are stored in the FitModel and can be accessed via s.FModel.params and s.FModel.errors.
If the FitModel is implemented correctly, s.peaks should also be updated to hold the values obtained from the fit. 
The available analysis functionality depends on the FitModel. A list of some basic properties provided by the superclass
"FitModel" is given as follows:
Attributes:
        self.params
        self.errors
        self.parnames
        self.ChiSquareRed
        self.MaxCorrelation
        self.CovMatrixPosDef
        self.Status
        self.pValue
Methods:
        getConfidenceInterval(self,x, cl = 0.68)
        getString(self)
        slope(self, x)
        evaluate(self, x)
        getString(self)


>> s.draw()


Plots the data together with the theory function, background, residuals, and a normal probability plot of the residuals.
If the residuals are normally distributed the latter should be a straight line.



2. Controlling the DDAVLL system with “controller.py”:


>> d = DDAVLL_Controller()


This object is used to coordinate the various hardware and software components of the lock. Before calling this function, make sure the global variables digilock_ip and box_port are set to the right destination.


A DDAVLL_Controller() controls a Toptica Digilock, a High Finesse Wavelength Meter and one “SignalBox” which is used to adjust the signal level of the probe and reference beam.


>> d.configureHardware()


Call this at the beginning to make sure all components are properly configured.


>> d.setWavelength(1014.916)


Use the wavelenght meter to move to a certain wavelength (infrared, in nm). 


>> d.calibrateBoxParams()


Find out the right values for the digipotis used in the two stages of the lock. Sets the variables “d.U1sig” (for the broad scan) and “d.U1flat” (for the lock signal). This can take some minutes.


>> d.recordNoise(100)


Concatenates 100 datasets from the digilock and shifts the time information with help of a time stamp. The return value is an array just like from “digilock.getScopeData()” - but 100 times as long. This can be used to measure the noise level.


>> d.recordSpecs(100,’time’)


Records 100 datasets from the digilock and saves them in the current working directory under the name “timeX.txt” where “X” is the index of the measurement. The files can be opened with “numpy.loadtxt(‘timeX.txt’)”. 


>> d.adjustSignal()


Uses wavelength meter, the doppler broadened spectrum and the doppler free spectrum to zoom in on the lock signal. Before calling this you should make sure that the doppler free signal is strong and that the SignalBox parameters are calibrated.


>> d.lockDL(-8)


Lock to -8 MHz from resonance. You can also set the PID parameters in the function call. 


3. Analyzing recorded data with “analize.py”:


This is a flexible tool to evaluate a larger number of spectra, e. g. recorded over time or with a changing parameter. As this file is not object oriented you have to make configurations before running it. In principle you have the following options:


global prop - The properties of the spectra to evaluate. Must be a list with methods from the “DFSpec” class returning a [value, error] pair.


global fitstyle - Set the style to use for the fit. Options in “spectrum.py”


global noise - Set a noiselevel here or, if noise=0, take the fitted noise from the first spectrum.


Of course modifications in the file are possible. It provides three functions


def analyzation(dateien):
…


Takes a list of files and analyzes them subsequently. They are fitted, the elements of “prop” are evaluated and the results stored to several lists, which are returned.


def saveLog(...):


Saves the results returned by “analyzation” in a .log file.


def loadLog(...):


Can be used to load a .log file created with “saveLog”.


The three functions are wrapped within a programm, that first provides a graphical interface to select files. If a .log file is selected, it is loaded and its content is plotted. If .txt or similar files containing spectral data are provided, they are analyzed and the results are saved in a log file. Lastly the results are plotted with matplotlib.
