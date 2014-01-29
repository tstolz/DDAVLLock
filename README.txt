DDAVLLock Code HowTo


1. Analyzing spectra with “spectrum.py”


>> from spectrum import *


This file contains the class “Spectrum” with subclasses “DFSpec” and “DBSpec”. These objects can contain transmission spectra, e.g. as recorded with digilock, and provide methods to analize them. 


>> s = Spectrum(data)


data must be a two dimensional array with x and y data. 


>> s.setStyle(‘voigt’)


Choose a style for the peaks in the spectrum. Standard types are “lorentz”, “gauss” and “voigt”. If “voigt” is chosen, there is a parameter s.vratio describing the ratio of gaussian to lorentzian width of the voigt profile. DFSpec also defines advanced types “simplelorentz”, “simplevoigt” (all peaks have the same width) and “fixedvoigt” (the gauss-lorentz-ratio is fixed). 


>> s.guessParams()


This command detects the peaks in the spectrum and determines their height, width and position. These values are converted to parameters of the theory function according to the style.


>> s.draw()


Plots the data together with the theory function and the noise underground. 


>> s.fit()


The spectrum is fitted with ROOT. Calls guessParams() before to obtain estimates for the fit.


>> s.noisefit1()


Subtract the theory curve from the datapoints and determines the RMS noise. Assigns the value to s.noise_level.


The fitted parameters & errors are contained within s.params and s.errors. Note that the errors have no real meaning unless the variable s.noise_level has been set to a reasonable value.


Order of s.params:
level, drift, curvature, height1, center1, width1, … , heightN, centerN, widhtN


level, drift and curvature define the signal underground (a quadratic fit). 


Specific features of DFSpec:


>> s.setAlpha(0.01)


Set an estimate for the angle between the beams (in rad). For the conversion to “s.vratio”, the global variable “wRatio” is used. 


>> s.lockPointY(freq=-8.9)
>> [0.123456789, 0.02]


Determine the lockpoint on the y-axis to lock to a frequency freq (in MHz). Returns an array with the lockpoint and its error estimate.


>> s.isValid()
>> True


Check if there are 6 peaks, the reduced chi-squared is ok and so on.


>> s.draw()


The x-axis is automatically translated into MHz.


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


The three functions are wrapped within a programm, that first provides a graphical interface to select files. If a .log file is selected, it is loaded and it content is plotted. If .txt or similar files containing spectral data are provided, they are analyzed and the results are saved in a log file. Lastly the results are plotted with matplotlib.