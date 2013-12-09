# -*- coding: utf-8 -*-
"""
Created on Mon Dec 09 16:24:53 2013

@author: Thomas Stolz
"""

from digilock import digilock
from WLM import WavelengthMeter
from findWL import findWavelength
from spectrum import Spectrum, DFSpec, DBSpec
from SignalBox import SignalBox

global digilock_ip

digilock_ip='localhost'
box_port='COM4'

class DDAVLL_Controller:
    def __init__(self):
        self.digilock=digilock(digilock_ip,60001)
        self.WLM=WavelengthMeter()
        self.box=SignalBox(box_port)
        
    def configureHardware(self):
        self.WLM.setContinuous()
        self.WLM.setExpMode()
        self.WLM.startMeasurement()
        
        #other configurations required for the lock
    
    def adjustSignal(self):
        #firstly the signal level has to be adjusted
        #we need just the reference signal for finding the broad peaks
        self.box.setPreGain(1)
        #... we need a function that detects level and range of the
        # signal and changes offset and gain until it is within the
        # accessible range of ~ (-1.5, +1.5) V
        
        #switch off the scan, otherwise findWavlength will fail
        self.digilock.switchscan(False)
        #got to 1014.91 nm which is approx. the center of the transition
        for i in xrange(3):
            if findWavelength(self.digilock,self.WLM,1014.91):
                continue
            else:
                if i==2:
                    return 'Error: frequency could not be adjusted'
        
        #turn on the scan        
        self.digilock.setscanfrequency(1)
        self.digilock.setscanrange(15)
        self.digilock.setscopetimescale(500)
        self.digilock.switchscan(True)
        
        #now get the spectrum
        for i in xrange(3):
            spec=DBSpec(self.digilock.getscopedata())
            spec.guessParams()
            if spec.isValid():
                continue
            else:
                if i==2:
                    return 'Error: Doppler broad spectrum invalid'
        
        #zoom in on the hg199 peak
        peak=spec.getHg199Peak()
        # peak is a list of parameters, [2] is the FWHM, [1] the center
        digilock.setscanrange(3*peak[2])
        digilock.setoffset(peak[1])
        
        #get spectrum again and see if it is ok
        
        #zoom in where the hg199 and hg201 peaks should be
        
        #adjust the pregain until the signal is within the range
        
        #then finetune the parameters pregain and offset until the
        #curvature becomes small enough 
        #-> try with guessed parameters first, otherwise time consuming
        
        #characterise the spectrum at last and center it
        #save the spectrum
        
        
        
        
        