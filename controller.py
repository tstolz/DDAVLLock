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
import numpy as np
import matplotlib.pyplot as plt
import time

global digilock_ip
global box_port
global maxvoltage
global secs
global U1sig
global U1flat

#time for a parameter change to affect the signal
secs=3
#maximum scan voltage for the digilock (absolute)
maxvoltage=30

U1sig=200
U1flat=66

digilock_ip='localhost'
box_port='COM4'

class DDAVLL_Controller:
    def __init__(self):
        self.digilock=digilock(digilock_ip,60001)
        self.WLM=WavelengthMeter()
        self.box=SignalBox(box_port)
        #idea: set different values for level and offset
        #store the behaviour of the signal
        #calculate the appropriate settings for the tasks
        #broadSpec and fineSpec        
        self.U1Gain_flat=U1flat
        self.U1Gain_sig=U1sig
        
    def configureHardware(self):
        self.WLM.setContinuous()
        self.WLM.setExpMode()
        self.WLM.startMeasurement()
        
        #other configurations required for the lock
        
    def recordNoise(self, num):
        start=time.time()
        noise=np.array([[],[],[],[]])
        for i in xrange(num):
            t=time.time()
            data=np.asarray(self.digilock.getscopedata())
            #shift the timescale so the values are with respect to 
            #the start of the measurement
            data[0]=data[0]+(t-start)
            noise=np.concatenate([noise,data],1)
        return noise
        
    def recordSpecs(self, num, name):
        start=time.time()
        for i in xrange(num):
            t=time.time()
            data=np.asarray(self.digilock.getscopedata())
            #shift the timescale so the values are with respect to 
            #the start of the measurement
            data[0]=data[0]+(t-start)
            np.savetxt('%s%d.txt' % (name,i), data)
            
    def setWavelength(self,wl):
        #switch off the scan, otherwise findWavlength will fail
        self.digilock.switchscan(False)
        #go to 1014.916 nm which is approx. the center of the 253.7 nm transition
        for i in xrange(3):
            if findWavelength(self.digilock,self.WLM,wl,threshold=0.001):
                return 1
            else:
                if i==2:
                    return 0

    def calibrateBoxParams(self):
        #go to an arbitrary point in the spectrum (at least one doppler-broadened
        #peak should be present) 
        self.digilock.setscanrange(5)
        self.digilock.setscanfrequency(2.5)
        self.digilock.setscopetimescale('200 ms')
        self.digilock.switchscan(True)
        self.box.setOutGain(0)
        time.sleep(secs)
        ranges=[]
        mins=[]
        gains=[]
        i=0
        # change U1Gain and save a quality parameter, that indicates the "flatness"
        # of the signal. the gain value with the flattest signal will be 
        # remembered. for obtainin a nice doppler broadened signal U2 should
        # dominate, so the signal baseline should be as negative as possible
        # while being well away from the zener-clipping region (~ -1.5 V)
        while True:
            self.box.setU1Gain(i)
            time.sleep(1)
            sig=self.digilock.getscopedata()
            sig=np.asarray([sig[3][50:],sig[1][50:]])
            if sig[1].max()>1.1:
                i+=10
                continue
            elif sig[1].min()<-1.1:
                break
            rng=abs(sig[1].mean()-np.median(sig[1]))
            ranges.append(rng)
            mins.append(sig[1].min())
            gains.append(i)
            if rng<0.01:
                i+=1
            elif rng<0.015:
                i+=5
            elif rng>0.05 and i > 120:
                i+=30
            else:
                i+=10
            if i>255:
                break
        ranges=np.array(ranges)
        mins=np.array(mins)
        self.U1Gain_sig=gains[mins.argmin()]
        print 'set signal gain: %f' % (self.U1Gain_sig)
        self.U1Gain_flat=gains[ranges.argmin()]
        print 'set flat gain: %f' % (self.U1Gain_flat)
            
    def adjustSignal(self):
        if not self.setWavelength(1014.916):
            return 'Error: Wavelength could not be adjusted'
        if abs(self.digilock.getoffset())>maxvoltage-10:
            return 'Error: Voltage out of range. Adjust laser wavelength!'
        
        #if no calibration has been done yet, do it now!
        if not self.U1Gain_flat or not self.U1Gain_sig:
            self.calibrateBoxParams()
        
        #turn on the scan - for the full spectrum only 0.5 Hz scan frequency
        #yield a good signal        
        self.digilock.setscanfrequency(0.5)
        self.digilock.setscanrange(20)
        self.digilock.setscopetimescale('1 s')
        self.digilock.switchscan(True)        
        
        #to find the broad spectrum, it is easiest to have no Gain
        self.box.setOutGain(0)
        #set the right U1Gain for the broad signal
        self.box.setU1Gain(self.U1Gain_sig)
        time.sleep(2*secs)
        
        #now we should be able to detect the peaks
        #signal[3] contains the scanning voltage
        #cut away first 50 points because of nasty edge effects
        signal=self.digilock.getscopedata()
        spec=DBSpec([signal[3][50:],signal[1][50:]])
        spec.guessParams()
        if not spec.isValid():
            return 'Error: invalid broad signal'
        #zoom in to 5V amplitude
        self.digilock.setscanrange(5)
        self.digilock.setoffset(spec.getHg199Peak())
        self.digilock.setscanfrequency(2.5)
        self.digilock.setscopetimescale('200 ms')
        time.sleep(secs)
        signal=self.digilock.getscopedata()
        #clip the signal because of nasty edge effect
        spec=Spectrum([signal[3][50:],signal[1][50:]])
        spec.guessParams(numPeaks=1)
        #if there is more than one peak something is wrong
        if not len(spec.peaks)==1:
            return 'Error: invalid signal at 5 V amplitude'
        #zoom in to 3V amplitude
        self.digilock.setoffset(spec.peaks[0][1])
        self.digilock.setscanrange(3)
        #check signal again
        signal=self.digilock.getscopedata()
        spec=Spectrum([signal[3][50:],signal[1][50:]])
        spec.guessParams(numPeaks=1)
        #if there is more than one peak something is wrong
        if not len(spec.peaks)==1:
            return 'Error: invalid signal at 3 V amplitude'
        #re-center the signal and zoom to the fine spectrum
        #position value based on experience (spec.peaks[0][1] is the center
        #and spec.peaks[0][2] is the width of the gaussian distribution)
        self.digilock.setoffset(spec.peaks[0][1]-spec.peaks[0][2]*0.15)
        self.digilock.setscanrange(0.5)
        
        #flatten the signal 
        self.box.setU1Gain(self.U1Gain_flat)
        
        #turn up the overall gain!
        self.box.setOutGain(255)
        time.sleep(secs)
        #check if we have a signal
        signal=self.digilock.getscopedata()
        spec=DFSpec([signal[3],signal[1]])
        #be prepared for narrow and thin peaks
        spec.guessParams(numPeaks=6)
        if not spec.isValid():
            return 'Error: no fine spectrum obtained'
        #zoom in further
        left=spec.peaks[0][1]
        right=spec.peaks[5][1]
        center=(left+right)/2
        rng=abs(left-right)*2
        self.digilock.setoffset(center)
        self.digilock.setscanrange(rng)
        
        #characterise the spectrum at last and center it
        #save the spectrum
        
        
        
        
        