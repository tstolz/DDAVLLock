# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 09:04:31 2013

@author: tom
"""

import digilock
import WLM
import time

global threshold
global max_offset
global test_offset
max_offset=30
test_offset=5
threshold=0.002

i=0
while True:
    
    data=D.getscopedata()
    np.savetxt('05_12_time%d.txt' % (i),data)
    i+=1
    time.sleep(1)

def findWavelength(digiLK, wavelengthMT, wavelength, factor=0):
    '''factor is the proportionality wavelength/offset'''
    #delta is the distance to the desired wavelength
    delta=wavelength-wavelengthMT.getWL()
    if abs(delta)<=threshold:
        print 'wavelength successfully set'
        return 1
    print '%f nm to go' % delta        
    offset=digiLK.getoffset()
    #calculate the needed offset where delta becomes zero
    #assuming linear dependence
    if not factor:
        #make sure we don't leave the range while testing
        if offset<=max_offset-test_offset:
            test=test_offset
        else:
            test=-test_offset
        offset+=test
        digiLK.setoffset(offset)
        time.sleep(0.5)
        new_delta=wavelength-wavelengthMT.getWL()
        factor=-(new_delta-delta)/test
        print 'calculated a factor of %f' % factor
        if factor==0:
            print 'hardware error (zero factor)'
            return 0
        delta=new_delta
    #print delta
    offset+=delta/factor
    #check if voltage out of range
    if abs(offset)>max_offset:
        print 'wavelength out of range'
        return 0
    #set the voltage and check the result
    digiLK.setoffset(offset)
    time.sleep(0.5)
    tmp=wavelengthMT.getWL()
    print 'set offset %f got wavelength %f' % (offset, tmp)
    new_delta=wavelength-tmp
    #print '%f nm to go' % new_delta
    if abs(new_delta)>=abs(delta):
        print 'failed (maybe threshold too low)'
        return 0
    #if the result is not good enough do the procedure again (recursive)
    if abs(new_delta)>threshold:
        #factor=factor*(new_delta-delta)/delta
        #print 'calculated a factor of %f' % factor                
        findWavelength(digiLK, wavelengthMT, wavelength, factor)
    else:
        print 'wavelength successfully set'        
        return 1
        
class DigiDummy:
    def __init__(self, laser):
        self.offset=0
        self.laser=laser
    def getoffset(self):
        return self.offset
    def setoffset(self,offset):
        self.laser.setWL(1014*(1+0.001*offset**3))
        self.offset=offset

class WLMDummy:
    def __init__(self):
        self.WL=280
    def setWL(self,WL):
        self.WL=WL
    def getWL(self):
        return self.WL

D=digilock.digilock('localhost',60001)
W=WLM.WavelengthMeter()
W.setContinuous()
W.setExpMode(True)
WD=WLMDummy()
DD=DigiDummy(WD)