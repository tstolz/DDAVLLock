# -*- coding: utf-8 -*-
"""
Created on Wed Aug 22 12:49:03 2012

@author: bernd
"""

import digilock
#import pid
import time


def lock(num_averages=10):
    """force the voltage offset (center of allowed SCC output
    volltage window) to follow the laserdrift averaging over
    the last num_averages measurements.
    
    """
    
 #   p=pid.PID(0.001,0.0005,0)
 #   p.setPoint(0)
 #   pidv=0
    initial_offset=digilock.getoffset()
    offsets=[initial_offset]*10
    while True:
        time.sleep(0.5)
        measurement_value=digilock.getoutputvoltage()
        #pidv = p.update(measurement_value)
        #offset=initial_offset+pidv
        offsets.append(measurement_value)
        if len(offsets) > num_averages:
            del offsets[1]
        offset=mean(offsets)
        digilock.setoffset(offset)
        time.sleep(0.2)
        digilock.getoffset()
        print str(measurement_value)+"; "+str(offset)

def mean(numberList):
    """return mean of list given as argument"""
    if len(numberList) == 0:
        return float('nan')
    floatNums = [float(x) for x in numberList]
    return sum(floatNums) / len(numberList)