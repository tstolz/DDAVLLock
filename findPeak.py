# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 14:57:32 2013

@author: universe
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import digilock
import WLM
import time

def broadSpecTheory(x, level, drift, c1, c2, c3, c4, c5, w1, w2, w3, w4, w5, d1, d2, d3, d4, d5):
    result = level + x*drift
    peaks = ((c1, w1, d1),(c2,w2,d2),(c3,w3,d3),(c4,w4,d4),(c5,w5,d5))
    for peak in peaks:
        result -= peak[2]*np.exp(-(x-peak[0])**2 / peak[1]**2)
        
    return result

def simpleTest(x, a1, a2, a3):
    a = (a1, a2, a3)
    result=a[0]*np.exp(-(x-a[1])**2/(a[2])**2)
    return result
    
def findPeaks(data):
    ''' will return position of the five peaks (moving average method)'''

def fitSpectrum(data):
    '''data[0] should contain x, data[1] y axis data'''
    x=np.array(data[0])
    y=np.array(data[1])
    pts=len(y)
    xrng=x[-1]
    level=y[int(0.05*pts)]
    drift=(y[int(0.95*pts)]-level)/(xrng*0.9)
    c=[]
    for i in xrange(5):
        c.append(0.1*xrng+(i+1)/5.*0.8*xrng)
    w=5*[0.04*xrng]
    d=[(y.max()-y.min())*0.75]*5
    print level, drift, c, w, d
    
    
        
data=np.load('testspec.npy')

