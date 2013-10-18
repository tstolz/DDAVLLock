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
    
def findPeaks(data, win_len=31):
    ''' will return position of the peaks (moving window method)'''
    x=np.array(data[0])
    y=np.array(data[1])
    pts=len(y)
    profile=[0]    
    #profile is like the discrete derivative of data
    for i in xrange(pts-1):
        profile.append((y[i+1]-y[i])/(x[i+1]-x[i]))
    #create a spectrum containing the sum of values 
    #within a window around the center point    
    win_width=int(win_len/2)
    enhanced=[0]*win_width
    for i in xrange(win_width,pts-win_width):
        window=np.array(profile[i-win_width:i+win_width])
        enhanced.append(window.sum())
    enhanced+=[0]*win_width
    enhanced=np.array(enhanced)
    #go through the spectrum from left to right
    #memorizing values above (maxcnt) and below (mincnt)
    #a threshold (thresh) until len(maxcnt) approx. len(mincnt)
    #-> peak
    #counts are deleted whenever values are ~0 for too long
    thresh=enhanced.max()*0.2
    mincnt=[]
    maxcnt=[]
    nocnt=0
    lastcnt=1
    PEAKS=[]
    for i in xrange(pts):
        val=enhanced[i]
        if val<-thresh:
            nocnt=0
            if lastcnt==1:
                mincnt=[]
                maxcnt=[]
            mincnt.append(i)
            lastcnt=-1
        elif val>thresh:
            nocnt=0
            lastcnt=1
            maxcnt.append(i)
        else:
            nocnt+=1
            if nocnt>len(mincnt) and nocnt>len(maxcnt):
                mincnt=[]
                maxcnt=[]
            elif len(mincnt)==0 or len(maxcnt)==0:
                pass
            elif abs(len(mincnt)-len(maxcnt))/len(mincnt) < 0.2:
                mins=np.array(mincnt)
                maxs=np.array(maxcnt)  
                index=int((mins.mean()+maxs.mean())/2)
                PEAKS.append([x[index],y[index]])
                mincnt=[]
                maxcnt=[]
                
    PEAKS=np.array(PEAKS)
    plt.subplot(121)
    plt.plot(x,enhanced)
    plt.subplot(122)
    plt.plot(x,y,PEAKS[:,0],PEAKS[:,1],'or')
    plt.show()
    return PEAKS
    

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

