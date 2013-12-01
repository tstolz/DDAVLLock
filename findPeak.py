# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 14:57:32 2013

@author: universe
"""

import numpy as np
import matplotlib.pyplot as plt

def findPeaks(data,win_len=31,show=False):
    ''' will return position of the peaks (moving window method) - 
    peaks are actually dips'''
    x=np.array(data[0])
    y=np.array(data[1])
    pts=len(data[0])
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
    if show:
        plt.subplot(121)
        plt.plot(x,enhanced)
        plt.subplot(122)
        plt.plot(x,y,PEAKS[:,0],PEAKS[:,1],'or')
        plt.show()
    return PEAKS

