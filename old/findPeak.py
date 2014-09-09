# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 14:57:32 2013

@author: universe
"""

import numpy as np
import matplotlib.pyplot as plt


def findPeaks(data,level,drift,win_len=31,show=False):
    ''' will return position of the peaks (zeros of first derivative) and their height and width.
    level and drift of the signal help identify the peaks.'''
    x=np.array(data[0])
    y=np.array(data[1])-(level+drift*x)
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
    #checking for values above and below 
    #a threshold (thresh) until there is a jump in the sign of these values.
    #find the zero points between the jumps, these are our peak-center-candidates
    thresh=enhanced.max()*0.1
    Ythresh=np.abs(x).max()*0.2
    lastcnt=0
    sgn=0
    PEAKS=[]
    for i in xrange(pts):
        val=enhanced[i]
        if abs(val)<thresh:
            continue
        sgntemp=np.sign(val)
        if sgn==sgntemp or sgntemp==0:
            lastcnt=i
            continue
        if not lastcnt:
            lastcnt=i
            sgn=sgntemp
            continue
        print sgntemp
        
        #determine the center of the peak
        c=lastcnt+np.abs(y[lastcnt:i]).argmax()
        #check if the peak has the right curvature for its sign and is not to flat
        if y[c] < Ythresh:
            continue
        elif not np.sign(y[c]) == sgn:
            continue
        
        #finally determine the width
        left=0
        right=0
        j=1
        while c-j>0 and c+j<pts:
            if abs(y[c-j])<abs(y[c])/2 and not left:
                left=c-j
            if abs(y[c+j])<abs(y[c])/2 and not right:
                right=c+j
            if left and right:
                break
            j+=1
        width=x[right]-x[left]
        
        PEAKS.append([y[c],x[c],width])        
        
        lastcnt=i
        sgn=sgntemp
        

                
    PEAKS=np.array(PEAKS)
    if show:
        plt.subplot(121)
        plt.plot(x,enhanced,'o')
        plt.subplot(122)
        plt.plot(x,y,PEAKS[:,0],PEAKS[:,1],'or')
        plt.show()
    return PEAKS

data=np.load('testspec.npy')
findPeaks(data,show=True)
