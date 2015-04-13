# -*- coding: utf-8 -*-
"""
Created on Mon Oct 27 09:31:38 2014

@author: Thomas Stolz
"""

import numpy as np
from ROOT import TGraph, TGraphErrors, TCanvas

#--------------------UTILITY FUNCTIONS--------------------

def enhancedDerivative(data,win_len=31):
    pts=len(data[0])
    x=np.array(data[0])
    y=np.array(data[1])
    d=[0]    
    #d is a discrete derivative of data
    for i in xrange(pts-1):
        if x[i+1]==x[i]:
            d.append(0)
            continue
        d.append((y[i+1]-y[i])/(x[i+1]-x[i]))
    #create a spectrum containing the sum of values 
    #within a window around the center point    
    half_width=int(win_len/2)
    enhanced=[0]*half_width
    for i in xrange(half_width,pts-half_width):
        window=np.array(d[i-half_width:i+half_width+1])
        enhanced.append(window.sum())
    enhanced+=[0]*half_width
    return np.array(enhanced)    
    
def removePeaks(data):
    pts=len(data[0])
    #need x and y as lists in order to use the .remove method
    x=list(data[0])
    y=list(data[1])
    #peaks will be recognized in the smoothed derivative spektrum
    enhanced=enhancedDerivative(data)
    #if the slope is more than 10 % of the maximum, this is a peak
    thresh=np.abs(enhanced).max()*0.1
    #after a peak we want to cut away a litle bit more, so we need to remember
    overcount=0
    #go through the spektrum in reverse order so that order is not messed up
    for i in xrange(pts-1,0,-1):
        #cut away point if the absolute derivative is too big
        #add 3 to overcount if positive
        if enhanced[i]>thresh:
            del(x[i])
            del(y[i])
            overcount+=3
        #subtract 3 from overcount if negative
        #after a (symmetric) peak overcount will be nearly balanced
        elif enhanced[i]<-thresh:
            del(x[i])
            del(y[i])
            overcount-=3
        #make sure the first 5% of the signal are kept
        # - if there is a peak the signal is not useful anyway
        elif i/float(pts) < 0.05:
            break
        #count down overcount after a peak depending on the asymmetry
        elif overcount:
            del(x[i])
            del(y[i])
            overcount=np.sign(overcount)*(abs(overcount)-1)
    y=np.array(y)
    #calculate median and median of deviations and then reject outliers
    md=np.median(y)
    res=np.abs(y-md)
    mdev=np.median(res)
    return np.array([x,y])[:,res<5*mdev]

def findPeaks(data,show=False,win_len=31):
    '''Returns [height, center, width] of all peaks by looking for 
    zeros of the first derivative. Expects a signal without background, i.e. 
    level at zero and no drifts.'''
    x=np.array(data[0])
    y=np.array(data[1])
    pts=len(data[0])
    enhanced=enhancedDerivative(data,win_len)
    #go through the spectrum from left to right
    #checking for values above and below 
    #a threshold (thresh) until there is a jump in the sign of these values.
    #find the zero points between the jumps, these are our peak-center-candidates
    thresh=np.abs(enhanced).max()*0.1 #threshold for peak recognition (derivative data)
    Ythresh=np.abs(y).max()*0.1 # threshold for peak recognition (signal)
    lastcnt=0 # 
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
        #if loop did not break sgn is different from sgntemp
        #---> jump in sign of derivative
        #---> there was a peak candidate
        #determine the center of the peak (maximum in the data)
        c=lastcnt+np.abs(y[lastcnt:i]).argmax()
        #check if the peak has the right curvature considering its sign 
        #and if it is not too flat
        if abs(y[c]) < Ythresh:
            lastcnt=i
            sgn=sgntemp
            continue
        elif not np.sign(y[c]) == sgn:
            lastcnt=i
            sgn=sgntemp            
            continue
        #finally determine the width
        left=0 #hold index of left half maximum
        right=0 #hold index of right half maximum
        j=1 #step away from maximum
        while c-j>0 and c+j<pts:
            if abs(y[c-j])<abs(y[c])/2 and not left:
                left=c-j
            if abs(y[c+j])<abs(y[c])/2 and not right:
                right=c+j
            if left and right:
                break
            j+=1
        width=x[right]-x[left] 
        #append the peak to results
        PEAKS.append([y[c],x[c],width])            
        lastcnt=i
        sgn=sgntemp      
    if len(PEAKS):
        PEAKS=np.array(PEAKS)
    else:
        print "WARNING: NO PEAKS FOUND!"
        PEAKS=np.array([[0,0,1]])
    if show:
        t1=create_TGraph([x,enhanced],0,0)
        t2=create_TGraph([x,y],0,0)
        t3=create_TGraph([PEAKS[:,1],PEAKS[:,0]],0,0)
        t3.SetMarkerSize(1)
        t3.SetMarkerColor(2)
        canvas=TCanvas('Peakfinder','Peakfinder')
        canvas.Divide(2,1)
        canvas.cd(1)
        t1.Draw("ap")
        canvas.cd(2)
        t2.Draw("ap")
        t3.Draw("psame")
        canvas.Show()
        raw_input("press any key")
        for obj in (t1,t2,t3,canvas):
            if obj:
                obj.IsA().Destructor(obj)
    return PEAKS

def create_TGraph(data,xerr,yerr):
    n=len(data[0])    
    #create float64 arrays for ROOT
    x=np.array(data[0],dtype=np.float64)
    y=np.array(data[1],dtype=np.float64)
    #if there are errors, have to create a TGraphErrors
    if yerr==0:
        T2gr=TGraph(n,x,y)
    else:
        xnoise=np.asarray([xerr]*n,dtype=np.float64)
        ynoise=np.asarray([yerr]*n,dtype=np.float64)
        T2gr=TGraphErrors(n,x,y,xnoise,ynoise)
    #appearance settings
    T2gr.SetLineColor( 4 )
    T2gr.SetLineWidth( 1 )
    T2gr.SetMarkerColor( 4 )
    T2gr.SetMarkerStyle( 21 )
    T2gr.SetMarkerSize(0.3)
    xaxis=T2gr.GetXaxis()
    xaxis.SetTitle( 'scan voltage [V]' )
    yaxis=T2gr.GetYaxis()
    yaxis.SetTitle( 'signal voltage [V]' )
    for ax in (xaxis, yaxis):
            ax.SetTitleSize(0.03)
            ax.SetLabelSize(0.025)
    return T2gr
