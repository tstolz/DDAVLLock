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
import ROOT
import dataprocessing

class broadSpectrum:
    def __init__(self, data):
        self.data=np.array(data)
        self.peaks=None
        self.pk_centers=[0]*5
        self.x_range=data[0].max()-data[0].min()
        self.points=len(data[0])
        self.ChiSquare=0
        self.params=None
    def broadSpecTheory(self, x, params):
        level, drift, w1, w2, w3, w4, w5, d1, d2, d3, d4, d5 = params
        c1, c2, c3, c4, c5  = self.pk_centers
        result = level + x[0]*drift
        peaks = ((c1, w1, d1),(c2,w2,d2),(c3,w3,d3),(c4,w4,d4),(c5,w5,d5))
        for peak in peaks:
            result -= peak[2]*np.exp(-(x[0]-peak[0])**2 / peak[1]**2)        
        return result
    def findPeaks(self,win_len=31,show=False):
        ''' will return position of the peaks (moving window method)'''
        x=np.array(self.data[0])
        y=np.array(self.data[1])
        pts=self.points
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
        self.peaks=PEAKS
        self.pk_centers=list(PEAKS[:,0].T)
        return self.peaks
    def guessLevel(self, data):
        centiV=np.array(data[1]*100,dtype=int)
        #shift centiV to positive values so we can use bincount
        centiVmin=centiV.min()
        centiV-=centiVmin
        freq=np.bincount(centiV)
        level=(freq.argmax()+centiVmin)*0.01
        return level
    def guessDrift(self):
        left=self.guessLevel(self.data[:,0:int(self.points/2)])
        right=self.guessLevel(self.data[:,int(self.points/2)+1:])
        drift=(right-left)/self.x_range * 3
        return drift
    def fitSpectrum(self, parameter=None, showgraphs=True):
        '''data[0] must contain x-, data[1] y-axis data, parameters a guess of the 
        17 parameters in broadSpecTheory'''
        data=self.data
        if not parameter:
            parameter=[]
            peaks=self.findPeaks()
            level=self.guessLevel(data)
            drift=self.guessDrift()
            rng=self.x_range
            #level is assumed to return the central signal level
            #we need the level on the left
            parameter+=list([level-drift*rng/2,drift])
            parameter+=list([rng*0.02]*5)
            parameter+=list(level-peaks[:,1].T)
            print 'guessed parameters:'
            print parameter        
    
        ParNames=["level","drift","width1","width2","width3","width4","width5","depth1","depth2","depth3","depth4","depth5"]
        gr=dataprocessing.create_TGraph(data)
        fit1=ROOT.TF1("fit1",self.broadSpecTheory, min(data[0]), max(data[0]), 12)
        #fit1.SetParameters(0.001,0.05,0.002,0.001,0.055,0.002,0.01)
        for i in range(len(ParNames)):
            fit1.SetParName(i, ParNames[i])
        for n in range(len(parameter)):
            fit1.SetParameter(n,parameter[n])
        #fit1.SetParLimits(1,0.045,0.05)
        #fit1.SetParLimits(4,0.052,0.06)
        if showgraphs:
            fit1.Draw()
            #print parameter
            raw_input("pre fit, press any key")
            print 'fitting...'
        gr.Fit(fit1, "", "", data[0].min(), data[0].max())
        pars=fit1.GetParameters()
        params=[]
        for i in xrange(len(parameter)):
            params.append((ParNames[i],pars[i]))
        #print params
        self.ChiSquare=fit1.GetChisquare()
        params.append(('ChiSquare',self.ChiSquare))
    #    #append SNR. Mean of amplitude1 and amplitude2 is used für Signal height
    #    ParNames.append("SNR")
    #    ParNames.append("Chisquare")
        gr.Draw()
        self.params=params
        return params
    
        
data=np.load('testspec.npy')
B=broadSpectrum(data)

