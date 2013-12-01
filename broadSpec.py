# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 14:57:32 2013

@author: universe
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
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
        '''params must have the form "level, drift, depth1, center1, width1, ... depthN, centerN, widthN'''
        #level, drift, w1, w2, w3, w4, w5, d1, d2, d3, d4, d5, c1, c2, c3, c4, c5 = params
        #c1, c2, c3, c4, c5  = self.pk_centers
        #params=list(params)
        level = params[0]
        drift = params[1]
        peaks = [[params[i],params[i+1],params[i+2]] for i in range(2,len(params),3)]
        result = level + x[0]*drift
        for peak in peaks:
            result += peak[0]*np.exp(-(x[0]-peak[1])**2 / (2*peak[2]**2))        
        return result
    def draw(self):
        if self.params:
            #params=np.array(np.array(self.params[:-1])[:,1],dtype=float)
            Fit=map(lambda x:self.broadSpecTheory(x,self.params),self.data.T[:,:2])
            plt.plot(self.data[0],self.data[1],'o',self.data[0],Fit,'r')
        else:
            plt.plot(self.data[0],self.data[1],'o')
        plt.show()
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
    def guessMeanLevel(self, data):
        centiV=np.array(data[1]*100,dtype=int)
        #shift centiV to positive values so we can use bincount
        centiVmin=centiV.min()
        centiV-=centiVmin
        freq=np.bincount(centiV)
        maxfreq=[]
        for i in xrange(int(0.02*len(freq))):
            argmax=freq.argmax()
            maxfreq.append(argmax)
            freq[argmax]=0
        mean=np.array(maxfreq).mean()
        level=(mean+centiVmin)*0.01
        return level
    def guessDrift(self):
        #guess level in the right and left half of the spectrum and 
        #divide by one third of the xrange
        left=self.guessMeanLevel(self.data[:,0:int(self.points/2)])
        right=self.guessMeanLevel(self.data[:,int(self.points/2)+1:])
        drift=(right-left)/self.x_range * 3
        return drift
    def guessParams(self):
        parameters=[]
        rng=self.x_range
        peaks=self.findPeaks()
        drift=self.guessDrift()
        level=self.guessMeanLevel(self.data)-drift*rng/2
        #level is assumed to return the central signal level
        #we need the level on the left
        parameters+=list([level,drift])
        for peak in peaks:
            parameters.append(-(level+drift*peak[0]-peak[1]))
            parameters.append(peak[0])
            parameters.append(rng*0.02)
#        parameters+=list([rng*0.02]*5)
#        parameters+=map(lambda i: level+drift*peaks[i,0]-peaks[i,1],xrange(len(peaks)))
#        parameters+=self.pk_centers
        return parameters
    def fitSpectrum(self, parameters=None, showgraphs=True):
        '''data[0] must contain x-, data[1] y-axis data, parameters a guess of the 
        parameters in broadSpecTheory'''
        data=np.array(self.data)
        if not parameters:
            parameters=self.guessParams()
            print 'guessed parameters:'
            print parameters        
        
        ParNames=["level","drift"]
        for i in xrange((len(parameters)-2)/3):
            ParNames.append("depth%d" % i)  
            ParNames.append("center%d" % i)
            ParNames.append("width%d" % i)
            
        print ParNames
        gr=dataprocessing.create_TGraph(data)
        fit1=ROOT.TF1("fit1","pol1+gaus(2)+gaus(5)+gaus(8)+gaus(11)+gaus(14)", min(data[0]), max(data[0]))
        for i in range(len(ParNames)):
            fit1.SetParName(i, ParNames[i])
        for n in range(len(parameters)):
            fit1.SetParameter(n,parameters[n])
        #fit1.SetParLimits(1,0.045,0.05)
        #fit1.SetParLimits(4,0.052,0.06)
        if showgraphs:
            fit1.Draw()
            #print parameter
            raw_input("pre fit, press any key")
            print 'fitting...'
        tm=time.time()
        gr.Fit(fit1, "V","", data[0].min(), data[0].max())
        print "Fitting time: %f s" % (time.time()-tm)
        self.params=[]
        pars=fit1.GetParameters()
        for i in xrange(len(parameters)):
            self.params.append(pars[i])
#        params=[]
#        for i in xrange(len(parameter)):
#            params.append((ParNames[i],pars[i]))
        self.ChiSquare=fit1.GetChisquare()
#        params.append(('ChiSquare',self.ChiSquare))
    #    #append SNR. Mean of amplitude1 and amplitude2 is used für Signal height
    #    ParNames.append("SNR")
    #    ParNames.append("Chisquare")
        #gr.Draw()
#        self.params=params
        if showgraphs:        
            self.draw()
        return True
    
        
data=np.load('testspec.npy')
B=broadSpectrum(data)

