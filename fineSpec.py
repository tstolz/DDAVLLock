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
import findPeak
import os

class fineSpec:
    def __init__(self, data):
        self.data=np.array(data)
        self.peaks=None
        self.pk_centers=[0]*6
        self.x_range=data[0].max()-data[0].min()
        self.points=len(data[0])
        self.ChiSquare=0
        self.params=None
        self.depth=[0]*6
    def draw(self):
        params=np.array(np.array(self.params[:-1])[:,1],dtype=float)
        Fit=map(lambda x:self.fineSpecTheory(x,params),self.data.T[:,:2])
        plt.plot(self.data[0],self.data[1],'o',self.data[0],Fit,'r')
        plt.show()
    def fineSpecTheory(self, x, params):
        level, drift, w1, w2, w3, w4, w5, w6= params
        c1, c2, c3, c4, c5, c6  = self.pk_centers
        d1, d2, d3, d4, d5, d6 = self.depth
        result = level + x[0]*drift
        peaks = ((c1, w1, d1),(c2,w2,d2),(c3,w3,d3),(c4,w4,d4),(c5,w5,d5),(c6,w6,d6))
        #in the fine spec the peaks are lorentzian         
        for peak in peaks:
            result -= peak[2]*(peak[1]/2)**2/((x[0]-peak[0])**2 + (peak[1]/2)**2)        
        return result
    def setPeaks(self,win_len=31,show=False):
        data=[self.data[0],-self.data[1]]
        multi=findPeak.findPeaks(data,win_len,show)
        multi=self.refinePeaks(multi,data,show)
        single=findPeak.findPeaks(self.data,win_len,show)
        self.refinePeaks(single,self.data,show)
        tmp=np.array(list([multi[0]])+ list([[single[0][0],-single[0][1]]])+ list(multi[1:]))        
        self.peaks=np.array([tmp.T[0],-tmp.T[1]]).T  
        self.pk_centers=list(self.peaks[:,0].T)
    def refinePeaks(self,peaks,data,show):
        width=int(len(data[0])*0.05)
        for i in xrange(len(peaks)):
            index=list(data[0]).index(peaks[i][0])
            area=data[1][index-width:index+width]
            m=area.min()
            shift=list(area).index(m)-width
            new=index+shift
            peaks[i]=[data[0][new],data[1][new]]
        if show:
            plt.plot(data[0],data[1],peaks[:,0],peaks[:,1],'or')
            plt.show()
        return peaks
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
        '''data[0] must contain x-, data[1] y-axis data, parameter a guess of the 
        17 parameters in broadSpecTheory'''
        data=self.data
        self.setPeaks()
        if not parameter:
            parameter=[]
            peaks=self.peaks
            level=self.guessLevel(data)
            drift=self.guessDrift()
            rng=self.x_range
            #level is assumed to return the central signal level
            #we need the level on the left
            parameter+=list([level-drift*rng/2,drift])
            parameter+=list([rng*0.01]*6)
            self.depth=list(level-peaks[:,1].T)
            print 'guessed parameters:'
            print parameter        
    
        ParNames=["level","drift","width1","width2","width3","width4","width5","width6"]
        gr=dataprocessing.create_TGraph(data)
        fit1=ROOT.TF1("fit1",self.fineSpecTheory, min(data[0]), max(data[0]), 8)
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
    
        
data=np.load('25_11_temp0.npy')
F=fineSpec(data)

def fit_all_files_in(path="C:\\Dokumente und Einstellungen\\universe\\Desktop\\Projects\\DAVLL\\Spektren\\25_11_2013", showgraphs=False):
    peaks=list()
    directory = os.path.join(path)
    for files in os.walk(directory):
        for file in files[2]:
            if file.endswith(".npy"):
                print 'Processing: '+file
                try:
                    data=np.load(path+'\\'+file)
                    F=fineSpec(data)
                    F.setPeaks()
                    pks=F.peaks
                    peaks.append(pks)
                    print pks
                except:
                    continue
    return peaks