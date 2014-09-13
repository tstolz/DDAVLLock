# -*- coding: utf-8 -*-
"""
Created on Wed Sep 10 13:14:10 2014

@author: Thomas Stolz
"""

from fitmodels import *
import numpy as np
import time
import matplotlib.pyplot as plt
from ROOT import TMath, TGraph, TGraphErrors, TF1, TCanvas, TVirtualFitter

class SpecFitter:
    def __init__(self, data, noiselevel=None, FModel=Hg199_204, BModel=Linear):  
        self.Data=None    
        self.Background=None
        self.Residuals=None
        self.BModel=None
        self.FModel=None
        self.setData(data,noiselevel,FModel,BModel)
        return
        
    def setData(self,data,noiselevel=None,FModel=None,BModel=None):
        Data=np.asarray(data)
        if Data.ndim==1:
            self.numpoints=len(Data)
            l=np.arange(self.numpoints)
            self.raw_data=[l,Data]
        elif Data.ndim==2:
            if not len(Data[0])==len(Data[1]):
                return 'invalid data'
            else:
                self.raw_data=[Data[0],Data[1]]
                self.numpoints=len(Data[0])
        else:
            return "invalid data"
        self.xmin=self.raw_data[0].min()
        self.xmax=self.raw_data[0].max()
        self.raw_background=removePeaks(self.raw_data)
        if not noiselevel:
            self.noise=np.std(self.raw_background[1])
        self._rootDestruct(self.Background)
        self.Background=create_TGraph(self.raw_background,0,0)
        self._rootDestruct(self.Data)
        self.Data=create_TGraph(self.raw_data,0,self.noise)
        self.setBModel(BModel, FModel)
    def setBModel(self, BModel=None, FModel=None):
        if not BModel:
            BModel=self.BModel.__class__
            self._rootDestruct(self.BModel.TF1)
        self.BModel=BModel(self)
        self._fit(self.Background, self.BModel)
        #subtract fitted background to get foreground
        self.raw_foreground=[self.raw_data[0], [(self.raw_data[1][i] - \
        self.BModel.evaluate(self.raw_data[0][i])) for i in xrange(self.numpoints)]]
        self.findPeaks(FModel)
    def findPeaks(self, FModel=None, show=False, win_len=31):
        self.peaks=findPeaks(self.raw_foreground, show, win_len)
        self.setFModel(FModel)
    def setFModel(self,FModel=None):
        if not FModel:
            FModel=self.FModel.__class__
            self._rootDestruct(self.FModel.TF1)
        self.FModel=FModel(self, self.BModel)
        self.calcResiduals()
    def calcResiduals(self):
        self.raw_residuals=[self.raw_data[0], [(self.raw_data[1][i] - \
        self.FModel.evaluate(self.raw_data[0][i])) for i in xrange(self.numpoints)]]
        self._rootDestruct(self.Residuals)
        self.Residuals=create_TGraph(self.raw_residuals,0,0)
    
    def _rootDestruct(self,obj):
        if obj:
            obj.IsA().Destructor(obj)
    def _fit(self, Graph, Model, opt="MEX0S"):
        tm=time.time()
        result=Graph.Fit(Model.TF1, opt, "", self.xmin, self.xmax)
        print "Fitting time: %f s" % (time.time()-tm)
        Model.processFitResults(result)
        print "Fit converged with reduced chi-squared: %f" % (Model.ChiSquareRed)

    def fit(self, opt="MEX0S"):
        'fit the FModel to Data'
        self._fit(self.Data, self.FModel, opt)
        self.calcResiduals()
    
    def fitWithFixedBackground(self, opt="MEX0S"):
        'fit the FModel to Data with background parameters fixed'
        self.FModel.fixBackground()
        self._fit(self.Data, self.FModel, opt)
        self.FModel.releaseBackground()
        self.calcResiduals()
    
    def draw(self):
        canvas=TCanvas("Fit", "Fit")
        canvas.Divide(1,2)
        canvas.cd(1)
        self.Data.Draw("ap")
        self.FModel.TF1.Draw("same")
        canvas.cd(2)
        self.Residuals.Draw("ap")
        raw_input("press any key")
        self._rootDestruct(canvas)
        
        

#--------------------UTILITY FUNCTIONS--------------------

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
    return np.array([x,y])

def enhancedDerivative(data,win_len=31):
    pts=len(data[0])
    x=np.array(data[0])
    y=np.array(data[1])
    profile=[0]    
    #profile is like the discrete derivative of data
    for i in xrange(pts-1):
        if x[i+1]==x[i]:
            profile.append(0)
            continue
        profile.append((y[i+1]-y[i])/(x[i+1]-x[i]))
    #create a spectrum containing the sum of values 
    #within a window around the center point    
    win_width=int(win_len/2)
    enhanced=[0]*win_width
    for i in xrange(win_width,pts-win_width):
        window=np.array(profile[i-win_width:i+win_width])
        enhanced.append(window.sum())
    enhanced+=[0]*win_width
    return np.array(enhanced)    

def findPeaks(data,show=False,win_len=31):
    ''' will return position of the peaks (zeros of first derivative) and their height and width.
    level and drift of the signal help identify the peaks.'''
    x=np.array(data[0])
    y=np.array(data[1])
    pts=len(data[0])
    enhanced=enhancedDerivative(data,win_len)
    #go through the spectrum from left to right
    #checking for values above and below 
    #a threshold (thresh) until there is a jump in the sign of these values.
    #find the zero points between the jumps, these are our peak-center-candidates
    thresh=np.abs(enhanced).max()*0.1
    Ythresh=np.abs(y).max()*0.1
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
        
        #determine the center of the peak
        c=lastcnt+np.abs(y[lastcnt:i]).argmax()
        #check if the peak has the right curvature considering its sign and is not too flat
        if abs(y[c]) < Ythresh:
            lastcnt=i
            sgn=sgntemp
            continue
        elif not np.sign(y[c]) == sgn:
            lastcnt=i
            sgn=sgntemp            
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
        

    if len(PEAKS):
        PEAKS=np.array(PEAKS)
    else:
        print "WARNING: NO PEAKS FOUND!"
        PEAKS=np.array([[0,0,1]])
    if show:
        plt.subplot(121)
        plt.plot(x,enhanced,'o')
        plt.subplot(122)
        plt.plot(x,y,PEAKS[:,1],PEAKS[:,0],'or')
        plt.show()
    return PEAKS

def create_TGraph(data,xerr,yerr):
    n=len(data[0])    
    x=np.asarray(data[0],dtype=np.float64)
    y=np.asarray(data[1],dtype=np.float64)
    if yerr==0:
        T2gr=TGraph(n,x,y)
    else:
        xnoise=np.asarray([xerr]*n,dtype=np.float64)
        ynoise=np.asarray([yerr]*n,dtype=np.float64)
        T2gr=TGraphErrors(n,x,y,xnoise,ynoise)
    T2gr.SetLineColor( 4 )
    T2gr.SetLineWidth( 1 )
    T2gr.SetMarkerColor( 4 )
    T2gr.SetMarkerStyle( 21 )
    T2gr.SetMarkerSize(0.3)
    T2gr.SetTitle( 'Fit' )
    T2gr.GetXaxis().SetTitle( 'Scan Voltage' )
    T2gr.GetYaxis().SetTitle( 'Signal Voltage' )
    return T2gr

#quickstart
d=np.loadtxt('timenew17.txt')
l=np.linspace(0,1,1000)
s=SpecFitter([d[3],d[1]], BModel=Linear)