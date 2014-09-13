# -*- coding: utf-8 -*-
"""
Created on Wed Sep 10 13:29:44 2014

@author: Thomas Stolz
"""
from ROOT import TMath, TGraph, TGraphErrors, TF1, TCanvas, TVirtualFitter
import numpy as np

class FitModel:
    def __init__(self, Fitter):
        self.Fitter=Fitter
        self.FitResult=None
        self.TF1=None
        self.parnames=None
        self.ChiSquareRed=None
        self.MaxCorrelation=None
        self.CovMatrixPosDef=False
        self.Status=None
        self.createTF1()
        self.params=np.zeros(self.getNumPars())
        self.errors=np.zeros(self.getNumPars())
        self.setParNames()
        self.setInitialValues()
        self.setAppearance()
        
    def createTF1(self):
        if self.TF1:
            self.TF1.IsA().Destructor(self.TF1)
        self.TF1=TF1(self.getName(), self.getString(), self.Fitter.xmin, self.Fitter.xmax)
        
    def getName(self):
        return self.__class__.__name__
    
    def getNumPars(self):
        return self.TF1.GetNpar()        
        
    def evaluate(self, x):
        return self.TF1.Eval(x)
        
    def slope(self, x):
        x=np.float64(x)
        return self.TF1.Derivative(x)        
    
    def getConfidenceInterval(self,x):
        x=np.array([x], dtype=np.float64)
        ci=np.zeros(0,dtype=np.float64)
        self.FitResult.GetConfidenceIntervals(1,1,1,x,ci)
        return ci[0]
        
    def getString(self):
        string= self.makeString()
        return string
    
    def setAppearance(self):
        # Subclass can override
        self.TF1.SetLineWidth(2)
        self.TF1.SetLineColor(2)
    
    def setInitialValues(self):
        #Subclass can implement
        pass
    
    def setParNames(self):
        # Subclass can implement
        pass        
        
    def makeString(self):
        # Subclass has to implement
        return ""
        
    def processFitResults(self, FitResult):
        self.FitResult=FitResult
        pars=FitResult.GetParams()
        errs=FitResult.GetErrors()
        self.params, self.errors= np.zeros((2,self.getNumPars()), dtype=np.float64)
        for i in xrange(self.getNumPars()):
            self.params[i]=pars[i]
            self.errors[i]=errs[i]
        self.ChiSquareRed=FitResult.Chi2()/self.Fitter.numpoints
        self.CovMatrixPosDef=(FitResult.CovMatrixStatus()==3)
        n=FitResult.NFreeParameters()
        #get the correlation matrix
        self.CorrelationMatrix=np.matrix([[abs(FitResult.Correlation(i,j)) for i in xrange(n)] for j in xrange(n)])
        #determine the biggest correlation coefficient
        self.MaxCorrelation=np.max(self.CorrelationMatrix-np.eye(n))
        self.Status=FitResult.Status()

class Linear(FitModel):
    def makeString(self):
        return "[0]+[1]*x"
    def setParNames(self):
        self.parnames=["offset","slope"]
        for i in xrange(self.getNumPars()):
            self.TF1.SetParName(i, self.parnames[i])

class Quadratic(FitModel):
    def makeString(self):
        return "[0]+[1]*x+[2]*x*x"
    def setParNames(self):
        self.parnames=["offset","slope","curvature"]
        for i in xrange(self.getNumPars()):
            self.TF1.SetParName(i, self.parnames[i])

class ForegroundModel(FitModel):
    def __init__(self, Fitter, BModel):
        self.BModel=BModel
        self.numBPars=self.BModel.getNumPars()
        FitModel.__init__(self,Fitter)
        self.numFPars=self.getNumPars()-self.numBPars
        for i in xrange(self.numFPars, self.getNumPars()):
            self.params[i]=self.BModel.params[i-self.numFPars]
    def getName(self):
        return self.__class__.__name__ + "_" + self.BModel.getName() 
    def getString(self):
        #referencing the BModel's ROOT name automatically sets its function and
        #parameters; indexing starts after the Foreground parameters; names
        #are not set
        string=self.makeString()
        string+="+"+self.BModel.getName()
        return string    
    def processFitResuls(self, FitResult):
        FitModel.processFitResults(self, FitResult)
        #only the foreground model's parameters are relevant
        self.paramsToPeaks(self.params[:self.getNumPars()-self.numBPars])
    def paramsToPeaks(self, params):
        #Subclass has to implement
        pass
    def fixBackground(self):
        for i in xrange(self.numFPars, self.getNumPars()):
            self.TF1.FixParameter(i, self.params[i])
    def releaseBackground(self):
        for i in xrange(self.numFPars, self.getNumPars()):
            self.TF1.ReleaseParameter(i) 
        
class Hg199_204(ForegroundModel):
    '''lorentzians with minimal amount of parameters:
    --> (center199, Delta199-204, zeeman199, zeeman204, width, height0, 
        height1, height2, height3, height4, height5)
    checks, if there are 6 peaks'''

    def setParNames(self):
        ParNames=["center199", "Delta199-204", "zeeman199", "zeeman204", "width"]
        for i in xrange(6):
                ParNames.append("height%d" % i)
        self.parnames=ParNames
        for i in xrange(11):
            self.TF1.SetParName(i, ParNames[i])
    def makeString(self):
        string=""
        center=["[0]-[2]","[0]","[0]+[2]","[0]+[1]-[3]","[0]+[1]","[0]+[1]+[3]"]
        for j in xrange(6):
            string+="+[%d]*TMath::Pi()*[4]*0.5*TMath::CauchyDist(x," % (5+j)
            string+=center[j]
            string+=",[4]*0.5)" 
        return string[1:]
    def paramsToPeaks(self, pars):
        peaks=np.zeros((6,3))
        peaks[1][1]=pars[0]
        peaks[4][1]=pars[0]+pars[1]
        peaks[0][1]=pars[0]-pars[2]
        peaks[2][1]=pars[0]+pars[2]
        peaks[3][1]=pars[0]+pars[1]-pars[3]
        peaks[5][1]=pars[0]+pars[1]+pars[3]
        for j in xrange(6):
            peaks[j][2]=pars[4]
            peaks[j][0]=pars[j+5]
        self.Fitter.peaks=peaks
    def setInitialValues(self):
        peaks=np.asarray(self.Fitter.peaks)
        #if there are less than 6 peaks, this is the wrong model
        if len(peaks)<6:
            return 'invalid data'
        #if there are more than 6 peaks, delete the smallest
        while len(peaks)>6:
            #obtain a list of the peak amplitudes by transposing peaks
            heights=np.abs(np.asarray(peaks).T[0])
            peaks=np.delete(peaks,heights.argmin(),0)
        #calculate mean full-width-half-maximum
        FWHM=np.mean(peaks,axis=0)[2]
        #create the initial list of parameters
        params=[peaks[1][1],peaks[4][1]-peaks[1][1],(peaks[2][1]-peaks[0][1])/2.,(peaks[5][1]-peaks[3][1])/2.,FWHM]
        for peak in peaks:
            params.append(peak[0])
        for i in xrange(11):
            self.params[i]=params[i]
            self.TF1.SetParameter(i, self.params[i])
    def getXLockPoint(self, ratio, D199_204=97):
        '''Get a point on the x-axis from 'ratio', which is the position
        relative to the Hg199 peak divided by the splitting between 199 and 204 
        (D199_204).'''
        return self.params[0]+ float(ratio)/D199_204 * self.params[1]
    def getYLockPoint(self, ratio, D199_204=97):
        '''Get the function value at the point returned by getXLockPoint(self,
        ratio, D199_204).'''
        return self.evaluate(self.getXLockPoint(ratio, D199_204))
    def getLockPointSlope(self, ratio, D199_204=97):
        '''Get the slope of the function at the point returned by 
        getXLockPoint(self, ratio, D199_204).'''
        return self.slope(self.getXLockPoint(ratio,D199_204))
    def getRatioTo199_204(self, x, D199_204=97):
        '''Get the ratio of an x-value normalized to the given D199_204.'''
        return float(x)/self.params[1]*D199_204
    