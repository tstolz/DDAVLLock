# -*- coding: utf-8 -*-
"""
Created on Wed Sep 10 13:14:10 2014

@author: Thomas Stolz
"""

from fitmodels import *
from spectrum_utility import *
import numpy as np
from scipy import stats
import time
from ROOT import TCanvas, gStyle
from ROOT.Math import GoFTest

class SpecFitter:
    def __init__(self, data, noiseRMS=None, numpeaks=None, FModel=Hg199_204, BModel=Linear):  
        self.Data=None    
        self.Background=None
        self.Residuals=None
        self.NormProbPlt=None
        self.NormProbPltChi2=0
        self.BModel=None
        self.FModel=None
        self.alive=True
        self.setData(data,noiseRMS,numpeaks,FModel,BModel)
        return
    def isAlive(self):
        return self.alive
    def kill(self):
        self.alive=False
    def setData(self,data,noiseRMS=None,numpeaks=None,FModel=None,BModel=None):
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
        self._rootDestruct(self.Background)
        self.Background=create_TGraph(self.raw_background,0,0)
        self.Background.SetTitle("Background")
        self._rootDestruct(self.Data)
        self.setBModel(BModel, noiseRMS, FModel, numpeaks)
    def setBModel(self, BModel=None, noiseRMS=None, FModel=None, numpeaks=None):
        if not BModel:
            BModel=self.BModel.__class__
            self._rootDestruct(self.BModel.TF1)
        self.BModel=BModel(self)
        self._fit(self.Background, self.BModel, opt="MQNS")
        #subtract fitted background to get foreground
        self.raw_foreground=[self.raw_data[0], [(self.raw_data[1][i] - \
        self.BModel.evaluate(self.raw_data[0][i])) for i in xrange(self.numpoints)]]
        if not noiseRMS:        
        #subtract the fitted background so we don't overestimate the noise
        #in case of a strong background
            self.noise=np.std([(self.raw_background[1][i] - \
                        self.BModel.evaluate(self.raw_background[0][i])) \
                        for i in xrange(len(self.raw_background[0]))])
        else:
            self.noise=noiseRMS    
        #now we can create the Data graph, no uncertainty for x-values
        self.Data=create_TGraph(self.raw_data,0,self.noise)
        self.Data.SetTitle("Spectrum")
        self.findPeaks(numpeaks, FModel)
    def findPeaks(self, numpeaks=None, FModel=None, show=False):
        if not numpeaks:
            win_len=33
        else:
            win_len=self.numpoints/(numpeaks*5)
        self.peaks=findPeaks(self.raw_foreground, show, win_len)
        if numpeaks:
            #problem if there are less than numpeaks
            if len(self.peaks)<numpeaks:
                self.kill()
            #if there are more than 6 peaks, delete the smallest
            while len(self.peaks)>numpeaks:
                #obtain a list of the peak amplitudes by transposing peaks
                heights=np.abs(np.asarray(self.peaks).T[0])
                self.peaks=np.delete(self.peaks,heights.argmin(),0)
        self.errors=np.zeros(np.shape(self.peaks))
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
        self.Residuals.SetTitle("Residuals")
        #ceate a normal probability plot
        self.raw_normprobplt=np.array(stats.probplot(self.raw_residuals[1],fit=False,plot=None))
        self._rootDestruct(self.NormProbPlt)
        self.NormProbPlt=create_TGraph(self.raw_normprobplt,0,0)
        self.NormProbPlt.GetYaxis().SetTitle('quantiles')
        self.NormProbPlt.GetXaxis().SetTitle('residuals')
        self.NormProbPlt.SetTitle("Normal Probability Plot")
        result=self.NormProbPlt.Fit("pol1","QS")    
        self.NormProbPltChi2=result.Chi2()
    def _rootDestruct(self,obj):
        #destructor method to safely dispose of ROOT objects
        if obj:
            obj.IsA().Destructor(obj)
    def _fit(self, Graph, Model, opt="MQS"):
        tm=time.time()
        result=Graph.Fit(Model.TF1, opt, "", self.xmin, self.xmax)
        if opt.find("Q")==(-1):
            print "Fitting time: %f s" % (time.time()-tm)
        Model.processFitResults(result)
    def fit(self, opt="MS"):
        'Fit FModel to Data'
        self._fit(self.Data, self.FModel, opt)
        self.calcResiduals()  
        if opt.find("Q")==(-1):
            print "AD-test p-value: %f" % (self.getADTestPValue())
    def fitWithFixedBackground(self, opt="MS"):
        'Fit FModel to Data with background parameters fixed'
        self.FModel.fixBackground()
        self._fit(self.Data, self.FModel, opt)
        self.FModel.releaseBackground()
        self.calcResiduals() 
    def getTheoryArray(self):
        'Calculate theory function values of FModel corresponding to x-data'
        xvals=self.raw_data[0]
        return np.array([xvals,[self.FModel.evaluate(x) for x in xvals]])
    def getADTestPValue(self):
        'Perform an Anderson-Darling-Test to check normal distribution of residuals'
        g = GoFTest(self.numpoints,np.array(self.raw_residuals[1]),GoFTest.kGaussian)
        pval=g.AndersonDarlingTest()
        return pval
    def getKSTestPValue(self):
        'Perform a Kolmogorov-Smirnov-Test to check normal distribution of residuals'
        g = GoFTest(self.numpoints,np.array(self.raw_residuals[1]),GoFTest.kGaussian)
        pval=g.KolmogorovSmirnovTest()
        return pval
    def getFourierSpectrum(self,data):
        '''This function returns a Fourier power spectrum assuming that data[0]
        contains evenly spaced time data in seconds'''        
        fft=np.abs(np.fft.rfft(data[0]))**2
        freq=np.fft.rfftfreq(len(data[1]),np.ptp(data[0])/(np.size(data)-1))
        return [freq, fft]
    def residualsWithin5Sigma(self):
        '''Checks if the residuals are all within 5 standard deviations of the
        estimated noiselevel'''
        return (np.sort(np.abs(self.raw_residuals[1]))[-1]<5*self.noise)
    def signalToNoise(self):
        '''Returns the estimated signal-to-noise ratio'''
        return self.peaks[:,0].max()/self.noise
    def save(self,fname):
        d=np.array([self.raw_data[0],self.raw_data[1],self.getTheoryArray()[1],self.raw_residuals[1]])
        np.savetxt(fname, d.T)
    def draw(self):
        '''Draw data with fit, residuals, normal probability plot and 
        background'''
        #display the p-value and chi-square in a box when drawing
        gStyle.SetOptFit(1111)
        canvas=TCanvas("Fit", "Fit")
        canvas.Divide(1,3)
        canvas.cd(1)
        self.Data.Draw("ap")
        #self.FModel.TF1.Draw("same")
        pad=canvas.cd(2)
        #gStyle.SetOptFit(0000)
        pad.Divide(2,1)
        pad.cd(1)
        self.Residuals.Draw("ap")
        pad.cd(2)
        self.NormProbPlt.Draw("ap")
        canvas.cd(3)
        self.Background.Draw("ap")
        self.BModel.TF1.Draw("same")
        raw_input("press any key")
        self._rootDestruct(canvas)
        self._rootDestruct(pad)



#quickstart
d=np.loadtxt('spectrum11.txt')
d=d[:,:]
s=SpecFitter([d[0],d[1]],numpeaks=6, FModel=Hg199_204, BModel=Linear)
s.fit()

#d=np.load('26_11_bigcell_mag1.24A.npy')
#d=d[:,200:600]
#s=SpecFitter([d[0],d[1]], noiseRMS=0.0058573036693445706, numpeaks=6, FModel=Hg199_204, BModel=Quadratic)
#s.fit()
#s.draw()
#xvals=[s.FModel.getRatioTo199_204(d[0][i]-s.FModel.params[0]) for i in xrange(len(d[0]))]
#t=SpecFitter([xvals,d[1]], noiseRMS=0.0058573036693445706, numpeaks=6, FModel=Hg199_204, BModel=Quadratic)
#t.fit()
#t.draw()

