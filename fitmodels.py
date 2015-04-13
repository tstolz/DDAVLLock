# -*- coding: utf-8 -*-
"""
Created on Wed Sep 10 13:29:44 2014

@author: Thomas Stolz
"""
from ROOT import TF1
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
        self.pValue=None
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
        self.TF1.SetNpx(1000)
        
    def getName(self):
        return self.__class__.__name__
    
    def getNumPars(self):
        return self.TF1.GetNpar()        
        
    def evaluate(self, x):
        return self.TF1.Eval(x)
        
    def slope(self, x):
        x=np.float64(x)
        return self.TF1.Derivative(x)        
    
    def getConfidenceInterval(self,x, cl = 0.68):
        x=np.array([x], dtype=np.float64)
        ci=np.zeros(1,dtype=np.float64)
        # get confidence interval from ROOT with 99 % cl
        self.FitResult.GetConfidenceIntervals(1,1,1,x,ci, cl)
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
        self.ChiSquareRed=FitResult.Chi2()/FitResult.Ndf()
        self.pValue=FitResult.Prob()
        self.CovMatrixPosDef=(FitResult.CovMatrixStatus()==3)
        n=FitResult.NFreeParameters()
        #get the correlation matrix
        self.CorrelationMatrix=np.matrix([[FitResult.Correlation(i,j) for i in xrange(n)] for j in xrange(n)])
        #determine the biggest correlation coefficient
        self.MaxCorrelation=np.max(np.abs(self.CorrelationMatrix)-np.eye(n))
        self.Status=FitResult.Status()
    
    def isValid(self):
        if not self.FitResult:
            return False
        elif self.MaxCorrelation>0.95 or not \
        self.CovMatrixPosDef or not (self.Status==0 or self.Status==4000):
            print "reduced chi-2: %f" % (self.ChiSquareRed)
            print "max. correlation: %f" % (self.MaxCorrelation)
            print "covariance matrix pos.def.: "+ str(self.CovMatrixPosDef)
            print "status code: %d" % (self.Status)            
            return False
        else:
            return True
            
#--- Background Models
class Constant(FitModel):
    def makeString(self):
        return "[0]"
    def setParNames(self):
        self.parnames=["offset"]
        for i in xrange(self.getNumPars()):
            self.TF1.SetParName(i, self.parnames[i])

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

#---Superclass for Foreground Models
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
    def processFitResults(self, FitResult):
        FitModel.processFitResults(self, FitResult)
        #only the foreground model's parameters are relevant
        self.paramsToPeaks(self.params[:self.getNumPars()-self.numBPars],self.errors[:self.getNumPars()-self.numBPars])
    def paramsToPeaks(self, params, errors):
        #Subclass has to implement
        pass
    def fixBackground(self):
        for i in xrange(self.numFPars, self.getNumPars()):
            self.TF1.FixParameter(i, self.params[i])
    def releaseBackground(self):
        for i in xrange(self.numFPars, self.getNumPars()):
            self.TF1.ReleaseParameter(i) 
        
class Hg199_204(ForegroundModel):
    '''lorentzians with minimal parametrization:
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
    def paramsToPeaks(self, pars, errs):
        peaks=np.zeros((6,3))
        errors=np.zeros((6,3))
        peaks[1][1]=pars[0]
        errors[1][1]=errs[0]
        peaks[4][1]=pars[0]+pars[1]
        errors[4][1]=np.sqrt(errs[0]**2+errs[1]**2)      
        peaks[0][1]=pars[0]-pars[2]
        errors[0][1]=np.sqrt(errs[0]**2+errs[2]**2)   
        peaks[2][1]=pars[0]+pars[2]
        errors[2][1]=np.sqrt(errs[0]**2+errs[2]**2)   
        peaks[3][1]=pars[0]+pars[1]-pars[3]
        errors[3][1]=np.sqrt(errs[0]**2+errs[1]**2+errs[3]**2)
        peaks[5][1]=pars[0]+pars[1]+pars[3]
        errors[5][1]=np.sqrt(errs[0]**2+errs[1]**2+errs[3]**2)
        for j in xrange(6):
            peaks[j][2]=pars[4]
            errors[j][2]=errs[4]
            peaks[j][0]=pars[j+5]
            errors[j][0]=errs[j+5]
        self.Fitter.peaks=peaks
        self.Fitter.errors=errors
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

class Gaussian(ForegroundModel):
    def makeString(self):
        string=""
        for i in xrange(len(self.Fitter.peaks)):
            string+="+gaus(%d)" % (i*3)
        return string[1:]
    def setInitialValues(self):
        params=[]
        for peak in self.Fitter.peaks:        
            params.append(peak[0])
            params.append(peak[1])
            params.append(peak[2]/(2*np.sqrt(2*np.log(2))))
        for i in xrange(len(params)):
            self.params[i]=params[i]
            self.TF1.SetParameter(i, self.params[i])
    def setParNames(self):
        self.parnames=[]
        for i in xrange(len(self.Fitter.peaks)):
            self.parnames+=["height%d" % (i), "center%d" % (i), "sigma%d" % (i)]
            for j in (0,1,2):
                self.TF1.SetParName(3*i+j, self.parnames[3*i+j])
    def paramsToPeaks(self, pars, errs):
        for i in xrange(self.numFPars/3):
            self.Fitter.peaks[i][0]=pars[i]
            self.Fitter.peaks[i][1]=pars[i+1]
            self.Fitter.peaks[i][2]=pars[i+2]*2*np.sqrt(2*np.log(2))
            self.Fitter.errors[i][0]=errs[i]
            self.Fitter.errors[i][1]=errs[i+1]
            self.Fitter.errors[i][2]=errs[i+2]*2*np.sqrt(2*np.log(2))            

class Lorentzian(ForegroundModel):      
    def makeString(self):
        string=""
        for i in xrange(len(self.Fitter.peaks)):
                string+="+[%d]*TMath::Pi()*[%d]*0.5*TMath::CauchyDist(x,[%d],[%d]*0.5)" % (i*3,i*3+2,i*3+1,i*3+2)
        return string[1:]
    def setInitialValues(self):
        params=np.ravel(self.Fitter.peaks)
        for i in xrange(len(params)):
            self.params[i]=params[i]
            self.TF1.SetParameter(i, self.params[i])
    def setParNames(self):
        self.parnames=[]
        for i in xrange(len(self.Fitter.peaks)):
            self.parnames+=["height%d" % (i), "center%d" % (i), "width%d" % (i)]
            for j in (0,1,2):
                self.TF1.SetParName(3*i+j, self.parnames[3*i+j])
    def paramsToPeaks(self, pars, errs):
        for i in xrange(self.numFPars/3):
            self.Fitter.peaks[i][0]=pars[3*i]
            self.Fitter.peaks[i][1]=pars[3*i+1]
            self.Fitter.peaks[i][2]=pars[3*i+2]
            self.Fitter.errors[i][0]=errs[3*i]
            self.Fitter.errors[i][1]=errs[3*i+1]
            self.Fitter.errors[i][2]=errs[3*i+2]           


class Hg199_204_Voigt(ForegroundModel):
    '''Voigt profiles with minimal parametrization:
    --> (center199, Delta199-204, zeeman199, zeeman204, widthL, widthG, height0, 
        height1, height2, height3, height4, height5)
    checks, if there are 6 peaks'''
    #DANGER: Fails if widthG is too small!!!
    def setParNames(self):
        ParNames=["center199", "Delta199-204", "zeeman199", "zeeman204", "widthL", "widthG"]
        for i in xrange(6):
                ParNames.append("height%d" % i)
        self.parnames=ParNames
        for i in xrange(12):
            self.TF1.SetParName(i, ParNames[i])
    def makeString(self):
        string=""
        center=["[0]-[2]","[0]","[0]+[2]","[0]+[1]-[3]","[0]+[1]","[0]+[1]+[3]"]
        for j in xrange(6):
            string+="+[%d]*sqrt(2*pi)*[5]*exp(-[4]*[4]/(8*[5]*[5]))/TMath::Erfc(sqrt(2)*[4]/(4*[5]))*TMath::Voigt(x-(" % (6+j)
            string+=center[j]
            string+="),[5],[4],5)"
        return string[1:]
    def paramsToPeaks(self, pars, errs):
        peaks=np.zeros((6,3))
        errors=np.zeros((6,3))
        peaks[1][1]=pars[0]
        errors[1][1]=errs[0]
        peaks[4][1]=pars[0]+pars[1]
        errors[4][1]=np.sqrt(errs[0]**2+errs[1]**2)      
        peaks[0][1]=pars[0]-pars[2]
        errors[0][1]=np.sqrt(errs[0]**2+errs[2]**2)   
        peaks[2][1]=pars[0]+pars[2]
        errors[2][1]=np.sqrt(errs[0]**2+errs[2]**2)   
        peaks[3][1]=pars[0]+pars[1]-pars[3]
        errors[3][1]=np.sqrt(errs[0]**2+errs[1]**2+errs[3]**2)
        peaks[5][1]=pars[0]+pars[1]+pars[3]
        errors[5][1]=np.sqrt(errs[0]**2+errs[1]**2+errs[3]**2)
        for j in xrange(6):
            peaks[j][2]=pars[4]
            errors[j][2]=errs[4]
            peaks[j][0]=pars[j+6]
            errors[j][0]=errs[j+6]
        self.Fitter.peaks=peaks
        self.Fitter.errors=errors
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
        #calculate mean full-width-half-maximum for alpha == 4
        FWHM=np.mean(peaks,axis=0)[2]
        #create the initial list of parameters
        #set vratio to 0.001 as a start
        params=[peaks[1][1],peaks[4][1]-peaks[1][1],(peaks[2][1]-peaks[0][1])/2.,(peaks[5][1]-peaks[3][1])/2.,FWHM,0.1*FWHM]
        for peak in peaks:
            params.append(peak[0])
        for i in xrange(12):
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

class Hg199_204_Voigt2(ForegroundModel):
    '''Voigt profiles with minimal amount of parameters:
    --> (center199, Delta199-204, zeeman199, zeeman204, width, vratio, height0, 
        height1, height2, height3, height4, height5)
    checks, if there are 6 peaks
    uses the approximate FWHM of the Voigt profile <width> and a shape parameter <vratio> '''
    #DANGER: Fails if vratio is too large!!!
    def setParNames(self):
        ParNames=["center199", "Delta199-204", "zeeman199", "zeeman204", "width", "vratio"]
        for i in xrange(6):
                ParNames.append("height%d" % i)
        self.parnames=ParNames
        for i in xrange(12):
            self.TF1.SetParName(i, ParNames[i])
    def makeString(self):
        string=""
        center=["[0]-[2]","[0]","[0]+[2]","[0]+[1]-[3]","[0]+[1]","[0]+[1]+[3]"]
        for j in xrange(6):
            string+="+[%d]*exp(-[5])*[4]*(0.5*sqrt(146.40920767581318*exp(-2*[5]) + 5.718885409907285)\
            - 1.37348977270357)/(14.24662272319341*exp(-2*[5]) - 0.1777807548824028)*exp(-exp(2*[5])/8)*\
            sqrt(pi*2)/TMath::Erfc(sqrt(2)*exp([5])/(4))*TMath::Voigt(x-(" % (6+j)
            string+=center[j]
            string+="),exp(-[5])*([4]*(0.5*sqrt(146.40920767581318*exp(-2*[5]) + 5.718885409907285) - \
            1.37348977270357))/(14.24662272319341*exp(-2*[5]) - 0.1777807548824028),\
            ([4]*(0.5*sqrt(146.40920767581318*exp(-2*[5]) + 5.718885409907285) - \
            1.37348977270357))/(14.24662272319341*exp(-2*[5]) - 0.1777807548824028),4)"
        return string[1:]
    def paramsToPeaks(self, pars, errs):
        #this is not correct yet
        peaks=np.zeros((6,3))
        errors=np.zeros((6,3))
        peaks[1][1]=pars[0]
        errors[1][1]=errs[0]
        peaks[4][1]=pars[0]+pars[1]
        errors[4][1]=np.sqrt(errs[0]**2+errs[1]**2)      
        peaks[0][1]=pars[0]-pars[2]
        errors[0][1]=np.sqrt(errs[0]**2+errs[2]**2)   
        peaks[2][1]=pars[0]+pars[2]
        errors[2][1]=np.sqrt(errs[0]**2+errs[2]**2)   
        peaks[3][1]=pars[0]+pars[1]-pars[3]
        errors[3][1]=np.sqrt(errs[0]**2+errs[1]**2+errs[3]**2)
        peaks[5][1]=pars[0]+pars[1]+pars[3]
        errors[5][1]=np.sqrt(errs[0]**2+errs[1]**2+errs[3]**2)
        for j in xrange(6):
            peaks[j][2]=pars[4]
            errors[j][2]=errs[4]
            peaks[j][0]=pars[j+6]
            errors[j][0]=errs[j+6]
        self.Fitter.peaks=peaks
        self.Fitter.errors=errors
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
        #calculate mean full-width-half-maximum for alpha == 4
        FWHM=np.mean(peaks,axis=0)[2]
        #create the initial list of parameters
        #set vratio to 0.001 as a start
        params=[peaks[1][1],peaks[4][1]-peaks[1][1],(peaks[2][1]-peaks[0][1])/2.,(peaks[5][1]-peaks[3][1])/2.,FWHM,1.5]
        for peak in peaks:
            params.append(peak[0])
        for i in xrange(12):
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

#Sample for creating new models

#class Gaussian(ForegroundModel):
#    def makeString(self):
#    
#    def setInitialValues(self):
#
#    def setParNames(self):
#        
#    def paramsToPeaks(self):
#        