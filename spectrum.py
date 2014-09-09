# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 14:57:32 2013

@author: Thomas Stolz
"""

from spectrumstyles import *
import numpy as np
from scipy.special import erf
import matplotlib.pyplot as plt
from pylab import get_current_fig_manager
import time
from ROOT import TMath, TGraph, TGraphErrors, TF1, TCanvas, TVirtualFitter

global wRatio
#ratio of the Doppler width to the convoluted peak width (FWHM)
#depends on the linewidth of the laser and saturation or pressure broadening
wRatio=50
#can be calculated from doppler width (430 MHz) and fit params! 

class Spectrum:
    def __init__(self, data, style="voigt"):
        self.data=np.asarray(data)
        self.peaks=None
        if not style in ["lorentz","gauss","voigt"]:
            print "unknown style, set to lorentz"
            self.style="lorentz"
        else:
            self.style=style
        self.x_range=self.data[0].max()-self.data[0].min()
        self.points=len(self.data[0])
        self.ChiSquareRed=None
        self.params=None
        self.errors=None
        self.noise_level=None
        #ratio of the gaussian width in the voigt profile
        self.vratio=0
        self.vratioErr=0
        #indicator if the spectrum has been fitted yet
        self.guessed=False
        self.fitted=False
        self.exact=False
    def setStyle(self,style):
        if not style in ["lorentz","gauss","voigt"]:
            print "unknown style"
        else:
            self.style=style
    def theory(self, x, params):
        '''return the theory value at x for a spectrum with parameters "params"

        params must have the form "level, drift, curvature, height1, center1, width1, ... heightN, centerN, widthN'''
        level = params[0]
        drift = params[1]
        curv = params[2]
        peaks = [[params[i],params[i+1],params[i+2]] for i in range(3,len(params),3)]
        result = level + x*drift + x**2 * curv
        for peak in peaks:
            if self.style=="gauss":
                result += peak[0]*TMath.Gaus(x,peak[1],peak[2])
            elif self.style=="lorentz":
                #ROOT defined the CauchyDist with HWHM and Voigt with FWHM
                #so we get a factor 0.5 in front of peak[2] (FWHM)
                result += peak[0]*TMath.Pi()*0.5*peak[2]*TMath.CauchyDist(x,peak[1],peak[2]*0.5)
            elif self.style=="voigt":
                #vratio is now the angle between the beams (after the theory...) 
                #92.4065 is the ratio between the gaussian and the lorentzian width                
                result += peak[0]*TMath.Voigt(x-peak[1],self.vratio*peak[2],peak[2])
        return result
    def noisefit1(self):
        '''subtract the fitted signal from data and calc noiselevel of the rest'''
        if not self.params==None:
            noise=self.data[1]-[self.theory(x,self.params) for x in self.data[0]]
            self.noise_level=np.sqrt(np.var(noise))
        else:
            self.noise_level=None
        return self.noise_level
    def noisefit2(self):
        """find the rms noise in the the first 1/100 of data"""
        #1/100 seems to be to small!! noise has low frequency
        noiseregion=self.data[1][:len(self.data[1])/100]
        sigma=np.sqrt(np.mean([noiseregion[i]**2 for i in range(len(noiseregion))])-(np.mean(noiseregion))**2)
        self.noise_level=sigma
        return sigma
    def draw(self):
        plt.ion()
        self.plot()
        plt.show()
    def plot(self):
        if not self.params==None:
            Fit=np.asarray([self.theory(x,self.params) for x in self.data[0]])
            plt.plot(self.data[0],self.data[1],'k,',self.data[0],Fit,'r-',self.data[0],self.data[1]-Fit,'g-')
        else:
            plt.plot(self.data[0],self.data[1],'k,')
    def peaksToParams(self):
        '''convert the previously obtained peak values to fit parameter values.
        As guessParams returns the absolute height and FWHM of the peaks, these
        values have to be converted for some line shapes.'''
        for peak in self.peaks:
            #height and width (FWHM) have to be adjusted to the form of the curve
            if self.style=='gauss':
                #not yet FWHM
                width=peak[2]/(2*np.sqrt(2*np.log(2)))
                height=peak[0]
            elif self.style=='lorentz' or self.vratio==0:
                width=peak[2]
                height=peak[0]
            elif self.style=='voigt':
                #the following are approximations taken from 
                #Zing-Qiao et al., adapted to our case
                width=2*peak[2]/(1.0692 + np.sqrt(0.86639 + 4*8 * np.log(2) * self.vratio**2))
                height=peak[0]*np.sqrt(2*np.pi)*self.vratio*width/(np.exp(0.5*(2*self.vratio)**(-2))*(1-erf(1/(2*np.sqrt(2)*self.vratio))))
            self.params.append(height)
            self.params.append(peak[1])
            self.params.append(width)
            
    def guessParams(self,numPeaks=6, show=False):
        '''guess all the parameters required for the fit. order: level, drift, curvature, (height, center, width) of all peaks'''
        self.params=[]

        noPeaks=removePeaks(self.data)

        curvature, drift, level = np.polyfit(*noPeaks,deg=2)
        self.params+=list([level,drift,curvature])
        
        #estimate the best window lenght for "enhancedDerivative()"
        #from the expected number of peaks
        #calculation based on experience
        if not numPeaks:
            win_len=31
        else:
            win_len=self.points/(numPeaks*5)
        self.peaks=findPeaks(self.data,level,drift,curvature,show=show,win_len=win_len)
        
        self.peaksToParams()
        self.guessed=True
        return self.params
    
    def exactFit(self, parameters=None, showgraphs=False):
        '''Fits the spectrum two times, first to get the correct
        parameters and determine the noiselevel, then again to 
        get correct error estimates based on the noiselevel. 
        
        The uncertainties are chosen in a way, that the reduced
        chi-squared becomes ~1. By this we assume the theory 
        function to be exact.'''
        self.fit(parameters,showgraphs)
        self.noisefit1()
        self.exact=True
        return self.fit(self.params, showgraphs)
        
    def createCmdString(self):
        '''Create the theory function string for the ROOT fit.'''
        if self.style=='gauss':
            cmdstring="pol2"
            for i in xrange(len(self.peaks)):
                cmdstring+="+gaus(%d)" % (3+i*3)
        elif self.style=='lorentz':
            cmdstring="pol2"
            for i in xrange(len(self.peaks)):
                cmdstring+="+[%d]*3.14*[%d]*0.5*TMath::CauchyDist(x,[%d],[%d]*0.5)" % (3+i*3,5+i*3,4+i*3,5+i*3)
        elif self.style=='voigt':
            cmdstring="pol2"
            for i in xrange(len(self.peaks)):
                cmdstring+="+[%d]*TMath::Voigt(x-[%d],[3]*[%d],[%d])" % (4+i*3,5+i*3,6+i*3,6+i*3)
        
        return cmdstring
        
    def createParLists(self,parameters):
        '''Create two lists with the fit parameters and their names,
        that are passed to the ROOT fit.'''
        ParList=list(parameters)
        
        if self.style=="voigt":
            ParList.insert(3,self.vratio)

        ParNames=["level","drift","curvature"]        
        
        if self.style=="voigt":
            ParNames.append("vratio")

        for i in xrange((len(parameters)-3)/3):
            ParNames.append("height%d" % i)  
            ParNames.append("center%d" % i)
            ParNames.append("width%d" % i)
        
        return ParList, ParNames
        
    def processFitResults(self, pars, errs):
        '''Fills self.params, self.peaks and self.errors with the fit results.
        For voigt profiles also the vratio has to be extracted and stored 
        in self.vratio. Similar things can happen in subclasses, where this
        function is overridden.'''
        self.params=[]
        self.errors=[]
        if self.style=='gauss' or self.style=='lorentz':
            numpars=3+3*len(self.peaks)
            for i in xrange(numpars):
                if i>2:
                    self.peaks[(i-3)/3][i%3]=pars[i]
                self.params.append(pars[i])
                
        elif self.style=='voigt':
            numpars=3+1+3*len(self.peaks)
            for i in xrange(numpars):
                if i>2:
                    if i==3:
                        self.vratio=pars[i]
                        continue
                    else:
                        j=i-1
                        self.peaks[(j-3)/3][j%3]=pars[i]
                self.params.append(pars[i])
            
        for i in xrange(numpars):
            if i==3 and self.style=="voigt":
                self.vratioErr=errs[i]
                continue
            self.errors.append(errs[i])
                    
    def fit(self, parameters=None, showgraphs=True, verbose=True, quiet=False):
        '''Fits the spectrum with initial values "parameters". If None are given,
        the parameters are guessed first. Showgraphs displays the guessed curve
        before the fit. Verbose sets the "v" flag when calling the ROOT fit. 
        If "quiet" is False, the fit time and reduced chi-squared are printed.'''
        data=np.asarray(self.data,dtype=np.float64)
        numpoints=len(data[0])
        x_range=max(data[0])-min(data[0])
        #the error estimate for the adc time values is half the stepwidth
        xerr=x_range/numpoints*0.5
        if parameters==None:
            parameters=self.guessParams()      
        
        # insert signal and style dependent string into the TF1
        # first create the string
        cmdstring=self.createCmdString()
        
        # create lists with initial parameters and their names
        ParList, ParNames = self.createParLists(parameters)
        print ParList
        
        gr=create_TGraph(data, xerr, self.noise_level)
        fit1=TF1("fit1",cmdstring, min(data[0]), max(data[0]))
        for i in range(len(ParNames)):
            fit1.SetParName(i, ParNames[i])
        for n in range(len(ParList)):
            fit1.SetParameter(n,ParList[n])
        
        fit1.SetLineWidth(2)
        fit1.SetLineColor(1)
        # we could set limits for the parameters here if we want        
        
        if showgraphs:
            canvas=TCanvas("Fit","Fit")
            gr.Draw("ap")
            fit1.Draw("same")
            #print initial curve
            raw_input("pre fit, press any key")
            print 'fitting...'
            
        tm=time.time()
        if verbose:
            opt="VMEX0S"
        else:
            opt="QMEX0S"
        # Tell root to do the fit. This can take some minutes.
        resultPtr = gr.Fit(fit1, opt,"AP", data[0].min(), data[0].max())
        
        if not quiet:
            print "Fitting time: %f s" % (time.time()-tm)
        
        pars=fit1.GetParameters()
        errs=fit1.GetParErrors()
        # store the obtained fit parameters according to the style
        self.processFitResults(pars,errs)
        
        self.ChiSquareRed=fit1.GetChisquare()/float(numpoints)
        if not quiet:
            print "Fit converged with reduced chi-squared: %f" % (self.ChiSquareRed)
            print "Status: %d" % (resultPtr.Status())
        print resultPtr.GetCorrelationMatrix().GetMatrixArray()

#       ROOT status codes for Minuit2 (issues if >0):
#       status = 1    : Covariance was made pos defined
#       status = 2    : Hesse is invalid
#       status = 3    : Edm is above max
#       status = 4    : Reached call limit
#       status = 5    : Any other failure
       
        self.fitted=True        
        if showgraphs:        
            fit1.Draw("same")
            canvas.Update()
            canvas.Show()
            raw_input("done, press any key")
            canvas.IsA().Destructor(canvas)
        
        # Release ROOT objects to avoid memory issues
        gr.IsA().Destructor(gr)
        fit1.IsA().Destructor(fit1)
        #resultPtr.IsA().Destructor(resultPtr)
        #this should return a fitresult
        return None


#------------------------Subclasses-------------------------------

class DFSpec(Spectrum):
    '''Subclass for the DDAVLL. Includes special styles, where the
    peak widths are fixed, as these should all be the same in the lock
    spectrum. Also utility functions to obtain properties of the spectrum,
    especially the lock point and its error.'''
    def __init__(self,data,style='slimlorentz',freq_diff=97, freq_err=0.18):
        '''freq_diff is the isotope shift in MHz, freq_err its fractional error'''        
        Spectrum.__init__(self, data)
        self.setStyle(style)
        self.freq_diff=freq_diff
        self.freq_err=freq_err
    def setStyle(self, style):
        '''"style" should be one of 
        
        simplevoigt: voigt profiles with equal peak widths
        simplelorentz: lorentz profiles with equal peak widths
        fixedvoigt: voigt profiles with a fixed gauss-lorentz-ratio
        slimlorentz: lorentzians with minimal amount of parameters 
        --> (center199, Delta199-204, zeeman199, zeeman204, width, height0, height1, height2, height3, height4, height5)
        
        General styles from superclass are also possible.'''
        
        if not style in ["simplevoigt","simplelorentz","fixedvoigt","slimlorentz"]:
            Spectrum.setStyle(self,style)
        else:
            self.style=style
            
    def setAlpha(self,alpha):
        '''Set the gauss-lorentz-ratio assuming an angle alpha between 
        the beams.'''
        self.vratio=alpha*wRatio
    def plot(self,background=True):
        '''Valid spectra are plotted with theory curve in red,
        data points in black and signal background in green. The x-axis
        is in MHz.'''
        if not self.isValid():
            Spectrum.plot(self)
            return
        Fit=np.asarray([self.theory(x,self.params) for x in self.data[0]])
        xaxis=self.getFrequency(self.data[0])
        if background:
            plt.plot(xaxis,self.data[1],'k,',xaxis,Fit,'r-',xaxis,self.data[1]-Fit,'g-')
        plt.plot(xaxis,self.data[1],'k,',xaxis,Fit,'r-')
        plt.xlabel('frequency [MHz]')
        plt.ylabel('diode signal [V]')
        plt.grid(True)
        get_current_fig_manager().window.raise_()
    def theory(self, x, params):
        '''return the theory value at x for a spectrum with parameters "params"

        params must have the form "level, drift, curvature, height1, center1, width1, ... heightN, centerN, widthN'''
        if not self.style in ["simplevoigt","simplelorentz","fixedvoigt","slimlorentz"]:
            return Spectrum.theory(self,x,params)
            
        level = params[0]
        drift = params[1]
        curv = params[2]
        
        if self.style=="slimlorentz":
            peaks=np.zeros((6,3))
            peaks[1][1]=params[3]
            peaks[4][1]=params[3]+params[4]
            peaks[0][1]=params[3]-params[5]
            peaks[2][1]=params[3]+params[5]
            peaks[3][1]=params[3]+params[4]-params[6]
            peaks[5][1]=params[3]+params[4]+params[6]
            for j in xrange(6):
                peaks[j][2]=params[7]
                peaks[j][0]=params[j+8]
        else:
            peaks = [[params[i],params[i+1],params[i+2]] for i in range(3,len(params),3)]
        
        result = level + x*drift + x**2 * curv
        for peak in peaks:
            if self.style=="simplevoigt":
                #this approximation of a voigt profile with fixed height
                #is good for vratio < 1.5 with less than 2 % error                   
                result += peak[0]*peak[2]/(0.637-0.255*self.vratio-0.998*self.vratio**2+1.78*self.vratio**3-1.15*self.vratio**4+0.265*self.vratio**5)*TMath.Voigt(x-peak[1],self.vratio*peak[2],peak[2])
            elif self.style=="simplelorentz" or self.style=='slimlorentz':
                #ROOT defined the CauchyDist with HWHM and Voigt with FWHM
                #so we get a factor 0.5 in front of peak[2] (FWHM)
                result += peak[0]*TMath.Pi()*0.5*peak[2]*TMath.CauchyDist(x,peak[1],peak[2]*0.5)        
        return result
    def getFrequency(self,x,pars=None):
        '''return the associated frequency relative to the Hg-199 
        transition in MHz for an x-value,
        derived from the distance between the isotope lines'''
        if type(pars)==type(None):
            pars=self.params
        if not self.isValid():
            print 'invalid spectrum'
            return 0
        if self.style=='slimlorentz':
            scale=self.freq_diff/pars[4]
            zero=pars[3]
        else:
            scale=self.freq_diff/(pars[16]-pars[7])
            zero=pars[7]
        return (x-zero)*scale
    def freq2volts(self,f,pars=None):
        '''convert detuning to voltage. inverts getFrequency'''
        if type(pars)==type(None):
            pars=self.params
        if not self.isValid():
            print 'invalid spectrum'
            return 0
        scale=self.freq_diff/(pars[16]-pars[7])
        return (f/scale+pars[7])
    def lockPointX(self,freq=-8.8):
        freq=float(freq)
        err=np.sqrt((1-freq/self.freq_diff)**2*self.errors[7]**2 + (freq/self.freq_diff)**2 * self.errors[16]**2 + (freq/self.freq_diff)**2 * (self.params[16]-self.params[7])**2 * self.freq_err**2)
        return [self.freq2volts(freq),err]
    def lockPointY(self,freq=-8.8):
        '''calculates the lock point (y-axis) for a frequency detuning "freq"
        and the associated error for a relative error in that frequency "freq_err" '''
        #lockPointY
        lPY=self.theory(self.freq2volts(freq),self.params)
        #the function for the error-calculation
        F=lambda x: self.theory(self.freq2volts(freq,pars=x),x)
        
        #calculate the error by quadratically adding the errors
        #caused by the parameter errors. Instead of partial derivatives
        #we directly evaluate the function. 
        squaredsum=0
        pars=np.asarray(self.params)
        errs=np.asarray(self.errors)
        num=len(pars)
        for i in xrange(num):
            #print F(pars+np.eye(num)[i]*errs)
            squaredsum+=max([(lPY-F(pars+np.eye(num)[i]*errs*j))**2 for j in (1,-1)])
        squaredsum+=max([(lPY-self.theory(self.freq2volts(freq*(1+self.freq_err*j)),self.params))**2 for j in (1,-1)])
        err=np.sqrt(squaredsum)
        return [lPY,err]
    def getSlope(self,x):
        '''calculate the slope of the theory curve at x'''
        dx=self.x_range/self.points*0.001
        dy=self.theory(x+dx,self.params)-self.theory(x,self.params)
        return dy/dx
    def lockPrecision(self,freq=-8.8):
        '''returns the systematic uncertainty of a possible lock at
        frequency "freq"'''
        yerr=self.lockPointY(freq)[1]
        freq_err=self.getFrequency(self.peaks[1][1]+abs(yerr/self.getSlope(self.lockPointX(freq)[0])))
        return freq_err
    def getFrequencyError(self, x):
        x1=self.params[16]
        x0=self.params[7]
        fq=self.freq_diff
        err=np.sqrt(self.freq_err**2 * ((x-x0)/(x1-x0))**2 + self.errors[7]**2 * (fq*(x-x1)/(x0-x1)**2)**2 + self.errors[16]**2 * (fq*(x-x0)/(x1-x0)**2)**2)
        return err
        #to be filled! composed of uncertainty of freq_diff and peak errors
    def getChiSquareRed(self):
        return [self.ChiSquareRed,0]
    def isValid(self):
        '''Checks, if the spectrum has the right number of peaks and, if
        it has been fitted, whether the fit was successful.'''
        if not self.guessed:
            return False
        if len(self.peaks)!=6:
            return False
        if self.exact and abs(self.ChiSquareRed-1) > 0.9:
            return False
        if self.fitted and self.ChiSquareRed>30:
            return False
        #further possibilities for the spectrum to be invalid
        return True
    def alpha(self):
        '''Give an estimate of the angle between the beams. Only
        reasonable if style is simplevoigt.'''
        wRatio=430/self.meanPeakWidth()[0]
        return [self.vratio/wRatio, self.vratioErr/wRatio]
    def peaksToParams(self):
        #if too many peaks remove the weakest until 6 peaks
        while len(self.peaks)>6:
            heights=np.abs(np.asarray(self.peaks).T[0])
            self.peaks=np.delete(self.peaks,heights.argmin(),0)
        
        if not self.style in ["simplevoigt","simplelorentz","fixedvoigt","slimlorentz"]:
            Spectrum.peaksToParams(self)
            return
        
        FWHM=[]
        for peak in self.peaks:
            FWHM.append(peak[2])
        FWHM=np.mean(FWHM)
        
        if self.style=='simplevoigt' or self.style=='fixedvoigt':
            #the following approximation for the width of voigt profiles 
            #is taken from Zing-Qiao et al., adapted to our case
            width=2*FWHM/(1.0692 + np.sqrt(0.86639 + 4*8 * np.log(2) * self.vratio**2))
        elif self.style=='simplelorentz' or self.style=='slimlorentz':
            width=FWHM
            
        if self.style=='slimlorentz':
            if len(self.peaks)!=6:
                print "Peak detection error: slimlorentz needs 6 peaks!"
            # center199, Delta199-204, zeeman199, zeeman204, width
            self.params+=[self.peaks[1][1],self.peaks[4][1]-self.peaks[1][1],(self.peaks[2][1]-self.peaks[0][1])/2.,(self.peaks[5][1]-self.peaks[3][1])/2.,width]
        
        for peak in self.peaks:
            self.params.append(peak[0])
            if self.style=='slimlorentz':
                continue
            self.params.append(peak[1])
            self.params.append(width)
        
    def createCmdString(self):
        if self.style=='simplevoigt':
            cmdstring="pol2"
            for i in xrange(len(self.peaks)):
                #this approximation of a voigt profile with fixed height
                #is good for vratio < 1.5 with less than 2 % error   
                #seems to be invalid!!!
                cmdstring+="+[%d]*[4]/(0.637-0.255*[3]-0.998*[3]*[3]+1.78*[3]*[3]*[3]-1.15*[3]*[3]*[3]*[3]+0.265*[3]*[3]*[3]*[3]*[3])*TMath::Voigt(x-[%d],[3]*[4],[4])" % (5+i*2,6+i*2)
        elif self.style=='simplelorentz':
            cmdstring="pol2"
            for i in xrange(len(self.peaks)):
                cmdstring+="+[%d]*3.14159*[3]*0.5*TMath::CauchyDist(x,[%d],[3]*0.5)" % (4+i*2,5+i*2)
        elif self.style=='slimlorentz':
            cmdstring="pol2"
            center=["[3]-[5]","[3]","[3]+[5]","[3]+[4]-[6]","[3]+[4]","[3]+[4]+[6]"]
            for i in xrange(len(self.peaks)):
                cmdstring+="+[%d]*TMath::Pi()*[7]*0.5*TMath::CauchyDist(x," % (8+i)
                cmdstring+=center[i]
                cmdstring+=",[7]*0.5)" 
        else:
            cmdstring=Spectrum.createCmdString(self)
        return cmdstring
        
    def createParLists(self,parameters):
        if self.style=='simplevoigt' or self.style=='simplelorentz':
            ParList=list(parameters)
            #calculate the mean width
            width=[]
            i=5
            while i<len(ParList):
                width.append(ParList.pop(i))
                i+=2
            width=np.mean(width)
            if self.style=='simplevoigt':
                ParList.insert(3,self.vratio)
                ParList.insert(4,width)
            elif self.style=='simplelorentz':
                ParList.insert(3,width)
                
            ParNames=["level","drift","curvature"]        
            
            if self.style=='simplevoigt':            
                ParNames.append("vratio")
            ParNames.append("width")
            for i in xrange((len(parameters)-3)/3):
                ParNames.append("height%d" % i)  
                ParNames.append("center%d" % i)   
                
        elif self.style=='slimlorentz':
            ParList=parameters
            ParNames=["level","drift","curvature", "center199", "Delta199-204", "zeeman199", "zeeman204", "width"]
            for i in xrange(6):
                ParNames.append("height%d" % i)
        else:
            ParList, ParNames=Spectrum.createParLists(self,parameters)
        
        return ParList, ParNames
        
    def processFitResults(self, pars, errs):
        if self.style=='simplevoigt':
            numpars=3+2+2*len(self.peaks)
            self.params=[]
            self.errors=[]
            for i in xrange(numpars):
                if i==3:
                    self.vratio=pars[i]
                    continue
                elif i==4:
                    width=pars[4]
                    continue
                elif i>4:
                    j=i-2
                    self.peaks[(j-3)/2][(j-3)%2]=pars[i]
                    if i%2==0:
                        self.peaks[(j-3)/2][2]=width
                        self.params.append(pars[i])
                        self.params.append(width)
                        continue
                self.params.append(pars[i])
            
            for i in xrange(numpars):
                if i==3:
                    self.vratioErr=errs[i]
                    continue
                elif i==4:
                    widthErr=errs[4]
                    continue
                self.errors.append(errs[i])
                if i>4 and i%2==0:
                    self.errors.append(widthErr)
                    
        elif self.style=='simplelorentz':
            numpars=3+1+2*len(self.peaks)
            for i in xrange(numpars):
                if i==3:
                    width=pars[i]
                    continue
                elif i>4:
                    j=i-1
                    self.peaks[(j-3)/2][(j-3)%2]=pars[i]
                    if i%2==1:
                        self.peaks[(j-3)/2][2]=width
                        self.params.append(pars[i])
                        self.params.append(width)
                        continue
                self.params.append(pars[i])
            
            for i in xrange(numpars):
                if i==3:
                    widthErr=errs[3]
                    continue
                self.errors.append(errs[i])
                if i>3 and i%2==1:
                    self.errors.append(widthErr)
                    
        elif self.style=='slimlorentz':
            self.params=[]
            self.errors=[]
            self.peaks[1][1]=pars[3]
            self.peaks[4][1]=pars[3]+pars[4]
            self.peaks[0][1]=pars[3]-pars[5]
            self.peaks[2][1]=pars[3]+pars[5]
            self.peaks[3][1]=pars[3]+pars[4]-pars[6]
            self.peaks[5][1]=pars[3]+pars[4]+pars[6]
            for j in xrange(6):
                self.peaks[j][2]=pars[7]
                self.peaks[j][0]=pars[j+8]
            for i in xrange(14):
                self.params.append(pars[i])
                self.errors.append(errs[i])
                    
        else:
            Spectrum.processFitResults(self, pars, errs)
    def meanPeakWidth(self):
        '''the average peak width in MHz. Error calc is not exact yet!'''
        if not self.isValid():
            print 'invalid spectrum'
            return [0,0]
        width=np.take(self.params,[5,8,11,14,17,20]).mean()
        width=self.getFrequency(self.peaks[1][1]+width)
        err=np.take(self.errors,[5,8,11,14,17,20]).max()
        err=self.getFrequency(self.peaks[1][1]+err)
        err+=self.getFrequencyError(self.peaks[1][1]+width)
        return [width, err]
    def leftZeemanSplitting(self):
        if not self.isValid():
            print 'invalid spectrum'
            return [0,0]
        #calculation of errors missing
        return [self.getFrequency(self.peaks[0][1]),0]
    def rightZeemanSplitting(self):
        if not self.isValid():
            print 'invalid spectrum'
            return 0
        #calculation of errors missing
        return [self.getFrequency(self.peaks[2][1]),0]
    def leftPD_ratio(self):
        if not self.isValid():
            print 'invalid spectrum'
            return [0,0]
        l=self.peaks[0][0]
        lE=self.errors[3]
        d=self.peaks[1][0]
        dE=self.errors[6]
        return [l/d,np.sqrt((lE/l)**2+(dE/d)**2)*abs(l/d)]
    def rightPD_ratio(self):
        if not self.isValid():
            print 'invalid spectrum'
            return [0,0]
        r=self.peaks[2][0]
        rE=self.errors[9]
        d=self.peaks[1][0]
        dE=self.errors[6]
        return [r/d,np.sqrt((rE/r)**2+(dE/d)**2)*abs(r/d)]
    def lineDistance(self):
        if not self.isValid():
            print 'invalid spectrum'
            return [0,0]
        l, lErr = self.dipCenter()
        r, rErr = self.refLineCenter()
        return [r-l, np.sqrt(lErr**2+rErr**2)]
    def symmetry199(self):
        if not self.isValid():
            print 'invalid spectrum'
            return [0,0]
        l = self.params[3]
        lErr = self.errors[3]
        r = self.params[9]
        rErr = self.errors[9]
        return [l/r, np.sqrt((lErr/l)**2+(rErr/r)**2)*abs(l/r)]
    def leftPeakCenter(self):
        return [self.params[4],self.errors[4]]
    def dipCenter(self):
        return [self.params[7],self.errors[7]]
    def rightPeakCenter(self):
        return [self.params[10],self.errors[10]]
    def leftRefCenter(self):
        return [self.params[13],self.errors[13]]
    def refLineCenter(self):
        return [self.params[16],self.errors[16]]
    def rightRefCenter(self):
        return [self.params[19],self.errors[19]] 
     
class DBSpec(Spectrum):
    def __init__(self,data):
        Spectrum.__init__(self, data, style='gauss')
    def isValid(self):
        if not len(self.peaks)==5:
            return False
        return True
    def getFrequency(self,x):
        '''return the associated frequency relative to the Hg-199 
        transition in MHz for an x-value,
        derived from the distance between the isotope lines'''
        # need to lookup frequency differences first        
        pass
    def getHg199Peak(self):
        return self.peaks[0][1]
    
#data=np.load('25_11_peakdipratio_mediumfocus.npy')
#S=DFSpec(data)

#--------------------UTILITY FUNCTIONS--------------------

def removePeaks(data):
    pts=len(data[0])
    y=np.array(data[1])
    enhanced=enhancedDerivative(data)
    thresh=np.abs(enhanced).max()*0.1
    overcount=0
    Yhist=[]
    lastY=0
    for i in xrange(pts):
        if enhanced[i]>thresh:
            y[i]=lastY
            overcount+=3
        elif enhanced[i]<-thresh:
            y[i]=lastY
            overcount-=3
        #make sure the last 5% of the signal are kept
        # - if there is a peak the signal is not useful anyway
        elif i/float(pts) > 0.95:
            break
        elif overcount:
            y[i]=lastY
            overcount=np.sign(overcount)*(abs(overcount)-1)
        else:
            Yhist.append(y[i])
            lastY=np.mean(Yhist[-pts/50:])
    return np.array([data[0],y])

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

def findPeaks(data,level,drift,curvature=0,show=False,win_len=31):
    ''' will return position of the peaks (zeros of first derivative) and their height and width.
    level and drift of the signal help identify the peaks.'''
    x=np.array(data[0])
    y=np.array(data[1])-(level+drift*x+curvature*x**2)
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
        #check if the peak has the right curvature for its sign and is not to flat
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
    T2gr.SetLineColor( 1 )
    T2gr.SetLineWidth( 1 )
    T2gr.SetMarkerColor( 4 )
    T2gr.SetMarkerStyle( 21 )
    T2gr.SetMarkerSize(0.3)
    T2gr.SetTitle( 'Fit' )
    T2gr.GetXaxis().SetTitle( 'Scan Voltage' )
    T2gr.GetYaxis().SetTitle( 'Signal Voltage' )
    return T2gr

# the next two functions are old and no longer used
def guessMeanLevel(data):
    '''this algorithm is based on the assumption, that the most frequent 
    values in the signal are belonging to the level'''
    #calculate the stretching factor so we obtain a good statistic in bincount
    #experience value: 100 for a range of 0 to 1        
    factor=100/(np.max(data[1])-np.min(data[1]))
    centiV=np.array(data[1]*factor,dtype=int)
    #shift centiV to positive values so we can use bincount
    centiVmin=centiV.min()
    centiV-=centiVmin
    freq=np.bincount(centiV)
    maxfreq=[]
    for i in xrange(7):
        argmax=freq.argmax()
        maxfreq.append(argmax)
        freq[argmax]=0
    mean=np.array(maxfreq).mean()
    level=(mean+centiVmin)/factor
    return level
    
def guessLevelDrift(data, steps=5):
    '''guess level in 'steps' parts of the data and do a linear fit'''
    centers=[]
    levels=[]
    points=len(data[0])        
    for i in xrange(steps):
        lowbound=i*int(points/steps)
        upbound=(i+1)*int(points/steps)
        region=data[:,lowbound:upbound]
        centers.append(np.mean(region[0]))
        levels.append(guessMeanLevel(region))
    drift, level= np.polyfit(centers,levels,1)
    return level,drift
    
#create a sample spectrum
d=np.loadtxt('timenew17.txt')
l=np.linspace(0,1,1000)
s=DFSpec([d[3],d[1]],'slimlorentz')
s.fit(verbose=True,showgraphs=True)