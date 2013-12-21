# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 14:57:32 2013

@author: universe
"""

import numpy as np
import matplotlib.pyplot as plt
import time
from ROOT import TMath, TGraphErrors, TF1

#voigt is not yet working, problem with amplitude!   
class Spectrum:
    def __init__(self, data, style="lorentz"):
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
        self.vratio=0.5
        #indicator if the spectrum has been fitted yet
        self.fitted=False
        self.exact=False
    def setStyle(self,style):
        if not style in ["lorentz","gauss","voigt"]:
            print "unknown style"
        else:
            self.style=style
    def theory(self, x, params):
        '''params must have the form "level, drift, curvature, height1, center1, width1, ... heightN, centerN, widthN'''
        level = params[0]
        drift = params[1]
        curv = params[2]
        peaks = [[params[i],params[i+1],params[i+2]] for i in range(3,len(params),3)]
        result = level + x*drift + x**2 * curv
        for peak in peaks:
            if self.style=="gauss":
                result += peak[0]*TMath.Gaus(x,peak[1],peak[2])
            elif self.style=="lorentz":
                result += peak[0]*TMath.CauchyDist(x,peak[1],peak[2])
            elif self.style=="voigt":
                result += peak[0]*TMath.Voigt(x-peak[1],self.vratio*peak[2],(1-self.vratio)*peak[2])
        return result
    def noisefit1(self):
        '''subtract the fitted signal from data and calc noiselevel of the rest'''
        if self.params:
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
        if self.params:
            #params=np.array(np.array(self.params[:-1])[:,1],dtype=float)
            Fit=[self.theory(x,self.params) for x in self.data[0]]
            plt.plot(self.data[0],self.data[1],'k,',self.data[0],Fit,'r-')
        else:
            plt.plot(self.data[0],self.data[1],'k,')
    def guessParams(self,numPeaks=None, show=False):
        '''guess all the parameters required for the fit. order: level, drift, curvature, (height, center, width) of all peaks'''
        parameters=[]
        #level, drift=self.guessLevelDrift()
        noPeaks=removePeaks(self.data)

        curvature, drift, level = np.polyfit(*noPeaks,deg=2)
        
        #estimate the best window lenght for "enhancedDerivative()"
        #from the expected number of peaks
        #calculation based on experience
        if not numPeaks:
            win_len=31
        else:
            win_len=self.points/(numPeaks*5)
        self.peaks=findPeaks(self.data,level,drift,curvature,show=show,win_len=win_len)
        #level is assumed to return the central signal level
        #we need the level on the left
        parameters+=list([level,drift,curvature])
        for peak in self.peaks:
            #the height and width (FWHM) has to be adjusted to the form of the curve
            if self.style=='gauss':
                width=peak[2]/(2*np.sqrt(2*np.log(2)))
                height=peak[0]
            elif self.style=='lorentz':
                width=peak[2]/2
                height=peak[0]*np.pi*width
            else:
                #to be filled for voigt
                height=peak[0]
                width=peak[2]
                
            parameters.append(height)
            parameters.append(peak[1])
            parameters.append(width)
#        parameters+=list([rng*0.02]*5)
#        parameters+=map(lambda i: level+drift*peaks[i,0]-peaks[i,1],xrange(len(peaks)))
#        parameters+=self.pk_centers
        self.params=parameters
        self.noisefit1()
        return parameters
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
    def fit(self, parameters=None, showgraphs=False, verbose=False):
        '''data[0] must contain x-, data[1] y-axis data, parameters a guess of the 
        parameters in broadSpecTheory'''
        data=np.asarray(self.data,dtype=np.float64)
        numpoints=len(data[0])
        x_range=max(data[0])-min(data[0])
        #the error estimate for the adc time values is half the stepwidth
        xerr=x_range/numpoints*0.5
        if not parameters:
            parameters=self.guessParams()      
            
        if self.style=="voigt":
            parameters.insert(3,self.vratio)

        ParNames=["level","drift","curvature"]        
        
        if self.style=="voigt":
            ParNames.append("vratio")
        for i in xrange((len(parameters)-2)/3):
            ParNames.append("height%d" % i)  
            ParNames.append("center%d" % i)
            ParNames.append("width%d" % i)
            
        gr=create_TGraph(data, xerr, self.noise_level)
        
        # insert signal dependent string into the TF1
        # first create the string
        cmdstring="pol2"
        
        for i in xrange(len(self.peaks)):
            if self.style=="gauss":
                cmdstring+="+gaus(%d)" % (3+i*3)
            elif self.style=="lorentz":
                cmdstring+="+[%d]*TMath::CauchyDist(x,[%d],[%d])" % (3+i*3,4+i*3,5+i*3)
            elif self.style=="voigt":
                cmdstring+="+[%d]*TMath::Voigt(x-[%d],[3]*[%d],(1-[3])*[%d])" % (4+i*3,5+i*3,6+i*3,6+i*3)
        fit1=TF1("fit1",cmdstring, min(data[0]), max(data[0]))
        for i in range(len(ParNames)):
            fit1.SetParName(i, ParNames[i])
        for n in range(len(parameters)):
            fit1.SetParameter(n,parameters[n])
        #fit1.SetParLimits(1,0.045,0.05)
        #fit1.SetParLimits(4,0.052,0.06)
        if showgraphs:
            gr.Draw()
            #print parameter
            raw_input("pre fit, press any key")
            print 'fitting...'
        tm=time.time()
        gr.Fit(fit1, "Q","", data[0].min(), data[0].max())
        
        if verbose:
            print "Fitting time: %f s" % (time.time()-tm)
        
        self.params=[]
        pars=fit1.GetParameters()
        for i in xrange(len(parameters)):
            if i==3 and self.style=="voigt":
                self.vratio=pars[i]
                continue
            self.params.append(pars[i])
            if i>2:
                self.peaks[(i-3)/3][i%3]=pars[i]
                
        errs=fit1.GetParErrors()
        self.errors=[]
        for i in xrange(len(parameters)):
            if i==3 and self.style=="voigt":
                continue
            self.errors.append(errs[i])
        
        self.ChiSquareRed=fit1.GetChisquare()/float(numpoints)
        if verbose:
            print "Fit converged with reduced chi-squared: %f" % (self.ChiSquareRed)
        self.fitted=True        
        if showgraphs:        
            fit1.Draw()
        return fit1

#------------------------Subclasses-------------------------------

class DFSpec(Spectrum):
    def getFrequency(self,x):
        '''return the associated frequency relative to the Hg-199 
        transition in MHz for an x-value,
        derived from the distance between the isotope lines'''
        #TO BE LOOKED UP! This is a rough estimate.
        if not self.isValid():
            print 'invalid spectrum'
            return 0
        freq_diff=70
        scale=freq_diff/(self.peaks[4][1]-self.peaks[1][1])
        return (x-self.peaks[1][1])*scale
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
    def getFrequencyError(self, x):
        return 0
        #to be filled! composed of uncertainty of freq_diff and peak errors
    def isValid(self):
        if len(self.peaks)!=6:
            return False
        if self.exact and abs(self.ChiSquareRed-1) > 0.3:
            return False
        if self.fitted and self.ChiSquareRed>30:
            return False
        #further possibilities for the spectrum to be invalid
        return True
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

def findDips(data,win_len=20,show=False):
    ''' will return position of the peaks (moving window method) - 
    peaks are actually dips'''
    x=np.array(data[0])
    y=np.array(data[1])
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
    return PEAKS

def create_TGraph(data,xerr,yerr):
    n=len(data[0])    
    x=np.asarray(data[0],dtype=np.float64)
    y=np.asarray(data[1],dtype=np.float64)
    xnoise=np.asarray([xerr]*n,dtype=np.float64)
    ynoise=np.asarray([yerr]*n,dtype=np.float64)
    T2gr=TGraphErrors(n,x,y,xnoise,ynoise)
    T2gr.SetLineColor( 2 )
    T2gr.SetLineWidth( 2 )
    T2gr.SetMarkerColor( 4 )
    T2gr.SetMarkerStyle( 21 )
    T2gr.SetMarkerSize(0.1)
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