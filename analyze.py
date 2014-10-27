from spectrum import SpecFitter
from fitmodels import Lorentzian, Hg199_204, Hg199_204_Voigt, Hg199_204_Voigt2, Constant, Linear, Quadratic
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['figure.subplot.hspace']=0.5
mpl.rcParams['backend']='TkAgg'
import Tkinter as tk
import tkFileDialog as tkf
import re
import os.path
import sys
from ROOT import TGraphErrors, TCanvas, TH1F
#from uncertainties import ufloat
#from uncertainties.umath import sqrt

#an example function, if the spectrum class doesn't contain a 
#needed method
def LockPoint(spec):
    return [spec.FModel.getYLockPoint(-8.8), spec.FModel.getConfidenceInterval(spec.FModel.getXLockPoint(-8.8))]

def LockPointMinusLevel(spec):
    return [spec.FModel.getYLockPoint(-8.8)-spec.BModel.params[0], spec.FModel.getConfidenceInterval(spec.FModel.getXLockPoint(-8.8))]

def yNoiseToMHz(spec):
    return spec.FModel.getRatioTo199_204(1)/spec.FModel.getLockPointSlope(-8.8)

def Chi2(spec):
    return spec.FModel.ChiSquareRed

def NormProbChi2(spec):
    return spec.NormProbPltChi2

def pValue(spec):
    return spec.FModel.pValue

def MaxCorrelation(spec):
    return spec.FModel.MaxCorrelation

def AndersonDarling(spec):
    return spec.getADTestPValue()

def KolmogorovSmirnov(spec):
    return spec.getKSTestPValue()

def Noise(spec):
    return spec.noise
    
def SignalToNoise(spec):
    return spec.signalToNoise()/spec.FModel.ChiSquareRed

#a list with the properties of interest. entries are from the DFSpec methods
#or other supported functions
#supported functions have to take one DFSpec argument and return either a value or
#a list of the form [value, error]
global prop
prop=[LockPoint]#LockPoint, yNoiseToMHz, Chi2, AndersonDarling, KolmogorovSmirnov, MaxCorrelation, Noise, SignalToNoise]

global fitmodel
global bmodel
fitmodel=Hg199_204
bmodel=Linear

global noise
#noise=0.025616899310231329
noise=None

#if multiple measurement points shall be averaged for one property point
#set this to something greater than 1
numPoints=0

#prevent plt from blocking
#plt.ion()

#------------------UTILITY FUNCTIONS---------------------------
def analyzation(dateien):
    #list with names of the fitted files
    files=[]
    
    #list with names of the evaluated parameters
    parNames=[]    
    
    #list with names of the evaluated properties
    propNames=[]
    for p in prop:
        propNames.append(p.__name__)    
    
    #list containing lists of every spectrums parameters
    params=[]
    #list containing lists of every spectrums parameters' errors
    paramErrs=[]
    #list with lists of every evaluated poperty's values
    propVals=[]
    #list containing errors of the evaluated properties
    propErrs=[]
    
    #list containing the Spectrum objects
    specs=[]
    
    for i in prop:
        propVals.append([])
        propErrs.append([])
    
    #index for deciding when to begin a new plot
    j=0
    
    for datei in dateien:
        print datei
        if datei.endswith('.npy'):
            data=np.load(datei)
        elif datei.endswith('.csv') or datei.endswith('.txt'):
            try:
                data=np.loadtxt(datei)
            except:
                print 'could not open file '+datei
                continue
        else:
            print 'unknown format '+ datei
            continue
        s=SpecFitter([data[0],data[1]], noise, numpeaks=6, FModel=fitmodel, BModel=bmodel)
        if not s.isAlive():
            print 'invalid spectrum, skipped: '+datei[-15:]
            continue
        s.fit()
        
#        xvals=[s.FModel.getRatioTo199_204(data[3][i]-s.FModel.params[0]) for i in xrange(len(data[3]))]        
#        
#        #rescale the x-axis to frequency
#        s=SpecFitter([xvals,data[1]], noise, numpeaks=6, FModel=fitmodel, BModel=bmodel)
#        if not s.isAlive():
#            print 'invalid spectrum, skipped: '+datei[-15:]
#            continue
#        s.fit()
    
        if not s.FModel.isValid(): #or not s.residualsWithin5Sigma() or not s.getADTestPValue()>0.05:
            print 'invalid spectrum, skipped: '+datei[-15:]
            continue
        
        files.append('...'+datei[-15:-4])
        specs.append(s)
        N=len(s.FModel.params)
        
        if j==0:
            for i in xrange(N):
                parNames.append(s.FModel.TF1.GetParName(i))
                params.append([])
                paramErrs.append([])
            s.draw()
                
        for i in xrange(N):
            params[i].append(s.FModel.params[i])
            paramErrs[i].append(s.FModel.errors[i])
        
        #now get the desired properties
        for i in xrange(len(prop)):
            #prop[i] is a function or class method taking s as an argument
            r=prop[i](s)
            #maybe prop[i] returns a tuple of value and error or only one value without error
            if type(r)==list or type(r) == tuple:
                propVals[i].append(r[0])
                propErrs[i].append(r[1])
            else:
                #there is no error for this value so set it to zero
                propVals[i].append(r)
                propErrs[i].append(0)
                
        #create a new figure for the spectra if the current is full
#        if j%20 == 0:
#            plt.figure(figsize=(20,20))
#        plt.subplot(5,4,j%20+1)
#        s.plot()
#        plt.draw()
        #s.draw()
        j+=1

    #get the folder of the files
    folder=os.path.dirname(dateien[0])+'/'
    filename=folder+'_'.join([str(x) for x in time.localtime()[:6]])+'.log'
    #write the results to file
    saveLog(files, parNames, propNames, params, paramErrs, propVals, propErrs, filename)

    return files, parNames, propNames, params, paramErrs, propVals, propErrs
    
def saveLog(files, parNames, propNames, params, paramErrs, propVals, propErrs, filename=None, listsep='\n', valsep=' '):
    ''' this function saves the data contained in the arguments. files and propnames must be lists of values,
    the other arguments lists of lists of values. listsep and valsep are the seperation characters for lists and values.
    if no filename is given, the date and time is taken.'''
    #header for the logfile.
    log='This is a logfile generated by analyze.py - lines containing the following: 1. filenames of the fitted spectra 2. names of the evaluated parameters 3. fitparameters of the spectra (one line for every parameter, in order of the filenames, this applies to the rest of the file) 4. errors of the fitparameters 5. evaluated properties 6. errors of the properties\n'
    #the lists are seperated by listsep, the values are seperated by valsep    
    for x in files, parNames, propNames:    
        log+=valsep.join([str(y) for y in x])
        log+=listsep
    
    for x in params, paramErrs, propVals, propErrs:
        for y in x:
            log+=valsep.join([str(z) for z in y])
            log+=listsep
    
    #remove the last listsep
    log=log[:-(len(listsep))]
    
    if not filename:
        filename='_'.join([str(x) for x in time.localtime()[:6]])+'.log'
    
    f=open(filename,'w')
    f.write(log)
    f.close()

def loadLog(filename, listsep='\n', valsep=' '):
    f=open(filename, 'r')
    header=f.readline()
    print header
    raw_data=f.read()
    lists=[x.split(valsep) for x in raw_data.split(listsep)]
    files, parNames, propNames = lists[:3]
    numprops=len(propNames)
    numparams=len(parNames)
    params=[[float(x) for x in l] for l in lists[3:numparams+3]]
    paramErrs=[[float(x) for x in l] for l in lists[numparams+3:2*numparams+3]]
    propVals=[[float(x) for x in l] for l in lists[-2*numprops:-numprops]]
    propErrs=[[float(x) for x in l] for l in lists[-numprops:]]
    
    return files, parNames, propNames, params, paramErrs, propVals, propErrs
    

#-----------------PROGRAM START------------------------

root=tk.Tk()
root.title("data analysis")

if sys.platform=='win32':
    dateien=unicode(tkf.askopenfilenames())
    dateien=re.findall('{.*?}',dateien)
    root.destroy()
     
    for i in xrange(len(dateien)):
        dateien[i]=dateien[i].replace('{','')
        dateien[i]=dateien[i].replace('}','')
else:
    dateien=list(tkf.askopenfilenames())

root.destroy()

print dateien

if len(dateien) == 1 and dateien[0].endswith('.log'):
    files, parNames, propNames, params, paramErrs, propVals, propErrs = loadLog(dateien[0])
else:
    files, parNames, propNames, params, paramErrs, propVals, propErrs = analyzation(dateien)


#this part is for the case multiple measurements have been performed
#with the same parameters -> set numpoints to this value
#if numPoints>1:
#    newpropVals=[]
#    newpropErrs=[]
#    #construct an array with values and errors using the uncertainties package
#    #valNerr=[[ufloat(propVals[i][j],propErrs[i][j]) for j in xrange(len(propVals[i]))] for i in xrange(len(propVals))]
#    valNerr=np.array(valNerr)
#    for line in valNerr:
#        valtemp=[]
#        errtemp=[]
#        for i in xrange(len(line)/numPoints):
#            #temporary, this should be checked for correct propagation of errors!        
#            l=line[numPoints*i:numPoints*(i+1)]            
#            mean=l.mean()
#            std=sqrt(l.var())
#            valtemp.append(mean.n)
#            errtemp.append(std.n)
#        newpropVals.append(valtemp)
#        newpropErrs.append(errtemp)
#    propVals=newpropVals
#    propErrs=newpropErrs


#create a canvas for the parameters and errors
j=0
canvas=[]
g=[]
for (names, p, pErrs) in (([parNames[0],parNames[1],parNames[11]], [params[0],params[1],params[11]], [paramErrs[0],paramErrs[1],paramErrs[11]]),(propNames, propVals, propErrs)):
    canvas.append(TCanvas("Analyzation%d"%(j),"Analyzation"))
    n=int(np.ceil(np.sqrt(len(p))))
    canvas[j].Divide(n,n)
    g.append([])
    #plot the parameters
    for i in xrange(len(p)):
        canvas[j].cd(i+1)
        canvas[j].SetGridx()
        canvas[j].SetGridy()
        n=len(p[i])
        x=np.arange(len(p[i]), dtype=np.float64)
        y=np.array(p[i], dtype=np.float64)
        xnoise=np.zeros(len(p[i]), dtype=np.float64)
        ynoise=np.array(pErrs[i], dtype=np.float64)
        g[j].append(TGraphErrors(n,x,y,xnoise,ynoise))
        g[j][i].SetMarkerStyle(4)
        g[j][i].SetMarkerSize(0.4)
        g[j][i].GetXaxis().SetLabelSize(0.03)
        g[j][i].GetXaxis().SetTitleSize(0.07)
        g[j][i].GetXaxis().SetTitle("measurement")
        g[j][i].GetXaxis().CenterTitle()
        g[j][i].GetXaxis().SetLabelFont(132)
        g[j][i].GetXaxis().SetTitleFont(132)
        g[j][i].GetXaxis().SetRangeUser(-1,len(p[i]))
        g[j][i].GetYaxis().SetLabelSize(0.03)
        g[j][i].GetYaxis().SetTitleSize(0.07)
        g[j][i].GetYaxis().SetTitle(names[i])
        g[j][i].GetYaxis().CenterTitle()
        g[j][i].GetYaxis().SetLabelFont(132)
        g[j][i].GetYaxis().SetTitleFont(132)        
        #set the function name as plot title
        g[j][i].SetTitle(names[i])
        #ignore error bars and be quiet
        g[j][i].Fit('pol1','WQ')
        g[j][i].Draw("ap")
    j+=1

raw_input("Done, press any key!")
new=TCanvas('new')
t, s = np.polyfit(range(len(propVals[0])),propVals[0],deg=1)
res=[propVals[0][i]-(s+t*i) for i in xrange(len(propVals[0]))]
np.savetxt('histogram_data.dat',np.transpose(res))
limit=np.std(res)*3
hist=TH1F('residuals','residuals',20,-limit,limit)
for i in xrange(len(res)):
    hist.Fill(res[i])
hist.Draw()

#for c in canvas:
#    if c:
#        c.IsA().Destructor(c)
