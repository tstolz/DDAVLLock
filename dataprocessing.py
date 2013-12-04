# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 11:45:32 2012

@author: bernd
"""
#todo: Noise estimation verbessern mit smoothen und dann noise ausrechnen
# fit quality parameter ausgeben damit man weiss ob alle fits sinnvoll sind ?DoNe?
import digilock
import math
import ROOT
from ROOT import gROOT, TCanvas, TGraph, gStyle, TMath
from array import array
import numpy

def calib(voltage):
    return (voltage * 1248.7367556174609)
def import_data(source):
    """return spectrum data imported from digilock or file
    give "digilock" as argument to import live data from digilock
    """
    data=[]
    if source=="digilock":
        imported=digilock.getscopedata()
        data.append(array("d",imported[2]))
        data.append(array("d",imported[1]))
    else:
        imported=import_spectradata_from_file(source)
        data.append(array("d", map(calib, imported[2])))
        data.append(array("d", imported[1]))
    return data
def guess_extrema(data,threshold_max=0.1,threshold_min=-0.1):
    """find all the extrema in a spectrum and return them as [maxima,minima]
    please note: thresholds are overwritten in the function!!!!
    """
    print data[1] 
    zero=mean([max(data[1]),min(data[1])])
    threshold_max=(max(data[1])-zero)/3+zero
    threshold_min=(min(data[1])-zero)/3+zero
    print threshold_max, threshold_min
    overthresh=numpy.where(numpy.array(data[1])>threshold_max)[0].tolist()
    underthresh=numpy.where(numpy.array(data[1])<threshold_min)[0].tolist()
    sep_overthresh=seperate_extrema(overthresh)
    sep_underthresh=seperate_extrema(underthresh)
    maxima=isolate_extrema(sep_overthresh, data)
    minima=isolate_extrema(sep_underthresh, data)
    return [maxima,minima]
def isolate_extrema(seperated_extrema_list,data):
    """find and return the extrema (index,x,y) of each list when given a list
    with one list of consecutive indices for each extremum
    !!!*****returns positive values for negative minima!!!!!****!!!
    """
    extrema=[]
    extremapos=[]
    extremax=[]
    extremay=[]
    for no_max in range(len(seperated_extrema_list)):
        extrema.append(max([abs(data[1][i]) for i in seperated_extrema_list[no_max]]))
    for no_max in range(len(seperated_extrema_list)):
        extremapos.append(map(abs,data[1]).index(extrema[no_max]))
    [extremax.append(data[0][i]) for i in extremapos]
    [extremay.append(data[1][i]) for i in extremapos]
#   for n in range(len(extrema)):
 #       extremay.append(seperated_extrema_list.find())    
    print extrema
    print extremapos
    print extremax
    extrempoint=zip(*[extremapos,extremax,extremay])
    return extrempoint
def seperate_extrema(outsidethresh):
    """find gaps in the list of indices where threshold is exceeded,
    and seperate the list at the positions of the gaps.
    
    """
    
    m=[]#separation marker
    for i in range(len(outsidethresh)-1):
        if outsidethresh[i+1]-outsidethresh[i]>5:
            m.append(i+1)
    extrema=partition(outsidethresh,m)
    return extrema
     
    
    
def picasso(data):
    #gROOT.Reset()
    c1 = TCanvas( 'c1', 'ScanData', 200, 10, 700, 500 )
    #c1.SetFillColor( 42 )
    c1.SetGrid()
    gROOT.SetStyle("Plain")
    gr, results, parnames=fit_errorsig(data)
    gr.Draw('AL')
    c1.Update()
    raw_input("euuft")


def pygaus( x, par ):
    arg1 = 0
    scale1 =0
    ddx = 0.01  
    
    if (par[2] != 0.0):
    #    print par
    #    print x
        arg1 = (x-float(par[1]))/float(par[2])
        scale1 = (ddx*0.39894228)/par[2]
        h1 = par[0]/(1+par[3])

        gauss = h1*scale1*math.exp(-0.5*arg1*arg1)
    else:
        gauss = 0.
    return gauss
def pygaus2(x,  par):
	""" definition of the guassian function
        par[0]=sigma
        par[1]=x offset
        par[2]=Amplitude
        """
	y = par[2] / ( math.sqrt(2. * math.pi) * par[0] )  * math.exp(-0.5 * ((x - par[1] + 0.0) / par[0])**2 )
	return y
def pylorentz(x, par):
    """ lorenzian function
        par[0]=gamma (fwhm)
        par[1]=x offset
        par[2]=amplitude
        The real voltage amplitude is 2*par[2]/(par[0]*Pi)
        """
    y=par[2] * par[0] / (2 * math.pi) / ((x-par[1])**2 + (par[0]/2)**2)
    return y
def errorsig(x,par):#when this function is changed, also change parnames in fit_errorsig
    #print list(par)[0]
    pars=list(par)
    y=pylorentz(x[0],pars[0:3])+pylorentz(x[0],pars[3:6])+pars[6]
    return y
def fit_errorsig(data, parameter=[0.001,0.05,0.002,0.001,0.055,-0.002,0.01], rangemin=0.03, rangemax=0.07, showgraphs=False):
    ParNames=["width1","x_offset1","amplitude1","width2","x_offset2","amplitude2","y offset"]
    gr=create_TGraph(data)
    fit1=ROOT.TF1("fit1",errorsig, min(data[0]), max(data[0]), 7)
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
    gr.Fit("fit1", "Q+", "LEGO", rangemin, rangemax)
    pars=fit1.GetParameters()
    params=[pars[0],pars[1],pars[2]*2/pars[0]/math.pi,pars[3],pars[4],pars[5]*2/pars[3]/math.pi,pars[6]]
    #print params
    params.append((abs(params[2])+abs(params[5]))/(2*noisefit(data)))
    params.append(fit1.GetChisquare())
    #append SNR. Mean of amplitude1 and amplitude2 is used für Signal height
    ParNames.append("SNR")
    ParNames.append("Chisquare")
    return gr,params,ParNames
    raw_input("after fit, press any key")
        
    #gr.Fit(fit1, "Q+", "LEGO", 0.03, 0.07)

    #res=list(fit1.GetParameters())
    #print res
    
#def fit_voigt(data):
#    gr=create_TGraph(data)
#    fitfunc=TMath
#    f1=ROOT.TF1("f1", "TMath.voigt(x)-TMath.voigt(x)",0.03,0.07)    
#    gr.Fit(f1,"Q+","LEGO", 0.03, 0.07)
def create_TGraph(data): #source="digilock" to import live data, and "filename" for file
#veraltet, wegen des datenimports
#    if source=="digilock":
#        times=array("d",digilock.getscopedata()[0])
#        data2=array("d",digilock.getscopedata()[1])
#    else:
#        times=array("d", import_spectradata_from_file(source)[0])
#        data2=array("d", import_spectradata_from_file(source)[1])
    #data1=array("d",[1,2,3,4,5,6])
    #data2=array("d",[2,1,2,3,2,1])
    T2gr=ROOT.TGraph(len(data[0]),data[0],data[1])
    T2gr.SetLineColor( 2 )
    T2gr.SetLineWidth( 2 )
    T2gr.SetMarkerColor( 4 )
    T2gr.SetMarkerStyle( 21 )
    T2gr.SetMarkerSize(0.1)
    T2gr.SetTitle( 'T2=[0]' )
    T2gr.GetXaxis().SetTitle( 'seconds' )
    T2gr.GetYaxis().SetTitle( 'Voltage' )
#    T2gr.Draw( 'AL' )
   #raw_input("press any key.")
    return T2gr
def import_spectradata_from_file(filename):
    file=open(filename, "r")
    scope_data=file.read()
   # del scope_data[1]
    print scope_data
    scope_data_list=scope_data.split(chr(13)+chr(10))
    print scope_data_list
    del scope_data_list[0]
    scope_data_splitted=map(digilock.ssplit,scope_data_list)
    del scope_data_splitted[-1]
    scope_data_splitted_float=[map(float,x) for x in scope_data_splitted]#[:-1]
    scope_data_transposed=zip(*scope_data_splitted_float)
    print scope_data_transposed
    return scope_data_transposed
def mean(numberList):
    if len(numberList) == 0:
        return float('nan')
 
    floatNums = [float(x) for x in numberList]
    return sum(floatNums) / len(numberList)
def partition(l, indexes):
    result, indexes = [], indexes+[len(l)]
    reduce(lambda x, y: result.append(l[x:y]) or y, indexes, 0)
    return result
def fit_all(data, showgraphs=False):
    fitresults=[]
    parameter=[1,0.05,0.002,1,0.055,-0.002,0.01]
    #data=dataprocessing.import_data(source)
    extrema=guess_extrema(data,0.1,-0.1)
    for i in range(len(extrema[0])):
        extrpos=[extrema[1][i][1],extrema[0][i][1]]        
        fitrangemin=3*min(extrpos)-2*max(extrpos)
        fitrangemax=3*max(extrpos)-2*min(extrpos)
        parameter[1]=extrema[0][i][1]
        parameter[2]=extrema[0][i][2]
        parameter[4]=extrema[1][i][1]
        parameter[5]=extrema[1][i][2]
       # print parameter
        gr, fitresult, parnames=fit_errorsig(data, parameter, fitrangemin, fitrangemax, showgraphs)
       # gr=result["gr"]
        if showgraphs:
            c1 = TCanvas( 'c1', 'ScanData', 200, 10, 700, 500 )
            c1.SetGrid()
            gROOT.SetStyle("Plain")
            gr.Draw('AL')
            c1.Update()
            raw_input("euuft")
        fitresults.append(fitresult)
    return fitresults, parnames
def noisefit(data):
    """find the rms noise in the the first 1/100 of data"""
    noiseregion=data[1][:len(data[1])/100]
    x=data[0][:len(data[1])/100]
    var=math.sqrt(mean([noiseregion[i]**2 for i in range(len(noiseregion))])-(mean(noiseregion))**2)
    print var
    return var
    
    