# -*- coding: utf-8 -*-
"""
Created on Wed Sep 12 12:27:39 2012

@author: bernd
"""

import dataprocessing
import digilock
from ROOT import gROOT, TCanvas, TGraph, gStyle, TMath

def fit_all(data):
    fitresults=[]
    parameter=[0.0005,0.05,0.002,0.0005,0.055,-0.002,0.01]
    #data=dataprocessing.import_data(source)
    extrema=dataprocessing.guess_extrema(data,0.2,-0.2)
    for i in range(len(extrema[0])):
        extrpos=[extrema[1][i][1],extrema[0][i][1]]        
        fitrangemin=3*min(extrpos)-2*max(extrpos)
        fitrangemax=3*max(extrpos)-2*min(extrpos)
        parameter[1]=extrema[0][i][1]
        parameter[2]=extrema[0][i][2]*0.0047
        parameter[4]=extrema[1][i][1]
        parameter[5]=extrema[1][i][2]*0.0047
        c1 = TCanvas( 'c1', 'ScanData', 200, 10, 700, 500 )
        c1.SetGrid()
        gROOT.SetStyle("Plain")
       # print parameter
        gr, fitresult=dataprocessing.fit_errorsig(data, parameter, fitrangemin, fitrangemax)
       # gr=result["gr"]
        gr.Draw('AL')
        c1.Update()
        raw_input("euuft")
        fitresults.append(fitresult)
    print str(fitresults)
        
        