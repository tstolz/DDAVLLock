# -*- coding: utf-8 -*-
"""
Created on Wed Sep 12 18:08:44 2012

@author: bernd

für jede messreihe einen eigenen ordner anlegen und den sich
verändernden parameter als dateinamen verwenden
"""

import dataprocessing
import os
from os.path import join, getsize
import ROOT
from array import array
from ROOT import gROOT, TCanvas, TGraph, gStyle, TMath

def fit_all_files_in(path="/home/bernd/Dropbox/work/Hg (co)magnetometer/laser setup/HgSpektren/CoilCurrent", showgraphs=False):
    parameters=list()
    results=list()
    directory = os.path.join(path)
    for files in os.walk(directory):
        for file in files[2]:
            if file.endswith(".logo") or file.endswith(".txt"):
                param=float(file[0:file.find(".",0,len(file))].replace("_","."))
                parameters.append(param)
                print path+"/"+file
                data=dataprocessing.import_data(path+"/"+file)
                result, parnames=dataprocessing.fit_all(data, showgraphs)
                results.append(result)
                print parameters, results
    return [parameters,results], parnames


def create_param_graphs(path="/home/bernd/Dropbox/work/Hg (co)magnetometer/laser setup/HgSpektren/CoilCurrent"):
    results, parnames=fit_all_files_in(path)
    graphs=list()
    for n in range(len(parnames)):
        for m in range(len(results[1][1])):
            y=[float(results[1][i][m][n]) for i in range(len(results[1]))]
            T2gr=ROOT.TGraph(len(results[0]), array("d", results[0]), array("d", y))
            T2gr.SetLineColor( 2 )
            T2gr.SetLineWidth( 2 )
            T2gr.SetMarkerColor( 4 )
            T2gr.SetMarkerStyle( 21 )
            T2gr.SetTitle( 'T2=[0]' )
            T2gr.GetXaxis().SetTitle( 'changed param' )
            T2gr.GetYaxis().SetTitle( 'fitted param' )
            save_param_graph(T2gr, path, parnames[n]+"_signal"+str(m))
            T2gr.SaveAs(path+"/"+parnames[n]+"_signal"+str(m)+".root")
            graphs.append([[n],results[0],y])
    return T2gr, parnames, graphs
    raw_input("yo!")
    
def get_zeeman_splitting(graphs, path="/home/bernd/Dropbox/work/Hg (co)magnetometer/laser setup/HgSpektren/CoilCurrent"):
    T2grs=[]
    for g in range(len(graphs)/7):
        gincr=g
        if len(graphs)/7 == 1:
            diff=[graphs[1][2][i]-graphs[4][2][i] for i in range(len(graphs[1][2]))]
        else:
            diff=[graphs[1*len(graphs)/7+gincr ][2][i]-graphs[4*len(graphs)/7+gincr][2][i] for i in range(len(graphs[1][2]))]
        T2gr=ROOT.TGraph(len(graphs[1][1]), array("d", graphs[1][1]), array("d", diff))
        T2gr.SetLineColor( 2 )
        T2gr.SetLineWidth( 2 )
        T2gr.SetMarkerColor( 4 )
        T2gr.SetMarkerStyle( 21 )
        T2gr.SetTitle( 'T2=[0]' )
        T2gr.GetXaxis().SetTitle( 'CoilCurrent' )
        T2gr.GetYaxis().SetTitle( 'ZeemanSplittint' )
        save_param_graph(T2gr, path+"signal"+str(g), "zeemanSplitt")
        T2gr.SaveAs(path+"/"+"Zeemansplit"+"signal"+str(g)+".root")
        T2grs.append(T2gr)
    return T2grs

def fit_zeeman_splitting(graphs, path="/home/bernd/Dropbox/work/Hg (co)magnetometer/laser setup/HgSpektren/CoilCurrent"):
    zeemangr=get_zeeman_splitting(graphs, path)
    for g in range(len(zeemangr)):
        gROOT.SetStyle("Plain")
        c1 = TCanvas( 'c1', 'ScanData', 200, 10, 700, 500 )
        c1.SetGrid()
        zeemangr[g].Draw('A*')
        splitting=ROOT.TF1("splitting", "[0]+[1]*x", 0,5)
        splitting.SetParameters(0, -12)
        zeemangr[g].Fit("splitting", "Q+", "LEGO", 0, 7)
        c1.Update()
        p=splitting.GetParameters()
        yoffs=p[0]
        slope=p[1]
        print "yoffset="+str(yoffs)+" , slope="+str(slope)
        c1.Print(path+"/"+"ZeemanplitFitted"+"_signal"+str(g)+".png")
        c1.SaveAs(path+"/"+"ZeemansplitFitted"+"_signal"+str(g)+".root")
        raw_input("euuft")
        
def signal_to_noise(graphs, path):
    return 1
    
    

def draw_param_graph(path="/home/bernd/Dropbox/work/Hg (co)magnetometer/laser setup/HgSpektren/CoilCurrent"):
    """void"""
    gROOT.SetStyle("Plain")
    gr, parnames=create_param_graphs(path)
    c1 = TCanvas( 'c1', 'ScanData', 200, 10, 700, 500 )
    c1.SetGrid()
    gr.Draw('A*')
    c1.Update()
    print parnames
    raw_input("euuft")

def save_param_graph(gr, path, parname):
    gROOT.SetStyle("Plain")
    c1 = TCanvas( 'c1', 'ScanData', 200, 10, 700, 500 )
    c1.SetGrid()
    gr.Draw('A*')
    c1.Update()
    c1.Print(path+"/"+parname+".png")
    


