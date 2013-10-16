# -*- coding: utf-8 -*-
"""
Created on Wed Sep 12 18:08:44 2012

@author: bernd

für jede messreihe einen eigenen ordner anlegen und den sich
verändernden parameter als dateinamen verwenden
"""
# todo: dictionary für die results anlegen
#log files benutzen und einlesen 
import dataprocessing
import os
from os.path import join, getsize
import ROOT
from array import array
from ROOT import gROOT, TCanvas, TGraph, gStyle, TMath

class analizeparams(object):
    def __init__(self, path, showgraphs=False):
        self.path=path
        self.showgraphs=showgraphs
        self.graphs=list()
        self.parnames=list()
        self.yaxisparnames=["width (MHz)", "x offset (MHz)", "signal amplitude (V)", "width (MHz)", "x offset (MHz)", "signal amplitude (V)", "SNR", "Chisquare"]
        #self.xaxisparnames=[]
        #self.data

    def fit_all_files_in(self):
        parameters=list()
        results=list()
        directory = os.path.join(self.path)
        for files in os.walk(directory):
            for file in files[2]:
                if file.endswith(".logo") or file.endswith(".txt"):
                    param=float(file[0:file.find(".",0,len(file))].replace("_","."))
                    parameters.append(param)
                    print self.path+"/"+file
                    data=dataprocessing.import_data(self.path+"/"+file)
                    result, parnames=dataprocessing.fit_all(data, self.showgraphs)
                    results.append(result)
                    print parameters, results
        return [parameters,results], parnames
    
    
    def create_param_graphs(self):
        """
        fit all the spectra in path and draw and save graphs showing the fitted 
        parameters over the measurement parameter that was changed in path
        """
        self.results, self.parnames=self.fit_all_files_in()
       # graphs=list()
        for n in range(len(self.results[0])):
            for m in range(len(self.results[1][n])):
                #y=[float(results[1][i][m][n]) for i in range(len(results[1]))]
                y=[]
                for i in range(len(self.results[1])):
                    try:
                        #results[1][i][m][n]
                        y.append(float(self.results[1][i][m][n]))
                    except:
                        pass
                T2gr=ROOT.TGraph(len(self.results[0]), array("d", self.results[0]), array("d", y))
                T2gr.SetLineColor( 2 )
                T2gr.SetLineWidth( 2 )
                T2gr.SetMarkerColor( 4 )
                T2gr.SetMarkerStyle( 21 )
                T2gr.SetTitle( 'T2=[0]' )
                T2gr.GetXaxis().SetTitle(self.path.split("/")[-1])
                T2gr.GetYaxis().SetTitle(self.yaxisparnames[n] )
                self.save_param_graph(T2gr, self.path, self.parnames[n]+"_signal"+str(m))
                self.save_param_graph_ASCII(self.results[0], y, self.parnames[n]+"_signal"+str(m), self.parnames)
                T2gr.SaveAs(self.path+"/"+self.parnames[n]+"_signal"+str(m)+".root")
                self.graphs.append([self.parnames[n],self.results[0],y])
        return T2gr, self.parnames, self.graphs
        raw_input("yo!")
        33
    def get_zeeman_splitting(self):
        T2grs=[]
        NumGrs=[]
 #       for g in range(len(self.graphs)/7):
 #           gincr=g
 #           if len(self.graphs)/9 == 1:
 #               diff=[self.graphs[1][2][i]-self.graphs[5][2][i] for i in range(len(self.graphs[1][2]))]
 #           else:
 #               diff=[self.graphs[1*len(self.graphs)/9+gincr ][2][i]-self.graphs[4*len(self.graphs)/9+gincr][2][i] for i in range(len(self.graphs[1][2]))]
        for g in range(len(self.graphs)):
            if self.graphs[g][0]=="x_offset1":
                xoff1=self.graphs[g][2]
                for g in range(len(self.graphs)): 
                    if self.graphs[g][0]=="x_offset2":
                        xoff2=self.graphs[g][2]
                        diff=[xoff1[i]-xoff2[i] for i in range(len(self.graphs[1][2]))]
                        T2gr=ROOT.TGraph(len(self.graphs[1][1]), array("d", self.graphs[1][1]), array("d", diff))
                        NumGr=[self.graphs[1][1],diff]
                        T2gr.SetLineColor( 2 )
                        T2gr.SetLineWidth( 2 )
                        T2gr.SetMarkerColor( 4 )
                        T2gr.SetMarkerStyle( 21 )
                        T2gr.SetTitle( 'T2=[0]' )
                        T2gr.GetXaxis().SetTitle( 'CoilCurrent' )
                        T2gr.GetYaxis().SetTitle( 'ZeemanSplittint' )
                        self.save_param_graph(T2gr, self.path, "ZeemanSplit_signal"+str(g))
                        T2gr.SaveAs(self.path+"/"+"Zeemansplit"+"signal"+str(g)+".root")
                        T2grs.append(T2gr)
                        NumGrs.append(NumGr)
                        self.save_param_graph_ASCII(self.graphs[1][1], diff, "Zeemansplit"+"signal"+str(g)+".ascii", "zeeman")
        return T2grs
    
    def fit_zeeman_splitting(self):
        zeemangr=self.get_zeeman_splitting()
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
            c1.Print(self.path+"/"+"ZeemanplitFitted"+"_signal"+str(g)+".png")
            c1.SaveAs(self.path+"/"+"ZeemansplitFitted"+"_signal"+str(g)+".root")
            raw_input("euuft")
            
    def signal_to_noise(self, graphs):
        pass        
        
    
    def draw_param_graph(self):
        """void"""
        gROOT.SetStyle("Plain")
        gr, parnames=self.create_param_graphs(self.path)
        c1 = TCanvas( 'c1', 'ScanData', 200, 10, 700, 500 )
        c1.SetGrid()
        gr.Draw('A*')
        c1.Update()
        print parnames
        raw_input("euuft")
    
    def save_param_graph(self, gr, filename, parname):
        gROOT.SetStyle("Plain")
        c1 = TCanvas( 'c1', 'ScanData', 200, 10, 700, 500 )
        c1.SetGrid()
        gr.Draw('A*')
        c1.Update()
        c1.Print(filename+"/"+parname+".png")
        #c1.Print(filename+"/"+parname+".png")
    
    def save_param_graph_ASCII(self, x, y, filename, parname):
        FILE=open(self.path+"/"+filename+".ascii","w")
        for i in range(0,len(y),1):
            FILE.write(str(x[i])+"    "+str(y[i])+"\n")
        FILE.close()


