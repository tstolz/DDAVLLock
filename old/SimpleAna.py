# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 11:43:19 2012

@author: bernd
"""
from analizeparamsclass import analizeparams
#import dataprocessing
execfile("analizeparamsclass.py")
path=input("please enter the path containing the to-be-evaluated files")
showgraphs=input("want to see the Graphs? Enter True or False")
if showgraphs=="True":
    showgraphs=True
else:
    showgraphs=False
a=analizeparams(path, showgraphs)
a.create_param_graphs()
zee=input("want create Zeeman Graphs? Enter True or False")
if zee == True:
    a.fit_zeeman_splitting()
else:
    print zee
#    pass