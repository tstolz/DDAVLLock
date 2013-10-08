# -*- coding: utf-8 -*-
#known issue: make sure that the socket is also closed in case of a crash
"""
Created on Tue Aug 21 18:07:36 2012

@author: bernd
"""
import sys
import re
import socket
import time
global ip 
global port
global piezochan#chan number of piezo output voltage in scope
ip = "10.155.59.248"
port=60001
piezochan=2
def simplequery(ip, port, command):
    """Send simple query to Digilock and return the requested value"""
    s=socket.create_connection((ip,port))
    s.sendall(command+chr(13)+chr(10))
    time.sleep(0.1)
    try:
        answer=recv_timeout2(s,1)
        #print answer
        indexi=answer.find(command[0:len(command)-1]+"=",0,len(answer))+len(command)-1
        indexe=answer.find(chr(13), indexi, len(answer))
        answercut=answer[indexi+1:indexe]
        #print answercut
        if "m" in answercut:
            answercut=answercut[0:len(answercut)-1]
            value=float(answercut)*0.001
        elif "u" in answercut:
            answercut=answercut[0:len(answercut)-1]
            value=float(answercut)*0.000001
        else:
            value=float(answercut)
        s.close()
    except:
        s.close()
        print "Unexpected error:", sys.exc_info()[0]
    return value
def simplecommand(ip, port, command):
    s=socket.create_connection((ip,port))
    s.sendall(command+chr(13)+chr(10))
    s.close()
def longquery(ip, port, command):
    s=socket.create_connection((ip,port))
    s.sendall(command+chr(13)+chr(10))
    time.sleep(0.1)
    try:#exeption handling: when something in the code between socket create and close happens, s.close() is executed anyway
        answer=recv_timeout(s,1)
        #print answer
        indexi=answer.find(command[0:len(command)-1]+"=",0,len(answer))+len(command)-1
        #print indexi
        #indexe=answer.find(chr(13), indexi, len(answer))
        answercut=answer[indexi+1:len(answer)]
        #print "this is the longquery output"
        #print answercut
        s.close()
    except:
        s.close()
        print "Unexpected error:", sys.exc_info()[0]
    return answercut
def recv_timeout(the_socket,timeout=2):
    the_socket.setblocking(0)
    total_data=[];data='';begin=time.time()
    while 1:
        #if you got some data, then break after wait sec
        #print total_data.find(">")
        #if total_data and total_data.find(">"+chr(13))!=-1:
        #    print total_data
        #    break
        if total_data and time.time()-begin>timeout:
            break
        #if you got no data at all, wait a little longer
        elif time.time()-begin>timeout*10:
            break
        try:
            data=the_socket.recv(8192)
            if data:
                total_data.append(data)
                begin=time.time()
            else:
                time.sleep(0.1)
        except:
            pass
    #total_data.split("/n")
    total_data=str(''.join(total_data))
    #total_data=total_data[0]
    #print "this is the recv_timeout output"
    print total_data
    #total_data=total_data[0]
    return total_data
def recv_timeout2(the_socket,timeout=2):
    the_socket.setblocking(0)
    total_data=[];data='';begin=time.time()
    while 1:
        #if you got some data, then break after wait sec
        #print total_data.find(">")
        #if total_data and total_data.find(">"+chr(13))!=-1:
        #    print total_data
        #    break
        if total_data:
            break
        #if you got no data at all, wait a little longer
        elif time.time()-begin>timeout*10:
            break
        try:
            data=the_socket.recv(8192)
            if data:
                total_data.append(data)
                begin=time.time()
            else:
                time.sleep(0.1)
        except:
            pass
    #total_data.split("/n")
    total_data=str(''.join(total_data))
    #total_data=total_data[0]
    #print "this is the recv_timeout output"
    #print total_data
    #total_data=total_data[0]
    return total_data
def getoutputvoltage():
    command="scope:ch"+str(piezochan)+":mean?"
    value=simplequery(ip,port,command)
    return value 
def setoffset(offsetvoltage):
    command="offset:value="+str(offsetvoltage)
    simplecommand(ip,port,command)
def getoffset():
    command="offset:value?"
    return simplequery(ip,port,command)
def setscanrange(ampl):
    command="scan:amplitude="+str(ampl)
    simplecommand(ip,port,command)
def setscanfrequency(freq):
    command="scan:frequency="+str(freq)
    simplecommand(ip,port,command)
def setaveraging(points):
    command="scope:smooth:number="+str(points)
    simplecommand(ip,port,command)
def setscopetimescale(time):
    command="scope:timescale="+str(time)
    simplecommand(ip,port,command)
def switchscan(bool):
    if bool:
        command="scan:enable=true"
    elif bool==False:
        command="scan:enable=false"
    else:
        print"boolean needed"
    return simplequery(ip,port,command)
def getscopedata():
    command="scope:graph?"
    scope_data=longquery(ip, port,command)
    #print scope_data
    #scope_data.split("/r/n")
    print "this is the getscopedataoutput"
    scope_data_list=scope_data.split(chr(13)+chr(10))
    print scope_data_list
    scope_data_splitted=map(ssplit,scope_data_list)
    del scope_data_splitted[-1]
    scope_data_splitted_float=[map(float,x) for x in scope_data_splitted]#[:-1]
    scope_data_transposed=zip(*scope_data_splitted_float)
    print scope_data_transposed
    return scope_data_transposed
def ssplit(string_to_split):
    return string_to_split.split(chr(9))


print "test"
