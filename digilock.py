# -*- coding: utf-8 -*-
#known issue: make sure that the socket is also closed in case of a crash
"""
Created on Tue Aug 21 18:07:36 2012

@author: bernd
"""
import sys
import socket
import time

global piezochan

piezochan=2

class digilock:
    def __init__(self, ip, port):
        self.ip=ip
        self.port=port
    def query(self, command):
        s=socket.create_connection((self.ip,self.port))
        try:#exeption handling: when something in the code between socket create and close happens, s.close() is executed anyway
            s.sendall(command+chr(13)+chr(10))
            time.sleep(0.1)
            answer=self.recv_timeout(s,1)
            indexi=answer.find(command[0:len(command)-1]+"=")+len(command)-1
            answercut=answer[indexi+1:len(answer)].strip(chr(13)+chr(10)+'>')
        except:
            print "Unexpected error:", sys.exc_info()[0]
        s.close()        
        return answercut
    def queryNUM(self, command):
        '''Send a query and interpret the result as numeric'''
        answercut=self.query(command)
        if "m" in answercut:
            answercut=answercut[0:len(answercut)-1]
            value=float(answercut)*0.001
        elif "u" in answercut:
            answercut=answercut[0:len(answercut)-1]
            value=float(answercut)*0.000001
        else:
            value=float(answercut)
        return value
    def recv_timeout(self, the_socket, timeout=2):
        the_socket.setblocking(0)
        total_data=[];data='';begin=time.time()
        while 1:
            #if you got some data, then break after wait sec
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
        total_data=str(''.join(total_data))
        print total_data
        return total_data
    def simplecommand(self, command):
        s=socket.create_connection((self.ip,self.port))
        s.sendall(command+chr(13)+chr(10))
        s.close()
    def getoutputvoltage(self):
        command="scope:ch"+str(piezochan)+":mean?"
        value=self.queryNUM(command)
        return value
    def setoffset(self, offsetvoltage):
        command="offset:value="+str(offsetvoltage)
        self.simplecommand(command)
    def getoffset(self):
        command="offset:value?"
        return self.queryNUM(command)
    def setscanrange(self,ampl):
        command="scan:amplitude="+str(ampl)
        self.simplecommand(command)
    def setscanfrequency(self,freq):
        command="scan:frequency="+str(freq)
        self.simplecommand(command)
    def setaveraging(self,points):
        command="scope:smooth:number="+str(points)
        self.simplecommand(command)
    def setscopetimescale(self,time):
        command="scope:timescale="+str(time)
        self.simplecommand(command)
    def switchscan(self,enable):
        if enable:
            command="scan:enable=true"
        elif enable==False:
            command="scan:enable=false"
        else:
            print"boolean needed"
        return self.query(command)
    def getscopedata(self):
        command="scope:graph?"
        scope_data=self.query(command)
        #print scope_data
        #scope_data.split("/r/n")
        print "this is the getscopedataoutput"
        scope_data_list=scope_data.split(chr(13)+chr(10))
        print scope_data_list
        scope_data_splitted=map(lambda x: x.split(chr(9)),scope_data_list)
        del scope_data_splitted[-1]
        scope_data_splitted_float=[map(float,x) for x in scope_data_splitted]#[:-1]
        scope_data_transposed=zip(*scope_data_splitted_float)
        print scope_data_transposed
        return scope_data_transposed

D=digilock('localhost',60001)
