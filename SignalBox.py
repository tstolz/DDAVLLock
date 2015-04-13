import time
import serial
import json
import matplotlib
import numpy as np
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import random

class SignalBox():
    '''Class for communication with a skript on an Arduino Uno, controlling two digipotis'''
    def __init__(self,devicePath):
        # timeout for waiting for the arduino's answer
        self.timeout=3
        # path to the port e.g. 'COM3'
        self.devicePath=devicePath
        # open the serial port
        self.port=serial.Serial(devicePath,115200,timeout=1)
        time.sleep(2)
    def __del__(self):
        self.port.close()
    def __communicateTaskWithJSON(self,task):
        'task is a dictionary containing keys "command" and (if command is write) "value"'
        # convert dictionary 'task' to json string
        string=json.dumps(task)
        # send string to the serial port
        self.port.write(string)
        # wait until the arduino answers or until timeout
        t=time.time()
        while self.port.inWaiting()==0 and not time.time()-t > self.timeout:
            pass
        # wait until the arduino has finished answering (python is MUCH faster than arduino skripts)
        time.sleep(0.1)
        # convert the answer to a dictionary using json
        try:
            result=json.loads(self.port.read(self.port.inWaiting()))
        except:
            result=dict(value=0,error=1,errorMessage="received bad string")
        return result        
    def reset(self):
        self.closePort()
        self.port=serial.Serial(self.devicePath,9600,timeout=1)
    def closePort(self):
        try:
            self.port.close()
        except: pass
    def setPreGain(self,value):
        task=dict(command="write1",value=value)
        return self.__communicateTaskWithJSON(task)
    def getPreGain(self):
        task=dict(command="read1")
        return self.__communicateTaskWithJSON(task)
    def increasePreGain(self):
        task=dict(command="inc1")
        return self.__communicateTaskWithJSON(task)
    def decreasePreGain(self):
        task=dict(command="dec1")
        return self.__communicateTaskWithJSON(task)
    def setGain(self,value):
        task=dict(command="write2",value=value)
        return self.__communicateTaskWithJSON(task)
    def getGain(self):
        task=dict(command="read2")
        return self.__communicateTaskWithJSON(task)
    def increaseGain(self):
        task=dict(command="inc2")
        return self.__communicateTaskWithJSON(task)
    def decreaseGain(self):
        task=dict(command="dec2")
        return self.__communicateTaskWithJSON(task)
    def _data4oszi(self):
        # generator function for oszi
        data=[]
        while True:
                read=self.port.read(self.port.inWaiting())
                # arduino sends data as integers separated by newlines
                read=read.splitlines()
                # the first and last entry could be incomplete
                for line in read[1:-1]:
                    try:
                        data.append(int(line))
                    except:
                        pass
                yield data
                data=[]
    def oszi(self,maxt,trigger=False):
        ''' "maxt" is the number of plottet datapoints '''
        try:
            self.port.write('{"command":"adcON"}')            
            fig = plt.figure()
            ax = fig.add_subplot(111)
            scope = Scope(ax, maxt)

            # pass a generator in to produce data for the update func
            ani = animation.FuncAnimation(fig, scope.update, self._data4oszi, fargs=[True], interval=0, blit=True)
            plt.show()
            
        except:
            # if the plot application is killed, make sure the arduino is told
            # to stop polling
            self.port.write('{"command":"adcOFF"}')
            time.sleep(2)
            # empty the serial buffer (contains the arduino's json-answer)
            self.port.read(self.port.inWaiting())
                

# Emulate an oscilloscope.  Requires the animation API introduced in
# matplotlib 1.0 SVN.

class Scope:
    def __init__(self, ax, maxt=5000):
        self.ax = ax
        self.maxt = maxt
        self.tdata = [0]
        self.ydata = [0]
        self.line = Line2D(self.tdata, self.ydata)
        self.ax.add_line(self.line)
        self.ax.set_ylim(-.1, 1026)
        self.ax.set_xlim(0, self.maxt)
        # maximum number of data points used for triggering
        self.trigzonedefault=int(self.maxt/2)
        self.trigzone=self.trigzonedefault
        self.triglevel=0

    def trigger(self, data):
        # works as a trigger for the oszilloscope.
        # triglevel is chosen as the maximum of the previous data set. 
        self.triglevel=np.array(self.ydata).max()
        # trigzone is the 'timeout' for triggering.
        self.trigzone-=len(data)
        if self.trigzone<0:
            zone=len(data)+self.trigzone
        else:
            zone=len(data)
        for i in xrange(zone):
            # if the new signal reaches the triglevel within some tolerance (5)
            # the plot is continued
            if abs(self.triglevel-data[i])<5:
                self.trigzone=self.trigzonedefault
                return data[i:]
        if self.trigzone<0:
            # if no match has been achieved within trigzone, reset the variable
            # and continue plotting
            self.trigzone=self.trigzonedefault
            return data[zone:]
        # no match but trigzone has not been reached -> wait for new data
        return None
    
    def update(self, data, trigger=False):
        # check if the end of the plotwindow has been reached
        if len(self.tdata)+len(data)>self.maxt:
            # fill remaining space
            self.ydata+=data[:self.maxt-len(self.tdata)]
            data=data[self.maxt-len(self.tdata):]
            # trigger the rest of data
            if trigger:
                data=self.trigger(data)
            # if trigger was successful, start new plot
            if data:
                self.ydata=data
        else:
            self.ydata+=data
                
        self.tdata=range(len(self.ydata))
        self.ax.figure.canvas.draw()
        self.line.set_data(self.tdata, self.ydata)
        return self.line,
