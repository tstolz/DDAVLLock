import socket
import time
import string
def simplequery(ip, port, command):
    time.sleep(1)
    s=socket.create_connection((ip,port))
    s.sendall(command+chr(13)+chr(10))
    time.sleep(1)
    answer=s.recv(1024)
    print answer
    indexi=answer.find("=",0,len(answer))
    indexe=answer.find(chr(13), indexi, len(answer))
    answercut=answer[indexi+1:indexe]
    if "m" in answercut:
        answercut=answercut[0:len(answercut)-1]
        value=float(answercut)*0.001
    else:
        value=float(answercut)
    print indexi
    print answer
    s.close()
    return value
print "test"
