import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
qDf=[]
pDf=[]

fq="/mnt/c/Users/YeLab/Desktop/pq/data/pqsim3_R1.fq"
fq="C:\\Users\\YeLab\\Desktop\\pq\\data\\pqsim3_R1.fq"
def q2p(base):
   return 10**((ord(base)-33)/(-10))
def getq(base):
    return ord(base)-33


linenum=0
for line in open(fq):
    linenum+=1
    if linenum>1000:
        break
    if linenum%4==0:
        lineQ=line[:-1]
        pDf.append(list(map(q2p,lineQ)))
        qDf.append(list(map(getq,lineQ)))

pmean=np.array(pDf).mean(axis=0)
qmean=np.array(qDf).mean(axis=0)

x=[ ((i+1)/(150))**3 for i in range(150)]
x2=[ ((i+1)/(150)) for i in range(150)]
x3=[ (1.1**(i+1))/(1.1**150) for i in range(150)]
x4=[ (1.02**(i+1))/(1.02**150) for i in range(150)]

plt.plot(x)
plt.plot(x2)
plt.plot(x3)
plt.plot(x4)

plt.show()

plt.plot(qmean)
plt.ylim([0,50])


plt.show()
plt.plot(pmean*pmean*pmean)
# plt.ylim([0,50])

plt.show()


    # print(line)
