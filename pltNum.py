import numpy as np

import json

import matplotlib.pyplot as plt
import glob
import re
#this script reads json files and plots photon and phonon numbers


groupNum=6

rowNum=0

inDir="./groupNew"+str(groupNum)+"/row"+str(rowNum)+"/"
photonJson=""
phononJson=""

for file in glob.glob(inDir+"./*.json"):
    # print(file)
    if "photon" in file:
        photonJson = file
    elif "phonon" in file:
        phononJson = file



with open(photonJson,"r") as fptr:
    photonData=json.load(fptr)

with open(phononJson,"r") as fptr:
    phononData=json.load(fptr)



photonNum=photonData["photonNum"]

timeAll=photonData["time"]

phononNum=phononData["phononNum"]


plt.figure()
tTot=np.max(timeAll)

#plot photon and phonon
plt.plot(timeAll,photonNum,color="blue",label="photon")
plt.plot(timeAll,phononNum,color="red",label="phonon")

xTicks=[0,1/4*tTot,2/4*tTot,3/4*tTot,tTot]
xTicks=[round(val,2) for val in xTicks]
plt.xticks(xTicks)

plt.xlabel("time")
plt.ylabel("number")
plt.legend(loc="upper left")
plt.savefig(inDir+"both.png")

plt.close()
#plt photon
plt.figure()
plt.plot(timeAll,photonNum,color="blue",label="photon")
plt.xticks(xTicks)
plt.xlabel("time")
plt.ylabel("photon number")
plt.legend(loc="upper left")
plt.savefig(inDir+"photon.png")
plt.close()

