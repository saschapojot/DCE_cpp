import numpy as np
from datetime import datetime
import json
import pandas as pd
import matplotlib.pyplot as plt
import glob
import re
import xml.etree.ElementTree as ET
#this script reads json files and plots photon and phonon numbers
#this script plots initial and final wavefunction

groupNum=6

rowNum=0

inDir="./groupNew"+str(groupNum)+"/row"+str(rowNum)+"/"
inParamFileName="./inParamsNew"+str(groupNum)+".csv"
dfstr=pd.read_csv(inParamFileName)
oneRow=dfstr.iloc[rowNum,:]

j1H=int(oneRow.loc["j1H"])
j2H=int(oneRow.loc["j2H"])

g0=float(oneRow.loc["g0"])
omegam=float(oneRow.loc["omegam"])
omegap=float(oneRow.loc["omegap"])
omegac=float(oneRow.loc["omegac"])
er=float(oneRow.loc["er"])#magnification

thetaCoef=float(oneRow.loc["thetaCoef"])

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
plt.title("$g_{0}=$"+str(g0)+", $e^{r}=$"+str(er)+", $\omega_{c}=$"+str(omegac))

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
plt.title("$g_{0}=$"+str(g0)+", $e^{r}=$"+str(er)+", $\omega_{c}=$"+str(omegac))
plt.xticks(xTicks)
plt.xlabel("time")
plt.ylabel("photon number")
plt.legend(loc="upper left")
plt.savefig(inDir+"photon.png")
plt.close()


#read initial wavefunction
initXML=""
finalXML=""
for file in glob.glob(inDir+"./*.xml"):
    if "init" in file:
        initXML=file
    if "final" in file:
        finalXML=file

matchN1=re.search(r"N1(\d+)N2",initXML)
if matchN1:
    N1=int(matchN1.group(1))

matchN2=re.search(r"N2(\d+)L1",initXML)
if matchN2:
    N2=int(matchN2.group(1))

matchL1=re.search(r"L1(\d+(\.\d+)?)L2",initXML)
if matchL1:
    L1=float(matchL1.group(1))

matchL2=re.search(r"L2(\d+(\.\d+)?)init",initXML)
if matchL2:
    L2=float(matchL2.group(1))


def xml2Array(xmlFileName):
    """

    :param xmlFileName: xml file containing the initial or final wavefunction
    :return:
    """
    tree = ET.parse(xmlFileName)
    root = tree.getroot()
    complex_numbers = []
    for item in root.findall(".//item"):
        real = float(item.find('real').text)
        imag = float(item.find('imag').text)
    # Create a complex number and append it to the list
        complex_number = complex(real, imag)
        complex_numbers.append(complex_number)
    complex_numbers=np.array(complex_numbers)
    outArr=complex_numbers.reshape((N1,N2))
    return np.abs(outArr)

tInitStart=datetime.now()
arrInit=xml2Array(initXML)
plt.figure()
plt.imshow(arrInit)
plt.xlabel("$x_{2}$")
plt.ylabel("$x_{1}$")
plt.title("init")
plt.colorbar()
plt.savefig(inDir+"init.png")
plt.close()

tInitEnd=datetime.now()
print("plot init time: ",tInitEnd-tInitStart)


tFinalStart=datetime.now()
arrFinal=xml2Array(finalXML)
plt.figure()
plt.imshow(arrFinal)
plt.xlabel("$x_{2}$")
plt.ylabel("$x_{1}$")
plt.title("init")
plt.colorbar()
plt.savefig(inDir+"final.png")
plt.close()

tFinalEnd=datetime.now()
print("plot final time: ",tFinalEnd-tFinalStart)