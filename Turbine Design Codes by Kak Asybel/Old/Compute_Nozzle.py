from math import degrees
from pydoc import doc
import numpy as np
import CoolProp
from CoolProp.CoolProp import PropsSI as Pr
from Compute_Rotor import*

def VariableNozz(NR,r4,Alpha4,b4,Ct4,rho4,mflow):
  global aOc,Betha2,Nn,r3,Ct3,Cm3,r2Or3,r2,s3,o3,Alpha3,t2Oc,t3Oc,tmaxOc,dOc
  aOc = 0.3               #Varopt
  Betha2=np.radians(21)   #Varopt
  Nn = NR + 1             #Varopt
  r3=r4*(1+(2*b4*np.sin(np.radians(Alpha4)))/r4)
  Ct3=Ct4*r4/r3
  Cm3=mflow/(2*np.pi*r3*b4*rho4)
  r2Or3=1.2       #Varopt
  r2=r2Or3*r3
  s3=2*np.pi*r3/Nn
  o3=s3*Cm3/(np.sqrt(Cm3**2+Ct3**2))
  Alpha3=np.arcsin(o3/s3)
  t2Oc = 0.03     #Varopt
  t3Oc = 0.015    #Varopt
  tmaxOc = 0.06   #Varopt
  dOc = 0.4       #Varopt
  return(Betha2,Nn,r3,Ct3,Cm3,r2,s3,o3,Alpha3,t2Oc,t3Oc,tmaxOc,dOc)
 
def GuessNozz():
  global bi,ci,ai
  bi=abs(np.random.uniform(0.1,0.5)*o3)
  ci=abs(np.random.uniform(0,1,1)*(r3-r2))  
  ai= ci*aOc
  return(ai,bi,ci)
    
def NozzCalc1(ai,bi,ci):
  #Variabel Desain
  global Y03,Tetha,Y2
  #Kalkulasi Y03
  X2=np.arctan(4*bi/(4*ai-ci))
  X3=np.arctan(4*bi/(3*ci-4*ai))
  Tetha=X2+X3
  Y2=Betha2+X2
  Y03=r2/r3*np.cos(Y2)
  return (Y03)
       
def NozzCalc2(ai,bi,ci):
  global a,b,c,Errora,Errorb,Errorc
  Y3=np.arccos(Y03)
  c=(r2-r3)/np.sin(Y3)
  a=c*aOc
  b=((np.sqrt(1+(4*np.tan(Tetha)**2)*(ai/ci-(ai/ci)**3-3/16)))-1)/(4*np.tan(Tetha))*ci
  Errora=(ai-a)
  Errorb=(bi-b)
  Errorc=(ci-c)
  return(a,b,c,Errora,Errorb,Errorc)  

def VariantSearchNozz():
    global SListA,SListB,SListC,ErrorList
    SListB=[]
    SListC=[]
    SListA=[]
    ErrorList=[]
    # Y03List=[]
    SListSize=1000
    for k in range(0,SListSize):
        Y03=2
        while Y03>1:
            bi=abs(np.random.uniform(0.1,0.8)*o3)
            ci=abs(np.random.uniform(0,5,1)*(r3-r2))  
            ai=ci*aOc
            Y03=NozzCalc1(ai,bi,ci)[0]
            if Y03<=1:
                break
        ListNozzCalc1=NozzCalc2(ai,bi,ci)
        Errora=ListNozzCalc1[3]
        Errorb=ListNozzCalc1[4]
        Errorc=ListNozzCalc1[5]
        Errorx=abs(Errora)+abs(Errorb)+abs(Errorc)
        SListB.append(bi)
        SListC.append(ci)
        SListA.append(ai)
        ErrorList.append(Errorx)
        # Y03.append(Y03)
        
        # print(SListA)
        # print()
        # print(SListB)
        # print()
        # print(SListC)
        # print()
        # print(ErrorList)
    return(SListA,SListB,SListC,ErrorList)
