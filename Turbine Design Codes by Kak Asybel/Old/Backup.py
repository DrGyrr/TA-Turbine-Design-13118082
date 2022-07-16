from Compute_Rotor import*
from Regress_Graph import*
from Compute_Nozzle import*
def Execution():
    global NR,r4,Alpha4,b4,Ct4,rho4,mflow,Reaction
    O=Compute(.2,1.1,'R245fa')
    NR,r4,Alpha4,b4,Ct4,rho4,mflow=O[2],O[3],O[4],O[5],O[6],O[7],O[8]
    print(NR,r4,Alpha4,b4,Ct4,rho4,mflow)
    VariableNozz(NR,r4,Alpha4,b4,Ct4,rho4,mflow)
    GuessNozz()
    NozzCalc1()
    Y03=NozzCalc1()[0]   
    LoopingNozz1()
    print(Y03)
    Errorabc=LoopingNozz1()
    LoopingNozz2(Errorabc)
    print()
    print(a,b,c)
Execution()



from math import degrees
import numpy as np
import CoolProp
from CoolProp.CoolProp import PropsSI as Pr
from Compute_Rotor import*
#Initial Variable

def VariableNozz(NR,r4,Alpha4,b4,Ct4,rho4,mflow):
  global aOc,Betha2,Nn,r3,Ct3,Cm3,r2Or3,r2,s3,o3,Alpha3
  aOc = 0.3
  Betha2=np.radians(21)
  Nn = NR + 1
  r3=r4*(1+(2*b4*np.sin(np.radians(Alpha4)))/r4)
  Ct3=Ct4*r4/r3
  Cm3=mflow/(2*np.pi*r3*b4*rho4)
  r2Or3=1.1
  r2=r2Or3*r3
  s3=2*np.pi*r3/Nn
  o3=s3*Cm3/(np.sqrt(Cm3**2+Ct3**2))
  Alpha3=np.arcsin(o3/s3)
  return(Betha2,Nn,r3,Ct3,Cm3,r2,s3,o3,Alpha3)
 
def GuessNozz():
  global bi,ci,ai
  bi=abs(np.random.uniform(0.1,0.5)*o3)
  ci=abs(np.random.uniform(0,1,1)*(r3-r2))  
  ai= ci*aOc
  return(ai,bi,ci)
    
def NozzCalc1(ai,bi,ci):
  #Variabel Desain
  global Y03,Tetha,Y2
  t2Oc = 0.03
  t3Oc = 0.015
  tmaxOc = 0.06
  dOc = 0.4
  #Kalkulasi Y03
  X2=np.arctan(4*bi/(4*ai-ci))
  X3=np.arctan(4*bi/(3*ci-4*ai))
  Tetha=X2+X3
  Y2=Betha2+X2
  Y03=r2/r3*np.cos(Y2)
  return (Y03)
       
def NozzCalc2():
  global a,b,c,Errora,Errorb,Errorc,Errorabc
  Y3=np.arccos(Y03)
  c=(r2-r3)/np.sin(Y3)
  a=ci*aOc
  b=((np.sqrt(1+(4*np.tan(Tetha)**2)*(ai/ci-(ai/ci)**3-3/16)))-1)/(4*np.tan(Tetha))*ci
  Errora=abs(ai-a)
  Errorb=abs(bi-b)
  Errorc=abs(ci-c)
  #Errorabc=Errora+Errorb+Errorc
  return(a,b,c,Errora,Errorb,Errorc)  

def LoopingNozz1():
    global Errorabc
    if Y03>1:
      GuessNozz()
      NozzCalc1()
    elif Y03<=1:
      NozzCalc2()
      Errorabc=Errora+Errorb+Errorc
    return(Errorabc)
      
def LoopingNozz2(Errorabc):
  while Errorabc > 10e-5:
    GuessNozz()
    NozzCalc1()
    LoopingNozz1()
    if Errorabc <= 10e-5:
      break

def NewVal(Errora,Errorb,Errorc,a,b,c):
  bi= b+Errorb
  ci= c+Errorc
  ai= a+Errora
  return(ai,bi,ci)

def ReturnVal():
  if Errorabc==False:
          GuessNozz()
  else:
      NewVal(Errorabc,Errora,Errorb,Errorc,a,b,c)


