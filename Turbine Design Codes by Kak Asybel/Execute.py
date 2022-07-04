from Compute_Rotor import*
from Regress_Graph import*
from Compute_Nozzle import*
from Compute_Volute import*
fluid='R245fa'

#Compute Rotor
O=Compute(.2,1.1,fluid)
NR,r4,Alpha4,b4,Ct4,rho4,mflow=O[2],O[3],O[4],O[5],O[6],O[7],O[8]

#Compute Nozzle
VariableNozz(NR,r4,Alpha4,b4,Ct4,rho4,mflow)
VarNozz=VariableNozz(NR,r4,Alpha4,b4,Ct4,rho4,mflow)
Y03=2
while Y03>1:
    GuessNozzList=GuessNozz()
    ai=GuessNozzList[0]
    bi=GuessNozzList[1]
    ci=GuessNozzList[2]
    Y03=NozzCalc1(ai,bi,ci)[0]
    if Y03<=1:
        break

#Looping Initial Value
ListNozzCalc1=NozzCalc2()
Errora=ListNozzCalc1[3]
Errorb=ListNozzCalc1[4]
Errorc=ListNozzCalc1[5]
Errorabc=abs(Errora)+abs(Errorb)+abs(Errorc)
a=ListNozzCalc1[0]
b=ListNozzCalc1[1]
c=ListNozzCalc1[2]

while Errorabc>=10e-5:
    NewValList=NewVal(Errora,Errorb,Errorc,a,b,c)
    ai=NewValList[0]
    bi=NewValList[1]
    ci=NewValList[2]
    Y03=NozzCalc1(ai,bi,ci)[0]
    while Y03>1:
        GuessNozzList=GuessNozz()
        ai=GuessNozzList[0]
        bi=GuessNozzList[1]
        ci=GuessNozzList[2]
        if Y03<=1:
            break
        Y03=NozzCalc1(ai,bi,ci)[0]
    if Errorabc<10e-5:
        break
    ListNozzCalc1=NozzCalc2()
    Errora=ListNozzCalc1[3]
    Errorb=ListNozzCalc1[4]
    Errorc=ListNozzCalc1[5]
    Errorabc=abs(Errora)+abs(Errorb)+abs(Errorc)
    a=ListNozzCalc1[0]
    b=ListNozzCalc1[1]
    c=ListNozzCalc1[2]
    print(Errorabc)
o3=VarNozz[7]
s3=VarNozz[6]
r2=VarNozz[5]
print(a,b,c,s3,o3,r2)

#Compute Volute
ComputeVol(fluid)