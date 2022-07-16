from Compute_Rotor import*
from Regress_Graph import*
from Compute_Nozzle import*
from Compute_Volute import*
NR,r4,Alpha4,b4,Ct4,rho4,mflow=22,0.0291,12.58,3.25064/1000,70.7,32.0081,0.3

#VariableNozzle
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
Errorx=abs(Errora)+abs(Errorb)+abs(Errorc)
print(ai,bi,ci,Errorx,Y03)
aOc=0.3
Betha2=np.radians(21)
while Errorx>=10e-5:
    bi=bi+Errorb/2
    ci=ci+Errorc/2
    ai=ci*aOc
    Y03=NozzCalc1(ai,bi,ci)[0]
    while Y03>1:
        GuessNozzList=GuessNozz()
        ai=GuessNozzList[0]
        bi=GuessNozzList[1]
        ci=GuessNozzList[2]
        Y03=NozzCalc1(ai,bi,ci)[0]  
        if Y03<=1:
            break
    ListNozzCalc1=NozzCalc2()
    Errora=ListNozzCalc1[3]
    Errorb=ListNozzCalc1[4]
    Errorc=ListNozzCalc1[5]
    Errorx=abs(Errora)+abs(Errorb)+abs(Errorc)
    if Errorx<10e-5:
        print(ai,bi,ci,Errorx)
        break

    print(Errora,Errorb,Errorc)
print(a,b,c,s3,o3,r2)