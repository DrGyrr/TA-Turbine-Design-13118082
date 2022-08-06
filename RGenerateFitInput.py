import matplotlib.pyplot as plt
from ifrturbinepackage.rotor import *
from ifrturbinepackage.definitions import *
from tqdm import tqdm

# Gather Data for specific cycle, gparams set, and rpm

# apa yang ingin di plot?
# => dijadikan input

cyclenum    = 10
gparamsetnum   = 1
rpmnum      = 4

MaxVar1 = 14            #Max flow Coeffx10 in axis
MinVar1 = 0.3             #Min flow Coeffx10 in axis
MaxVar2 = 30            #Max Work Coeffx10 in axis
MinVar2 = 2             #Min Work Coeffx10 in axis
ndata1  = 10
ndata2  = 10
z_datapos = []
x_datapos = []
y_datapos = []
dvar2   = (MaxVar2-MinVar2)/(ndata2-1)
dvar1   = (MaxVar1-MinVar1)/(ndata1-1)
dataindex=0
countneg =0
countValEr =0
ndatatocalculate = ndata1*ndata2

with tqdm (total=ndatatocalculate,desc='Computed Data',unit='pts',) as pbar:
    for i in range(0,ndata2):


            Var2        = MinVar2 + i*dvar2
            for j in range (0,ndata1):
                try:
                    Var1    = MinVar1 + j*dvar2
                    flowcoeff10 = Var1
                    workcoeff10 = Var2
                    result  = ComputeR3(flowcoeff10,workcoeff10,cyclenum,gparamsetnum,rpmnum)
                    ztoplot = result['efficiency']['Effts']
                    if ztoplot>=0:
                        z_datapos.append(ztoplot)
                        x_datapos.append(Var1)
                        y_datapos.append(Var2)
                    if ztoplot<=0:
                        countneg += 1
                    pbar.update(1)
                except ValueError:
                    pbar.update(1)
                    countValEr += 1
                    continue
                # pbar.set_postfix
                # pbar.display(msg='\ntest', pos=None)
                # print(f"\n Negative Data : {countneg}",end="\r")
                # print(f"Value Error Data: {countValEr}",end="\r")
                
                
pbar.close()
print(f"Negative data amount: {countneg}")
print(f"Value Error data amount: {countValEr}")
print("3D scatter plot is opened in new window")

fig = plt.figure(figsize=(10,10))
ax  = plt.axes(projection = "3d")
ax.scatter(x_datapos,y_datapos,z_datapos,lw=1,c= z_datapos,cmap="viridis")
ax.set_title("Efficiency vs Flow Coeff and Work Coeff Fluida 1")
ax.set_xlabel("Flow Coefficient")
ax.set_ylabel("Work Coefficient")
ax.set_zlabel("Efficiency (%)")
plt.show()

dosave = input("Would you like to save genrated data and the scatter plot? [Y/N]")