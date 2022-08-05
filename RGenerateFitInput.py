import matplotlib.pyplot as plt
from ifrturbinepackage.rotor import *
from ifrturbinepackage.definitions import *

# Gather Data for specific cycle, gparams set, and rpm

# apa yang ingin di plot?
# => dijadikan input

cyclenum    = 10
gparamsetnum   = 1
rpmnum      = 4

MaxVar1 = 11            #Max flow Coeffx10 in axis
MinVar1 = 0.3             #Min flow Coeffx10 in axis
MaxVar2 = 11            #Max Work Coeffx10 in axis
MinVar2 = 5             #Min Work Coeffx10 in axis
ndata1  = 5
ndata2  = 5
z_datapos = []
x_datapos = []
y_datapos = []
dvar2   = (MaxVar2-MinVar2)/(ndata2-1)
dvar1   = (MaxVar1-MinVar1)/(ndata1-1)
dataindex=0
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

            except ValueError:
                continue

fig = plt.figure(figsize=(10,10))
ax  = plt.axes(projection = "3d")
ax.scatter(x_datapos,y_datapos,z_datapos,lw=1,c= z_datapos,cmap="viridis")
ax.set_title("Efficiency vs Flow Coeff and Work Coeff Fluida 1")
ax.set_xlabel("Flow Coefficient")
ax.set_ylabel("Work Coefficient")
ax.set_zlabel("Efficiency (%)")
plt.show()
