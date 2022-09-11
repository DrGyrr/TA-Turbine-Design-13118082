
# %% Imports
import matplotlib.pyplot as plt
from ifrturbinepackage.rotor import *
from ifrturbinepackage.definitions import *
from tqdm import tqdm
import os


# %% 

''' ComputeR4 is used
data generating mode:
[1] -> variasikan flow & work coefficient
[2] -> variasikan gparamset


'''


genmode = 3

# %%


cyclenum    = 10
gparamsetnum   = 1
rpmnum      = 4

MaxVar1 = 14            #Max flow Coeffx10 in axis
MinVar1 = 0.7           #Min flow Coeffx10 in axis
MaxVar2 = 27            #Max Work Coeffx10 in axis
MinVar2 = 2             #Min Work Coeffx10 in axis
ndata1  = 70
ndata2  = 70
z_datapos = []
x_datapos = []
y_datapos = []
dvar2   = (MaxVar2-MinVar2)/(ndata2-1)
dvar1   = (MaxVar1-MinVar1)/(ndata1-1)
dataindex   = 0
countneg    = 0
countValEr  = 0
countuse    = 0
ndatatocalculate = ndata1*ndata2

with tqdm (total=ndatatocalculate,desc='Computed',unit=' data') as pbar:
    for i in range(0,ndata2):

            
            Var2        = MinVar2 + i*dvar2
            for j in range (0,ndata1):
                try:
                    Var1    = MinVar1 + j*dvar2
                    flowcoeff10 = Var1
                    workcoeff10 = Var2
                    result  = ComputeR4(flowcoeff10,workcoeff10,cyclenum,gparamsetnum)
                    ztoplot = result['efficiency']['Effts']
                    if ztoplot>=0:
                        z_datapos.append(ztoplot)
                        x_datapos.append(Var1)
                        y_datapos.append(Var2)
                        countuse += 1
                    if ztoplot<=0:
                        countneg += 1
                    pbar.update(1)
                except ValueError:
                    pbar.update(1)
                    countValEr += 1
                    continue

                print(f"\n\x1b[2KNegative data amount: {countneg} ({countneg/ndatatocalculate*100:.1f}%)")
                print(f"\x1b[2KValue Error data amount: {countValEr} ({countValEr/ndatatocalculate*100:.1f}%)")
                print(f"\x1b[2KUseful data amount: {countuse} ({countuse/ndatatocalculate*100:.1f}%)",end="\r\033[A\033[A\033[A")
                
                
                
pbar.close()

print("\n\n\n->> 3D scatter plot is opened in new window. After finished viewing it and set the angle to save, close to proceed ")

fig = plt.figure(figsize=(12,12))
ax  = plt.axes(projection = "3d")
ax.scatter(x_datapos,y_datapos,z_datapos,lw=1,c= z_datapos,cmap="viridis")
ax.set_title("Efficiency vs Flow Coeff and Work Coeff Fluida 1")
ax.set_xlabel("Flow Coefficient x10")
ax.set_ylabel("Work Coefficient x10")
ax.set_zlabel("Efficiency ,%")



saveddir  = os.path.join(ROOT_DIR,"Outputs","FitInputData")
savedfigdir = os.path.join(saveddir,"Figures")
fig.savefig(os.path.join(savedfigdir,'temp.png'),dpi=300,transparent=True)
fig.savefig(os.path.join(savedfigdir,'temp.jpg'),dpi=300)
plt.show()
# print(ROOT_DIR,)


while True:
    dosave = str(input("Would you like to save generated data and the scatter plot? [Y/N] : "))
    if dosave.upper() == 'Y' or dosave.upper() == 'YES':
        # prevsaveddata = os.listdir(os.path.join(ROOT_DIR,"Outputs","FitInputData"))
        # for filepath in fileslist:
        #     timestamp   = time.strftime( '%m/%d/%Y :: %H:%M:%S', time.gmtime(os.path.getmtime(filepath)) )
        #     print(timestamp,'-->',filepath)
        print(f"Previously saved FitInput files:")
        print(f"\nModified date & time{' '*8}File name")
        print(f"{'-'*22}{' '*6}{'-'*15}")
        listcontent(saveddir)
        savenumber  = input("\nInsert save number to assign to generated data: ")
        FitInputdict = {
            '10xFlowCoeff'              : x_datapos,
            '10xWorkCoeff'              : y_datapos,
            'Total-Static_Efficiency'   : z_datapos
        }
        df = pd.DataFrame(FitInputdict)
        datasavedas = os.path.join(saveddir,f"{savenumber}.csv")
        plotsavedas = os.path.join(savedfigdir,f"{savenumber}.png")
        df.to_csv(datasavedas)
        fig.savefig(plotsavedas,dpi=300,transparent=True)
        print(f"Generating FitInput Data done.")
        print(f"Generated data saved as {savenumber}.csv")
        print(f"Generated plot saved as {savenumber}.png and {savenumber}.jpg")
        print("Plot in default position is also saved under name of temp.png and temp.jpg")
        print(f"Now you can proceed to use generated data for curve fitting use")
        # print("OK saved")
        # print(dosave.upper())
        break
    elif dosave.upper() == 'N' or dosave.upper() == 'NO':
        break
    else :
        print("Come on, wrong input. '_'")

print("Program finished. You can close this now '_'")

