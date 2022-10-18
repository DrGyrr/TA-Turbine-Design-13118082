from ifrturbinepackage.definitions  import *
from ifrturbinepackage.inputs       import *
from ifrturbinepackage.rotor        import *
from ROptimize                      import *
import matplotlib.pyplot as plt

cyclenum = 10
# cyclenum = int(input("which cycle #?: "))

optresR = optimizeR(cyclenum)
# askproceed = str(input("proceed? [Y/N]"))
# if askproceed.upper() == 'Y' or askproceed.upper() == 'YES':
R2DContour  = Gen2DContour(indict=optresR,z5=0,dataamount=70)
Rangles     = BladeAngles(indict=R2DContour,ns=4)
Meri        = meridional(Rangles)
Vect        = VectorComp(Meri)
print(f"this is tb4 tb5 {R2DContour['tb4']},{R2DContour['tb5']} ")
R3DCoord    = Coord3D(Vect,R2DContour['tb4'],R2DContour['tb5'])

df3DRotor   = pd.DataFrame(R3DCoord)
shroudpos = ['XcorSP','YcorSP','ZcorSP']
shroudneg = ['XcorSN','YcorSN','ZcorSN']
hubpos = ['XcorHP','YcorHP','ZcorHP']
hubneg = ['XcorHN','YcorHN','ZcorHN']
# df3DRotor.to_csv(savetoas,index=False)
# for whatfile in [shroudpos,shroudneg,hubpos,hubneg]:
savetoas    = os.path.join(ROOT_DIR,"Outputs",f"{cyclenum}_shroudpos_3D Coordinate.txt")
df3DRotor.to_csv(savetoas,columns=shroudpos,index=None,header=None,sep=' ',mode='w')
savetoas    = os.path.join(ROOT_DIR,"Outputs",f"{cyclenum}_shroudneg_3D Coordinate.txt")
df3DRotor.to_csv(savetoas,columns=shroudneg,index=None,header=None,sep=' ',mode='w')
savetoas    = os.path.join(ROOT_DIR,"Outputs",f"{cyclenum}_hubpos_3D Coordinate.txt")
df3DRotor.to_csv(savetoas,columns=hubpos,index=None,header=None,sep=' ',mode='w')
savetoas    = os.path.join(ROOT_DIR,"Outputs",f"{cyclenum}_hubneg_3D Coordinate.txt")
df3DRotor.to_csv(savetoas,columns=hubneg,index=None,header=None,sep=' ',mode='w')
# print(optresR)
# print(R2DContour)
# plt.show(plotkeun[4])