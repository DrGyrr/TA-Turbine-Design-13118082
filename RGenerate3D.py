from ifrturbinepackage.rotor import *
from ifrturbinepackage.definitions import *


savenumber    = 12



tenflow_coeff = 1.5
tenwork_coeff = 15
cyclenum      = 10
gparamsetnum  = 2
rpmnum        = 4

result = ComputeR3(tenflow_coeff,tenwork_coeff,cyclenum,gparamsetnum,rpmnum)
# n   = 4
Z5  = 0
# Quasi    = QuasiNorm(result['proceed'],n,Z5)
# Zreg     = Zrregs(Quasi)

Contour2D   = Gen2DContour(result['proceed'],Z5,dataamount=150) 
Angles      = BladeAngles(Contour2D,4)
Mer         = meridional (Angles)
Vect        = VectorComp(Mer)
CoordR3D    = Coord3D(Vect,result['geometry']['tb4'],result['geometry']['tb5'])

df3drotor   = pd.DataFrame(CoordR3D)
savetoas    = os.path.join(ROOT_DIR,"Outputs",f"rotor{savenumber}coordinate.csv")
shroudpos = [XcorSP,YcorSP,ZcorSP]
shroudneg = [XcorSN,YcorSN,ZcorSN]
hubpos = [XcorHP,YcorHP,ZcorHP]
hubneg = [XcorHN,YcorHN,ZcorHN]
df3drotor.to_csv(savetoas,index=False)
for whatfile in [shroudpos,shroudneg,hubpos,hubneg]:
    savetoas    = os.path.join(ROOT_DIR,"Outputs",f"{savenumber}_{whatfile}_3D Coordinate.txt")
    df3drotor.to_csv(savetoas,columns=whatfile,index=None,header=None,sep=' ',mode='w')


