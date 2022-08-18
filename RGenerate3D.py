from ifrturbinepackage.rotor import *
from ifrturbinepackage.definitions import *


savenumber    = 11



tenflow_coeff = 1.5
tenwork_coeff = 15
cyclenum      = 10
gparamsetnum  = 6
rpmnum        = 4

result = ComputeR3(tenflow_coeff,tenwork_coeff,cyclenum,gparamsetnum,rpmnum)
n   = 4
Z5  = 0
Quasi    = QuasiNorm(result['proceed'],n,Z5)
Zreg     = Zrregs(Quasi)
Mer      = meridional (Zreg)
Vect     = VectorComp(Mer)
CoordR3D = Coord3D(Vect,result['geometry']['tb4'],result['geometry']['tb5'])

df3drotor   = pd.DataFrame(CoordR3D)
savetoas    = os.path.join(ROOT_DIR,"Outputs",f"rotor{savenumber}coordinate.csv")
df3drotor.to_csv(savetoas,index=False)

