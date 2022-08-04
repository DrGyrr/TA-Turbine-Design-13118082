from ifrturbinepackage.rotor import *
from ifrturbinepackage.definitions import *

tenflow_coeff = 2
tenwork_coeff = 10
cyclenum      = 10
gparamsetnum  = 1
rpmnum        = 4

savenumber    = 1


result = ComputeR3(tenflow_coeff,tenwork_coeff,cyclenum,gparamsetnum,rpmnum)

n   = 4
Z5  = 0
Quasi    = QuasiNorm(result['proceed'],n,Z5)
Zreg     = Zrregs(Quasi)
Mer      = meridional (Zreg)
Vect     = VectorComp(Mer)
CoordR3D = Coord3D(Vect,result['geometry']['tb4'],result['geometry']['tb5'])
# Coorddict= dict()
# Coorddict= {
#     'xshroud':CoordR3D['XcorS']
#     'yshroud':CoordR3D['Ycors']
#     'zshroud':CoordR3D['ZcorS']
#     'xhub'   :CoordR3D['XcorH']
#     'yhub'   :CoordR3D['YcorH']
#     'zhub'   :CoordR3D['ZcorH']
# }
df3drotor   = pd.DataFrame(CoordR3D)
savetoas    = os.path.join(ROOT_DIR,"Outputs",f"rotor{savenumber}coordinate.csv")
df3drotor.to_csv(savetoas,index=False)
# print("lets go Mer",type(Mer),len(Mer))
# print("lets go quasi",type(Quasi),len(Quasi))
# print("lets go Zreg", type(Zreg),len(Zreg))
# print(f"done, save in {savetoas}")
# # for i in ('Z','r','Zh','rh','m','Ash','Bsh','Csh','Dsh','Esh','Fsh'):
# print(Zreg['tethah'])
