# import glob
# import os
# import time
# from ifrturbinepackage.definitions import *
# dir_name = ROOT_DIR

# dir_name = os.path.join(ROOT_DIR,"Outputs","FitInputData")
# # Get list of all files only in the given directory
# list_of_files = filter( lambda x: os.path.isfile(os.path.join(dir_name, x)),
#                         os.listdir(dir_name) )
# # Sort list of files based on last modification time in ascending order
# list_of_files = sorted( list_of_files,
#                         key = lambda x: os.path.getmtime(os.path.join(dir_name, x))
#                         )
# # Iterate over sorted list of files and print file path 
# # along with last modification time of file 
# for file_name in list_of_files:
#     file_path = os.path.join(dir_name, file_name)
#     timestamp_str = time.strftime(  '%m/%d/%Y || %H:%M:%S',
#                                 time.gmtime(os.path.getmtime(file_path))) 
#     print(timestamp_str, ' -->', file_name)  
# ////
# from ifrturbinepackage.definitions import *
# from ifrturbinepackage.inputs import *
# from ifrturbinepackage.rotor import *

# result=ComputeR3(1.53,15,10,6,7)
# print(result['efficiency'])
# print(result['geometry'])
# print(result['velocity'])
# print(result['thermo'])
# ////

# from CoolProp.CoolProp import PropsSI as Props
# from ifrturbinepackage.rotor import *
# fluid='R32'
# res=Props('S','P',2.15*10**6,'T',350,fluid)
# res2=Phase('P',2.15*10**6,'T',350,fluid)
# print(f"S={res}:{res2}")
# satpv=getsatpvfors(2.15*10**6,res,fluid)
# print(f"found at P={satpv}")
# print(Props('Q','P',satpv['P'],'S',res,fluid))

# # pmin=Props('P_min',fluid)
# # pmax=Props('Pcrit',fluid)
# print(pmin,pmax)

# x=3500.5
# y=4.7
# print(x//y)
# print(Props('D','T',425.22,'P',0.82567*10**6,'R32'))


''' Preview 2D Contour for n=(2,9)'''
# %%
from ifrturbinepackage.definitions import *
from ifrturbinepackage.inputs import *
from ifrturbinepackage.rotor import *

# %%
sol=ComputeR4(1.2,15,10,1)
print(sol)



# %%
