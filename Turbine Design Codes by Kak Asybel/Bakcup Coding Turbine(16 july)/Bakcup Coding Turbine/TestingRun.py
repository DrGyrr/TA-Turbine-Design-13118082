import numpy as np
mhi=[]
mhi.append(0)
dmh=[1,2,3,4,5,6,7,8,9.10]
for i in range(1,len(dmh)):
    mhi.append(mhi[i-1]+dmh)