from tqdm import tqdm
x=0
for i in tqdm(range(10)):
    x=x**(x**x)
    print(i)