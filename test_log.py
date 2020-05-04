import numpy as np
   
frames = 10**5
speed_mul = 100
size_mul = 1.2

for i in np.arange(round(frames/speed_mul)):
    if np.log10(i) < 1:
        speed_mul = 10
    speed_mul np.log10(i)= np.log10(i)
    print (round(frames/speed_mul))