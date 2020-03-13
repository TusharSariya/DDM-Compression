from PIL import Image
import numpy as np
from matplotlib import pyplot as plt
from sys import getsizeof

#read data
DMMData = np.loadtxt("DMMRaw2.txt", dtype='i', delimiter=',')

#plt.imshow(DMMData, interpolation='nearest')
#plt.savefig('DMMImage.png')
