import numpy as np # for maths 
import matplotlib as mpl # for plotting 
import matplotlib.pyplot as plt

from tqdm import tqdm # tqdm is a package that lets you make progress bars to see how a loop is going

import os 

# configure notebook for plotting
#matplotlib inline

mpl.style.use('seaborn-colorblind') # colourblind-friendly colour scheme

# subsequent lines default plot settings
mpl.rcParams['image.origin'] = 'lower'
mpl.rcParams['figure.figsize']=(8.0,6.0)   
mpl.rcParams['font.size']=16              
mpl.rcParams['savefig.dpi']= 300             

import warnings
warnings.filterwarnings('ignore')

## GENERATE FAKE DATA
x = np.linspace(0,2.5,100) # Generates 100 points evenly spaced between values 0 and 2.5
noise = 2*np.random.randn(len(x))
y = 12*x -5  + noise
plt.plot(x,12*x-5,'-k',label='True')
plt.plot(x,y,'.',label='Data')
plt.legend()
plt.show()

A = np.vander(x,2) # the Vandermonde matrix of order N is the matrix of polynomials of an input vector 1, x, x**2, etc
b, residuals, rank, s = np.linalg.lstsq(A,y)
print('True parameters: 12, -5. Recovered parameters: %.2f, %.2f' % (b[0],b[1]))
reconstructed = A @ b # @ is shorthand for matrix multiplication in python
plt.plot(x,y,'.',label='Data')
plt.plot(x,reconstructed,'-r',label='Reconstructed')
plt.plot(x,12*x-5,'-k',label='True')
plt.legend()
plt.show()