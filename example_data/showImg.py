import sys
import numpy as np
import matplotlib.pyplot as plt

if __name__=="__main__":
    """A quick and dirty viewer for ascii
    data in a matrix
    """
  
    img = np.loadtxt(sys.argv[1])
    plt.imshow(img,interpolation="none")
    plt.show()
  

