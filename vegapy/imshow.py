import numpy as np
import matplotlib.pyplot as plt

def imshow(image, scale='lin', **kwargs):
    if scale == 'log':
        try:
            image = np.log10(image)
        except TypeError as e:
            image = np.log10(image.value)
    elif scale == 'ln':
        image = np.log(image)
    if 'imshow' in kwargs:
        plt.imshow(image, **kwargs['imshow'])
    else:
        plt.imshow(image)
    if 'axis' in kwargs:
        plt.axis(**kwargs['axis'])
    if 'colorbar' in kwargs:
        plt.colorbar(**kwargs['colorbar'])
    plt.show()
