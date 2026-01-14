from matplotlib import pyplot as plt 
import numpy as np 

def plot_cutouts(images, zp, pxscale, imsize=2):

    fig, axs = plt.subplots(1, len(images), figsize=(imsize*len(images), imsize))
    for ax, img in zip(axs, images):
        if img.dtype in [int, bool]:
            ax.imshow(img, origin='lower', cmap='gray', vmin=0, vmax=np.max(img)*0.8)
        else:
            ax.imshow(-2.5*np.log10(np.abs(img)/pxscale**2)+zp, origin='lower', cmap='gray', vmin=19, vmax=27) 
    
    for ax in axs:
        ax.axis('off')
    plt.subplots_adjust(wspace=0, hspace=0, left=0, right=1, top=1, bottom=0)
    return axs