# this script explores the scale size dependence on L shell
import matplotlib.pyplot as plt
from datetime import datetime
import scale_size
import os

# Set L shell parameters
Lrange = range(3, 8)

for (lL, uL) in zip(Lrange[:-1], Lrange[1:]):
    # Catalog dir and name
    cDir = '/home/mike/research/ac6-microburst-scale-sizes/data/flash_catalogues'
    cName = 'flash_catalogue_v2_sorted.txt'
    # Normalization dir and name
    nDir = '/home/mike/research/ac6-microburst-scale-sizes/data/norm'
    nName = 'ac6_norm_{}_L_{}.txt'.format(lL, uL)
    
    ss = scale_size.ScaleSizeDist(os.path.join(cDir, cName), os.path.join(nDir, nName), 
                    fltDict={'burstType':'flash', 'Lm_OPQ':[lL,uL]})

    # Plot data
    fig, ax = plt.subplots(2, sharex=True)
    ss.plot_hist(ax=ax[0])
    ax[0].set_xlim(0, 100)
    ss.plot_norm(bx=ax[1])
    ax[0].set_title('Flash Scale Sizes | {} < L < {}'.format(lL, uL))
    ax[1].set_title('Normalization')

    plt.savefig(('/home/mike/Dropbox/0_grad_work/ac6-flashes-curtains/'
                'plots/{}_scale_sizes/'
                'flash_scale_size_{}_L_{}.png').format(datetime.now().date(), lL, uL))
