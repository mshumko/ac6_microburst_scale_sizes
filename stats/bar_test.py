# Test code to make the global distribution of X using a dial plot.

import matplotlib.pyplot as plt 
import matplotlib.patches
import matplotlib
import numpy as np

r = np.arange(3, 8)
N = 24
width = 2*np.pi/N
theta = width*np.arange(0, N)

rr, tt = np.meshgrid(r, theta)

rr = rr.flatten()
tt = tt.flatten()

ax = plt.subplot(111, polar=True)
ax.set_theta_zero_location("S")
bars = ax.bar(tt, np.ones_like(tt), width, rr)

for (i, bar) in enumerate(bars):
    bar.set_facecolor(plt.cm.plasma(5*i+1))

cmap = matplotlib.cm.plasma
norm = matplotlib.colors.Normalize(vmin=0, vmax=10)

cb1 = matplotlib.colorbar.ColorbarBase(ax, cmap=cmap,
                                norm=norm,
                                orientation='horizontal')




# Draw dark side of Earth
patch = matplotlib.patches.Rectangle((3*np.pi/2,0), np.pi, 1, color='k')
ax.add_artist(patch)
plt.plot(np.linspace(0, 2*np.pi, N), np.ones(N), 'k')

# # Plot the Earth
# angle = 0
# center=(0,0)
# radius=1
# colors = ('w', 'k')
# # theta1, theta2 = angle, angle + 180
# theta1, theta2 = angle, angle + 2*np.pi
# w1 = Wedge(center, radius, theta1, theta2, fc=colors[0])
# w2 = Wedge(center, radius, theta2, theta1, fc=colors[1])
# for wedge in [w1, w2]:
#     ax.add_artist(wedge)

plt.show()