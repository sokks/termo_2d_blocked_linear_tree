import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
import sys
from scipy.interpolate import griddata
from matplotlib.colors import ListedColormap


cmap_name = 'coolwarm'
filename = 'data/temp/000000.out.txt'
plot_dat = 'temp'
max_lvl = 9

max_lvl   = int(sys.argv[1])
filename  = sys.argv[2]
imgfile   = sys.argv[3]
cmap_name = sys.argv[4]
plot_dat  = sys.argv[5]

cells = np.loadtxt(filename, delimiter=',')
# print(cells)

L = cells[:,0]
X = cells[:,1]
Y = cells[:,2]
B = cells[:,3]
P = cells[:,4]
Z = cells[:,5]



# define grid
xi = np.linspace(0,2**max_lvl,2**max_lvl)
yi = np.linspace(0,2**max_lvl,2**max_lvl)


cmp = cm.get_cmap(cmap_name)
tcks = []
cbar_label = ''
v_min = 0
v_max = 10
ttl = ''

zi = griddata((Y, X), Z, (xi[None,:], yi[:,None]), method='linear')
if plot_dat == 'temp':
    cbar_label = 'temperature'
    v_min = 0
    v_max = np.round(np.max(Z))
    tcks = np.linspace(v_min, v_max, num=11)
    # tcks = [i for i in range(int(v_min), int(v_max) + 1)]
    plt.xticks(np.linspace(0, 2**max_lvl, num=6), ['{:1.1f}'.format(x) for x in np.linspace(0, 1, num=6)])
    plt.yticks(np.linspace(0, 2**max_lvl, num=6), ['{:1.1f}'.format(x) for x in np.linspace(0, 1, num=6)])
    if 'base_grid' in filename:
        ttl = 'Начальное распределение температуры'

if plot_dat == 'lvls':
    zi = griddata((Y, X), L, (xi[None,:], yi[:,None]), method='nearest')
    viridis = cm.get_cmap(cmap_name, 12)
    v_min = L.min()
    v_max = L.max()
    newcolors = viridis(np.linspace(0, 1, int(v_max - v_min) + 1))
    cmp = ListedColormap(newcolors)
    tcks = [i for i in range(int(v_min), int(v_max) + 1)]
    cbar_label = 'grid level'


elif plot_dat == 'procs':
    n_procs = int(sys.argv[6])
    zi = griddata((Y, X), P, (xi[None,:], yi[:,None]), method='nearest')
    viridis = cm.get_cmap(cmap_name, 100)
    v_min = P.min()
    v_max = P.max()
    newcolors = viridis(np.linspace(0, 1, int(v_max - v_min) + 1))
    cmp = ListedColormap(newcolors)
    tcks = [i for i in range(n_procs)]
    cbar_label = 'processor number'


elif plot_dat == 'blks':
    zi = griddata((Y, X), B, (xi[None,:], yi[:,None]), method='nearest')
    viridis = cm.get_cmap(cmap_name, 100)
    v_min = B.min()
    v_max = B.max()
    newcolors = viridis(np.linspace(0, 1, int(v_max - v_min) + 1))
    cmp = ListedColormap(newcolors)
    # tcks = [i for i in range(int(v_min), int(v_max) + 1)]
    cbar_label = 'blocks'


pl = plt.pcolor(xi, yi, zi, cmap=cmp, vmin=v_min, vmax=v_max, clip_on=False)
# plt.axis('off')
pl.axes.set_frame_on(False)
cbar = None
print(tcks)
if plot_dat == 'temp' and tcks.any() or type(tcks) == list and tcks:
    cbar = plt.colorbar(pl, ticks=tcks)
else:
    cbar = plt.colorbar(pl)
if cbar_label != '':
    cbar.ax.get_yaxis().labelpad = 20
    cbar.ax.set_ylabel(cbar_label, rotation=270)
# plt.show()
# plt.axes([0,0,1,1], frameon=False)

if ttl:
    plt.title(ttl)

plt.savefig(imgfile, transparent=True)


# plt.scatter(Y, X, c=Z, cmap=plt.get_cmap(cmap_name))
# plt.tripcolor(X, Y, L, cmap="RdBu_r")
# ax2.tricontour(x, y, z, levels=14, linewidths=0.5, colors='k')
# plt.show()

# fig, ax = plt.subplots()
# p = ax.pcolor(X, Y, Z, cmap=cm.RdBu, vmin=abs(Z).min(), vmax=abs(Z).max())
# cb = fig.colorbar(p)
# plt.show()

