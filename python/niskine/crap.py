# library & dataset
from matplotlib import pyplot as plt
import numpy as np
 
#colormap='RdBu_r'

colormap = plt.cm.get_cmap('RdBu_r')
color_sn = colormap(0)
color_in = colormap(51)
color_wn = colormap(102)
color_wp = colormap(153)
color_ip = colormap(204)
color_sp = colormap(255)

# create data
x = np.random.rand(15)
y = x+np.random.rand(15)
z = x+np.random.rand(15)
z=z*z
 
x=np.linspace(0,1,2)

# Use it with a call in cmap
#plt.scatter(x, y, s=z*2000, c=1, cmap=colormap, alpha=0.4, edgecolors="grey", linewidth=2)
#plt.scatter(x, 2.*y, s=z*2000, c=2, cmap=colormap, alpha=0.4, edgecolors="grey", linewidth=2)
#plt.scatter(x, 3.*y, s=z*2000, c=4, cmap=colormap, alpha=0.4, edgecolors="grey", linewidth=2)

plt.plot(x,1*x,color=color_sn, linewidth=2)
plt.plot(x,2*x,color=color_in, linewidth=2)
plt.plot(x,3*x,color=color_wn, linewidth=2)
plt.plot(x,4*x,color=color_wp, linewidth=2)
plt.plot(x,5*x,color=color_ip, linewidth=2)
plt.plot(x,6*x,color=color_sp, linewidth=2)
 
plt.show()

# You can reverse it:
#plt.scatter(x, y, s=z*2000, c=x, cmap="BuPu_r", alpha=0.4, edgecolors="grey", linewidth=2)
 
# OTHER: viridis / inferno / plasma / magma
#plt.scatter(x, y, s=z*2000, c=x, cmap="plasma", alpha=0.4, edgecolors="grey", linewidth=2)
