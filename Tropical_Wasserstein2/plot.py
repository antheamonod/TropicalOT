import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np
import csv

directory="./data"

#--------------------------------------------------
#   Getting n1 and nt
#--------------------------------------------------

with open("{}/parameters.csv".format(directory)) as F:
    csvReader = csv.reader(F)
    for i in csvReader:
        n1 = int(i[0])
        n2 = int(i[1])
        nt = int(i[2])

#--------------------------------------------------
#   Getting Rho Data
#--------------------------------------------------

with open("{}/rho.csv".format(directory)) as F:
    for i in F:
        rho = list(map(lambda x: float(x), i[:-1].split(",")))

rho = np.array(rho).reshape(nt,n1,n2)



#--------------------------------------------------
#   Create animation
#--------------------------------------------------

# First set up the figure, the axis, and the plot element we want to animate
fig, ax = plt.subplots(1,1,figsize=(5,4))
cax = ax.imshow(rho[0], cmap='jet')
fig.colorbar(cax)
plt.axis('off')

# animation function.  This is called sequentially
def animate(n):
    # fig.clear()
    cax.set_array(np.flipud(rho[n]))
    cax.set_clim(np.min(rho[n]), np.max(rho[n]))
    plt.suptitle("t={:.4}\n{}".format(n/(nt-1),np.sum(rho[n])))
    return cax, 

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, 
                               frames=nt, interval=10, blit=True)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
anim.save("video.mp4", fps=5)