import matplotlib
matplotlib.use("TKAgg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import user


def plot_3d(data, name, color):

    X = np.linspace(0, 5, np.shape(data)[0])
    Y = np.linspace(0, 5, np.shape(data)[1])
    x, y = np.meshgrid(X, Y)

    mappable = plt.cm.ScalarMappable(cmap=color)
    mappable.set_array(data)
    #mappable.set_clim(0.5,1.1)

    fig = plt.figure(figsize=(10, 4))
    ax1 = fig.add_subplot(122, projection="3d")
    ax1.plot_surface(x, y, data, cmap=mappable.cmap, norm=mappable.norm,
                     linewidth=0, antialiased=False)
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    ax1.set_zlabel("$"+name+"$"+" $m/s$")
    ax1.set_xlim(np.min(X), np.max(X))
    ax1.set_ylim(np.min(Y), np.max(Y))
    #ax1.set_zlim(0.5,1.1)

    ax2 = fig.add_subplot(121)
    ax2.imshow(data, cmap=mappable.cmap, norm=mappable.norm,
               extent=(np.min(x),
                       np.max(x),
                       np.min(y),
                       np.max(y)), interpolation = "none")
    plt.savefig(name, dpi=1000)
    plt.colorbar(mappable)
    plt.tight_layout()


def Animation(cube, name, color, time):
    fig = plt.figure()
    # x = np.arange(np.shape(cube)[0])
    # y = np.arange(np.shape(cube)[0])

    mappable = plt.cm.ScalarMappable(cmap=color)
    mappable.set_array(cube[-1])
    # mappable.set_clim(0,1)
    ims = []
    for i in range(len(cube)):
        im = plt.imshow(cube[i], cmap=mappable.cmap,
                        interpolation="gaussian", animated=True)
        im1 = plt.text(user.Nx, user.Ny + 5, "time = {:.3f} s".format(time[i]))
        ims.append([im, im1])
        #plt.savefig("imagen_"+str(i))
    ani = animation.ArtistAnimation(fig, ims, interval=50, blit=True,
                                    repeat_delay=1000)
    #plt.streamplot()
    plt.colorbar(mappable)
    plt.tight_layout()
    plt.xlabel("$x$")
    plt.ylabel("$y$")
    plt.title("2D field solution")
    ani.save(name+".gif", dpi=80, writer="imagemagick")
    # plt.show()


