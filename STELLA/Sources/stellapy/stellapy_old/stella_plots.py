# -*- coding: utf-8 -*-

## some coding for commonly used plots

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

# setup some plot defaults
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size=30)
rcParams.update({'figure.autolayout': True})

def plot_1d(x,y,xlab,title='',ylab=''):

    fig = plt.figure(figsize=(12,8))
    plt.plot(x,y)
    plt.xlabel(xlab)
    if len(ylab) > 0:
        plt.ylabel(ylab)
    if len(title) > 0:
        plt.title(title)
    return fig

def logyplot_1d(x,y,xlab,title='',ylab=''):

    fig = plt.figure(figsize=(12,8))
    plt.semilogy(x,y)
    plt.xlabel(xlab)
    if len(ylab) > 0:
        plt.ylabel(ylab)
    if len(title) > 0:
        plt.title(title)
    return fig

def plot_2d(z,xin,yin,zmin,zmax,xlab='',ylab='',title='',cmp='RdBu'):

    fig = plt.figure(figsize=(12,8))
    x,y = np.meshgrid(xin,yin)
    plt.imshow(z, cmap=cmp, vmin=zmin, vmax=zmax,
               extent=[x.min(),x.max(),y.min(),y.max()],
               interpolation='nearest', origin='lower', aspect='auto')
    plt.axis([x.min(), x.max(), y.min(), y.max()])
    plt.colorbar()
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(title)
    return fig

def movie_2d(z,xin,yin,zmin,zmax,nframes,outfile,xlab='',ylab='',title='',step=1,cmp='RdBu'):

    from matplotlib import animation

    fig = plt.figure(figsize=(12,8))
    x,y = np.meshgrid(xin,yin)
    im = plt.imshow(z[0,:,:], cmap=cmp, vmin=zmin[0], vmax=zmax[0],
               extent=[x.min(),x.max(),y.min(),y.max()],
               interpolation='nearest', origin='lower', aspect='auto')
    plt.colorbar()
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(title)
    
    ims = []
    ims.append([im])

    for i in range(1,nframes,step):
        im = plt.imshow(z[i,:,:], cmap=cmp, vmin=zmin[i], vmax=zmax[i],
               extent=[x.min(),x.max(),y.min(),y.max()],
               interpolation='nearest', origin='lower', aspect='auto')
        ims.append([im])

    ani = animation.ArtistAnimation(fig,ims,interval=50,blit=True)
    ani.save(outfile)

def movie_1d(x,y,xmin,xmax,ymin,ymax,nframes):

    from matplotlib import animation

    fig = plt.figure(figsize=(12,8))
    ax=plt.axes(xlim=(xmin,xmax),ylim=(ymin,ymax))
    line, = ax.plot([],[],lw=2)

    def init():
        line.set_data([],[])
        return line,

    def animate(i):
        line.set_data(x,y[i,:])
        return line,

    anim=animation.FuncAnimation(fig, animate, init_func=init, 
                                 frames=nframes, interval=20)

    anim.save('test.mp4')
