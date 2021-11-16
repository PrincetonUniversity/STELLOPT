def load_plotbox2y(x_data=None, y1_data=None, y2_data=None,\
         x_label='x', y_label='y', key1=None, key2=None, x_range=None,\
         y_range=None, fig_size=(8.5, 7.5), wp=0, ax=None, mkt="o", mkc='red', title='',\
         hline1=None, hline2=None, vshadow=None, ls1=1, ls2=2, sciform=True):

    # Create a new figure if no existing axis was give 
    if not ax:
        f, ax = plt.subplots(figsize=fig_size)

    # Change the appearance of the figure
    ax.grid(color='grey', linestyle='-', linewidth=0.3)  
    ax.legend(loc=1,labelspacing=0.0)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_xlim(x_range)
    ax.set_ylim(y_range)
    ax.set_title(title)
    
    if hline1 != None:
        ax.axhline(y=hline1, linestyle='-', linewidth=4, color='black')
    if hline2 != None:
        ax.axhline(y=hline2, linestyle='-', linewidth=4, color='gray')
    if vshadow != None:
        xmin, xmax, ymin, ymax = ax.axis()
        ax.barh(ymin, vshadow[1]-vshadow[0], height=(ymax-ymin), \
                  left=vshadow[0], facecolor='yellow', alpha=0.3)
    if wp == 0:
        ax.plot(x_data, y1_data, ls(1)[0], color=ls(ls1)[2], linewidth=ls(ls1)[3], label=key1)
        ax.plot(x_data, y2_data, ls(2)[0], dashes=ls(ls2)[1], color=ls(ls2)[2], linewidth=ls(ls2)[3], label=key2)
    elif wp == 1:
        ax.plot(x_data, y2_data, ls(1)[0], color=ls(ls2)[2], linewidth=ls(ls2)[3], label=key2,\
                marker=mkt, markersize=10, markeredgewidth=2.5, markeredgecolor=ls(ls2)[2], markerfacecolor="white")
        ax.plot(x_data, y1_data, ls(1)[0], color=ls(ls1)[2], linewidth=ls(ls1)[3], label=key1,\
                marker=mkt, markersize=10, markeredgewidth=2.5, markeredgecolor=ls(ls1)[2], markerfacecolor="white")
        
    ax.yaxis.set_minor_locator(AutoMinorLocator(10))
    ax.xaxis.set_minor_locator(AutoMinorLocator(10))
    ax.legend(loc='best',labelspacing=0.0, prop={'size':2})
    if sciform: ax.ticklabel_format(style='sci', scilimits=(0,0))
    if sciform == False: ax.ticklabel_format(useOffset=False)
    return ax
