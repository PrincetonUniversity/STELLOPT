from numpy import *
from tkinter import *
from tkinter.filedialog import *
from matplotlib import *
from os import listdir
from scipy.io import netcdf
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import math
from matplotlib.ticker import AutoMinorLocator
from scipy import interpolate

#sys.path.append('/home/antonio/PROGRAMAS/stella/stellapy')

#import stella_diag
#from stella_diag import *
#import stella_read
#from stella_read import *


############################################ Some utils ##################################################################
def format1(value):
    return "%.3e" % value
def format2(value):
    return "%14.6e" % value
def format3(value):
    return "%4.2f" % value
def format4(value):
    return "%6.2f" % value
def format6(value):
    return "%7.3f" % value
def format5(value):
    return "%.5e" % value
def format7(value):
    return "%22.3f" % value
def format8(value):
    return "%04d" % value
def format9(value):
    return "%7.5f" % value
# Some utils ended
##########################################################################################################################

################################################ PLOT FUNCTIONS ##########################################################

##########################################################################################################################

#Plot several curves in the same graph

def pl2d(xrange=None, yrange=None, xlabel=None, ylabel=None,\
         fig_size=(8.5, 7.5), title=None, ax=None):
        
    if not ax:
        f, ax = plt.subplots(figsize=fig_size)
    ax.xaxis.grid(color='grey', linestyle='-', linewidth=0.3)
    ax.yaxis.grid(color='grey', linestyle='-', linewidth=0.3)  
    ax.set_xlabel(xlabel, fontsize=40)
    ax.set_ylabel(ylabel, fontsize=40)
    ax.set_xlim(xrange)
    ax.set_ylim(yrange)
    ax.ticklabel_format(style='sci', scilimits=(0,0))
    ax.set_title(title, fontsize=40)

    return ax


#########################################################################################################################

#To make subplots

def plxy(x1data=None, y1data=None, \
         xlabel='x', ylabel='y', key1=None, xrange=None,\
         yrange=None,fig_size=(8.5, 7.5), title=None, ax=None):

    ax.grid(color='grey', linestyle='-', linewidth=0.3)
    ax.plot(x1data, y1data, linestyle='-', color='red', linewidth=3.5, label=key1)
    ax.set_xlabel(xlabel, fontsize=25)
    ax.set_ylabel(ylabel, fontsize=25)
    ax.set_xlim(xrange)
    ax.set_ylim(yrange)
    ax.ticklabel_format(style='sci', scilimits=(0,0))

    return ax

########################################################################################################################

############################################### USEFUL FUNCTIONS #######################################################

########################################################################################################################

################################ To create a directory (if it does not exist) ##########################################

def createFolder(directory):
        if not os.path.exists(directory):
            os.mkdir(directory)

########################################## To quit the extension #######################################################

#To quit a directory extension

def separate(case):
    dis=case.split(".")
    dim=size(dis)
    element=[]
    use=dis[0]
            #To avoid . problems
    for i in arange(1,dim-1):
        element = dis[i]
        use=use + "." +element
    return use

########################################################################################################################

def normali(real,imaginary):
    r=size(real)
    l=size(imaginary)
    if r==l:
        absolu=empty((r))
        for i in arange(0,size(imaginary)):
	        absolu[i]=sqrt((imaginary[i]*imaginary[i])+(real[i]*real[i]))
        return absolu
    else:
        print("Vectors must have the same length")


########################################################################################################################

def factordivide(v,f):
    s=size(v)   
    divi=empty((s))
    for i in arange(0,s):
        divi[i]=v[i]/f
    return divi


#######################################################################################################################

def factormult(v,f):
    s=size(v)
    mult=empty((s))
    for i in arange(0,s):
        mult[i]=v[i]*f
    return mult    

######################################################################################################################

def vecpot(v,p):
    s=size((v))
    pot=ones((s))
    for i in arange(0,s):
        for l in arange(0,p):
            pot[i]=pot[i]*v[i]
    return pot

#######################################################################################################################

def vecabs(v):
    s=size(v)
    ab=empty((s))
    for i in arange(0,s):
        ab[i]=abs(v[i])
    return ab
############################################### GUI FUNCTIONS #########################################################

############################################ Information Window #######################################################

def extra_win(geom,textinf):
        Window2 = Toplevel() 
        Window2.title("Information")
        Window2.geometry(geom)
        Inform=Label(Window2,text=textinf, width=40, height=5).place(x=40,y=20)
        Exit=Button(Window2,text='Accept',command=Window2.destroy).place(x=170,y=160)

#######################################################################################################################
 
#############################################  STELLA FUNCTIONS  #####################################################

#################################### omega(t) plot #############################################################

def omega_t_plot(file_list):
 
        a=  file_list
        case=[]
        colors=[]
        col=['r','g','b','y','m','y','k']  
        times=[]
        Im_max=[]
        Im_min=[]
        Re_max=[]
        Re_min=[]
        xlab  = "$t$"
        ylab  = "$\\omega a/v_{\mathrm{th},\mathrm{i}}$"
        tit1  = '$\\Re(\\omega)/v_{\mathrm{th},\mathrm{i}}$'
        tit2  = '$\\Im(\\omega)/v_{\mathrm{th},\mathrm{i}}$'    
 
        for i in arange(0,100):
            colors=colors+col           
        for i in arange(0,size(a)):
           element=separate(a[i])+'.omega' 
           case=case+[element] 


        for l in arange(0,size(a)):
            dat=loadtxt(case[l],dtype='float')
            time_grap=dat[:,0]
            times=times+[max(time_grap)]
            imaginary_grap=dat[:,4]
            Im_max=Im_max+[max(imaginary_grap)]
            Im_min=Im_min+[min(imaginary_grap)]
            real_grap=dat[:,3]
            Re_max=Re_max+[max(real_grap)]
            Re_min=Re_min+[min(real_grap)]

        ax1 = pl2d(xrange=[min(time_grap)-0.2,max(times)],yrange=[min(Re_min)-0.1*max(vecabs(Re_max)),max(Re_max)+0.1*max(vecabs(Re_max))],xlabel=xlab, ylabel=ylab, fig_size=(8.5, 7.5), title=tit1)
        ax2 = pl2d(xrange=[min(time_grap)-0.2,max(times)],yrange=[min(Im_min)-0.1*max(vecabs(Im_max)),max(Im_max)+0.1*max(vecabs(Im_max))],xlabel=xlab, ylabel=ylab, fig_size=(8.5, 7.5), title=tit2)

        for l in arange(0,size(a)):
            dat=loadtxt(case[l],dtype='float')
            imaginary=dat[:,4]
            real=dat[:,3]
            time=dat[:,0]
            ky=str(dat[0,1])
            kx=str(dat[0,2])   
            key="$k_y \\rho$=" + ky + "   ," + "   $k_x \\rho=$" + kx  
            
            
            ax1.plot(time, real, linestyle='-', marker="o",\
                markersize=5, color=colors[l], markerfacecolor="white", linewidth=4, label=key)
            ax2.plot(time, imaginary, linestyle='-', marker="s",\
                markersize=5, color=colors[l], markerfacecolor="white", linewidth=4, label=key)
            ax1.legend(loc='best',labelspacing=0.0, prop={'size':24})      
            ax2.legend(loc='best',labelspacing=0.0, prop={'size':24}) 


    



############################################### omega(k) plot #########################################################

def omega_k_plot(file_ext_free):

       case=file_ext_free+'.omega'
       data=loadtxt(case,dtype='float')
       ky=data[0,1]
       im=data[:,4]
       re=data[:,3]
       real=[]
       imaginary=[]
       si=size(im)
       sr=size(re)
       for i in arange(0,sr):
            if math.isnan(re[i]) != TRUE:
                element=re[i]
                real=real+[element]
       for j in arange(1,si):
            if math.isnan(im[j]) != TRUE:
                element=im[j]
                imaginary=imaginary+[element]
       realpoint=real[-1]
       impoint=imaginary[-1]
       limits=append(impoint,realpoint)
       alim=abs(limits)
       xlab="$k_y\\rho$"
       ylab="$\\omega a/v_{\mathrm{th},\mathrm{i}}$"
       key1= '$\\Re(\\omega)/v_{\mathrm{th},\mathrm{i}}$'
       key2= '$\\Im(\\omega)/v_{\mathrm{th},\mathrm{i}}$'

       ax = pl2d(xrange=[ky-0.5,ky+0.5], yrange=[min(limits)-0.2*max(alim),max(limits)+0.2*max(alim)], xlabel=xlab, ylabel=ylab,\
           fig_size=(8.5, 7.5))
       ax.plot(ky, realpoint, marker="o",\
           markersize=12, color="black", markerfacecolor="black", label=key1)
       ax.plot(ky, impoint, marker="o",\
           markersize=12, color="blue", markerfacecolor="blue", label=key2)
       ax.legend(loc='best',labelspacing=0.0, prop={'size':24})      

##############################################  multi omega(k) plot ###################################################

def multi_omega_k_plot(file_list):
            

        a=[]
        case=[]
        data=empty((size(file_list),7))
        om=empty((size(file_list),7)) 
        num=[]                        
        for j in arange(0,size(file_list)):
           element=separate(file_list[j])+'.omega'
           case=case+[element]
        for i in arange(0,size(file_list)):
            dat=loadtxt(case[i],dtype='float')
            #om[i,:]=dat[-1]
            for n in arange(0,size(dat)/7):
                    row=dat[int(n)]
                    if math.isnan(row[3]) == FALSE or math.isnan(row[4]) == FALSE:                 
                        om[i,:]=row


        imaginary=abs(om[:,4])
        real=abs(om[:,3])
        ky=om[:,1]
        kx=str(om[0,2])  
        xlab  = "$k_{y}\\rho_i$"
        y1lab = "$a|\\Re(\\omega)|/v_{\mathrm{th},\mathrm{i}}$"
        y2lab = "$a|\\Im(\\omega)|/v_{\mathrm{th},\mathrm{i}}$"
        tit   = "$\\omega(k_{y})$ "+"for "+"$k_{x}\\rho_{i}=$"+kx          
        ax1 = pl2d(xrange=[0,max(ky)+0.1], yrange=[min(real)-0.1*max(abs(real)),max(real)+0.1*max(abs(real))], xlabel=xlab, ylabel=y1lab,\
                    fig_size=(8.5, 7.5), title=tit)
        ax1.plot(ky, real, linestyle='-', marker='o', markersize=5, color="purple", markerfacecolor="white", linewidth=4)

        ax2 = pl2d(xrange=[0,max(ky)+0.1], yrange=[min(imaginary)-0.1*max(abs(imaginary)),max(imaginary)+0.1*max(abs(imaginary))], xlabel=xlab, ylabel=y2lab,\
                    fig_size=(8.5, 7.5), title=tit)
        ax2.plot(ky, imaginary, linestyle='-', marker='o', markersize=5, color="green", markerfacecolor="white", linewidth=4)  


################################################### Potential plot  #####################################################

def potential_plot(file_exten_free):

        case=file_exten_free+'.final_fields'
        data=loadtxt(case,dtype='float')
        ky=str(data[0,2])
        kx=str(data[0,3])       
        imaginary=data[:,5]
        real=data[:,4]
        zed=factormult(data[:,0],3.282)
        limits=append(real,imaginary)
        alim=abs(limits)
 
        tit   = "$k_y \\rho$=" + ky + "   ," + "$   k_x \\rho=$" + kx
        xlab  = '$\\zeta$ ($\\cdot \\pi$)'
        ylab  = '$\\tilde{\\varphi}(\\zeta)$'
  
        key1='$\\Re(\\tilde{\\varphi})(\\zeta)$'
        key2='$\\Im(\\tilde{\\varphi})(\\zeta)$'
    
        ax = pl2d(xrange=[min(zed),max(zed)], yrange=[min(limits)-0.2*max(alim),max(limits)+0.2*max(alim)], xlabel=xlab, ylabel=ylab,\
            fig_size=(8.5, 7.5), title=tit)
        ax.plot(zed, real, linestyle='-', marker="o",\
            markersize=5, color="black", markerfacecolor="white", linewidth=4, label=key1)
        ax.plot(zed, imaginary, linestyle='-', marker="s",\
            markersize=5, color="blue", markerfacecolor="white", linewidth=4, label=key2)
        ax.legend(loc='best',labelspacing=0.0, prop={'size':24})

################################################ normalized potential plot ###############################################

def potential_norm_plot(file_list):

        a=  file_list
        case=[]
        colors=[]
        col=['r','g','b','y','m','y','k']  
        max_pot_norm=[]
        xlab  = '$\\zeta$'
        ylab  = '$\\tilde{\\varphi}(\\zeta)/\\tilde{\\varphi}(\\zeta=0)$'
        tit = 'Normalized Potential'  
 
        for i in arange(0,100):
           colors=colors+col           
        for i in arange(0,size(a)):
           element=separate(a[i])+'.final_fields' 
           case=case+[element] 
        for n in arange(0,size(a)):
            data=loadtxt(case[n],dtype='float')
            im=data[:,5]
            re=data[:,4]
            zeta=factormult(data[:,0],3.282) 
            nz0_graph=int(size(zeta)/2)
            potential_graph=normali(re,im)
            p0_g=potential_graph[nz0_graph]
            potential_norm_graph=factordivide(potential_graph,p0_g)
            max_pot_norm=max_pot_norm+[max(potential_norm_graph)]

        ax = pl2d(xrange=[-10,10], yrange=[0,max(max_pot_norm)], xlabel=xlab, ylabel=ylab,fig_size=(8.5, 7.5),title=tit)

        for l in arange(0,size(a)):
            dat=loadtxt(case[l],dtype='float')
            ky=str(dat[0,2])
            kx=str(dat[0,3])
            imaginary=dat[:,5]
            real=dat[:,4]
            zed=factormult(dat[:,0],3.282) 
            nz0=int(size(zed)/2)
            potential=normali(real,imaginary)
            p0=potential[nz0]
            potential_norm=factordivide(potential,p0)
 
            key   = "$k_y \\rho$=" + ky + "   ," + "$   k_x \\rho=$" + kx
            ax.plot(zed, potential_norm, linestyle='-', color=colors[l], linewidth=3,label=key)
            ax.legend(loc='best',labelspacing=0.0, prop={'size':24})

################################################# geo together plot ######################################################

def geo_tg(file_free_extend):

        out=file_free_extend+'.vmec_geo'
        data=loadtxt(out,dtype='float',skiprows=4)
        zeta      = data[:,1]
        bmag      = data[:,2]
        gradpar   = data[:,3]
        gbdrift   = data[:,9]
        gbdrift0  = data[:,10]
        cvdrift   = data[:,11]
        cvdrift0  = data[:,12]
        gs2      = data[:,4]
        gds21     = data[:,5]
        gds22     = data[:,6]
    
        labbmag     = '$B/B_{\mathrm{ref}}$'
        labgradpar  = '$L_{\\mathrm{ref}}\\nabla_{\|} z$'
        labgbdrift  = '$2B_{\\mathrm{ref}}L_{\mathrm{ref}}\\mathbf{B}\\times\\nabla B\\cdot\\nabla y/B^3$'
        labgbdrift0 = '$\\hat{s}2B_{\mathrm{ref}}L_{\mathrm{ref}}\\mathbf{B}\\times\\nabla B\\cdot\\nabla x/B^3$'
        labcdrift   = '$2B_{\mathrm{ref}}L_{\mathrm{ref}}\\mathbf{B}\\times\mathbf{\kappa}\\cdot\\nabla y/B^2$'
        labcdrift0  = '$\\hat{s}2B_{\mathrm{ref}}L_{\mathrm{ref}}\\mathbf{B}\\times\\kappa\\cdot\\nabla x/B^2$'
        labgs2   = '$|\\nabla y|^2$'
        labgds21 = '$\hat{s}^2\\nabla x\\cdot\\nabla y$'
        labgds22 = '$\hat{s}^2|\\nabla x|^2$' 
    
        fig1 = plt.figure(figsize=(18, 12))
        fig1.subplots_adjust(hspace=0.6, wspace=0.6)
    
        ax1 = fig1.add_subplot(231)
        ax2 = fig1.add_subplot(232)
        ax3 = fig1.add_subplot(233)
        ax4 = fig1.add_subplot(234)
        ax5 = fig1.add_subplot(235)
        ax6 = fig1.add_subplot(236)

        plxy(zeta, bmag, xrange=[min(zeta),max(zeta)], yrange=[min(bmag)-0.2*max(abs(bmag)),max(bmag)+0.2*max(abs(bmag))], xlabel='$\\zeta$', ylabel=labbmag,      title='\\texttt{bmag}'   , ax=ax1)
        plxy(zeta, gradpar, xrange=[min(zeta),max(zeta)], yrange=[min(gradpar)-0.2*max(abs(gradpar)),max(gradpar)+0.2*max(abs(gradpar))], xlabel='$\\zeta$', ylabel=labgradpar,   title='\\texttt{gradpar}', ax=ax2)
        plxy(zeta, gbdrift, xrange=[min(zeta),max(zeta)], yrange=[min(gbdrift)-0.2*max(abs(gbdrift)),max(gbdrift)+0.2*max(abs(gbdrift))], xlabel='$\\zeta$', ylabel=labgbdrift,   title='\\texttt{gbdrift}', ax=ax3)
        plxy(zeta, gbdrift0, xrange=[min(zeta),max(zeta)], yrange=[min(gbdrift0)-0.2*max(abs(gbdrift0)),max(gbdrift0)+0.2*max(abs(gbdrift0))], xlabel='$\\zeta$', ylabel=labgbdrift0,  title='\\texttt{gbdrift0}',ax=ax4)
        plxy(zeta, cvdrift, xrange=[min(zeta),max(zeta)], yrange=[min(cvdrift)-0.2*max(abs(cvdrift)),max(cvdrift)+0.2*max(abs(cvdrift))], xlabel='$\\zeta$', ylabel=labcdrift,    title='\\texttt{cdrift}',  ax=ax5)
        plxy(zeta, cvdrift0, xrange=[min(zeta),max(zeta)], yrange=[min(cvdrift0)-0.2*max(abs(cvdrift0)),max(cvdrift0)+0.2*max(abs(cvdrift0))], xlabel='$\\zeta$', ylabel=labcdrift0,   title='\\texttt{cdrift0}', ax=ax6)    

        fig2 = plt.figure(figsize=(18, 5))
        fig2.subplots_adjust(wspace=0.6)
    
        ax7 = fig2.add_subplot(131)
        ax8 = fig2.add_subplot(132)
        ax9 = fig2.add_subplot(133)
   
        plxy(zeta, gs2, xrange=[min(zeta),max(zeta)], yrange=[min(gs2)-0.2*max(abs(gs2)),max(gs2)+0.2*max(abs(gs2))], xlabel='$\\zeta$', ylabel=labgs2,      title='\\texttt{gs2}'   , ax=ax7)
        plxy(zeta, gds21, xrange=[min(zeta),max(zeta)], yrange=[min(gds21)-0.2*max(abs(gds21)),max(gds21)+0.2*max(abs(gds21))], xlabel='$\\zeta$', ylabel=labgds21,    title='\\texttt{gds21}',   ax=ax8)
        plxy(zeta, gds22, xrange=[min(zeta),max(zeta)], yrange=[min(gds22)-0.2*max(abs(gds22)),max(gds22)+0.2*max(abs(gds22))], xlabel='$\\zeta$', ylabel=labgds22,    title='\\texttt{gds22}',   ax=ax9)

###################################################  geo split plot #####################################################

def geo_sp(file_free_exten):

        out=file_free_exten+'.vmec_geo'
        data=loadtxt(out,dtype='float',skiprows=4)
        zeta      = data[:,1]
        bmag      = data[:,2]
        gradpar   = data[:,3]
        gbdrift   = data[:,9]
        gbdrift0  = data[:,10]
        cvdrift   = data[:,11]
        cvdrift0  = data[:,12]         
        gs2      = data[:,4]
        gds21     = data[:,5]
        gds22     = data[:,6]
    
        labbmag     = '$B/B_{\mathrm{ref}}$'
        labgradpar  = '$L_{\\mathrm{ref}}\\nabla_{\|} z$'
        labgbdrift  = '$2B_{\\mathrm{ref}}L_{\mathrm{ref}}\\mathbf{B}\\times\\nabla B\\cdot\\nabla y/B^3$'
        labgbdrift0 = '$\\hat{s}2B_{\mathrm{ref}}L_{\mathrm{ref}}\\mathbf{B}\\times\\nabla B\\cdot\\nabla x/B^3$'
        labcdrift   = '$2B_{\mathrm{ref}}L_{\mathrm{ref}}\\mathbf{B}\\times\mathbf{\kappa}\\cdot\\nabla y/B^2$'
        labcdrift0  = '$\\hat{s}2B_{\mathrm{ref}}L_{\mathrm{ref}}\\mathbf{B}\\times\\kappa\\cdot\\nabla x/B^2$'
        labgs2   = '$|\\nabla y|^2$'
        labgds21 = '$\hat{s}^2\\nabla x\\cdot\\nabla y$'
        labgds22 = '$\hat{s}^2|\\nabla x|^2$' 

        ax1 = pl2d(xrange=[min(zeta),max(zeta)], yrange=[min(bmag)-0.2*max(abs(bmag)),max(bmag)+0.2*max(abs(bmag))], xlabel='$\\zeta$', ylabel=labbmag,      title='\\texttt{bmag}')
        ax2 = pl2d(xrange=[min(zeta),max(zeta)], yrange=[min(gradpar)-0.2*max(abs(gradpar)),max(gradpar)+0.2*max(abs(gradpar))],xlabel='$\\zeta$', ylabel=labgradpar,   title='\\texttt{gradpar}')
        ax3 = pl2d(xrange=[min(zeta),max(zeta)], yrange=[min(gbdrift)-0.2*max(abs(gbdrift)),max(gbdrift)+0.2*max(abs(gbdrift))],xlabel='$\\zeta$', ylabel=labgbdrift,   title='\\texttt{gbdrift}')
        ax4 = pl2d(xrange=[min(zeta),max(zeta)], yrange=[min(gbdrift0)-0.2*max(abs(gbdrift0)),max(gbdrift0)+0.2*max(abs(gbdrift0))],xlabel='$\\zeta$', ylabel=labgbdrift0,  title='\\texttt{gbdrift0}')
        ax5 = pl2d(xrange=[min(zeta),max(zeta)], yrange=[min(cvdrift)-0.2*max(abs(cvdrift)),max(cvdrift)+0.2*max(abs(cvdrift))],xlabel='$\\zeta$', ylabel=labcdrift,    title='\\texttt{cdrift}')
        ax6 = pl2d(xrange=[min(zeta),max(zeta)], yrange=[min(cvdrift0)-0.2*max(abs(cvdrift0)),max(cvdrift0)+0.2*max(abs(cvdrift0))],xlabel='$\\zeta$', ylabel=labcdrift0,   title='\\texttt{cdrift0}')
        ax7 = pl2d(xrange=[min(zeta),max(zeta)], yrange=[min(gs2)-0.2*max(abs(gs2)),max(gs2)+0.2*max(abs(gs2))],xlabel='$\\zeta$', ylabel=labgs2,   title='\\texttt{gs2}')
        ax8 = pl2d(xrange=[min(zeta),max(zeta)], yrange=[min(gds21)-0.2*max(abs(gds21)),max(gds21)+0.2*max(abs(gds21))],xlabel='$\\zeta$', ylabel=labgds21,   title='\\texttt{gds21}')
        ax9 = pl2d(xrange=[min(zeta),max(zeta)], yrange=[min(gds22)-0.2*max(abs(gds22)),max(gds22)+0.2*max(abs(gds22))],xlabel='$\\zeta$', ylabel=labgds22,   title='\\texttt{gds22}')

        
        ax1.plot(zeta, bmag, color='r', linewidth=4)
        ax2.plot(zeta, gradpar, color='r', linewidth=4)
        ax3.plot(zeta, gbdrift, color='r', linewidth=4)
        ax4.plot(zeta, gbdrift0, color='r', linewidth=4)
        ax5.plot(zeta, cvdrift, color='r', linewidth=4)
        ax6.plot(zeta, cvdrift0, color='r', linewidth=4)    
        ax7.plot(zeta,gs2, color='r', linewidth=4)
        ax8.plot(zeta,gds21, color='r', linewidth=4)
        ax9.plot(zeta,gds22, color='r', linewidth=4)

################################################## Field line plot ######################################################



def field_line(file_list,nm):
        
        a=  file_list 
        case=[]
        colors=[]
        col=['r','g','b','y','m','y','k']       
        for i in arange(0,100):
            colors=colors+col
        ax = pl2d(xrange=[-pi/nm,pi/nm], yrange=[-pi,pi], xlabel='$\\zeta$', ylabel='$\\theta$', fig_size=(8.5, 7.5))            
        for i in arange(0,size(a)):
           element=separate(a[i])+'.vmec_geo' 
           case=case+[element]
        for l in arange(0,size(a)):
            dat=loadtxt(case[l],dtype='float',skiprows=4, usecols=1)
            data0=loadtxt(case[l],dtype='float',usecols=(1,2,3,4))
            q=data0[0,0]
            qu=str(q)
            i="{0:.4f}".format(1/q)
            iot=str(i)
            zed = factordivide(dat,10)
            zed_use_theta=loadtxt(case[l],dtype='float',skiprows=4, usecols=1)
            #for p in arange(0,size(dat)):                                
            #    zed_use_theta[p]=dat[p]
            for n in arange(1,100):
                    for i in arange(0,size(dat)):
                        if dat[i] <= pi/nm:
                            dat[i]=dat[i]
                        else:
                            dat[i]=dat[i]-2*pi/nm
                    for i in arange(0,size(dat)):
                        if dat[i] >= -pi/nm:
                            dat[i]=dat[i]
                        else:
                            dat[i]=dat[i]+2*pi/nm
            theta=factordivide(zed_use_theta,q)
            for n in arange(1,100):
                    for i in arange(0,size(theta)):
                        if theta[i] <= pi:
                            theta[i]=theta[i]
                        else:
                            theta[i]=theta[i]-2*pi
                    for i in arange(0,size(theta)):
                        if theta[i] >= -pi:
                            theta[i]=theta[i]
                        else:
                           theta[i]=theta[i]+2*pi
            #runs=a[l].split('/')
            #run=runs[size(runs)-2]+'/'+runs[size(runs)-1]
            #ax.plot(,theta,colors[l]+'o',markersize=3,label='$'+run+'$')
            ax.plot(dat,theta,marker='o',color=colors[l],markersize=3,label='case '+str(l+1)+ ' q= '+qu, linewidth=2)
            ax.legend(loc='best',labelspacing=0.0, prop={'size':24})
            

########################################################################################################################################################

def z_value(file_free_exten,quant,zed, PRINT = False ,PLOT = False):
        out=file_free_exten+'.vmec_geo'

        data=loadtxt(out,dtype='float',skiprows=4)
        data0=loadtxt(out,dtype='float',usecols=(0,1,2,3,4))

        q         = data0[0,1]
        shat      = data0[0,2]
        zeta      = data[:,1]
        bmag      = data[:,2]
        gradpar   = data[:,3]
        gbdrift   = data[:,9]
        gbdrift0  = data[:,10]
        cvdrift   = data[:,11]
        cvdrift0  = data[:,12]         
        gs2       = data[:,4]
        gds21     = data[:,5]
        gds22     = data[:,6]
        theta=factordivide(zeta,q)
        
        m=min(zeta)
        M=max(zeta)
        zeta_interp=linspace(m,M,1000000)

        ax = pl2d(xrange=[m,M], yrange=[min(theta),max(theta)], xlabel='$\\zeta(\cdot\pi)$', ylabel='$\\theta$', fig_size=(8.5, 7.5))
        ax.plot(zeta,theta)

        if quant=='bmag':
           f=interp1d(zeta,bmag,kind='cubic')
           ax.plot(zeta,bmag,color='r',linewidth=5,label='bmag')
           ax.plot(zeta_interp,f(zeta_interp),color='g',linewidth=1,label='bmag interpol')
        if quant=='gradpar':
           f=interp1d(zeta,gradpar,kind='cubic')
           ax.plot(zeta,gradpar,color='r',linewidth=5,label='gradpar')
           ax.plot(zeta_interp,f(zeta_interp),color='g',linewidth=1,label='gradpar interpol')
        if quant=='gbdrift':
           f=interp1d(zeta,gbdrift,kind='cubic')
           ax.plot(zeta,gbdrift,color='r',linewidth=5,label='gbdrift')
           ax.plot(zeta_interp,f(zeta_interp),color='g',linewidth=1,label='gbdrift interpol')
        if quant=='gbdrift0':
           f=interp1d(zeta,gbdrift0,kind='cubic')
           ax.plot(zeta,gbdrift0,color='r',linewidth=5,label='gbdrift0')
           ax.plot(zeta_interp,f(zeta_interp),color='g',linewidth=1,label='gbdrift0 interpol')
        if quant=='cvdrift':
           f=interp1d(zeta,cvdrift,kind='cubic')
           ax.plot(zeta,cvdrift,color='r',linewidth=5,label='cvdrift')
           ax.plot(zeta_interp,f(zeta_interp),color='g',linewidth=1,label='cvdrift interpol')
        if quant=='gvdrift0':
           f=interp1d(zeta,cvdrift0,kind='cubic')
           ax.plot(zeta,cvdrift0,color='r',linewidth=5,label='cvdrift0')
           ax.plot(zeta_interp,f(zeta_interp),color='g',linewidth=1,label='cvdrift0 interpol')
        if quant=='gds2':
           f=interp1d(zeta,gs2,kind='cubic')
           ax.plot(zeta,gs2,color='r',linewidth=5,label='gds2')
           ax.plot(zeta_interp,f(zeta_interp),color='g',linewidth=1,label='gds2 interpol')
        if quant=='gds21':
           f=interp1d(zeta,gds21,kind='cubic')
           ax.plot(zeta,gds21,color='r',linewidth=5,label='gds21')
           ax.plot(zeta_interp,f(zeta_interp),color='g',linewidth=1,label='gds21 interpol')
        if quant=='gds22':
           f=interp1d(zeta,gds22,kind='cubic')
           ax.plot(zeta,gds22,color='r',linewidth=5,label='gds22')
           ax.plot(zeta_interp,f(zeta_interp),color='g',linewidth=1,label='gds22 interpol')

        
        ax.legend(loc='best',labelspacing=0.0, prop={'size':24})
        
        theta_show=zed/q
        theta_s=str(theta_show)
        value_int=f(zed)
        value_int_s=str(value_int)

        if PRINT == True:
           print('theta for the $\\zeta$ value is ' + theta_s)
           print(quant + ' for the $\\zeta$ value is ' + value_int_s)

        if PLOT == True:                
           plt.show()

        
        return theta_show, value_int, q, shat

################################################################ LOAD QUANTITIES ##############################################

def load_case(file_free_exten):

        out=file_free_exten+'.vmec_geo'

        data=loadtxt(out,dtype='float',skiprows=4)
        data0=loadtxt(out,dtype='float',usecols=(0,1,2,3,4))

        q         = data0[0,1]
        shat      = data0[0,2]
        zeta      = data[:,1]
        bmag      = data[:,2]
        gradpar   = data[:,3]
        gbdrift   = data[:,9]
        gbdrift0  = data[:,10]
        cvdrift   = data[:,11]
        cvdrift0  = data[:,12]         
        gs2       = data[:,4]
        gds21     = data[:,5]
        gds22     = data[:,6]
        theta=factordivide(zeta,q)

        return q,shat,zeta,bmag,gradpar,gbdrift,gbdrift0,cvdrift,cvdrift0,gs2,gds21,gds22,theta



######################################################### ky read ############################################################# 

def ky(filen):
    # get ky grid
    out=separate(filen)+'.out.nc'
    ncfile    = netcdf.netcdf_file(out,'r')
    ky        = numpy.copy(ncfile.variables['ky'][:])
    naky      = ncfile.dimensions['ky']
    return ky, naky

######################################################### kx read ############################################################

def kx(filen):
    # get kx grid
    # this is the index of the first negative value of kx
    # note stella orders kx as (0, dkx, ..., kx_max, -kx_max, -kx_max+dkx, ..., -dkx)
    out=separate(filen)+'.out.nc'    
    ncfile    = netcdf.netcdf_file(out,'r')
    kx_stella = numpy.copy(ncfile.variables['kx'][:])
    nakx      = ncfile.dimensions['kx']
    nakx_mid  = nakx//2+1
    kx        = numpy.concatenate((kx_stella[nakx_mid:],kx_stella[:nakx_mid]))
    return kx, nakx, nakx_mid

######################################################### t read #############################################################

def time(filen):
    # get time grid
    out=separate(filen)+'.out.nc'
    ncfile    = netcdf.netcdf_file(out,'r')
    time      = numpy.copy(ncfile.variables['t'][:])
    ntime     = time.size
    return time, ntime

############################################### read omega if more than one k #################################################

def omega(case, last=False, view=False, yrange=None, xrange=None):
    # Note that:
    # time   ky   kx   Re[om]   Im[om]   Re[omavg]  Im[omavg]
    datafile = separate(case)+'.omega'
    dat=loadtxt(datafile,dtype='float')

    n_ky     = ky(case)[1]
    n_kx     = kx(case)[1]
    n_time   = int(size(dat[:,0])/n_ky) # Including t=0
    #n_time=time(case)[1]

    data     = loadtxt(datafile, dtype='float').reshape(n_time, n_kx, n_ky, 7)

    if view:
        pl2y(ky(case)[0], y1data=abs(data[n_time - 2, 0,:,3]), y2data=data[n_time - 2, 0,:,4],\
             xlabel='$k_y\\rho_{i}$', ylabel='$\\omega a/v_{\mathrm{th},\mathrm{i}}$',\
             key1='$a\\Re(\\omega)/v_{\mathrm{th},\mathrm{i}}$',\
             key2='$a\\Im(\\omega)/v_{\mathrm{th},\mathrm{i}}$', xrange=xrange,\
             yrange=yrange,fig_size=(8.5, 7.5), wp=1, ax=None, mkt="o", mkc='red',\
             title="$k_x="+str(kx(case)[0][0])+"$",\
             hline1=None, hline2=None, vshadow=None, ls1=1, ls2=2)

        if n_kx == 1 and n_ky > 1:
            # The evolution of omega(ky, t) can be represented
            #t      = time(case)[0][1:]
            kyrhoi = ky(case)[0]
            ky_re  = data[:, 0,:,3]
            ky_im  = data[:, 0,:,4]
            [X,Y]  = meshgrid(t, kyrhoi)
            Z      = ky_im.T

            surf(X, Y, Z, xlabel='$t$', ylabel='$k_y\\rho_{i}$',\
                 zlabel='$a\\Im(\\omega)/v_{\mathrm{th},\mathrm{i}}$',\
                 title="$k_x="+str(kx(case)[0][0])+"$", zrange=yrange, ax=None)

    plt.show()
            
    if last:
        return data[n_time-2, :, :, :]
    else:
        return data

########################################################### LnTLnN ############################################################

def ktkn(filen, process=False, PLOT=False):



    
    nLn = 9
    nLT = 9
    nky = 16
    run0 = 1
    irun=run0
    d=filen.split('/')
    direct=d[0]

    for p in arange(1,size(d)-2):
        direct=direct+'/'+d[p]
    outfile=direct+'/'+'ktkn.dat'
    
    omega_r     = empty((nLn, nLT, nky), dtype='float')
    omega_i     = empty((nLn, nLT, nky), dtype='float')
    omega_r_max = empty((nLn, nLT),     dtype='float')
    omega_i_max = empty((nLn, nLT),     dtype='float')    
    mLn         = empty((nLn, nLT),      dtype='float')
    mLT         = empty((nLn, nLT),      dtype='float')


    if process==True:

        f = open(outfile, 'w')
        f.write('# (0) run   (1) i_Ln   (2) i_LT   (3) nprim   (4) tprim   (5) Max(Re(omega))   (6) Max(Im(omega))' + '\n')
        print('\n'+ "Writing file :", outfile + '\n')
    
        for i in arange (0, nLn):  
  
            for j in arange (0, nLT): 
                
                case=direct+'/grads_00'+str(irun)+'/input.out.nc'
                ncfile = netcdf.netcdf_file(case,'r')
        
                nprim = numpy.copy(ncfile.variables['fprim'][0])
                tprim = numpy.copy(ncfile.variables['tprim'][0])
           

                mLn[i, j]         = nprim
                mLT[i, j]         = tprim
                
                omega_r[i, j, :]  = omega(filen, last=True)[0, :, 3]
                omega_i[i, j, :]  = omega(filen, last=True)[0, :, 4]
                omega_r_max[i, j] = max(omega(filen, last=True)[0, :, 3])
                omega_i_max[i, j] = max(omega(filen, last=True)[0, :, 4])

                #print(case, format8(i), format8(j), format2(nprim), format2(tprim),\
                #          format2(omega_r[i, j, :].max()), format2(omega_i[i, j, :].max()))
                
                f.write(case+'\t'+str(format8(i))+'\t'+str(format8(j))+'\t'+\
                        str(format2(nprim))+'\t'+str(format2(tprim))+'\t'+\
                        str(format2(max(omega_r[i, j, :])))+'\t'+\
                        str(format2(max(omega_i[i, j, :])))+'\n')
                          
                irun=irun+1  
                print(irun)
                #print(omega_r_max) 

        print('\n'+ "End of the writting"+'\n')

        f.close()
    
        

    if PLOT==True:

        print("Reading existing file : ", outfile)

        data = loadtxt(outfile, dtype='float', usecols=[3,4,5,6])
        mLn  = data[:,0].reshape((nLn, nLT))
        mLT  = data[:,1].reshape((nLn, nLT))
        omega_r_m = data[:,2].reshape((nLn, nLT))
        omega_i_m = data[:,3].reshape((nLn, nLT))      

        fig, ax = plt.subplots(figsize=(11.0, 8.4))
        col = ax.pcolormesh(mLn[:,0], mLT[0,:], omega_r_m.T, vmin=0, vmax=1.8)
        ax.set_xlabel('$a/L_{n_i}$')
        ax.set_ylabel('$a/L_{T_i}$')
        ax.xaxis.set_minor_locator(AutoMinorLocator(10))
        ax.yaxis.set_minor_locator(AutoMinorLocator(10))
        cbar=plt.colorbar(col)
        cbar.set_label('$\\Im({\\omega})^{\\textrm{max}}$')
        cbar.set_label('$a\\Im({\\omega})^{max}/v_{\\mathrm{th,i}}$')
        ax.set_aspect(1.0)
        title = '$k_{x}=0$; $0\\le k_{y}\le 10$'
        ax.set_title(title, y=1.02)
        ax.grid(color='grey', linestyle='-', linewidth=0.3)
                
                #        print mLn[:,0]
                #        print mLT[0,:]
                #        print omega_i_m
                #        mLT  = data[0,:].reshape(nLn, nLT)
                
                #    cmap(xdata=mLn[:,0], ydata=mLT[0,:], zdata=omega_r_m,\
                #         xlabel='$-a/L_{n_i}$', ylabel='$-a/L_{T_i}$', zlabel='$\\Im(\\omega)/v_{t}$')
                
        plt.show()


 ################################################### EUTERPE READ ################################################
    

def interpol(vector_x, vector_y, value_x, der=0):

    tck     = interpolate.splrep(vector_x, vector_y, s=0)
    value_y = interpolate.splev(value_x,tck,der=der)

    return value_y

def euprof(fileprof=None, svalue=None):
    # This function reads a profile in the format of EUTERPE
    # and converts it to another with the parameters
    # that stella needs.
    profdata=loadtxt(fileprof, dtype='float')
    if shape(profdata)[1] == 9:
        S,DTI,TI,DTE,TE,DNI,NI,DNE,NE = arange(0,9)

    dsdrho = 2*sqrt(profdata[:,S])

    s       =  profdata[:,S]
    nine    =  profdata[:,NI]/profdata[:,NE]
    tite    =  profdata[:,TI]/profdata[:,TE]
    tprim_i = -profdata[:,DTI]*dsdrho
    tprim_e = -profdata[:,DTE]*dsdrho
    fprim_i = -profdata[:,DNI]*dsdrho
    fprim_e = -profdata[:,DNE]*dsdrho
    t_i     =  profdata[:,TI]/1000.  # stella uses keV for temp
    t_e     =  profdata[:,TE]/1000.  # stella uses keV for temp
    dens_e  =  profdata[:,NE]/1.0E19 # stella uses keV for temp
    dens_i  =  profdata[:,NI]/1.0E19 # stella uses keV for temp

    if svalue != None:
        print("&vmec_parameters")
        print("torflux="+str(format9(svalue))+'\n')
        print("&parameters")
        print("nine="  +str(format9(interpol(profdata[:,S],nine,svalue))))
        print("tite="  +str(format9(interpol(profdata[:,S],tite,svalue)))+'\n')
        print("&species_parameters")
        print("dens="  +str(format9(interpol(profdata[:,S],dens_i,svalue))))
        print("temp="  +str(format9(interpol(profdata[:,S],t_i,   svalue))))
        print("tprim=" +str(format9(interpol(profdata[:,S],tprim_i[:],svalue))))
        print("fprim=" +str(format9(interpol(profdata[:,S],fprim_i[:],svalue)))+'\n')
        print("\n")
        print("&species_parameters_2")        
        print("dens="  +str(format9(interpol(profdata[:,S],dens_e,svalue))))
        print("temp="  +str(format9(interpol(profdata[:,S],t_e,   svalue))))
        print("tprim=" +str(format9(interpol(profdata[:,S],tprim_e,svalue))))
        print("fprim=" +str(format9(interpol(profdata[:,S],fprim_e,svalue))))

    return s, nine, tite, dens_i, t_i, tprim_i, fprim_i, dens_e, t_e, tprim_e, fprim_e 








