import numpy as np
from tkinter import *
from tkinter.filedialog import *
from matplotlib import *
#from turtle import *
from os import listdir
from scipy.io import netcdf
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import math
from math import*
from matplotlib.ticker import AutoMinorLocator
from matplotlib.pyplot import fill_between
from scipy import interpolate
from matplotlib import pylab, mlab, cm
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import LinearLocator, FixedLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.mplot3d import axes3d, Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from numpy import *
import stella_read
import stella_diag
import aux_functions
from aux_functions import *
from stella_read import *
from stella_diag import *





#################################### omega(t) plot #############################################################

def omega_t_plot(file_list,split_op=False):
    
    if split_op==False:
    
        n_ky=ky(file_list[0])[1]
        if n_ky==1:
 
            a=  file_list
            case=[]
            colors=color()  
            times=[]
            Im_max=[]
            Im_min=[]
            Re_max=[]
            Re_min=[]
            xlab  = "$t$"
            ylab  = "$a\\omega /v_{\mathrm{th},\mathrm{i}}$"
            tit1  = '$\\Re(\\omega)/v_{\mathrm{th},\mathrm{i}}$'
            tit2  = '$\\Im(\\omega)/v_{\mathrm{th},\mathrm{i}}$'    
 
         
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
                ky_value=str(dat[0,1])
                kx_value=str(dat[0,2])   
                key="$k_y \\rho$=" + ky_value + "   ," + "   $k_x \\rho=$" + kx_value
                
                ax2.plot(time, imaginary, linestyle='-', marker="s",\
                    markersize=5, color=colors[l], markerfacecolor="white", linewidth=4, label=key)                                  
                ax1.plot(time, real, linestyle='-', marker="o",\
                    markersize=5, color=colors[l], markerfacecolor="white", linewidth=4, label=key)
                ax1.legend(loc='best',labelspacing=0.0, prop={'size':24})      
                ax2.legend(loc='best',labelspacing=0.0, prop={'size':24}) 

        else:

            xlab  = "$t$"
            ylab  = '$|\\Im(\\omega)/v_{\mathrm{th},\mathrm{i}}|$'
            y2lab  = '$|\\Re(\\omega)/v_{\mathrm{th},\mathrm{i}}|$'

            colors=color()
            case=separate(file_list[0])
            n_ky = ky(file_list[0])[1]   
            caseom=case+'.omega'
            s=(len(open(caseom,'r').readlines()))/(n_ky+1)
            t=empty((int(s-1)))
            im=empty((int(s-1)))
            re=empty((int(s-1)))
            data=loadtxt(caseom,dtype='float',skiprows=1)
            ky_value = ky(file_list[0])[0]      
            t_test=empty((int(s-1)))
            im_test=empty((int(s-1)))
            re_test=empty((int(s-1)))
            for j in arange(0,s-2):
                i=int(j)
                a=int(j*n_ky)
                casel=data[a,:]
                t_test[i]=casel[0]
                im_test[i]=abs(casel[4])
                re_test[i]=abs(casel[3])        
        
            ax1 = pl2d(xrange=[min(t_test),max(t_test)],xlabel=xlab, ylabel=ylab, fig_size=(8.5, 7.5))
            ax2 = pl2d(xrange=[min(t_test),max(t_test)],yrange=[0,2],xlabel=xlab, ylabel=y2lab, fig_size=(8.5, 7.5))
        
            for p in arange(0,n_ky):
                for k in arange(0,s-2):
                    i=int(k)
                    a=int(k*n_ky+p)
                    casel=data[a,:]
                    t[i]=casel[0]
                    im[i]=abs(casel[4])
                    re[i]=abs(casel[3])
                ax1.plot(t, im, linestyle='-', marker="o",markersize=5, color=colors[p], markerfacecolor="white", linewidth=4, label='$k_y=$'+format3(ky_value[p]))
                ax2.plot(t, re, linestyle='-', marker="o",markersize=5, color=colors[p], markerfacecolor="white", linewidth=4, label='$k_y=$'+format3(ky_value[p]))
                ax1.legend(loc='best',labelspacing=0.0, prop={'size':20})
                ax2.legend(loc='best',labelspacing=0.0, prop={'size':20})

    if split_op==True:
        n_ky=ky(file_list[0])[1]
        if n_ky==1:
 
            a=  file_list
            case=[]
            colors=color()  
            times=[]
            Im_max=[]
            Im_min=[]
            Re_max=[]
            Re_min=[]
            xlab  = "$t$"
            y1lab  = '$\\Re(\\omega)/v_{\mathrm{th},\mathrm{i}}$' 
            y2lab  = '$\\Im(\\omega)/v_{\mathrm{th},\mathrm{i}}$'    
 
         
            for i in arange(0,size(a)):
               element=separate(a[i])+'.omega' 
               case=case+[element] 
   
            for l in arange(0,size(a)):
                dat=loadtxt(case[l],dtype='float')
                imaginary=dat[:,4]
                real=dat[:,3]
                time=dat[:,0]
                Re_min=min(real)
                Re_max=max(real)
                Im_min=min(imaginary)
                Im_max=max(imaginary)
                ky_value=str(dat[0,1])
                kx_value=str(dat[0,2])   
                tit="$k_y \\rho$=" + ky_value + "   ," + "   $k_x \\rho=$" + kx_value  
                
                ax1 = pl2d(xrange=[min(time),max(time)],yrange=[Re_min,Re_max],xlabel=xlab, ylabel=y1lab, fig_size=(8.5, 7.5), title=tit)
                ax2 = pl2d(xrange=[min(time),max(time)],yrange=[Im_min,Im_max],xlabel=xlab, ylabel=y2lab, fig_size=(8.5, 7.5), title=tit)   
                ax2.plot(time, imaginary, linestyle='-', marker="s",\
                    markersize=5, color='b', markerfacecolor="white", linewidth=4)            
                ax1.plot(time, real, linestyle='-', marker="o",\
                    markersize=5, color='r', markerfacecolor="white", linewidth=4)


        else:

            xlab  = "$t$"
            y1lab  = '$|\\Im(\\omega)/v_{\mathrm{th},\mathrm{i}}|$'
            y2lab  = '$|\\Re(\\omega)/v_{\mathrm{th},\mathrm{i}}|$'

            colors=color()
            case=separate(file_list[0])
            n_ky = ky(file_list[0])[1]   
            caseom=case+'.omega'
            s=(len(open(caseom,'r').readlines()))/(n_ky+1)
            t=empty((int(s-1)))
            im=empty((int(s-1)))
            re=empty((int(s-1)))
            data=loadtxt(caseom,dtype='float',skiprows=1)
            ky_value = ky(file_list[0])[0]            
        
            for p in arange(0,n_ky):
                for k in arange(0,s-2):
                    i=int(k)
                    a=int(k*n_ky+p)
                    casel=data[a,:]
                    t[i]=casel[0]
                    im[i]=abs(casel[4])
                    re[i]=abs(casel[3])
                ax1 = pl2d(xrange=[min(t),max(t)],xlabel=xlab, ylabel=y1lab, fig_size=(8.5, 7.5), title='$k_y=$'+format3(ky_value[p]))
                ax2 = pl2d(xrange=[min(t),max(t)],yrange=[0,2],xlabel=xlab, ylabel=y2lab, fig_size=(8.5, 7.5), title='$k_y=$'+format3(ky_value[p]))
                ax1.plot(t, im, linestyle='-', marker="o",markersize=5, color='b', markerfacecolor="white", linewidth=4)
                ax2.plot(t, re, linestyle='-', marker="o",markersize=5, color='r', markerfacecolor="white", linewidth=4)
    



############################################### omega(k) plot #########################################################

def omega_k_plot(file_list,n=False):
       if size(file_list)==1:
            n_ky = ky(file_list[0])[1]
            n_kx = kx(file_list[0])[1]

            if n_kx==1 and n_ky==1:
                case=separate(file_list[0])+'.omega'
                data=loadtxt(case,dtype='float')
                ky_value=ky(file_list[0])[0]
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

                xlab="$k_y\\rho_i$"
                ylab="$a\\omega /v_{\mathrm{th},\mathrm{i}}$"
                key1= '$\\Re(\\omega)/v_{\mathrm{th},\mathrm{i}}$'
                key2= '$\\Im(\\omega)/v_{\mathrm{th},\mathrm{i}}$'

                ax = pl2d(xrange=[ky_value-0.5,ky_value+0.5], yrange=[min(limits)-0.2*max(alim),max(limits)+0.2*max(alim)], xlabel=xlab, ylabel=ylab, fig_size=(8.5, 7.5))
                ax.plot(ky_value, realpoint, marker="o", markersize=12, color="black", markerfacecolor="black", label=key1)
                ax.plot(ky_value, impoint, marker="o", markersize=12, color="blue", markerfacecolor="blue", label=key2)
                ax.legend(loc='best',labelspacing=0.0, prop={'size':24})


            if n_kx==1 and n_ky>1:
                max_abs_pot=empty(n_ky)
                casel=empty((n_ky,4))               
                if n==False:
                    n=0
                case=separate(file_list[0])+'.omega'
                s=len(open(case,'r').readlines())
                last=(s-n_ky-1)-n*(n_ky+1)
                data=loadtxt(case,dtype='float',skiprows=last)
                realin=abs(data[:,3])
                imaginaryin=abs(data[:,4])
                tin=data[:,0]
                kx_value=data[0,2]
                kyin=data[:,1]
                for m in arange(0,size(tin)):
                    if tin[m]==min(tin):
                        casel[m,0]=tin[m]
                        casel[m,1]=kyin[m]
                        casel[m,3]=imaginaryin[m]
                        casel[m,2]=realin[m]
                real=casel[:,2]
                imaginary=casel[:,3]
                t=casel[:,0]
                ky_value=casel[:,1]
                imaginary_value=zeros(size(imaginary))
                real_value=zeros(size(real))

                for l in arange(0,n_ky): #Bucle para calcular el potencial.

                        phi = phi_vs_t(file_list[0])
                        nt, nzed, nkx, nky, nphi = shape(phi)
    
                        phi_last = phi[nt-1,:,0,l,:]
                        phi_re   = phi_last[:,0]
                        phi_im   = phi_last[:,1]
                        abs_pot=max(abs(append(phi_im,phi_re)))
                        max_abs_pot[l]=abs_pot


                for n in arange(0,n_ky):  #La idea esq si el potencial es pequeño el growth rate sea NaN
                    if max_abs_pot[n]>100:
                        imaginary_value[n]=imaginary[n]
                        real_value[n]=real[n]
                    else:
                        imaginary_value[n]=NaN
                        real_value[n]=NaN

                xlab  = "$k_{y}\\rho_i$"
                y1lab = "$a|\\Re(\\omega)|/v_{\mathrm{th},\mathrm{i}}$"
                y2lab = "$a|\\Im(\\omega)|/v_{\mathrm{th},\mathrm{i}}$"
                tit   = "$\\omega(k_{y}\\rho_i)$ for $k_x\\rho_i$="+format3(kx_value) 
                
                ax1 = pl2d(xrange=[min(ky_value)-0.1,max(ky_value)+0.1], yrange=[0,max(real)+0.2*max(real)], xlabel=xlab, ylabel=y1lab,\
                         fig_size=(8.5, 7.5), title=tit)
                ax1.plot(ky_value, real_value, linestyle='-', marker='o', markersize=5, color="purple", markerfacecolor="white", linewidth=4)

                ax2 = pl2d(xrange=[min(ky_value)-0.1,max(ky_value)+0.1], yrange=[0,max(imaginary)+0.2*max(imaginary)], xlabel=xlab, ylabel=y2lab,\
                        fig_size=(8.5, 7.5), title=tit)
                ax2.plot(ky_value, imaginary_value, linestyle='-', marker='o', markersize=5, color="green", markerfacecolor="white", linewidth=4) 

       else:          
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
            ky_value=om[:,1]
            kx_value=str(om[0,2])  
            xlab  = "$k_{y}\\rho_i$"
            y1lab = "$a|\\Re(\\omega)|/v_{\mathrm{th},\mathrm{i}}$"
            y2lab = "$a|\\Im(\\omega)|/v_{\mathrm{th},\mathrm{i}}$"
            tit   = "$\\omega(k_{y}\\rho_i)$ "+"for "+"$k_{x}\\rho_{i}=$"+kx_value          
            ax1 = pl2d(xrange=[min(ky_value)-0.1,max(ky_value)+0.1], yrange=[min(real)-0.1*max(abs(real)),max(real)+0.1*max(abs(real))], xlabel=xlab, ylabel=y1lab,\
                        fig_size=(8.5, 7.5), title=tit)
            ax1.plot(ky_value, real, linestyle='-', marker='o', markersize=5, color="purple", markerfacecolor="white", linewidth=4)

            ax2 = pl2d(xrange=[min(ky_value)-0.1,max(ky_value)+0.1], yrange=[min(imaginary)-0.1*max(abs(imaginary)),max(imaginary)+0.1*max(abs(imaginary))], xlabel=xlab, ylabel=y2lab,\
                        fig_size=(8.5, 7.5), title=tit)
            ax2.plot(ky_value, imaginary, linestyle='-', marker='o', markersize=5, color="green", markerfacecolor="white", linewidth=4)     


################################################### Potential plot  #####################################################

def potential_plot(file_list):

    nky=ky(file_list[0])[1]

    if nky==1:
        if size(file_list)==1:

            case=separate(file_list[0])+'.final_fields'
            data=loadtxt(case,dtype='float')
            ky_value=format3(ky(file_list[0])[0])
            kx_value=format3(kx(file_list[0])[0])       
            imaginary=data[:,5]
            real=data[:,4]
            zeta = zed(file_list[0])[0]
            limits=append(real,imaginary)
            alim=abs(limits)
 
            tit   = "$k_y \\rho$=" + ky_value + "   ," + "$   k_x \\rho=$" + kx_value
            xlab  = '$\\zeta$'
            ylab  = '$\\tilde{\\varphi}(\\zeta)$'
  
            key1='$\\Re(\\tilde{\\varphi})(\\zeta)$'
            key2='$\\Im(\\tilde{\\varphi})(\\zeta)$'
    
            ax = pl2d(xrange=[min(zeta),max(zeta)], yrange=[min(limits)-0.2*max(alim),max(limits)+0.2*max(alim)], xlabel=xlab, ylabel=ylab,\
                fig_size=(8.5, 7.5), title=tit)
            ax.plot(zeta, real, linestyle='-', marker="o",\
                markersize=5, color="black", markerfacecolor="white", linewidth=4, label=key1)
            ax.plot(zeta, imaginary, linestyle='-', marker="s",\
                markersize=5, color="blue", markerfacecolor="white", linewidth=4, label=key2)
            ax.legend(loc='best',labelspacing=0.0, prop={'size':24})
        else:
            for i in arange(0,size(file_list)):
                case=separate(file_list[i])+'.final_fields'
                data=loadtxt(case,dtype='float')
                ky_value=format3(ky(file_list[i])[0])
                kx_value=format3(kx(file_list[i])[0])       
                imaginary=data[:,5]
                real=data[:,4]
                zeta = zed(file_list[i])[0]
                limits=append(real,imaginary)
                alim=abs(limits)
                tit   = "$k_y \\rho$=" + ky_value + "   ," + "$   k_x \\rho=$" + kx_value
                xlab  = '$\\zeta$'
                ylab  = '$\\tilde{\\varphi}(\\zeta)$'
                key1='$\\Re(\\tilde{\\varphi})(\\zeta)$'
                key2='$\\Im(\\tilde{\\varphi})(\\zeta)$'
                ax = pl2d(xrange=[min(zeta),max(zeta)], yrange=[min(limits)-0.2*max(alim),max(limits)+0.2*max(alim)], xlabel=xlab, ylabel=ylab,\
                    fig_size=(8.5, 7.5), title=tit)
                ax.plot(zeta, real, linestyle='-', marker="o",\
                    markersize=5, color="black", markerfacecolor="white", linewidth=4, label=key1)
                ax.plot(zeta, imaginary, linestyle='-', marker="s",\
                    markersize=5, color="blue", markerfacecolor="white", linewidth=4, label=key2)
            ax.legend(loc='best',labelspacing=0.0, prop={'size':24})

                
    else:
            key1='$\\Re(\\tilde{\\varphi})(\\zeta)$'
            key2='$\\Im(\\tilde{\\varphi})(\\zeta)$'
            xlab  = '$\\zeta$'
            ylab  = '$\\tilde{\\varphi}(\\zeta)$'
            vzed, nzed, iz0          = zed(file_list[0])[0], zed(file_list[0])[1], zed(file_list[0])[2]
            xlims = [vzed.min(), vzed.max()]            
            for m in arange(0,nky):

                phi                      = phi_vs_t(file_list[0])
                nt, nzed, nkx, nky, nphi = shape(phi)
    
                phi_last = phi[nt-1,:,0,m,:]
                phi_re   = phi_last[:,0]
                phi_im   = phi_last[:,1]
                ky_value=format3(ky(file_list[0])[0][m])
                #kx_value=format3(kx(a[0])[0][m])
                limits=append(phi_im,phi_re)
                alim=abs(limits)
                tit='$k_y\\rho_i=$  '+ky_value
                ax=pl2d(xrange=[min(vzed),max(vzed)], yrange=[min(limits)-0.2*max(alim),max(limits)+0.2*max(alim)], xlabel=xlab, ylabel=ylab, fig_size=(8.5, 7.5), title=tit)
                ax.plot(vzed,phi_re,linestyle='-', marker="o", markersize=5, color="black", markerfacecolor="white", linewidth=4, label=key1)
                ax.plot(vzed,phi_im,linestyle='-', marker="s", markersize=5, color="blue", markerfacecolor="white", linewidth=4, label=key2)
                ax.legend(loc='best',labelspacing=0.0, prop={'size':18})


#################################################################################################
def potential_norm_plot(a,split_op=False):

    if split_op==False:

        case=[]
        case_n=[]
        colors=color() 
        max_pot_norm=[]
        xlab  = '$\\zeta$'
        ylab  = '$\\tilde{\\varphi}(\\zeta)/\\mathrm{max}(\\tilde{\\varphi}(\\zeta))$'
        tit = 'Normalized Potential'  
        nky=ky(a[0])[1]
        
        if nky==1:
 
            if size(a)>1:
                for i in arange(0,size(a)):
                    element=separate(a[i])+'.final_fields' 
                    case=case+[element] 
                    element_n=separate(a[i])+'.out.nc'
                    case_n=case_n+[element_n]

                for n in arange(0,size(a)):
                    ncfile = netcdf.netcdf_file(case_n[n],'r')
                    zeta = np.copy(ncfile.variables['zed'][:])
                    data=loadtxt(case[n],dtype='float')
                    im=data[:,5]
                    re=data[:,4]
                    potential_graph=normali(re,im)
                    p0_g=max(potential_graph)
                    potential_norm_graph=factordivide(potential_graph,p0_g)

                ax = pl2d(xrange=[min(zeta),max(zeta)], yrange=[0,1], xlabel=xlab, ylabel=ylab,fig_size=(8.5, 7.5),title=tit)

                for l in arange(0,size(a)):
                    dat=loadtxt(case[l],dtype='float')
                    ky_value=ky(a[l])[0]
                    kx_value=kx(a[l])[0]
                    imaginary=dat[:,5]
                    real=dat[:,4]
                    ncfile = netcdf.netcdf_file(case_n[n],'r')
                    zet = np.copy(ncfile.variables['zed'][:]) 
                    potential=normali(real,imaginary)
                    p0=max(potential)
                    potential_norm=factordivide(potential,p0)
    
                    key   = "$k_y \\rho$=" + format3(ky_value) + "   ," + "$   k_x \\rho=$" + format3(kx_value)
                    ax.plot(zet, potential_norm, linestyle='-', color=colors[l], linewidth=3, label=key)
                    ax.legend(loc='best',labelspacing=0.0, prop={'size':24})
            else:
                element=separate(a[0])+'.final_fields' 
                case=case+[element] 
                element_n=separate(a[0])+'.out.nc'
                case_n=case_n+[element_n]
                ncfile = netcdf.netcdf_file(str(case_n[0]),'r')
                zeta = np.copy(ncfile.variables['zed'][:])
                data=loadtxt(str(case[0]),dtype='float')
                ky_value=ky(a[0])[0]
                kx_value=kx(a[0])[0]
                im=data[:,5]
                re=data[:,4]
                potential_graph=normali(re,im)
                p0_g=max(potential_graph)
                potential_norm_graph=factordivide(potential_graph,p0_g)
                ax = pl2d(xrange=[min(zeta),max(zeta)], yrange=[0,1], xlabel=xlab, ylabel=ylab,fig_size=(8.5, 7.5),title=tit)
                key   = "$k_y \\rho$=" + format3(ky_value) + "   ," + "$   k_x \\rho=$" + format3(kx_value)
                ax.plot(zeta, potential_norm_graph, linestyle='-', color='b', linewidth=3,label=key)
                ax.legend(loc='best',labelspacing=0.0, prop={'size':24})
        else:
            vzed, nzed, iz0          = zed(a[0])[0], zed(a[0])[1], zed(a[0])[2]
            xlims = [vzed.min(), vzed.max()]            
            ax1=pl2d(xrange=xlims,xlabel='$\\zeta$',yrange=[0,1],ylabel='$|\\tilde{\\varphi}(\\zeta)|/\\mathrm{max}(|\\tilde{\\varphi}(\\zeta)|)$',fig_size=(8.5,7.5),title='Normalized potential')
            for m in arange(0,nky):
                phi                      = phi_vs_t(a[0])
                nt, nzed, nkx, nky, nphi = shape(phi)
    
                phi_last = phi[nt-1,:,0,m,:]
                phi_re   = phi_last[:,0]
                phi_im   = phi_last[:,1]
                phi_0=max(max(abs(phi_re)),max(abs(phi_im)))
                phi_re_plot=factordivide(phi_re,phi_0)
                phi_im_plot=factordivide(phi_im,phi_0)
                ky_value=format3(ky(a[0])[0][m])
                npr=nprim(a[0])[0]
                tpr=tprim(a[0])[0]
                normal=normali(phi_re,phi_im)
                norm_0=max(normal)
                normal_plot=factordivide(normal,norm_0)
                key='$k_i\\rho_i=$  '+ky_value
                ax1.plot(vzed,normal_plot,label=key)
                ax1.legend(loc='best',labelspacing=0.0, prop={'size':18})

    if split_op==True:

        case=[]
        case_n=[]
        colors=color() 
        max_pot_norm=[]
        xlab  = '$\\zeta$'
        ylab  = '$\\tilde{\\varphi}(\\zeta)/\\mathrm{max}(\\tilde{\\varphi}(\\zeta))$'  
        nky=ky(a[0])[1]
        
        if nky==1:
 
            if size(a)>1:
                for i in arange(0,size(a)):
                    element=separate(a[i])+'.final_fields' 
                    case=case+[element] 
                    element_n=separate(a[i])+'.out.nc'
                    case_n=case_n+[element_n]

                for n in arange(0,size(a)):
                    ncfile = netcdf.netcdf_file(case_n[n],'r')
                    zeta = np.copy(ncfile.variables['zed'][:])
                    data=loadtxt(case[n],dtype='float')
                    imaginary=data[:,5]
                    real=data[:,4]
                    ky_value=ky(a[n])[0]
                    kx_value=kx(a[n])[0] 
                    potential=normali(real,imaginary)
                    p0=max(potential)
                    potential_norm=factordivide(potential,p0)
                    
                    title   = "$k_y \\rho$=" + format3(ky_value) + "   ," + "$   k_x \\rho=$" + format3(kx_value)
                    ax = pl2d(xrange=[min(zeta),max(zeta)], yrange=[0,1], xlabel=xlab, ylabel=ylab,fig_size=(8.5, 7.5))    
                    ax.plot(zeta, potential_norm, linestyle='-', color='b', linewidth=3)
            else:
                element=separate(a[0])+'.final_fields' 
                case=case+[element] 
                element_n=separate(a[0])+'.out.nc'
                case_n=case_n+[element_n]
                ncfile = netcdf.netcdf_file(str(case_n[0]),'r')
                zeta = np.copy(ncfile.variables['zed'][:])
                data=loadtxt(str(case[0]),dtype='float')
                ky_value=ky(a[0])[0]
                kx_value=kx(a[0])[0]
                im=data[:,5]
                re=data[:,4]
                potential_graph=normali(re,im)
                p0_g=max(potential_graph)
                potential_norm_graph=factordivide(potential_graph,p0_g)

                tit   = "$k_y \\rho$=" + format3(ky_value) + "   ," + "$   k_x \\rho=$" + format3(kx_value)
                ax = pl2d(xrange=[min(zeta),max(zeta)], yrange=[0,1], xlabel=xlab, ylabel=ylab,fig_size=(8.5, 7.5),title=tit)
                ax.plot(zeta, potential_norm_graph, linestyle='-', color='b', linewidth=3)
        else:
            vzed, nzed, iz0          = zed(a[0])[0], zed(a[0])[1], zed(a[0])[2]
            xlims = [vzed.min(), vzed.max()]            
            
            for m in arange(0,nky):
                phi                      = phi_vs_t(a[0])
                nt, nzed, nkx, nky, nphi = shape(phi)
    
                phi_last = phi[nt-1,:,0,m,:]
                phi_re   = phi_last[:,0]
                phi_im   = phi_last[:,1]
                phi_0=max(max(abs(phi_re)),max(abs(phi_im)))
                phi_re_plot=factordivide(phi_re,phi_0)
                phi_im_plot=factordivide(phi_im,phi_0)
                ky_value=format3(ky(a[0])[0][m])
                npr=nprim(a[0])[0]
                tpr=tprim(a[0])[0]
                normal=normali(phi_re,phi_im)
                norm_0=max(normal)
                normal_plot=factordivide(normal,norm_0)
                tit='$k_i\\rho_i=$  '+ky_value
                ax1=pl2d(xrange=xlims,xlabel='$\\zeta$',yrange=[0,1],ylabel='$|\\tilde{\\varphi}(\\zeta)|/\\mathrm{max}(|\\tilde{\\varphi}(\\zeta)|)$',fig_size=(8.5,7.5),title=tit)
                ax1.plot(vzed,normal_plot,color='b')

##########################################################################################################################





################################################# geo together plot ######################################################

def geo_tg(file_list):  
        case=separate(file_list[0])
        out=case+'.vmec_geo'
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

def geo_sp(file_list):

        case=separate(file_list[0])
        out=case+'.vmec_geo'
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
        colors=color()     

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

    if n_kx == 1 and n_ky > 1:
            # The evolution of omega(ky, t) can be represented
            t      = time(case)[0][1:]
            kyrhoi = ky(case)[0]
            ky_re  = data[:, 0,:,3]
            ky_im  = data[:, 0,:,4]
            [X,Y]  = meshgrid(t, kyrhoi)
            Z      = ky_im.T

    plt.show()
            
    if last:
        return data[n_time-2, :, :, :]
    else:
        return data


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


############################################growth rate (t) more than 1 ky#####################################
def gro_time(cases,n):
    
    xlab  = "$t$"
    ylab  = '$|\\Im(\\omega)/v_{\mathrm{th},\mathrm{i}}|$'
    y2lab  = '$|\\Re(\\omega)/v_{\mathrm{th},\mathrm{i}}|$'

    colors=color()

    case=separate(cases)
    n_ky = ky(cases)[1]   
    caseom=case+'.omega'
    s=(len(open(caseom,'r').readlines()))/(n_ky+1)
    t=empty((int(s-1)))
    im=empty((int(s-1)))
    re=empty((int(s-1)))
    data=loadtxt(caseom,dtype='float',skiprows=1)
    
    if size(n)==1:
        for j in arange(0,s-2):
            i=int(j)
            a=int(j*n_ky+n)
            casel=data[a,:]
            t[i]=casel[0]
            im[i]=abs(casel[4])
            re[i]=abs(casel[3])
        ax1 = pl2d(xrange=[min(t),max(t)],yrange=[min(im)-0.1*max(im),max(im)+0.1*max(im)],xlabel=xlab, ylabel=ylab, fig_size=(8.5, 7.5))
        ax1.plot(t, im, linestyle='-', marker="o",markersize=5, color=colors[j], markerfacecolor="white", linewidth=4)
        ax2 = pl2d(xrange=[min(t),max(t)],yrange=[min(re)-0.1*max(re),max(re)+0.1*max(re)],xlabel=xlab, ylabel=y2lab, fig_size=(8.5, 7.5))
        ax2.plot(t, re, linestyle='-', marker="o",markersize=5, color=colors[j], markerfacecolor="white", linewidth=4)

        
    if size(n)>1:  
        ky_n = ky(cases)[0]      
        t_test=empty((int(s-1)))
        im_test=empty((int(s-1)))
        re_test=empty((int(s-1)))
        for j in arange(0,s-2):
            i=int(j)
            a=int(j*n_ky+n[0])
            casel=data[a,:]
            t_test[i]=casel[0]
            im_test[i]=abs(casel[4])
            re_test[i]=abs(casel[3])        
        
        ax1 = pl2d(xrange=[min(t_test),max(t_test)],xlabel=xlab, ylabel=ylab, fig_size=(8.5, 7.5))
        ax2 = pl2d(xrange=[min(t_test),max(t_test)],yrange=[0,2],xlabel=xlab, ylabel=y2lab, fig_size=(8.5, 7.5))
    
        for p in arange(0,size(n)):
            for k in arange(0,s-2):
                i=int(k)
                a=int(k*n_ky+n[p])
                casel=data[a,:]
                t[i]=casel[0]
                im[i]=abs(casel[4])
                re[i]=abs(casel[3])
            pos=n[p]
            ax1.plot(t, im, linestyle='-', marker="o",markersize=5, color=colors[p], markerfacecolor="white", linewidth=4, label='$k_y=$'+format3(ky_n[pos]))
            ax1.legend(loc='best',labelspacing=0.0, prop={'size':20})
            ax2.plot(t, re, linestyle='-', marker="o",markersize=5, color=colors[p], markerfacecolor="white", linewidth=4, label='$k_y=$'+format3(ky_n[pos]))
            ax2.legend(loc='best',labelspacing=0.0, prop={'size':20})



####################################COLISIONS#####################################
def coulog(species, mass, charge, n, T):
#def coulog(species=[1,1], mass=[1,1], charge=[1,1], n=[1e19, 1.0e19], T=[1000., 1000.]):
    #
    # species==1 mean ions. species==2 mean electrons
    # Coulomb log for a pair of colliding species with densities
    # n=[na,nb] (in 1x10^19) and T=[Ta, Tb] (in eV). Function dens and temp provide that data.
    # ee collisions (NRL plasma formulary 2016)
    # conversion factor 1.0D6 needed to transform density to cgs
    
    spa, spb = species
    ma , mb  = mass
    za , zb  = charge
    na , nb  = n
    ta , tb  = T
 
    
        
    if spa == 2 and spb == 2:
        # e-e collisions
        coulog = 23.5 - log(sqrt(na/1.0E6)*(ta**(-5.0/4.0)))-\
                 sqrt( 1.0E-5 + ((log(ta)-2.0)**2) / 16.0)
    elif spa == 2 and spb != 2:
        # e-i and e-z collisions
        coulog = 24.0 - log(sqrt(na/1.0E6)/ta)
    elif spa != 2 and spb == 2:
        # i-e and z-e collisions
        coulog = 24.0 - log(sqrt(nb/1.0E6)/tb)
    elif spa != 2 and spb !=2:
        # i-i, i-z, z-i, z-z
#        coulog = 23.0 - LOG( (za*zb*(mua+mub)/(mua*tb+mub*ta))* &
#                             SQRT(na/1.0D6*(za**2)/ta + nb/1.0D6*(zb**2)/tb) )
        coulog = 23.0 - log( (za*zb*(ma+mb)/(ma*tb+mb*ta))*\
                             sqrt(na/1.0E6*(za**2)/ta + nb/1.0E6*(zb**2)/tb) )

    return coulog
def nuab_fac(species, mass, charge, n, T):
#def nuab_fac(species=[1,1], mass=[1,1], charge=[1,1], n=[1e19, 1.0e19], T=[1000, 1000]):
    # Calculates constant factor of the collision frequencies
    # for any two species a and b, according to
    # nu_ab = nb * (q1*q2)**2 * coulog / (4 * pi * eps0**2 *(m1**2)*(vtha**3))
    #
    # This function requires temperature (ta, tb) in eV, 
    # densities (na, nb) in m^-3 and the charge below (za, zb) in units of "e".

    e=1.6022e-19
    cloc =299792458
    eps_0=1e7/(4*pi*cloc*cloc)
    m_p=1.672621637e-27

    ma , mb  = mass
    za , zb  = charge
    vtha     = sqrt( 2.0 * T[0] * e / (ma*m_p))
    #vtha     = sqrt(T[0] * e / (ma*m_p))
    fac_0    = (e**2 / (eps_0 * m_p))**2 / ( 4 * pi )
    fac_nuab = fac_0 * n[1] * za**2 * zb**2 * coulog(species, mass, charge, n, T) / (ma**2 * vtha**3)
    
    return fac_nuab

def colfactor_d(xa=1.0, xb=1.0, plot=0):
    # Calculates (phi(xb)-G(xb))/xa**3
    # xa := v/vtha. ("b" is the target species)
    # xb := v/vthb. ("a" is the colliding species)
    # v  :  velocity of the colliding particle.
    # erf(xb) is phi and  erf(xb)/(2*xb**2) ) - 1.0 / (sqrt(pi) * xb) * exp(-xb**2) is G

    G=erf(xb)/(2*xb**2) - 1.0 / (sqrt(pi) * xb) * exp(-xb**2)
    factor = erf(xb) -G

    return factor / xa**3.

def colfactor_s(mass, T):
#def colfactor_s(xa=1.0, xb=1.0, mass=[1,1], T=[1000,1000]):
    xa=1
    xb=1
    ma ,mb= mass
    Ta, Tb =T
    G=erf(xb)/(2*xb**2) - 1.0 / (sqrt(pi) * xb) * exp(-xb**2)
    fac=2*Ta/Tb*(1+mb/ma)
    factor=G*fac/xa
    
    return factor

def colfactor_pararell(xa=1.0, xb=1.0):
    G=erf(xb)/(2*xb**2) - 1.0 / (sqrt(pi) * xb) * exp(-xb**2)
    factor=2*G/(xa**3)
    
    return factor

def nu_d_ab(species, mass, charge, n, T):
#def nu_d_ab(species=[1,1], mass=[1,1], charge=[1,1], n=[1e19, 1.0e19], T=[1000, 1000]):
    
    fac_nuab=nuab_fac(species, mass, charge, n, T)
    colfac_d=colfactor_d()
    nu_d_ab=fac_nuab*colfac_d

    return nu_d_ab


def nu_s_ab(species, mass, charge, n, T):
#def nu_s_ab(species=[1,1], mass=[1,1], charge=[1,1], n=[1e19, 1.0e19], T=[1000, 1000]):
    
    fac_nuab=nuab_fac(species, mass, charge, n, T)
    colfac_s=colfactor_s(mass,T)
    nu_s_ab=fac_nuab*colfac_s

    return nu_s_ab


def nu_par_ab(species, mass, charge, n, T):
#def nu_par_ab(species=[1,1], mass=[1,1], charge=[1,1], n=[1e19, 1.0e19], T=[1000, 1000]):
    
    fac_nuab=nuab_fac(species, mass, charge, n, T)
    colfac_par=colfactor_pararell()
    nu_par_ab=fac_nuab*colfac_par

    return nu_par_ab

########################################### Scalar Quant read ##################################
def vmec_quant(case):
        case_free=separate(case)
        out=case_free+'.vmec_geo'

        data=loadtxt(out,dtype='float',usecols=(0,1,2,3,4))

        rho= data[0,0]
        q=data[0,1]
        shat= data[0,2]
        aref=data[0,3]
        Bref=data[0,4]
        iota=1/q


        return q,iota,rho,shat,aref,Bref 

def phys_quant(mass,T,dens,charge,B):
    e=1.6022e-19
    cloc =299792458
    eps_0=1e7/(4*pi*cloc*cloc)
    m_p=1.672621637e-27
    vth=empty((size(T)))
    lar=empty((size(T)))
    freq=empty((size(T)))
    deb=empty((size(T)))
    mfp_def=empty((size(T)))
    for i in arange(0,size(T)):    
        vth[i]  = sqrt( 2.0 * T[i] * e / (mass[i]*m_p))
    for n in arange(0,size(T)):
        lar[n]=(mass[n]*m_p*vth[n])/(abs(charge[n]*e)*B)
        freq[n]=(e*B)/(mass[n]*m_p)
        deb[n]=sqrt(eps_0/e)*sqrt(T[n]/dens[n])

    return lar, freq, deb



##########################################

    
    
    





