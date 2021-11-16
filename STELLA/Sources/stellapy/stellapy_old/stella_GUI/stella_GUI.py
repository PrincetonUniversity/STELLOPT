
import os
import matplotlib.pyplot as plt
import stella_plots
import post_processing
import functions


from os import listdir
from numpy import *
from tkinter import *
from tkinter.filedialog import *
from matplotlib import *
from scipy.io import netcdf
from array import array
from functions import *
from stella_plots import *
from post_processing import *





                                          #FUNCTION DEFINITIONS#






#########################################  FUNCTION TO PLOT OPTIONS ######################################################
def options():
    ofile=dire.get()
    mult=che.get()
    multiple=mult.split(",,")
    directory=separate(ofile)
    oa=a.get()
    ob=b.get()
    oc=c.get()
    od=d.get()
    op=p.get()
    oo=o.get()
    ogeot=geo1.get()
    ogeos=geo2.get()

    
    if oa==1: #omega(t) plot
        omega_t_plot(multiple)

        


    if ob==1: #omega(k) plot       
        if size(multiple)==1:  
            omega_k_plot(directory)      
        else:
            multi_omega_k_plot(multiple)



    if oc==1: #Potential plot
        if size(multiple)==1:
            potential_plot(directory)
        else:
            extra_win("400x200",'You have chosen more than one file')

    
 
    if op==1: #Normalied potential plot        
            potential_norm_plot(multiple)




    if od==1 and ogeot==1: #geo plot (all graphs together)

        if size(multiple)==1:
            geo_tg(directory)            
        else:
            extra_win("400x200",'You have chosen more than one file')


    if od==1 and ogeos==1: #geo plot (split graphs)

        if size(multiple)==1:
            geo_sp(directory)
        else:
            extra_win("400x200",'You have chosen more than one file')


    if oo==1: #omega(t) plot
    
        field_line(multiple,10)
        

    plt.show()
   
                      




    #########################################  BROWSER FUNCTION ######################################################

def browse():
    directory=askopenfilenames(filetypes = (("in files","*.in"),("all files","*.*")))
    if size(directory)==1:
        direct=directory[0]
        lis=direct.split(".")
        s=size(lis)
        extension=lis[s-1]
        data=open(direct)
        txt=data.read()
        out.set(txt)
        if direct.endswith(".in"):
            dire.set(direct)
            che.set(direct)

        else:
            Window2 = Toplevel() 
            Window2.title("Information")
            Window2.geometry("300x200")
            Inform=Label(Window2,text='The file selected is a .' + extension + ' file \n \n Please browse for a .in file', width=30, height=5).place(x=40,y=20)
            Exit=Button(Window2,text='Accept',command=Window2.destroy).place(x=120,y=160)
    else:
        lis=directory[0].split("/")
        direct=""
        for i in arange(1,size(lis)-1):
            element=lis[i]
            direct=direct+"/"+element
        dire.set("Multiple files selected in "+ direct )
        director=directory[0]
        for n in arange(1,size(directory)):
            elemento=directory[n]
            if elemento.endswith(".in"):
                director=director + ",," +elemento   
            else:
                Window2 = Toplevel() 
                Window2.title("Information")
                Window2.geometry("600x200")
                Inform=Label(Window2,text='You have selected '+ elemento +' file \n \n Please, select only .in files' , width=75, height=5).place(x=3,y=20)
                Exit=Button(Window2,text='Accept',command=Window2.destroy).place(x=270,y=160)
        che.set(director)

  ############################################## EUTERPE BROWSER #########################################################

def browse2():
    eutdir=askopenfilenames(filetypes = (("at files","*.dat"),("all files","*.*")))
    eu.set(eutdir[0])
    
def eusearch():
    eutdir=eu.get()
    svalue=s.get()   
    euval=euprof(eutdir,svalue)
    inset='                                                '+'\n \n'+ 'vmec_parameters \n \n torflux='+str(format9(svalue))+'\n \n&parameters'+'\n \n'+' nine='+str(format9(interpol(euval[0],euval[1],svalue)))+'\n'+' tite='+str(format9(interpol(euval[0],euval[2],svalue)))+'\n \n&species_parameters'+'\n \n'+' dens='+str(format9(interpol(euval[0],euval[3],svalue)))+'\n'+' temp='+str(format9(interpol(euval[0],euval[4],svalue)))+'\n'+' tprim='+str(format9(interpol(euval[0],euval[5],svalue)))+'\n'+' fprim='+str(format9(interpol(euval[0],euval[6],svalue)))+'\n \n&species_parameters_2'+'\n \n'+' dens='+str(format9(interpol(euval[0],euval[7],svalue)))+'\n'+' temp='+str(format9(interpol(euval[0],euval[8],svalue)))+'\n'+' tprim='+str(format9(interpol(euval[0],euval[9],svalue)))+'\n'+' fprim='+str(format9(interpol(euval[0],euval[10],svalue)))



    inputset.set(inset)
        
  ############################################## INFORMATION FUNCTION #####################################################

def about_me():
    Inform_window=Toplevel()
    Inform_window.title("GUI Information")
    Inform_window.geometry("400x170")
    Information=Label(Inform_window,text="\n \n Version 0.3 \n \n This GUI is being created by Antonio González Jerez \n \n email: gugui95@gmail.com \n \n email: Antonio.Gonzalez@externos.ciemat.es" ).pack()




###################################### FUNCTION TO SHOW OUTPUTS FILES #####################################################
def outputs():    
    ok=dire.get()
    option=ext.get()
    if ok !="Browse your input directory.in":
        exten=option.split(".")
        extension=exten[1]
        filex=separate(ok)
        use=filex+"."+extension
        os.system('emacs '+use)

#        datainput=open(use)
#        inputtxt=datainput.read()
#        out.set(inputtxt)



################################################ SYMMETRY ######################################################

def symm():
    ofile=dire.get()
    infi=separate(ofile)
    mult=che.get()
    multiple=mult.split(",,")
    if size(multiple)==1:        
        infile=infi+'.out.nc'   
        pp_symmetry(infile)
    else:
            extra_win("400x200",'You have chosen more than one file')        
    plt.show()


################################################ GVMUS ######################################################

def gvmus():
    ext=gve.get()
    ofile=dire.get()
    infi=separate(ofile)
    mult=che.get()
    multiple=mult.split(",,")
    infi_s=infi.split('/')
    element=""
    directory=infi_s[0]
 
    if size(multiple)==1:
        for i in arange(1,size(infi_s)-1):
            element=infi_s[i]
            directory=directory + '/' + element
        outname=infi_s[size(infi_s)-1]
        d=directory + '/' + 'post_processing'
        createFolder(d)      
        infile=infi+'.out.nc' 
        f=d+'/'+outname
        pp_gvmus_video(infile,f,ext)
        video=outname+'_gvmus'+ext
        extra_win("400x200",'video saved in post_processing in'+'\n'+d)
        os.chdir(d)
        os.system('mpv '+video)

    else:
        extra_win("400x200",'You have chosen more than one file')  
          


################################################ GZVS ######################################################
def gzvs():
    ext=gze.get()
    ofile=dire.get()
    infi=separate(ofile)
    mult=che.get()
    multiple=mult.split(",,")
    infi_s=infi.split('/')
    element=""
    directory=infi_s[0]
 
    if size(multiple)==1:
        for i in arange(1,size(infi_s)-1):
            element=infi_s[i]
            directory=directory + '/' + element
        outname=infi_s[size(infi_s)-1]
        d=directory + '/'+ 'post_processing'
        createFolder(d)      
        infile=infi+'.out.nc' 
        f=d+'/'+outname
        pp_gzvs_video(infile,f,ext)
        video=outname+'_gzvs'+ext
        extra_win("400x200",'video saved in post_processing in'+'\n'+d)
        os.chdir(d)
        os.system('mpv '+video)

    else:
        extra_win("400x200",'You have chosen more than one file')        


  ##########################################################################################################################
  ##########################################################################################################################
  ##########################################################################################################################
  ##########################################################################################################################
                     

                                         #Window Creation#


Window=Tk()
Window.geometry("1900x1060")
Window.title("Stella GUI")

                      #Variable definition

a=IntVar()
b=IntVar()
c=IntVar()
d=IntVar()
p=IntVar()
o=IntVar()
s=DoubleVar()
geo1=IntVar()
geo2=IntVar()
dire=StringVar()
eu=StringVar()
out=StringVar()
inputset=StringVar()
runsdir=StringVar()
ext=StringVar()
gve=StringVar()
gze=StringVar()
directory=[]
che=StringVar()         


                     #Menu creation

topmenu=Menu(Window)
Window.config(menu=topmenu)

menuhelp=Menu(topmenu)
menuhelp.add_command(label="About",command=about_me)

menuexit=Menu(topmenu)
menuexit.add_command(label="Exit",command=Window.destroy)

topmenu.add_cascade(label="Help",menu=menuhelp)
topmenu.add_cascade(label="Exit",menu=menuexit)




                 #Checkbutton options creation


cha=Checkbutton(Window, text="\u03c9(t)", variable=a,  onvalue=1, offvalue=0).place(x=700, y=320) 
chb=Checkbutton(Window, text="\u03c9(k)", variable=b,  onvalue=1, offvalue=0).place(x=700, y=420)
chc=Checkbutton(Window, text="\u03c6(\u03b6)", variable=c,  onvalue=1, offvalue=0).place(x=700, y=520) 
chp=Checkbutton(Window, text="|\u03c6(\u03b6)|/|\u03c6(0)|", variable=p,  onvalue=1, offvalue=0).place(x=700, y=620)
chextra=Checkbutton(Window, text="", variable=d,  onvalue=1, offvalue=0).place(x=700, y=820) 

geo=  Menubutton (Window, text="geo", relief=RAISED )
geo.place(x=725,y=820)
geo.menu  =  Menu ( geo)
geo["menu"]  =  geo.menu
geo1.set(1)  #To set together as elected
geo.menu.add_checkbutton ( label="Plot together", variable=geo1, onvalue=1, offvalue=0)
geo.menu.add_checkbutton ( label="Plot Split", variable=geo2, onvalue=1, offvalue=0)

cho=Checkbutton(Window, text="field line (under maintenance)", variable=o,  onvalue=1, offvalue=0).place(x=700, y=720)


                #Plot button creation

run=Button(Window, text="APPLY", command=options, relief="raised").place(x=700, y=220) 

                #File box creation

dire.set("Browse your input directory.in")
direct=Entry(Window,textvariable=dire,relief="sunken", width=60, state='readonly').place(x=100, y=100)

                    #Browser creation

search=Button(Window,text="Browse", relief="raised",command=browse).place(x=650, y=100)

                  #2nd browser creation

eu.set("Browse your Euterpe file")
eut=Entry(Window,textvariable=eu,relief="sunken", width=30, state='readonly').place(x=100, y=350)
s_lab=Label(Window,text='s=',width=4).place(x=460,y=350)
s_elec=Entry(Window,textvariable=s,width=5).place(x=500,y=350)
search2=Button(Window,text="Browse", relief="raised",command=browse2).place(x=350, y=345)
view=Button(Window,text="VIEW", relief="raised",command=eusearch).place(x=560, y=345)


                   #input box creation


inputset.set("\n INPUT PARAMETERS")
frame1=Frame(Window,highlightbackground="green", highlightcolor="green", highlightthickness=1,width=300, height=450)
frame1.place(x=100, y=450)
canvas=Canvas(frame1,width=300, height=450,bg='white')
scrollbar1 = Scrollbar(frame1, orient='vertical', command=canvas.yview)
scrollbar1.pack(side='right',fill='y')
scrollbar2 = Scrollbar(frame1, orient='horizontal', command=canvas.xview)
scrollbar2.pack(side='bottom',fill='x')
canvas.pack()
frame2 = Frame(canvas)
canvas.configure(yscrollcommand=scrollbar1.set,xscrollcommand=scrollbar2.set,)
canvas.create_window(0, 0, window=frame2, anchor='nw')
linput=Label(frame2,bg="white",textvariable=inputset).pack()

                    #Output extension label creation

ext.set(".in")
labinput=OptionMenu(Window, ext,".in",".final_fields",".fluxes",".omega",".out").place(x=307,y=169)
Accept=Button(Window,text='OPEN FILE WITH EXTENSION',command=outputs).place(x=98, y=170)


                    #Python prompt creation

termf=Frame(Window,bg="black", highlightbackground="white", highlightcolor="black", highlightthickness=5,width=750, height=750)
termf.place(x=1080, y=200)
wid = termf.winfo_id()
os.system('xterm -into %d -geometry 400x500 -sb &' % wid)

                    #Trial button

SYM=Button(Window,text='Symmetry',command=symm).place(x=850, y=370)
GVMUS=Button(Window,text='gvmus video',command=gvmus).place(x=850, y=470)
gve.set(".gif")
gvmusext=OptionMenu(Window, gve,".gif",".mp4").place(x=980, y=470)
GZVS=Button(Window,text='gzvs video',command=gzvs).place(x=850, y=570)
gze.set(".gif")
gzvsext=OptionMenu(Window, gze,".gif",".mp4").place(x=980, y=570)



      
                    #Window close

Window.mainloop()








