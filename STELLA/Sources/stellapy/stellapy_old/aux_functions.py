import numpy
from numpy import *
################################ To create a directory (if it does not exist) ##########################################

def createFolder(directory):
        if not os.path.exists(directory):
            os.mkdir(directory)

########################################## To quit the extension #######################################################

#To quit a directory extension

def separate(case):
    dis=case.split(".in")
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
 
def color():
    colors=[]
    col=['r','g','b','y','m','c','k','gray','pink','chocolate','gold','slateblue','turquoise','darkseagreen','orange','tan','fuchsia','sienna','crimson','coral']
    for i in arange(0,100):
        colors=colors+col
    return colors
