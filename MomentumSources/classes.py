'''
classes.py contains functions and classes used within the plotting codes
'''
import gzip
import numpy as np
import matplotlib.dates as md
import matplotlib.patches as patches
import scipy.interpolate as sci
import datetime as dt

def gridpoints(x,N):
    ''' 
    xi = gridpoints(x,N) is a function that returns interpolated values along a grid
    INPUTS: x is a numpy array 
            N is an integer for the desired number of gridpoints
    OUTPUT: xi is a numpy array and is the interpolated over x
    '''
    xi = np.linspace(x.min(),x.max(),N)
    return xi

def twodinterpolate(x,y,xi,yi,z):
    '''
    twodinterpolate(x,y,xi,yi,z) interpolates a function in 2d using a cubic
    splines interpolation.
    INPUTS: x is a 1d numpy array and is the x-coordinate
            y is a 1d numpy array and is the ordinate
            xi is a 1d numpy array and is the interpolated x-coordinate
            yi is a 1d numpy array and is the interpolated ordinate
            z is a 2d numpy array and the array to be interpolated
    OUTPUT: zi is a 2d numpy array and is the interpolated array
    '''
    zi = sci.griddata((x,y),z,(xi[None,:],yi[:,None]),method='cubic')
    return zi

def verticaldraw(reversal):
    '''
    verticaldraw is a function that returns a string value
    INPUT: reversal is a string that takes two possible string values
            'runFlipNorth' is the South-to-North IMF configuration 
            'runFlipSouth' is the North-to-South IMF configuration 
    OUTPUT: timedraw is a string for the time of IMF polarity reversal
    '''
    if (reversal == 'runFlipNorth'):
        timedraw = '04:00'
    if (reversal == 'runFlipSouth'):
        timedraw = '08:00'
    return timedraw

def ordinate(cut):
    '''
    p = ordinate(cut) is a function that returns an integer depending on the plane 
    of data for the 2D BATS-R-US output (Z = 0 or Y = 0) this is used for reading 
    the data
    INPUT: cut is a string and can only take two possible values
           -'z' for the Z=0 cut
           -'y' for the Y=0 cut
    OUTPUT: p is an integer that is used to access the BATS-R-US output data
    '''
    if (cut == 'z'):
        p = 1
    if (cut == 'y'):
        p = 2
    return p

def timestring(k):
    '''
    hr,mn = timestring(k) returns the hour, minute value for a given k
    INPUT: k is an integer. This is used within for loops in the code to access
    the BATS-R-US files. The filenames contain hr and minutes.
    OUTPUT: hr,mn are integers and stand for hour and minute respectively
    '''
    hr = k//60             #updated for python 3 to return integer
    if (hr < 10):
        hr = '0'+str(hr)
    else:
        hr = str(hr)
    mn = k%60
    if (mn < 10):
        mn = '0'+str(mn)
    else:
        mn = str(mn)
    return hr,mn

def satellite_format(reversal):
    xfmt = md.DateFormatter('%H:%M')
    if (reversal == 'runFlipNorth'):
        timeticks = [md.date2num(dt.datetime.strptime('02:00','%H:%M')), \
                     md.date2num(dt.datetime.strptime('04:00','%H:%M')),\
                     md.date2num(dt.datetime.strptime('06:00','%H:%M')), \
                     md.date2num(dt.datetime.strptime('08:00','%H:%M')),\
                     md.date2num(dt.datetime.strptime('10:00','%H:%M')), \
                     md.date2num(dt.datetime.strptime('12:00','%H:%M'))]
    if (reversal == 'runFlipSouth'):
        timeticks = [md.date2num(dt.datetime.strptime('06:00','%H:%M')), \
                     md.date2num(dt.datetime.strptime('07:00','%H:%M')), \
                      md.date2num(dt.datetime.strptime('08:00','%H:%M')),\
                      md.date2num(dt.datetime.strptime('09:00','%H:%M')),\
                     md.date2num(dt.datetime.strptime('10:00','%H:%M')), \
                     md.date2num(dt.datetime.strptime('11:00','%H:%M')), \
                             md.date2num(dt.datetime.strptime('12:00','%H:%M'))]
    return(xfmt,timeticks)


def satellite_format_hour(reversal):
    xfmt = md.DateFormatter('%H')
    if (reversal == 'runFlipNorth'):
        timeticks = [md.date2num(dt.datetime.strptime('02','%H')), \
                     md.date2num(dt.datetime.strptime('04','%H')),\
                     md.date2num(dt.datetime.strptime('06','%H')), \
                     md.date2num(dt.datetime.strptime('08','%H')), \
                     md.date2num(dt.datetime.strptime('10','%H')), \
                     md.date2num(dt.datetime.strptime('12','%H'))]
    if (reversal == 'runFlipSouth'):
        timeticks = [md.date2num(dt.datetime.strptime('06','%H')), \
                     md.date2num(dt.datetime.strptime('07','%H')), \
                      md.date2num(dt.datetime.strptime('08','%H')),\
                      md.date2num(dt.datetime.strptime('09','%H')),\
                     md.date2num(dt.datetime.strptime('10','%H')), \
                     md.date2num(dt.datetime.strptime('11','%H')), \
                             md.date2num(dt.datetime.strptime('12','%H'))]
    return(xfmt,timeticks)


def readsolarwind(inputfile):
    dataline = 3 #line where data is read (starting from zero)
    #declare empty arrays
    bx,by,bz = [],[],[]
    vx,vy,vz= [],[],[]
    timestamp = []
    density,temperature = [],[]
    with open(inputfile,"r") as infile:
        i = 0
        for line in infile:
            if (i <= 3):
                i = i +1
            else:
                columns = line.rsplit()
                hour,minute,sec = columns[3],columns[4],columns[5]
                temptime = dt.datetime.strptime(hour+':'+minute+':'+sec,'%H:%M:%S')
                timestamp.append(temptime)
                bx.append(float(columns[7]))
                by.append(float(columns[8]))
                bz.append(float(columns[9]))
                vx.append(float(columns[10]))
                vy.append(float(columns[11]))
                vz.append(float(columns[12]))                                 
                density.append(float(columns[13]))
                temperature.append(float(columns[14]))

                i = i +1
        #sample every 60 seconds
        j = 1
        timestamp = md.date2num(timestamp[::j])
        bx,by,bz = np.array(bx[::j]),np.array(by[::j]),np.array(bz[::j])
        vx,vy,vz = np.array(vx[::j]),np.array(vy[::j]),np.array(vz[::j])
        density,temperature = np.array(density[::j]),np.array(temperature[::j])
    return(timestamp,bx,by,bz,vx,vy,vz,density,temperature)

def readBats(infile):
    '''
    readBats reads the tecplot formatted BATS-R-US file
    INPUT: infile is the filename (string)
    OUTPUT: -data is a 2d numpy array containing BATS-R-US simulation 
             information
            -connectivity is the connectivity list and describes which lines
             of the codes form the vertices of a cube. It is a 2d numpy array.
            -timetick is the simulation time (string)
    '''
    i=0
    data = []
    connectivity = []
    for line in infile: #read line by line the input file
        columns = line.rsplit()
        temp = line.rsplit()
        if (i == 0): #this case is meant to read the simulation time
            times=columns[5].replace('"',"")
            times = times.split('.',1)[0]
            timetick = times
        if (i >22): #read after header information from TecPlot file
            if (len(temp) >4):
                #columns = list(map(float,columns))
                data.append(columns)
            if (len(temp)==4):
                #columns = list(map(int,columns))
                connectivity.append(columns)
        i=i+1
    data = np.array(data)
    data = data.astype(float)
    connectivity = np.array(connectivity)
    connectivity = connectivity.astype(int)
    return(data,connectivity,timetick)



    
def readgeopause(infile):
    i=0
    data = np.loadtxt(infile,skiprows=1)
    data = np.array(data)
    row,column=data.shape
    xcoord,ycoord = np.zeros(row),np.zeros(row)
    xcoord[0:row]=data[0:row,0]
    ycoord[0:row]=data[0:row,1]
    return(xcoord,ycoord)

def readlastclosed(infile):
    '''
    readlastclosed reads the file containing the locus of points
    forming the last closed field lines.
    INPUT: infile is a string with the filename
    OUTPUT: xday a 1d numpy array containing the x dayside values
            yday is a 1d numpy array containing the dayside ordinate values
            xnight is a 1d numpy array containing the x nightside values
            ynight is a 1d numpy array containing the nightside ordinate values
    '''
    i=0
    data = []
    xnight,ynight,xday,yday = [],[],[],[]
    for line in infile: #read line by line the input file
        columns = line.rsplit()
        if (i > 2):
            if (len(columns)==4):
                xnight.append(float(columns[0]))
                ynight.append(float(columns[1]))
                xday.append(float(columns[2]))
                yday.append(float(columns[3]))
            if (len(columns)==2):
                xnight.append(float(columns[0]))
                ynight.append(float(columns[1]))
        i=i+1
    return(xnight,ynight,xday,yday)

def geopause_daynight(xcoord,ycoord): #separate the data into dayside and nightside
    '''
    geopause_daynight split the data into a dayside component and nightside component
    INPUT: xcoord is a 1d numpy array containing the x values
           ycoord is a 1d numpy array with the ordinate values
    OUTPUT: datadayx is a 1d numpy array containing the x dayside values
            datadayy is a 1d numpy array containing the dayside ordinate values
            datanightx is a 1d numpy array containing the x nightside values
            datanighty is a 1d numpy array containing the nightside ordinate values
    '''
    datadayx,datanightx = [],[]
    datadayy,datanighty = [],[]
    row = len(xcoord)
    '''
    The for loop splits the data according to the x value.
    If x is positive then it is along the dayside.
    If x is negative then it is along the nightside.
    Then, the code will append the value to the appropriate array.
    '''
    for i in range(0,row):
        x = xcoord[i]
        y = ycoord[i]
        if (x >= 0.):
            datadayx.append(x)
            datadayy.append(y)
        if (x < 0.):
            datanightx.append(x)
            datanighty.append(y)
    datadayx = np.array(datadayx)
    datadayy = np.array(datadayy)
    datanightx = np.array(datanightx)
    datanighty = np.array(datanighty)
    return datadayx,datadayy,datanightx,datanighty

def geopause_daynight_sunearth(xday,yday,xnight,ynight):
    '''
    geopause_daynight_sunearth takes the coordinates of the geopause
    and extracts them along the Sun-Earth line and finds the closest geopause
    to Earth. For the 2D
    BATS-R-US output, this corresponds to finding where Y=0 with the Z=0 cut
    or Z=0 for the Y=0 cut.
    INPUT: xday is a numpy array and the x coordinates of the dayside geopause
           yday is a numpy array and the ordinate of the dayside geopause
           xnight is a numpy array and the x coordinates of the nightside geopause
           ynight is a numpy array and the ordinates of the nightside geopause
    '''
    xtemp,xtemp2 =[],[]
    ytemp,ytemp2 =[],[]
    daysize = len(xday)
    nightsize = len(xnight)
    for i in range(0,daysize):
        if (yday[i] == 0.):
            if (xday[i] != 0.):
                xtemp.append(xday[i])
    for i in range(0,nightsize):
        if (ynight[i] == 0.):
            if (xnight[i] != 0.):
                xtemp2.append(xnight[i])
          
    xtemp = np.array(xtemp)
    xtemp2 = np.array(xtemp2)
    xdayearth = np.min(xtemp)
    xnightearth = np.max(xtemp2)
    return xdayearth, xnightearth

def geopause_daynight_sunearthall(xday,yday,xnight,ynight,tvert):
    '''
    geopause_daynight_sunearth takes the coordinates of the geopause
    and extracts them along the Sun-Earth line.
    BATS-R-US output, this corresponds to finding where Y=0 with the Z=0 cut
    or Z=0 for the Y=0 cut.
    INPUT: xday is a numpy array and the x coordinates of the dayside geopause
           yday is a numpy array and the ordinate of the dayside geopause
           xnight is a numpy array and the x coordinates of the nightside geopause
           ynight is a numpy array and the ordinates of the nightside geopause
    '''
    xtemp,xtemp2 =[],[]
    ytemp,ytemp2 =[],[]
    tdaytemp,tnighttemp = [],[]
    daysize = len(xday)
    nightsize = len(xnight)
    for i in range(0,daysize):
        if (yday[i] >= -0. and yday[i] <= 0.):
            if (xday[i] != 0.):
                xtemp.append(xday[i])
                tdaytemp.append(tvert)
    for i in range(0,nightsize):
        if (ynight[i] >= -0. and ynight[i] <=0.):
            if (xnight[i] != 0.):
                xtemp2.append(xnight[i])
                tnighttemp.append(tvert)
    xtemp = np.array(xtemp)
    xtemp2 = np.array(xtemp2)
    xdayearth = xtemp
    xnightearth = xtemp2
    tday = tdaytemp
    tnight = tnighttemp
    return xdayearth, xnightearth, tday,tnight

def numbertoletter(i,j): #map the array number to a letter for the figure
    if (i==0 and j==0):
        letter='(a)'
    if (i==0 and j==1):
        letter='(b)'
    if(i==1 and j==0):
        letter='(c)'
    if(i==1 and j==1):
        letter='(d)'
    if(i==2 and j==0):
        letter='(e)'
    if(i==2 and j==1):
        letter='(f)'
    return letter

def numbertoimfletter(j):
    number = j
    if (number==0):
        letter='(a)'
    if (number==1):
        letter='(b)'
    return letter
    
def ylabel(cut):
    if (cut == 'y'):
        ordinatelabel = r'Z [$R_E$]'
    if (cut == 'z'):
        ordinatelabel = r'Y [$R_E$]'
    return(ordinatelabel)
def IMFtitle(reversal):
    if (reversal=='runFlipNorth'):
        IMF = 'Northern IMF Configuration'
    if (reversal=='runFlipSouth'):
        IMF = 'Southern IMF Configuration'
    return(IMF)


class PlotProperties:

    def __init__(self,name,xmin,xmax,ymin,ymax,compressed,cut,mode):
        self.name=name             #filename
        self.xmin=xmin
        self.xmax=xmax
        self.ymin=ymin
        self.ymax=ymax
        self.compressed=compressed #whether file is compressed or not
        self.cut=cut               #cut is either Z=0 or Y=0
        self.mode=mode

    def get_data(self): 
        '''
        get_data extracts the BATS-R-US data
        OUTPUT: data is the raw BATS-R-US data returned as a numpy array
                connectivity is the connnectivity list returned as a numpy array
                             and designates the cell members
                timetick extracts the time in the file
        '''
        if (self.compressed == True): #for compressed data
            with gzip.open(self.name,"rt") as infile:
                data,connectivity,timetick =  readBats(infile)
        if (self.compressed == False):
            with open(self.name,"r") as infile:
                data,connectivity,timetick =  readBats(infile)
        return data,connectivity,timetick

    def data_filter(self,data):
        '''
        data_filter filters the BATS-R-US data to only the data within the domain of interest
        '''
        datanew = []
        row,column=data.shape
        p = ordinate(self.cut)
        for i in range(0,row):
            x = data[i,0]
            y = data[i,p]
            if (self.mode == '2d'): #The following are meant to extract a subset of data
                within_domain = x >=self.xmin and x <= self.xmax and y >= self.ymin and y <= self.ymax
            if (within_domain):
                datanew.append(data[i,0:column])
        datanew = np.array(datanew)
        return datanew

    
    def data_daynight(self,data): #separate the data into dayside and nightside
        dataday,datanight = [],[]    
        row,column=data.shape
        for i in range(0,row):
            x = data[i,0]
            if (x >= 0.):
                dataday.append(data[i,0:column])
            if (x <= 0.):
                datanight.append(data[i,0:column])
        dataday = np.array(dataday)
        datanight = np.array(datanight)
        return dataday,datanight

    def data_readable(self,data): #put data into different arrays
        '''
        data_readable processes the BATS-R-US data into readable formats
        OUTPUT: - xyz[a,b] is the position of the cell vertex and a numpy array.
                  a is the array index representing the cell vertex.
                  b is the component of the position vector.
                  b = 0 = x
                  b = 1 = y
                  b = 2 = z
                - bfield[a,b] is the magnetic field vector and a numpy array.
                  a is the array index representing the cell vertex.
                  b is the array index representing the component of the vector.
                  b takes values from 0 to 2
                - J[a,b] is the current density and is represented by a numpy array.
                  a is the  array index representing cell vertex. 
                  b is the componeffnt of the vector. 
                - stateFluid[a,b,c] is a numpy array containing the mass density,
                  velocity, and thermal pressure of the fluid. 
                  The a index represents the cell vertex
                  The b index takes values 
                  from 0 to 4 and represent the mass density(0), x-velocity(1),
                  y-velocity(2), z-velocity(3), and thermal pressure (4).
                  The c index represents the fluid. For two fluid
                  simulations, the index will take integer values from 0 to 1.
                  For three fluid simulations, the values will take values from 
                  0 to 2. The 0 value will always represent the solar wind. While
                  the other values represent the ionospheric species.
                - nFluid is an integer and represents the number of fluids
        '''
        row,column=data.shape
        nFluid=(column-9)//5 #9 corresponds to x,y,z,bx,by,bz,jx,jy,jz
        xyz,bfield,J = np.zeros((row,3)),np.zeros((row,3)),np.zeros((row,3))
        stateFluid=np.zeros((row,5,nFluid))
        xyz[0:row,0]= data[0:row,0] # 0 stands for x
        xyz[0:row,1]= data[0:row,1] # 1 stands for y
        xyz[0:row,2]= data[0:row,2] # 2 stands for z
        bfield[0:row,0]= data[0:row,7]
        bfield[0:row,1]= data[0:row,8]
        bfield[0:row,2]= data[0:row,9]
        J[0:row,0]= data[0:row,column-3]
        J[0:row,1]= data[0:row,column-2]
        J[0:row,2]= data[0:row,column-1]
        for i in np.arange(0,nFluid):
            if (i == 0):
                stateFluid[0:row,0,i]=data[0:row,3] #rho [amu/cm^3]
                stateFluid[0:row,1,i]=data[0:row,4] #vx
                stateFluid[0:row,2,i]=data[0:row,5] #vy
                stateFluid[0:row,3,i]=data[0:row,6] #vz
                stateFluid[0:row,4,i]=data[0:row,10]#p
            if (i != 0):
                stateFluid[0:row,0,i]=data[0:row,6+5*i] #rho [amu/cm^3]
                stateFluid[0:row,1,i]=data[0:row,7+5*i] #vx
                stateFluid[0:row,2,i]=data[0:row,8+5*i] #vy
                stateFluid[0:row,3,i]=data[0:row,9+5*i] #vz
                stateFluid[0:row,4,i]=data[0:row,10+5*i] #p
                #stateFluid[0:row-1,0:4,i]=data[0:row-1,6+5*i:10+5*i]
        return(xyz,bfield,J,stateFluid,nFluid)
    
    def data_readableSI(self,data): #put data into different arrays
        '''
        data_readableSI processes the BATS-R-US data into readable formats
        OUTPUT: - xyz[a,b] is the position of the cell vertex and a numpy array.
                  a is the array index representing the cell vertex.
                  b is the component of the position vector.
                  b = 0 = x
                  b = 1 = y
                  b = 2 = z
                - bfield[a,b] is the magnetic field vector and a numpy array.
                  a is the array index representing the cell vertex.
                  b is the array index representing the component of the vector.
                  b takes values from 0 to 2
                - J[a,b] is the current density and is represented by a numpy array.
                  a is the  array index representing cell vertex. 
                  b is the componeffnt of the vector. 
                - stateFluid[a,b,c] is a numpy array containing the mass density,
                  velocity, and thermal pressure of the fluid. 
                  The a index represents the cell vertex
                  The b index takes values 
                  from 0 to 4 and represent the mass density(0), x-velocity(1),
                  y-velocity(2), z-velocity(3), and thermal pressure (4).
                  The c index represents the fluid. For two fluid
                  simulations, the index will take integer values from 0 to 1.
                  For three fluid simulations, the values will take values from 
                  0 to 2. The 0 value will always represent the solar wind. While
                  the other values represent the ionospheric species.
                - nFluid is an integer and represents the number of fluids
        '''
        row,column=data.shape
        nFluid=int((column-9)/5) #9 corresponds to x,y,z,bx,by,bz,jx,jy,jz
        xyz,bfield,J = np.zeros((row,3)),np.zeros((row,3)),np.zeros((row,3))
        stateFluid=np.zeros((row,5,nFluid))
        qe = 1.602176621e-19
        amutokg = 1.660539e-27
        invcm3toinvm3 = 1.e6
        nPatoPa=1.e-9
        muAtoA=1.e-6
        kmtom = 1.e3
        Re =6.3781e6
        xyz[0:row,0]= data[0:row,0]*Re # 0 stands for x
        xyz[0:row,1]= data[0:row,1]*Re # 1 stands for y
        xyz[0:row,2]= data[0:row,2]*Re # 2 stands for z
        bfield[0:row,0]= data[0:row,7]*1.e-9
        bfield[0:row,1]= data[0:row,8]*1.e-9
        bfield[0:row,2]= data[0:row,9]*1.e-9
        J[0:row,0]= data[0:row,column-3]*muAtoA
        J[0:row,1]= data[0:row,column-2]*muAtoA
        J[0:row,2]= data[0:row,column-1]*muAtoA
        for i in np.arange(0,nFluid):
            if (i == 0):
                stateFluid[0:row,0,i]=data[0:row,3]*amutokg*invcm3toinvm3 #rho [kg/m^3]
                stateFluid[0:row,1,i]=data[0:row,4]*kmtom #vx
                stateFluid[0:row,2,i]=data[0:row,5]*kmtom #vy
                stateFluid[0:row,3,i]=data[0:row,6]*kmtom #vz
                stateFluid[0:row,4,i]=data[0:row,10]*nPatoPa#p
            if (i != 0):
                stateFluid[0:row,0,i]=data[0:row,6+5*i]*amutokg*invcm3toinvm3 #rho [kg/m^3]
                stateFluid[0:row,1,i]=data[0:row,7+5*i]*kmtom #vx
                stateFluid[0:row,2,i]=data[0:row,8+5*i]*kmtom #vy
                stateFluid[0:row,3,i]=data[0:row,9+5*i]*kmtom #vz
                stateFluid[0:row,4,i]=data[0:row,10+5*i]*nPatoPa #p
                #stateFluid[0:row-1,0:4,i]=data[0:row-1,6+5*i:10+5*i]
        return(xyz,bfield,J,stateFluid,nFluid)

def draw_earth(ax):
    '''
    draw_earth draws Earth with a dayside and nightside
    INPUT: ax is an Axes class
    '''
    earth =  patches.Wedge(
        (0.0, 0.0),     # (x,y)
        2.5,            # radius
        0,             # theta1 (in degrees)
        360,             # theta2
        facecolor='#808080',zorder=3)
    ax.add_patch(earth)
    dayside =  patches.Wedge(
        (0.0, 0.0),     # (x,y)
        1,            # radius
        270,             # theta1 (in degrees)
        90,             # theta2
        facecolor='w',zorder=4)
    ax.add_patch(dayside)
    nightside =  patches.Wedge(
        (0.0, 0.0),     # (x,y)
        1,            # radius
        90,             # theta1 (in degrees)
        270,             # theta2
        facecolor='k',zorder=4)
    ax.add_patch(nightside)
    return

