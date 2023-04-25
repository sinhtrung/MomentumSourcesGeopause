import sys
import os
import glob
import gzip
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.dates as md
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
import scipy.interpolate as sci
from matplotlib.colors import LogNorm
import matplotlib.colors as colors
import datetime as dt
from matplotlib.figure import Figure
from matplotlib import rc
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import cProfile
import multiprocessing
from classes import *



#geopause component draws a single geopause type onto a panel component
def geopausecomponent(ax,i,j,k,reversal,quantity,cut,Title,xmin,xmax,zmin,zmax):
    home=os.getcwd()+'/'
    file1,file2,file3,file4=[],[],[],[]
    hr,mn=timestring(k)
    time = 't'+hr+mn
    root1=home+'Data/Geopause/'+'SwIono/'+reversal+'/'
    root2=home+'Data/Geopause/'+'SwHIonoO/'+reversal+'/'
    root3=home+'Data/Geopause/'+'SwIonoO/'+reversal+'/'
    root4=home+'Data/Geopause/'+'SwIonoO28amu/'+reversal+'/'
    file1.extend(glob.glob(root1+quantity+cut+'*'+time+'*'))
    file2.extend(glob.glob(root2+quantity+cut+'*'+time+'*'))
    file3.extend(glob.glob(root3+quantity+cut+'*'+time+'*'))
    file4.extend(glob.glob(root4+quantity+cut+'*'+time+'*'))
    #letter = numbertoletter(i,j)
    filepath = file1[0]
    print(filepath)
    with open(filepath,'r') as infile:
        x1,y1=readgeopause(infile)
    filepath = file2[0]
    with open(filepath,'r') as infile:
        x2,y2=readgeopause(infile)
    filepath = file3[0]
    with open(filepath,'r') as infile:
        x3,y3=readgeopause(infile)
    filepath = file4[0]
    with open(filepath,'r') as infile:
        x4,y4=readgeopause(infile)
    l1=ax.scatter(x1,y1,s=1,color='k',label=Title[0])
    l2=ax.scatter(x2,y2,s=1,color='b',label=Title[1])
    l3=ax.scatter(x3,y3,s=1,color='r',label=Title[2])
    l4=ax.scatter(x4,y4,s=1,color='g',label=Title[3])
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(zmin,zmax)
    draw_earth(ax)
    ax.set_xlabel(r'X [$R_E$]')
    ordinatelabel = ylabel(cut)
   
    if (j == 1):
        ax.yaxis.tick_right()
        #ax.set_yticklabels([])
    ax.set_ylabel(ordinatelabel)
    IMF = IMFtitle(reversal)
    if (quantity == 'pressure'):
        qlabel = 'Pressure Geopause for the'
    if (quantity == 'geopause'):
        qlabel = 'Number Density Geopause for the'
    if (quantity == 'rho'):
        qlabel = 'Mass Density Geopause for the'
    #atitle = qlabel+' '+IMF
    '''
    if (i == 0): #placed for animation otherwise remove for general purpose
        ax.set_title(IMF,fontsize=8)
'''
    ax.set_aspect('equal')
    ax.grid()
    #ax.text(0.87,0.05,letter,transform=ax.transAxes,fontsize=8)
    return l1,l2,l3,l4
        
def lastclosedcomponent(ax,j,k,reversal,quantity,cut,Title,xmin,xmax,zmin,zmax):
    home=os.getcwd()+'/'
    root1=home+'Data/Geopause/'+'SwIono/'+reversal+'/'
    root2=home+'Data/Geopause/'+'SwHIonoO/'+reversal+'/'
    root3=home+'Data/Geopause/'+'SwIonoO/'+reversal+'/'
    root4=home+'Data/Geopause/'+'SwIonoO28amu/'+reversal+'/'
    file1,file2,file3,file4=[],[],[],[]
    hr,mn=timestring(k)
    time = 't'+hr+mn
    letter = numbertoimfletter(j)
    file1.extend(glob.glob(root1+quantity+cut+'*'+time+'*'))
    file2.extend(glob.glob(root2+quantity+cut+'*'+time+'*'))
    file3.extend(glob.glob(root3+quantity+cut+'*'+time+'*'))
    file4.extend(glob.glob(root4+quantity+cut+'*'+time+'*'))
    filepath = file1[0]
    print(filepath)
    with open(filepath,'r') as infile:
        xnight1,ynight1,xday1,yday1=readlastclosed(infile)
    filepath = file2[0]
    with open(filepath,'r') as infile:
        xnight2,ynight2,xday2,yday2=readlastclosed(infile)
    filepath = file3[0]
    with open(filepath,'r') as infile:
        xnight3,ynight3,xday3,yday3=readlastclosed(infile)
    filepath = file4[0]
    with open(filepath,'r') as infile:
        xnight4,ynight4,xday4,yday4=readlastclosed(infile)
    draw_earth(ax)
    ax.plot(xnight1,ynight1,color='k',label=Title[0],zorder=0)
    ax.plot(xnight2,ynight2,color='b',label=Title[1],zorder=0)
    ax.plot(xnight3,ynight3,color='r',label=Title[2],zorder=0)
    ax.plot(xnight4,ynight4,color='g',label=Title[3],zorder=0)
    l1,=ax.plot(xday1,yday1,color='k',zorder=0)
    l2,=ax.plot(xday2,yday2,color='b',zorder=0)
    l3,=ax.plot(xday3,yday3,color='r',zorder=0)
    l4,=ax.plot(xday4,yday4,color='g',zorder=0)
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(zmin,zmax)
  
    ordinatelabel = ylabel(cut)
    ax.set_xlabel(r'X [$R_E$]')
    if (j==0):
        ax.set_ylabel(ordinatelabel)
    if (j == 1):
        ax.set_yticklabels([])
    IMF = IMFtitle(reversal)
    ax.set_title(IMF,fontsize=8)
    ax.set_aspect('equal')
    ax.grid()
    ax.text(0.87,0.05,letter,transform=ax.transAxes,fontsize=8)
    return(l1,l2,l3,l4)


def allcomponent(ax,j,k,GM,reversal,cut,Title,xmin,xmax,zmin,zmax):
    home=os.getcwd()+'/'
    Legends = ['Density','Mass Density','Pressure','Magnetopause']
    minimum = 720
    plotnumber =1+12*60
    cut = 'y'
    root=home+'Data/Geopause/'+GM+'/'+reversal+'/'
    file1,file2,file3,file4 = [],[],[],[]
    hr,mn=timestring(k)
    time = 't'+hr+mn
    file1.extend(glob.glob(root+'geopause'+cut+'*'+time+'*'))
    file2.extend(glob.glob(root+'rho'+cut+'*'+time+'*'))
    file3.extend(glob.glob(root+'pressure'+cut+'*'+time+'*'))
    file4.extend(glob.glob(root+'lastclosed'+cut+'*'+time+'*'))
    filepath = file1[0]
    print(filepath)
    with open(filepath,'r') as infile:
        x1,y1=readgeopause(infile)
    filepath = file2[0]
    with open(filepath,'r') as infile:
        x2,y2=readgeopause(infile)
    filepath = file3[0]
    with open(filepath,'r') as infile:
        x3,y3=readgeopause(infile)
    filepath = file4[0]
    with open(filepath,'r') as infile:
        xnight,ynight,xday,yday=readlastclosed(infile)
    l1=ax.scatter(x1,y1,s=1,color='sienna',label=Legends[0])
    l2=ax.scatter(x2,y2,s=1,color='gold',label=Legends[1])
    l3=ax.scatter(x3,y3,s=1,color='k',label=Legends[2])
    l4,=ax.plot(xnight,ynight,color='dodgerblue',label=Legends[3],zorder=0)
    ax.plot(xday,yday,color='dodgerblue',zorder=0)
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(zmin,zmax)
    draw_earth(ax)
    ordinatelabel = ylabel(cut)
    ax.set_xlabel(r'X [$R_E$]')
    if (j==0):
        ax.set_ylabel(ordinatelabel)
    if (j == 1):
        ax.set_yticklabels([])
    IMF = IMFtitle(reversal)
    ax.set_title(IMF,fontsize=8)
    ax.set_aspect('equal')
    ax.grid()
    letter = numbertoimfletter(j)
    ax.text(0.87,0.05,letter,transform=ax.transAxes,fontsize=10)
    return(l1,l2,l3,l4)

def plotpanelgeopause(): #draw the geopauses on 2 x 2 grid
    home=os.getcwd()+'/'
    reversal = ['runFlipNorth','runFlipSouth']
    GM = ['SwIono','SwHIonoO','SwIonoO','SwIonoO28amu']
    Title=['Sw H + Iono H\n Case 1','Sw H + Iono O\n Case 2','Sw H + Iono (H + O)\n Case 3','Sw H + Iono (H + O) at 28 amu/cc \n Case 4']
    minimum = 720
    plotnumber =1+12*60
    mode = '2d'
    cut = ['y','z']
    rc('font',family='serif')
    rc('text', usetex=True)
    xmin,xmax = -60.,10.
    ymin,ymax = -20.,20.
    #for i in range(0,4):
   
   
    quantity=['pressure','geopause','rho']
    plottingname = ['pressurepanel.png','densitypanel.png','rhopanel.png']
    for m in range(0,3):
        fig,axs = plt.subplots(2,2)
        plotquantity=quantity[m]
        for i in range(0,2):
            for j in range(0,2):
            #file location step
            #read data from file step
            #plotting step
            #save file
                 for k in range(minimum,plotnumber):
                     hr,mn=timestring(k)
                     l1,l2,l3,l4=geopausecomponent(axs[j,i],j,i,k,reversal[i],plotquantity,cut[j],Title,xmin,xmax,ymin,ymax)
        plotpath = plottingname[m]
        if (plotquantity == 'pressure'):
            Suptitle = 'Pressure Geopause at'+' '+hr+':'+mn+' UT'
        if (plotquantity == 'geopause'):
            Suptitle = 'Number Density Geopause at'+' '+hr+':'+mn+' UT'
        if (plotquantity == 'rho'):
            Suptitle = 'Mass Density Geopause at'+' '+hr+':'+mn+' UT'
        leg = plt.figlegend((l1,l2,l3,l4),Title,(0.02,0.0),markerscale=3,ncol=4,fontsize=7,frameon=False)
        for t in leg.texts:
            t.set_multialignment('center')
        fig.suptitle(Suptitle)
        fig.savefig(plotpath,dpi=300)
        print(plotpath)
        plt.close(fig)
    
    return


def plotpanelmagnetopause(): #draw the geopauses on 2 x 2 grid
    home=os.getcwd()+'/'
    reversal = ['runFlipNorth','runFlipSouth']
    GM = ['SwIono','SwHIonoO','SwIonoO','SwIonoO28amu']
    Title=['Sw H + Iono H\n Case 1','Sw H + Iono O\n Case 2','Sw H + Iono (H + O)\n Case 3','Sw H + Iono (H + O) at 28 amu/cc \n Case 4']
    minimum = 720
    plotnumber =1+12*60
    mode = '2d'
    cut = ['y','z']
    rc('font',family='serif')
    rc('text', usetex=True)
    xmin,xmax = -40.,10.
    ymin,ymax = -20.,20.
    #for i in range(0,4):
    fig,axs = plt.subplots(1,2)
  
    quantity='lastclosed'
    for j in range(0,2):
    #file location step
    #read data from file step
    #plotting step
    #save file
         for k in range(minimum,plotnumber):
             hr,mn=timestring(k)
             l1,l2,l3,l4= lastclosedcomponent(axs[j],j,k,reversal[j],quantity,cut[0],Title,xmin,xmax,ymin,ymax)
    plotpath = 'lastclosed.png'
    if (quantity == 'lastclosed'):
        Suptitle = 'Last Closed Field Line at'+' '+hr+':'+mn+' UT'
    leg=plt.figlegend((l1,l2,l3,l4),Title,(0.,0.09),ncol=4,fontsize=7,frameon=False)
    for t in leg.texts:
        t.set_multialignment('center')
    fig.suptitle(Suptitle,x=0.5,y=0.78)
    fig.savefig(plotpath,dpi=300,bbox_inches = 'tight',
                pad_inches = 0.1)
    plt.close(fig)
    
    return
def plotpanelall(): #draw the geopauses on 1 x 2 grid
    home=os.getcwd()+'/'
    reversal = ['runFlipNorth','runFlipSouth']
    GM = ['SwIono','SwHIonoO','SwIonoO','SwIonoO28amu']
    Title=['Sw H + Iono H','Sw H + Iono O','Sw H + Iono (H + O)','Sw H + Iono (H + O) at 28 amu/cc']
    Legends = ['Number Density Geopause','Mass Density Geopause','Pressure Geopause','Magnetopause']
    minimum = 720
    plotnumber =1+12*60
    mode = '2d'
    cut = ['y','z']
    rc('font',family='serif')
    rc('text', usetex=True)
    xmin,xmax = -60.,10.
    ymin,ymax = -20.,20.
    for i in range(0,4):
        fig,axs = plt.subplots(1,2)
        for j in range(0,2):
        #file location step
        #read data from file step
        #plotting step
        #save file
             for k in range(minimum,plotnumber):
                 hr,mn=timestring(k)
                 l1,l2,l3,l4=allcomponent(axs[j],j,k,GM[i],reversal[j],cut,Title,xmin,xmax,ymin,ymax)
        plotpath = 'all'+GM[i]+'.png'
        Suptitle = Title[i]+' Geopauses at'+' '+hr+':'+mn+' UT'
        plt.figlegend((l1,l2,l3,l4),Legends,(0.01,0.15),markerscale=3,ncol=4,fontsize=7,frameon=False)
        fig.suptitle(Suptitle,x=0.5,y=0.72)
        fig.savefig(plotpath,dpi=300,bbox_inches = 'tight',
                pad_inches = 0.1)
        plt.close(fig)
    
    return

def plotloop(i): #plot density geopause
    home=os.getcwd()+'/'
    reversal = ['runFlipNorth','runFlipSouth']
    GM = ['SwIono','SwHIonoO','SwIonoO','SwIonoO28amu']
    Title=['Sw H + Iono H','Sw H + Iono O','Sw H + Iono H + O','Sw H + Iono H + O at 28 amu/cc']
    fluidcontent = ['HH','HO','HOO','HOO']
    plotnumber =1+12*60
    mode = '2d'
    cut = 'z'
    zmin,zmax = -11.,11.
    rc('font',family='serif')
    rc('text', usetex=True)
    xmin,xmax = -40.,20.
    ymin,ymax = -20.,20.
    #for i in range(0,4):
    for j in range(0,2):
        run =GM[i]+'/'+reversal[j]+'/'
        datalocation = home+run
        plotlocation = home+'Plots/'+run
        if not os.path.exists(plotlocation):
            os.makedirs(plotlocation)
            print('making a plot directory')
        for k in range(0,plotnumber):
            filelist = []
            hr,mn=timestring(k)
            time = 't00'+hr+mn
            if not os.path.exists(datalocation):
                sys.exit('path does not exist')
            filelist.extend(glob.glob(datalocation+cut+'*'+time+'*'))
            print(i,hr, mn, GM[i], reversal[j])
            filepath = filelist[0]
            filename = PlotProperties(filepath,xmin,xmax,ymin,ymax,True,cut,mode)
            data,connectivity,timetick = filename.get_data()
            data = filename.data_filter(data)
            xyz,bfield,J,stateFluid,nFluid=filename.data_readable(data)
            N = 200
            xi = np.linspace(xyz[:,0].min(),xyz[:,0].max(),N)
            if (cut == 'z'):
                p = 1
            if (cut == 'y'):
                p = 2
            yi = np.linspace(xyz[:,p].min(),xyz[:,p].max(),N)
            if (GM[i]=='SwIono'):
                swi = stateFluid[:,0,0]
                iono= stateFluid[:,0,1]
                frac = swi/(iono+swi)
            if (GM[i]=='SwHIonoO'):
                swi = stateFluid[:,0,0]
                nO = stateFluid[:,0,1]/16.
                ionoO= nO
                frac = swi/(ionoO+swi)
            if (GM[i]=='SwIonoO' or GM[i] =='SwIonoO28amu'):
                swi = stateFluid[:,0,0]
                iono= stateFluid[:,0,1]
                nO = stateFluid[:,0,2]/16.
                ionoO= nO
                frac = swi/(iono+ionoO+swi)
            frac = sci.griddata((xyz[:,0],xyz[:,p]),frac,(xi[None,:],yi[:,None]),method='cubic')
            fig,ax = plt.subplots()
            cbar=ax.pcolormesh(xi,yi,frac,zorder=0,vmin=0.,vmax=1.,cmap = cm.get_cmap('seismic') )
            draw_earth(ax)
            ax.set_xlabel(r'X [$R_E$]')
            ax.set_ylabel(r'Y [$R_E$]')
            ax.set_xlim(xmin,xmax)
            ax.set_ylim(ymin,ymax)
            cb = fig.colorbar(cbar)
            cb.set_label(r'$n_{sw}/n$')
            ax.set_title(Title[i])
            fig.text(0.85,0.05,'T='+hr+':'+mn+':'+'00')
            plotpath = plotlocation+'abundance'+cut+'t'+hr+mn+'00'+'.png'
            fig.savefig(plotpath,dpi=300)
            plt.close(fig)
    return

def writegeopauseloop(i):
    home=os.getcwd()+'/'
    reversal = ['runFlipNorth','runFlipSouth']
    GM = ['SwIono','SwHIonoO','SwIonoO','SwIonoO28amu']
    Title=['Sw H + Iono H','Sw H + Iono O','Sw H + Iono H + O','Sw H + Iono H + O at 28 amu/cc']
    plotnumber =1+12*60
    mode = '2d'
    cut = 'z'
    eps = 0.01
    #for i in range(0,4):
    for j in range(0,2):
        run =GM[i]+'/'+reversal[j]+'/'
        datalocation = home+run
        plotlocation = home+'Data/Geopause/'+run
        if not os.path.exists(plotlocation):
            os.makedirs(plotlocation)
            print('making a plot directory')
        for k in range(0,plotnumber):
            filelist = []
            hr,mn=timestring(k)
            time = 't00'+hr+mn
            if not os.path.exists(datalocation):
                sys.exit('path does not exist')
            filelist.extend(glob.glob(datalocation+cut+'*'+time+'*'))
            print(hr, mn, GM[i], reversal[j])
            filepath = filelist[0]
            filename = PlotProperties(filepath,-40.0,20.0,-20.0,20.0,True,cut,mode)
            data,connectivity,timetick = filename.get_data()
            xyz,bfield,J,stateFluid,nFluid=filename.data_readable(data)
            xcoord,ycoord=[],[]
            xrho,yrho=[],[]
            xp,yp =[],[]
            xbeta,ybeta=[],[]
            xbetastar,ybetastar=[],[]
            xmp,ymp = [],[]
            bmag2 = bfield[:,0]**2+bfield[:,1]**2+bfield[:,2]**2
            bmag2 = bmag2*(1.e-9)**2
            mu0 = 4.*np.pi*1.e-7
            pmag = bmag2/(2.*mu0)
            mpro = 1.6726219e-27
            p = ordinate(cut)
            if (GM[i]=='SwIono'):
                swi = stateFluid[:,0,0]
                iono= stateFluid[:,0,1]
                densgeopause = swi/(iono+swi)
                rhogeopause = swi/(iono+swi)
                pramsw = stateFluid[:,0,0]*mpro*(1.e2)**3
                pramsw = pramsw*(stateFluid[:,1,0]**2+stateFluid[:,2,0]**2
                                 +stateFluid[:,3,0]**2)*(1.e3)**2
                pramio = stateFluid[:,0,1]*mpro*(1.e2)**3
                pramio = pramio*(stateFluid[:,1,1]**2+stateFluid[:,2,1]**2
                                 +stateFluid[:,3,1]**2)*(1.e3)**2
                pramtot = pramsw+pramio
                psw = stateFluid[:,4,0]
                pio = stateFluid[:,4,1]
                ptot = psw+pio
                pgeopause = psw/ptot
                beta = ptot*1.e-9/pmag
                betastar = beta+pramtot/pmag
                mp = pramsw/pmag
            if (GM[i]=='SwHIonoO'):
                swi = stateFluid[:,0,0]
                rhoO = stateFluid[:,0,1]
                nO= rhoO/16.
                densgeopause = swi/(nO+swi)
                rhogeopause = swi/(rhoO+swi)
                pramsw = stateFluid[:,0,0]*mpro*(1.e2)**3
                pramsw = pramsw*(stateFluid[:,1,0]**2+stateFluid[:,2,0]**2
                                 +stateFluid[:,3,0]**2)*(1.e3)**2
                pramO = stateFluid[:,0,1]*mpro*(1.e2)**3
                pramO= pramO*(stateFluid[:,1,1]**2+stateFluid[:,2,1]**2
                                 +stateFluid[:,3,1]**2)*(1.e3)**2
                pramtot = pramsw+pramO
                psw = stateFluid[:,4,0]
                pO = stateFluid[:,4,1]
                ptot = psw+pO
                pgeopause = psw/ptot
                beta = ptot*1.e-9/pmag
                betastar = beta+pramtot/pmag
                mp = pramsw/pmag
            if (GM[i]=='SwIonoO' or GM[i] =='SwIonoO28amu'):
                swi = stateFluid[:,0,0]
                iono= stateFluid[:,0,1]
                rhoO=stateFluid[:,0,2]
                nO = rhoO/16.
                densgeopause = swi/(iono+nO+swi)
                rhogeopause=swi/(iono+rhoO+swi)
                pramsw = stateFluid[:,0,0]*mpro*(1.e2)**3
                pramsw = pramsw*(stateFluid[:,1,0]**2+stateFluid[:,2,0]**2
                                 +stateFluid[:,3,0]**2)*(1.e3)**2
                pramio = stateFluid[:,0,1]*mpro*(1.e2)**3
                pramio = pramio*(stateFluid[:,1,1]**2+stateFluid[:,2,1]**2
                                 +stateFluid[:,3,1]**2)*(1.e3)**2
                pramO = stateFluid[:,0,2]*mpro*(1.e2)**3
                pramO= pramO*(stateFluid[:,1,2]**2+stateFluid[:,2,2]**2
                                 +stateFluid[:,3,2]**2)*(1.e3)**2
                pramtot = pramsw+pramio+pramO
                psw = stateFluid[:,4,0]
                pio = stateFluid[:,4,1]
                pO = stateFluid[:,4,2]
                ptot = psw+pio+pO
                pgeopause = psw/ptot
                beta = ptot*1.e-9/pmag
                betastar = beta+pramtot/pmag
                mp = pramsw/pmag
            row,column = connectivity.shape
            eps = 0.5
            one = 1.
            '''subtract the connectivity
            index by 1 due to how python indexing works'''
            '''
            slope formula: y = (y1-y0)/(x1-x0)*(x-x0)+y0
            '''
            for m in np.arange(0,row):  
                p0 = connectivity[m,0]-1 
                p1 = connectivity[m,1]-1
                p2 = connectivity[m,2]-1
                p3 = connectivity[m,3]-1
                if((densgeopause[p0]-eps)*(densgeopause[p1]-eps)<0):
                    xtemp = eps*(xyz[p1,0]-xyz[p0,0])+densgeopause[p1]*xyz[p0,0]-densgeopause[p0]*xyz[p1,0]
                    xtemp = xtemp/(densgeopause[p1]-densgeopause[p0])
                    xcoord.append(xtemp)
                    ycoord.append(xyz[p0,p])
                if((densgeopause[p1]-eps)*(densgeopause[p2]-eps)<0):
                    ytemp = eps*(xyz[p2,p]-xyz[p1,p])+densgeopause[p2]*xyz[p1,p]-densgeopause[p1]*xyz[p2,p]
                    ytemp = ytemp/(densgeopause[p2]-densgeopause[p1])
                    ycoord.append(ytemp)
                    xcoord.append(xyz[p1,0])
                if((densgeopause[p3]-eps)*(densgeopause[p2]-eps)<0):
                    xtemp = eps*(xyz[p2,0]-xyz[p3,0])+densgeopause[p2]*xyz[p3,0]-densgeopause[p3]*xyz[p2,0]
                    xtemp = xtemp/(densgeopause[p2]-densgeopause[p3])
                    xcoord.append(xtemp)
                    ycoord.append(xyz[p2,p])
                if((densgeopause[p0]-eps)*(densgeopause[p3]-eps)<0):
                    ytemp = eps*(xyz[p3,p]-xyz[p0,p])+densgeopause[p3]*xyz[p0,p]-densgeopause[p0]*xyz[p3,p]
                    ytemp = ytemp/(densgeopause[p3]-densgeopause[p0])
                    ycoord.append(ytemp)
                    xcoord.append(xyz[p0,0])
                    
                if((rhogeopause[p0]-eps)*(rhogeopause[p1]-eps)<0):
                    xtemp = eps*(xyz[p1,0]-xyz[p0,0])+rhogeopause[p1]*xyz[p0,0]-rhogeopause[p0]*xyz[p1,0]
                    xtemp = xtemp/(rhogeopause[p1]-rhogeopause[p0])
                    xrho.append(xtemp)
                    yrho.append(xyz[p0,p])
                if((rhogeopause[p1]-eps)*(rhogeopause[p2]-eps)<0):
                    ytemp = eps*(xyz[p2,p]-xyz[p1,p])+rhogeopause[p2]*xyz[p1,p]-rhogeopause[p1]*xyz[p2,p]
                    ytemp = ytemp/(rhogeopause[p2]-rhogeopause[p1])
                    yrho.append(ytemp)
                    xrho.append(xyz[p1,0])
                if((rhogeopause[p3]-eps)*(rhogeopause[p2]-eps)<0):
                    xtemp = eps*(xyz[p2,0]-xyz[p3,0])+rhogeopause[p2]*xyz[p3,0]-rhogeopause[p3]*xyz[p2,0]
                    xtemp = xtemp/(rhogeopause[p2]-rhogeopause[p3])
                    xrho.append(xtemp)
                    yrho.append(xyz[p2,p])
                if((rhogeopause[p0]-eps)*(rhogeopause[p3]-eps)<0):
                    ytemp = eps*(xyz[p3,p]-xyz[p0,p])+rhogeopause[p3]*xyz[p0,p]-rhogeopause[p0]*xyz[p3,p]
                    ytemp = ytemp/(rhogeopause[p3]-rhogeopause[p0])
                    yrho.append(ytemp)
                    xrho.append(xyz[p0,0])

                if((pgeopause[p0]-eps)*(pgeopause[p1]-eps)<0):
                    xtemp = eps*(xyz[p1,0]-xyz[p0,0])+pgeopause[p1]*xyz[p0,0]-pgeopause[p0]*xyz[p1,0]
                    xtemp = xtemp/(pgeopause[p1]-pgeopause[p0])
                    xp.append(xtemp)
                    yp.append(xyz[p0,p])
                if((pgeopause[p1]-eps)*(pgeopause[p2]-eps)<0):
                    ytemp = eps*(xyz[p2,p]-xyz[p1,p])+pgeopause[p2]*xyz[p1,p]-pgeopause[p1]*xyz[p2,p]
                    ytemp = ytemp/(pgeopause[p2]-pgeopause[p1])
                    yp.append(ytemp)
                    xp.append(xyz[p1,0])
                if((pgeopause[p3]-eps)*(pgeopause[p2]-eps)<0):
                    xtemp = eps*(xyz[p2,0]-xyz[p3,0])+pgeopause[p2]*xyz[p3,0]-pgeopause[p3]*xyz[p2,0]
                    xtemp = xtemp/(pgeopause[p2]-pgeopause[p3])
                    xp.append(xtemp)
                    yp.append(xyz[p2,p])
                if((pgeopause[p0]-eps)*(pgeopause[p3]-eps)<0):
                    ytemp = eps*(xyz[p3,p]-xyz[p0,p])+pgeopause[p3]*xyz[p0,p]-pgeopause[p0]*xyz[p3,p]
                    ytemp = ytemp/(pgeopause[p3]-pgeopause[p0])
                    yp.append(ytemp)
                    xp.append(xyz[p0,0])

                if((beta[p0]-one)*(beta[p1]-one)<0):
                    xtemp = one*(xyz[p1,0]-xyz[p0,0])+beta[p1]*xyz[p0,0]-beta[p0]*xyz[p1,0]
                    xtemp = xtemp/(beta[p1]-beta[p0])
                    xbeta.append(xtemp)
                    ybeta.append(xyz[p0,p])
                if((beta[p1]-one)*(beta[p2]-one)<0):
                    ytemp = one*(xyz[p2,p]-xyz[p1,p])+beta[p2]*xyz[p1,p]-beta[p1]*xyz[p2,p]
                    ytemp = ytemp/(beta[p2]-beta[p1])
                    ybeta.append(ytemp)
                    xbeta.append(xyz[p1,0])
                if((beta[p3]-one)*(beta[p2]-one)<0):
                    xtemp = one*(xyz[p2,0]-xyz[p3,0])+beta[p2]*xyz[p3,0]-beta[p3]*xyz[p2,0]
                    xtemp = xtemp/(beta[p2]-beta[p3])
                    xbeta.append(xtemp)
                    ybeta.append(xyz[p2,p])
                if((beta[p0]-one)*(beta[p3]-one)<0):
                    ytemp = one*(xyz[p3,p]-xyz[p0,p])+beta[p3]*xyz[p0,p]-beta[p0]*xyz[p3,p]
                    ytemp = ytemp/(beta[p3]-beta[p0])
                    ybeta.append(ytemp)
                    xbeta.append(xyz[p0,0])
                    
                if((betastar[p0]-one)*(betastar[p1]-one)<0):
                    xtemp = one*(xyz[p1,0]-xyz[p0,0])+betastar[p1]*xyz[p0,0]-betastar[p0]*xyz[p1,0]
                    xtemp = xtemp/(betastar[p1]-betastar[p0])
                    xbetastar.append(xtemp)
                    ybetastar.append(xyz[p0,p])
                if((betastar[p1]-one)*(betastar[p2]-one)<0):
                    ytemp = one*(xyz[p2,p]-xyz[p1,p])+betastar[p2]*xyz[p1,p]-betastar[p1]*xyz[p2,p]
                    ytemp = ytemp/(betastar[p2]-betastar[p1])
                    ybetastar.append(ytemp)
                    xbetastar.append(xyz[p1,0])
                if((betastar[p3]-one)*(betastar[p2]-one)<0):
                    xtemp = one*(xyz[p2,0]-xyz[p3,0])+betastar[p2]*xyz[p3,0]-betastar[p3]*xyz[p2,0]
                    xtemp = xtemp/(betastar[p2]-betastar[p3])
                    xbetastar.append(xtemp)
                    ybetastar.append(xyz[p2,p])
                if((betastar[p0]-one)*(betastar[p3]-one)<0):
                    ytemp = one*(xyz[p3,p]-xyz[p0,p])+betastar[p3]*xyz[p0,p]-betastar[p0]*xyz[p3,p]
                    ytemp = ytemp/(betastar[p3]-betastar[p0])
                    ybetastar.append(ytemp)
                    xbetastar.append(xyz[p0,0])

                if((mp[p0]-one)*(mp[p1]-one)<0):
                    xtemp = one*(xyz[p1,0]-xyz[p0,0])+mp[p1]*xyz[p0,0]-mp[p0]*xyz[p1,0]
                    xtemp = xtemp/(mp[p1]-mp[p0])
                    xmp.append(xtemp)
                    ymp.append(xyz[p0,p])
                if((mp[p1]-one)*(mp[p2]-one)<0):
                    ytemp = one*(xyz[p2,p]-xyz[p1,p])+mp[p2]*xyz[p1,p]-mp[p1]*xyz[p2,p]
                    ytemp = ytemp/(mp[p2]-mp[p1])
                    ymp.append(ytemp)
                    xmp.append(xyz[p1,0])
                if((mp[p3]-one)*(mp[p2]-one)<0):
                    xtemp = one*(xyz[p2,0]-xyz[p3,0])+mp[p2]*xyz[p3,0]-mp[p3]*xyz[p2,0]
                    xtemp = xtemp/(mp[p2]-mp[p3])
                    xmp.append(xtemp)
                    ymp.append(xyz[p2,p])
                if((mp[p0]-one)*(mp[p3]-one)<0):
                    ytemp = one*(xyz[p3,p]-xyz[p0,p])+mp[p3]*xyz[p0,p]-mp[p0]*xyz[p3,p]
                    ytemp = ytemp/(mp[p3]-mp[p0])
                    ymp.append(ytemp)
                    xmp.append(xyz[p0,0])


                    
            fname1 = plotlocation+'geopause'+cut+'t'+hr+mn+'00'+'.txt'
            fname2 = plotlocation+'rho'+cut+'t'+hr+mn+'00'+'.txt'
            fname3 = plotlocation+'pressure'+cut+'t'+hr+mn+'00'+'.txt'
            fname4 = plotlocation+'beta'+cut+'t'+hr+mn+'00'+'.txt'
            fname5 = plotlocation+'betastar'+cut+'t'+hr+mn+'00'+'.txt'
            fname6 = plotlocation+'magnetopause'+cut+'t'+hr+mn+'00'+'.txt'
            np.savetxt(fname1,np.c_[xcoord,ycoord],header='density',fmt='%.6e')
            np.savetxt(fname2,np.c_[xrho,yrho],header='rho',fmt='%.6e')
            np.savetxt(fname3,np.c_[xp,yp],header='pressure',fmt='%.6e')
            np.savetxt(fname4,np.c_[xbeta,ybeta],header='beta',fmt='%.6e')
            np.savetxt(fname5,np.c_[xbetastar,ybetastar],header='betastar',fmt='%.6e')
            np.savetxt(fname6,np.c_[xmp,ymp],header='magnetopause',fmt='%.6e')
    return


#plots all the different boundaries
def plotall(i):
    home=os.getcwd()+'/'
    reversal = ['runFlipNorth','runFlipSouth']
    GM = ['SwIono','SwHIonoO','SwIonoO','SwIonoO28amu']
    Title=['Sw H + Iono H','Sw H + Iono O','Sw H + Iono H + O','Sw H + Iono H + O at 28 amu/cc']
    Legends = [r'n',r'$\rho$',r'p',r'$\beta$',r'$\beta$*','magnetopause']
    minimum = 0
    plotnumber =1+12*60
    mode = '2d'
    cut = 'y'
    xmin,xmax = -80,20
    zmin,zmax = -20.,20.
    rc('font',family='serif')
    rc('text', usetex=True)
    for j in range(0,2):
        root=home+'Data/Geopause/'+GM[i]+'/'+reversal[j]+'/'
        plotlocation = home+'Plots/'+'combined/'+GM[i]+'/'+reversal[j]+'/'
        if not os.path.exists(plotlocation):
            os.makedirs(plotlocation)
            print('making a plot directory')
        for k in range(minimum,plotnumber):
            file1,file2,file3,file4,file5,file6 = [],[],[],[],[],[]
            hr,mn=timestring(k)
            time = 't'+hr+mn
            if not os.path.exists(plotlocation):
                sys.exit('path does not exist')
            file1.extend(glob.glob(root+'geopause'+cut+'*'+time+'*'))
            file2.extend(glob.glob(root+'rho'+cut+'*'+time+'*'))
            file3.extend(glob.glob(root+'pressure'+cut+'*'+time+'*'))
            file4.extend(glob.glob(root+'beta'+cut+'*'+time+'*'))
            file5.extend(glob.glob(root+'betastar'+cut+'*'+time+'*'))
            file6.extend(glob.glob(root+'magnetopause'+cut+'*'+time+'*'))
            print(time,hr, mn, reversal[j], root)
            filepath = file1[0]
            print(filepath)
            with open(filepath,'r') as infile:
                x1,y1=readgeopause(infile)
            filepath = file2[0]
            with open(filepath,'r') as infile:
                x2,y2=readgeopause(infile)
            filepath = file3[0]
            with open(filepath,'r') as infile:
                x3,y3=readgeopause(infile)
            filepath = file4[0]
            with open(filepath,'r') as infile:
                x4,y4=readgeopause(infile)
            filepath = file5[0]
            with open(filepath,'r') as infile:
                x5,y5=readgeopause(infile)
            filepath = file6[0]
            with open(filepath,'r') as infile:
                x6,y6=readgeopause(infile)
            fig,ax = plt.subplots()
            ax.scatter(x1,y1,s=1,color='k',label=Legends[0])
            ax.scatter(x2,y2,s=1,color='b',label=Legends[1])
            ax.scatter(x3,y3,s=1,color='r',label=Legends[2])
            #ax.scatter(x4,y4,s=1,color='g',label=Legends[3])
            #ax.scatter(x5,y5,s=1,color='m',label=Legends[4])
            #ax.scatter(x6,y6,s=1,color='gray',label=Legends[5])
            ax.set_xlim(xmin,xmax)
            ax.set_ylim(zmin,zmax)
            draw_earth(ax)
            ax.set_xlabel(r'X [$R_E$]')
            ax.set_ylabel(r'Z [$R_E$]')
            ax.set_title(Title[0])
            ax.legend(loc='lower left')
            fig.text(0.85,0.03,'T='+hr+':'+mn+':'+'00')
            plotpath = plotlocation+'all'+cut+'t'+hr+mn+'00'+'.png'
            fig.savefig(plotpath,dpi=300)
            plt.close(fig)
    return

#plots geopause for rho, density, pressure, and last closed field line
def plotpaper(i):
    home=os.getcwd()+'/'
    reversal = ['runFlipNorth','runFlipSouth']
    GM = ['SwIono','SwHIonoO','SwIonoO','SwIonoO28amu']
    Title=['Sw H + Iono H','Sw H + Iono O','Sw H + Iono H + O','Sw H + Iono H + O at 28 amu/cc']
    Legends = [r'n',r'$\rho$',r'p',r'Magnetopause']
    minimum = 720
    plotnumber =1+12*60 #
    mode = '2d'
    cut = 'y'
    xmin,xmax = -140,20
    zmin,zmax = -20.,20.
    rc('font',family='serif')
    rc('text', usetex=True)
    for j in range(0,2):
        root=home+'Data/Geopause/'+GM[i]+'/'+reversal[j]+'/'
        plotlocation = home+'Plots/'+'paper/'+GM[i]+'/'+reversal[j]+'/'
        if not os.path.exists(plotlocation):
            os.makedirs(plotlocation)
            print('making a plot directory')
        for k in range(minimum,plotnumber):
            file1,file2,file3,file4 = [],[],[],[]
            hr,mn=timestring(k)
            time = 't'+hr+mn
            if not os.path.exists(plotlocation):
                sys.exit('path does not exist')
            file1.extend(glob.glob(root+'geopause'+cut+'*'+time+'*'))
            file2.extend(glob.glob(root+'rho'+cut+'*'+time+'*'))
            file3.extend(glob.glob(root+'pressure'+cut+'*'+time+'*'))
            file4.extend(glob.glob(root+'lastclosed'+cut+'*'+time+'*'))
            print(time,hr, mn, reversal[j], root)
            filepath = file1[0]
            print(filepath)
            with open(filepath,'r') as infile:
                x1,y1=readgeopause(infile)
            filepath = file2[0]
            with open(filepath,'r') as infile:
                x2,y2=readgeopause(infile)
            filepath = file3[0]
            with open(filepath,'r') as infile:
                x3,y3=readgeopause(infile)
            filepath = file4[0]
            with open(filepath,'r') as infile:
                xnight,ynight,xday,yday=readlastclosed(infile)
            fig,ax = plt.subplots()
            ax.scatter(x1,y1,s=1,color='k',label=Legends[0])
            ax.scatter(x2,y2,s=1,color='b',label=Legends[1])
            ax.scatter(x3,y3,s=1,color='r',label=Legends[2])
            ax.plot(xnight,ynight,color='g',label=Legends[3],zorder=0)
            ax.plot(xday,yday,color='g',zorder=0)
            #ax.scatter(x5,y5,s=1,color='m',label=Legends[4])
            #ax.scatter(x6,y6,s=1,color='gray',label=Legends[5])
            ax.set_xlim(xmin,xmax)
            ax.set_ylim(zmin,zmax)
            draw_earth(ax)
            ax.set_xlabel(r'X [$R_E$]')
            ax.set_ylabel(r'Z [$R_E$]')
            ax.set_title(Title[i])
            ax.legend(loc='lower left')
            fig.text(0.85,0.03,'T='+hr+':'+mn+':'+'00')
            ax.set_aspect('equal')
            plotpath = plotlocation+'all'+cut+'t'+hr+mn+'00'+'.png'
            fig.savefig(plotpath,dpi=300)
            plt.close(fig)
    return


#plotgeopause draws the geopause in the meridional or equatorial plane
def plotgeopause(i):
    home=os.getcwd()+'/'
    reversal = ['runFlipNorth','runFlipSouth']
    GM = ['SwIono','SwHIonoO','SwIonoO','SwIonoO28amu']
    Title=['Sw H + Iono H','Sw H + Iono O','Sw H + Iono H + O','Sw H + Iono H + O at 28 amu/cc']
    minimum = 720
    plotnumber =1+12*60
    mode = '2d'
    cut = 'z'
    quantity = 'lastclosed'
    xmin,xmax = -140,20
    zmin,zmax = -20.,20.
    rc('font',family='serif')
    rc('text', usetex=True)
    root1=home+'Data/Geopause/'+'SwIono/'+reversal[i]+'/'
    root2=home+'Data/Geopause/'+'SwHIonoO/'+reversal[i]+'/'
    root3=home+'Data/Geopause/'+'SwIonoO/'+reversal[i]+'/'
    root4=home+'Data/Geopause/'+'SwIonoO28amu/'+reversal[i]+'/'
    plotlocation = home+'paper/'+quantity+'/'+reversal[i]+'/'
    print(root1)
    if not os.path.exists(plotlocation):
        os.makedirs(plotlocation)
        print('making a plot directory')
    for k in range(minimum,plotnumber):
        file1 = []
        file2 = []
        file3 = []
        file4 = []
        hr = k/60
        if (hr < 10):
            hr = '0'+str(hr)
        else:
            hr = str(hr)
        mn = k%60
        if (mn < 10):
            mn = '0'+str(mn)
        else:
            mn = str(mn)
        time = 't'+hr+mn
        if not os.path.exists(plotlocation):
            sys.exit('path does not exist')
        file1.extend(glob.glob(root1+quantity+cut+'*'+time+'*'))
        file2.extend(glob.glob(root2+quantity+cut+'*'+time+'*'))
        file3.extend(glob.glob(root3+quantity+cut+'*'+time+'*'))
        file4.extend(glob.glob(root4+quantity+cut+'*'+time+'*'))
        print(time,hr, mn, reversal[i])
        filepath = file1[0]
        print(filepath)
        with open(filepath,'r') as infile:
            x1,y1=readgeopause(infile)
        filepath = file2[0]
        with open(filepath,'r') as infile:
            x2,y2=readgeopause(infile)
        filepath = file3[0]
        with open(filepath,'r') as infile:
            x3,y3=readgeopause(infile)
        filepath = file4[0]
        with open(filepath,'r') as infile:
            x4,y4=readgeopause(infile)
        fig,ax = plt.subplots()
        ax.scatter(x1,y1,s=1,color='k',label=Title[0])
        ax.scatter(x2,y2,s=1,color='b',label=Title[1])
        ax.scatter(x3,y3,s=1,color='r',label=Title[2])
        ax.scatter(x4,y4,s=1,color='g',label=Title[3])
        ax.set_xlim(xmin,xmax)
        ax.set_ylim(zmin,zmax)
        draw_earth(ax)
        ax.set_xlabel(r'X [$R_E$]')
        ax.set_ylabel(r'Z [$R_E$]')
        ax.set_title(r'Last Closed Field Line')
        ax.legend(loc='lower left')
        fig.text(0.85,0.03,'T='+hr+':'+mn)
        plotpath = plotlocation+quantity+cut+'t'+hr+mn+'00'+'.png'
        fig.savefig(plotpath,dpi=300)
        plt.close(fig)
    return

def plotlastclosed(i):
    '''
    plotlastclosed is a function that draws the last closed field line in the meridional plane
    
    '''
    home=os.getcwd()+'/'
    reversal = ['runFlipNorth','runFlipSouth']
    GM = ['SwIono','SwHIonoO','SwIonoO','SwIonoO28amu']
    Title=['Sw H + Iono H','Sw H + Iono O','Sw H + Iono H + O','Sw H + Iono H + O at 28 amu/cc']
    minimum = 720
    plotnumber =1+12*60
    mode = '2d'
    cut = 'y'
    quantity = 'lastclosed'
    xmin,xmax = -40,10
    zmin,zmax = -20.,20.
    rc('font',family='serif')
    rc('text', usetex=True)
    root1=home+'Data/Geopause/'+'SwIono/'+reversal[i]+'/'
    root2=home+'Data/Geopause/'+'SwHIonoO/'+reversal[i]+'/'
    root3=home+'Data/Geopause/'+'SwIonoO/'+reversal[i]+'/'
    root4=home+'Data/Geopause/'+'SwIonoO28amu/'+reversal[i]+'/'
    plotlocation = home+'Plots/paper/'+quantity+'/'+reversal[i]+'/'
    print(plotlocation)
    if not os.path.exists(plotlocation):
        os.makedirs(plotlocation)
        print('making a plot directory')
    for k in range(minimum,plotnumber):
        file1 = []
        file2 = []
        file3 = []
        file4 = []
        hr,mn=timestring(k)
        time = 't'+hr+mn
        if not os.path.exists(plotlocation):
            sys.exit('path does not exist')
        file1.extend(glob.glob(root1+quantity+cut+'*'+time+'*'))
        file2.extend(glob.glob(root2+quantity+cut+'*'+time+'*'))
        file3.extend(glob.glob(root3+quantity+cut+'*'+time+'*'))
        file4.extend(glob.glob(root4+quantity+cut+'*'+time+'*'))
        print(time,hr, mn, reversal[i])
        filepath = file1[0]
        print(filepath)
        with open(filepath,'r') as infile:
            xnight1,ynight1,xday1,yday1=readlastclosed(infile)
        filepath = file2[0]
        with open(filepath,'r') as infile:
            xnight2,ynight2,xday2,yday2=readlastclosed(infile)
        filepath = file3[0]
        with open(filepath,'r') as infile:
            xnight3,ynight3,xday3,yday3=readlastclosed(infile)
        filepath = file4[0]
        with open(filepath,'r') as infile:
            xnight4,ynight4,xday4,yday4=readlastclosed(infile)
        fig,ax = plt.subplots()
        draw_earth(ax)
        ax.plot(xnight1,ynight1,color='k',label=Title[0],zorder=0)
        ax.plot(xnight2,ynight2,color='b',label=Title[1],zorder=0)
        ax.plot(xnight3,ynight3,color='r',label=Title[2],zorder=0)
        ax.plot(xnight4,ynight4,color='g',label=Title[3],zorder=0)
        ax.plot(xday1,yday1,color='k',zorder=0)
        ax.plot(xday2,yday2,color='b',zorder=0)
        ax.plot(xday3,yday3,color='r',zorder=0)
        ax.plot(xday4,yday4,color='g',zorder=0)
        ax.set_xlim(xmin,xmax)
        ax.set_ylim(zmin,zmax)
        
        ax.set_xlabel(r'X [$R_E$]')
        ax.set_ylabel(r'Z [$R_E$]')
        ax.set_title(r'Last Closed Field Line')
        ax.legend(loc='lower left')
     
        plotpath = plotlocation+quantity+cut+'t'+hr+mn+'00'+'.png'
        ax.set_aspect('equal')
        fig.savefig(plotpath,dpi=300)
        plt.close(fig)
    return



def plotlinegeopause(i):
    home=os.getcwd()+'/'
    Geopause = ['Density Geopause','Mass Density Geopause','Pressure Geopause']
    reversal = ['runFlipNorth','runFlipSouth']
    IMF = ['South-to-North IMF','North-to-South IMF']
    GM = ['SwIono','SwHIonoO','SwIonoO','SwIonoO28amu']
    Title=['Sw H + Iono H','Sw H + Iono O','Sw H + Iono H + O','Sw H + Iono H + O at 28 amu/cc']
    plotnumber =1+12*60
    mode = '2d'
    cut = 'y'
    quantity = ['geopause','rho','pressure']
    xmin,xmax = -80,20
    zmin,zmax = -20.,20.
    rc('font',family='serif')
    rc('text', usetex=True)
    root1=home+'Data/Geopause/'+'SwIono/'+reversal[i]+'/'
    root2=home+'Data/Geopause/'+'SwHIonoO/'+reversal[i]+'/'
    root3=home+'Data/Geopause/'+'SwIonoO/'+reversal[i]+'/'
    root4=home+'Data/Geopause/'+'SwIonoO28amu/'+reversal[i]+'/'
    plotlocation = home+'Plots'+'sunearth'+'/'+reversal[i]+'/'
   
    timeflip = verticaldraw(reversal[i])
    if (reversal[i] == 'runFlipNorth'):
        tmin = '02:00'
    if (reversal[i] == 'runFlipSouth'):
        tmin = '06:00'
    if not os.path.exists(plotlocation):
        os.makedirs(plotlocation)
        print('making a plot directory')
    for j in range(0,len(quantity)):
        time1,time2,time3,time4 =[],[],[],[]
        xday1,xnight1,yday1,ynight1=[],[],[],[]
        xday2,xnight2,yday2,ynight2=[],[],[],[]
        xday3,xnight3,yday3,ynight3=[],[],[],[]
        xday4,xnight4,yday4,ynight4=[],[],[],[]
        for k in range(0,plotnumber):
            file1,file2,file3,file4 = [],[],[],[]
            hr,mn=timestring(k)
            time = 't'+hr+mn
            if not os.path.exists(plotlocation):
                sys.exit('path does not exist')
            file1.extend(glob.glob(root1+quantity[j]+cut+'*'+time+'*'))
            file2.extend(glob.glob(root2+quantity[j]+cut+'*'+time+'*'))
            file3.extend(glob.glob(root3+quantity[j]+cut+'*'+time+'*'))
            file4.extend(glob.glob(root4+quantity[j]+cut+'*'+time+'*'))
            print('Reading file at time:', hr+':'+mn, reversal[i])
            filepath = file1[0]
            with open(filepath,'r') as infile:
                x1,y1=readgeopause(infile)
            filepath = file2[0]
            with open(filepath,'r') as infile:
                x2,y2=readgeopause(infile)
            filepath = file3[0]
            with open(filepath,'r') as infile:
                x3,y3=readgeopause(infile)
            filepath = file4[0]
            with open(filepath,'r') as infile:
                x4,y4=readgeopause(infile)
            xd1,yd1,xn1,yn1=geopause_daynight(x1,y1)
            xd2,yd2,xn2,yn2=geopause_daynight(x2,y2)
            xd3,yd3,xn3,yn3=geopause_daynight(x3,y3)
            xd4,yd4,xn4,yn4=geopause_daynight(x4,y4)

            tvert = hr+":"+mn
            tdraw =  md.date2num(dt.datetime.strptime(tvert,'%H:%M'))
            xtemp1,xtemp2=geopause_daynight_sunearth(xd1,yd1,xn1,yn1)
            xday1.append(xtemp1)
            xnight1.append(xtemp2)
            xtemp1,xtemp2=geopause_daynight_sunearth(xd2,yd2,xn2,yn2)
            xday2.append(xtemp1)
            xnight2.append(xtemp2)
            xtemp1,xtemp2=geopause_daynight_sunearth(xd3,yd3,xn3,yn3)
            xday3.append(xtemp1)
            xnight3.append(xtemp2)
            xtemp1,xtemp2=geopause_daynight_sunearth(xd4,yd4,xn4,yn4)
            xday4.append(xtemp1)
            xnight4.append(xtemp2)
            time1.append(tdraw)
        time2 = time1
        time3 = time1
        time4 = time1
        fig,ax = plt.subplots()
        ax.plot(time1,xday1,color='k',label=Title[0])
        ax.plot(time2,xday2,color='b',label=Title[1])
        ax.plot(time3,xday3,color='r',label=Title[2])
        ax.plot(time4,xday4,color='g',label=Title[3])
        xfmt,timeticks=satellite_format(reversal[i])
        tdraw =  md.date2num(dt.datetime.strptime(timeflip,'%H:%M'))
        ax.axvline(tdraw,color='k',linestyle = 'dashed')
        ax.set_xticks(timeticks)
        ax.set_xlim([md.date2num(dt.datetime.strptime(tmin,'%H:%M')),md.date2num(dt.datetime.strptime('12:00','%H:%M'))])
        ax.xaxis.set_major_formatter(xfmt)
        #ax.set_ylim(3,10)
        ax.set_ylabel(r'X [$R_E$]')
        ax.set_title('Dayside '+Geopause[j]+' for '+IMF[i])
        ax.legend()
        ax.grid()
        plotpath = plotlocation+quantity[j]+'linedayside'+reversal[i]+'.png'
        fig.savefig(plotpath,dpi=300)
        plt.close(fig)
        fig,ax = plt.subplots()
        ax.plot(time1,xnight1,color='k',label=Title[0])
        ax.plot(time2,xnight2,color='b',label=Title[1])
        ax.plot(time3,xnight3,color='r',label=Title[2])
        ax.plot(time4,xnight4,color='g',label=Title[3])
        xfmt,timeticks=satellite_format(reversal[i])
        tdraw =  md.date2num(dt.datetime.strptime(timeflip,'%H:%M'))
        ax.axvline(tdraw,color='k',linestyle = 'dashed')
        ax.set_xticks(timeticks)
        ax.set_xlim([md.date2num(dt.datetime.strptime(tmin,'%H:%M')),md.date2num(dt.datetime.strptime('12:00','%H:%M'))])
        ax.xaxis.set_major_formatter(xfmt)
        #ax.set_ylim(3,10)
        ax.set_xlabel(r'Time [UT]')
        ax.set_ylabel(r'X [$R_E$]')
        ax.set_title('Nightside '+Geopause[j]+' for '+IMF[i])
        ax.legend()
        ax.grid()
        plotpath = plotlocation+quantity[j]+'linenightside'+reversal[i]+'.png'
        print(plotpath)
        fig.savefig(plotpath,dpi=300)
        plt.close(fig)
    return

def writelinegeopause(i):
    home=os.getcwd()+'/'
    Geopause = ['Density Geopause','Mass Density Geopause','Pressure Geopause']
    reversal = ['runFlipNorth','runFlipSouth']
    IMF = ['South-to-North IMF','North-to-South IMF']
    GM = ['SwIono','SwHIonoO','SwIonoO','SwIonoO28amu']
    Title=['Sw H + Iono H','Sw H + Iono O','Sw H + Iono H + O','Sw H + Iono H + O at 28 amu/cc']
    plotnumber =1+12*60
    mode = '2d'
    cut = 'y'
    quantity = ['geopause','rho','pressure']
    xmin,xmax = -80,20
    zmin,zmax = -20.,20.
    rc('font',family='serif')
    rc('text', usetex=True)
    root=[home+'Data/Geopause/'+'SwIono/'+reversal[i]+'/', \
          home+'Data/Geopause/'+'SwHIonoO/'+reversal[i]+'/', \
          home+'Data/Geopause/'+'SwIonoO/'+reversal[i]+'/', \
          home+'Data/Geopause/'+'SwIonoO28amu/'+reversal[i]+'/']
    datalocation = home+'Data/LineGeopause/'
    if not os.path.exists(datalocation):
        os.makedirs(datalocation)
        print('making a plot directory')
    for j in range(0,len(quantity)):
        if not os.path.exists(datalocation):
            sys.exit('path does not exist')
        for k in range(0,len(GM)):
            filelocation = datalocation+GM[k]
            fname1 = datalocation+GM[k]+'/'+'linedayside'+quantity[j]+reversal[i]+'.txt'
            fname2 = datalocation+GM[k]+'/'+'linenightside'+quantity[j]+reversal[i]+'.txt'
            if not os.path.exists(filelocation):
                os.makedirs(filelocation)
            print(root[k])
            writelinegeopauseloop(plotnumber,root[k],quantity[j],cut,fname1,fname2)
    return

def writelinegeopauseloop(plotnumber,root,quantity,cut,filenameday,filenamenight):
    xday,xnight = [],[]
    time1 = []
    for i in range(0,plotnumber):
        inputfile = []
        hr,mn=timestring(i)
        time = 't'+hr+mn
        inputfile.extend(glob.glob(root+quantity+cut+'*'+time+'*'))
        filepath = inputfile[0]
        with open(filepath,'r') as infile:
            x,y=readgeopause(infile)
        xd,yd,xn,yn=geopause_daynight(x,y)
        tvert = hr+":"+mn
        xtemp1,xtemp2=geopause_daynight_sunearth(xd,yd,xn,yn)
        xday.append(xtemp1)
        xnight.append(xtemp2)
        time1.append(tvert)
    np.savetxt(filenameday,np.c_[time1,xday],header='time, x', fmt="%s")
    np.savetxt(filenamenight,np.c_[time1,xnight],header='time, x', fmt="%s")
    return

def readlinegeopause(filename):
    a,Xgeo=np.loadtxt(filename,dtype = {'names':('time','X'),'formats':('S6','f')},unpack=True)
    tdraw =  [md.date2num(dt.datetime.strptime(date,'%H:%M')) for date in a]
    return(tdraw,Xgeo)

def plotlinegeopausecomponent(ax,filename,colors,legends):
    tgeo,Xgeo = readlinegeopause(filename)
    l1, = ax.plot(tgeo,Xgeo,color=colors,label=legends)
    return(l1)

def plotlinegeopauseaxes(ax,IMF,j,variable):
    '''
    plotlinegeopause labels the axes, draws the grids,
    writes the title
    '''
    xfmt,timeticks=satellite_format(IMF)
    timeflip = verticaldraw(IMF)
    tdraw =  md.date2num(dt.datetime.strptime(timeflip,'%H:%M'))
    if (IMF == 'runFlipNorth'):
        tmin = '02:00'
        tmax = '12:00'
        solarwind = 'South-to-North IMF'
    if (IMF == 'runFlipSouth'):
        tmin = '06:00'
        tmax = '12:00'
        solarwind = 'North-to-South IMF'
    if (variable == 'geopause'):
        ymin,ymax = -75,-2
        Title = 'Number Density Geopause'
    if (variable == 'rho'):
        ymin,ymax = -170,-2
        Title = 'Mass Density Geopause'
    if (variable == 'pressure'):
        ymin,ymax = -180,0
        Title = 'Pressure Geopause'
    ax.axvline(tdraw,color='k',linestyle = 'dashed')
    ax.set_xticks(timeticks)
    ax.set_xlim([md.date2num(dt.datetime.strptime(tmin,'%H:%M')),md.date2num(dt.datetime.strptime(tmax,'%H:%M'))])
    letter = numbertoimfletter(j)
    ax.text(0.9,0.5,letter,transform=ax.transAxes,fontsize=10)
    ax.set_ylim(ymin,ymax)
    ax.xaxis.set_major_formatter(xfmt)
    if(j==0):
        ax.set_ylabel(r'X [$R_E$]')
    ax.set_xlabel(r'Time [UT]')
    ax.set_title(Title+' for '+solarwind)
    ax.grid()
    return

def plotlinecombinedgeopauseaxes(ax,IMF,j,GM):
    '''
    plotlinegeopause labels the axes, draws the grids,
    writes the title
    '''
    xfmt,timeticks=satellite_format(IMF)
    timeflip = verticaldraw(IMF)
    tdraw =  md.date2num(dt.datetime.strptime(timeflip,'%H:%M'))
    if (IMF == 'runFlipNorth'):
        tmin = '02:00'
        tmax = '12:00'
        solarwind = 'South-to-North IMF'
    if (IMF == 'runFlipSouth'):
        tmin = '06:00'
        tmax = '12:00'
        solarwind = 'North-to-South IMF'
    if (GM == 'SwIono'):
        ymin,ymax = -50,0
        Title = 'Sw H + Iono H Geopauses'
    if (GM == 'SwHIonoO'):
        ymin,ymax = -180,0
        Title = 'Sw H + Iono O Geopauses'
    if (GM == 'SwIonoO'):
        ymin,ymax = -170,0
        Title = 'Sw H + Iono H + O Geopauses'
    if (GM == 'SwIonoO28amu'):
        ymin,ymax = -50,0
        Title = 'Sw H + Iono H + O at 28 amu/cc Geopauses'
    ax.axvline(tdraw,color='k',linestyle = 'dashed')
    ax.set_xticks(timeticks)
    ax.set_xlim([md.date2num(dt.datetime.strptime(tmin,'%H:%M')),md.date2num(dt.datetime.strptime(tmax,'%H:%M'))])
    ax.set_ylim(ymin,ymax)
    ax.xaxis.set_major_formatter(xfmt)
    letter = numbertoimfletter(j)
    ax.text(0.9,0.5,letter,transform=ax.transAxes,fontsize=10)
    if(j==0):
        ax.set_ylabel(r'X [$R_E$]')
    ax.set_xlabel(r'Time [UT]')
    ax.set_title(Title+' for '+solarwind)
    ax.grid()
    return

def plotlinegeopauseall():
    GM = ['SwIono','SwHIonoO','SwIonoO','SwIonoO28amu']
    colorsGM = ['k','b','r','g']
    Title=['Sw H + Iono H','Sw H + Iono O','Sw H + Iono (H + O)','Sw H + Iono (H + O) at 28 amu/cc']
    reversal = ['runFlipNorth','runFlipSouth']
    quantity = ['geopause','rho','pressure']
    filenamenight = ['Data/LineGeopause/'+outflow+'/linenightside'+variable+IMF+'.txt' for IMF in reversal for variable in quantity for outflow in GM]
    filenamenight = np.array(filenamenight)
    filenameday = ['Data/LineGeopause/'+outflow+'/linedayside'+variable+IMF+'.txt' for IMF in reversal for variable in quantity for outflow in GM]
    filenameday = np.array(filenameday)

    '''convert to 3d array
     index 1 is IMF orientation
    index 2 is the geopause type
    index 3 is the outflow type
    '''
    fday = filenameday.reshape((len(reversal),len(quantity),len(GM))) 
    fnight = filenamenight.reshape((len(reversal),len(quantity),len(GM)))
    rc('font',family='serif')
    rc('text', usetex=True)
    for j in range(0,len(quantity)):
        l = []
        fig,ax = plt.subplots(1,2,figsize=(14,6))
        for k in range(0,len(reversal)):
            for i in range(0,len(GM)):
                lc = plotlinegeopausecomponent(ax[k],fnight[k,j,i],colorsGM[i],Title[i])
                l.append(lc)
            plotlinegeopauseaxes(ax[k],reversal[k],k,quantity[j])
        plotpath = 'Paper2/sunearth'+quantity[j]+'.png'
        plt.figlegend((l[0],l[1],l[2],l[3]),Title,(0.2,0.0),ncol=4,fontsize=10,frameon=False)
        plt.subplots_adjust(left=0.07,right=0.95)
        fig.savefig(plotpath,dpi=300,pad_inches = 0.1)
        plt.close(fig)
    return

def plotlinegeopausecombined():
    GM = ['SwIono','SwHIonoO','SwIonoO','SwIonoO28amu']
    Geopause = ['Density Geopause','Mass Density Geopause','Pressure Geopause']
    Title=['Sw H + Iono H','Sw H + Iono O','Sw H + Iono H + O','Sw H + Iono H + O at 28 amu/cc']
    reversal = ['runFlipNorth','runFlipSouth']
    quantity = ['geopause','rho','pressure']
    colorsquant= ['sienna','gold','k']
    filenamenight = ['Data/LineGeopause/'+outflow+'/linenightside'+variable+IMF+'.txt' for IMF in reversal for variable in quantity for outflow in GM]
    filenamenight = np.array(filenamenight)
    filenameday = ['Data/LineGeopause/'+outflow+'/linedayside'+variable+IMF+'.txt' for IMF in reversal for variable in quantity for outflow in GM]
    filenameday = np.array(filenameday)

    '''convert to 3d array
     index 1 is IMF orientation
    index 2 is the geopause type
    index 3 is the outflow type
    '''
    fday = filenameday.reshape((len(reversal),len(quantity),len(GM))) 
    fnight = filenamenight.reshape((len(reversal),len(quantity),len(GM)))
    rc('font',family='serif')
    rc('text', usetex=True)   
    for i in range(0,len(GM)):
        l = []
        fig,ax = plt.subplots(1,2,figsize=(14,6))
        for k in range(0,len(reversal)):
            for j in range(0,len(quantity)):
                lc=plotlinegeopausecomponent(ax[k],fnight[k,j,i],colorsquant[j],Geopause[j])
                l.append(lc)
            plotlinecombinedgeopauseaxes(ax[k],reversal[k],k,GM[i])
        plotpath = 'Paper2/sunearthcombined'+GM[i]+'.png'
        plt.figlegend((l[0],l[1],l[2]),Geopause,(0.3,0.0),ncol=3,fontsize=10,frameon=False)
        plt.subplots_adjust(left=0.07,right=0.95)
        fig.savefig(plotpath,dpi=300,pad_inches = 0.1)
        plt.close(fig)
    return


def drawgeopause(ax,k,GM,reversal,cut):
    home=os.getcwd()+'/'
    Legends = ['Density Geopause','Mass Density Geopause','Pressure Geopause']
    #Colors=['sienna','gold','k']
    #Colors = ['green','violet','gold']
    Colors = ['green','fuchsia','black']
    quantity = ['geopause','rho','pressure']
    root=home+'Data/Geopause/'+GM+'/'+reversal+'/'
    hr,mn=timestring(k)
    time = 't'+hr+mn
    line = []
    for i in range(0,len(quantity)):
        file1 = []
        file1.extend(glob.glob(root+quantity[i]+cut+'*'+time+'*'))
        with open(file1[0],'r') as infile:
            x1,y1=readgeopause(infile)
        l = ax.scatter(x1,y1,s=1,color=Colors[i],label=Legends[i])
        line.append(l)
    return(line)



def artificialfrictionsw(k,GM,reversal,cut,mode):
    '''
    calculate the artificial friction term for the solar wind
    '''
    home=os.getcwd()+'/'
    run =GM+'/'+reversal+'/'
    datalocation = home+run
    hr,mn=timestring(k)
    time = 't00'+hr+mn
    if not os.path.exists(datalocation):
        sys.exit('path does not exist')
    filelist=[]
    filelist.extend(glob.glob(datalocation+cut+'*'+time+'*'))
    filepath = filelist[0]
    filename = PlotProperties(filepath,-40.0,20.0,-20.0,20.0,True,cut,mode)
    data,connectivity,timetick = filename.get_data()
    xyz,bfield,J,stateFluid,nFluid=filename.data_readable(data)
    tau = 1000.     #cutoff time scale
    uc = 100.*1.e3  #100 km/s to m/s cutoff velocity
    alpha = 2
    amutokg = 1.660539e-27
    invcm3toinvm3 = 1.e6
    friction = np.zeros((len(xyz),3))
    vdiff = np.zeros((len(xyz),3))
    if (GM=='SwIono'):
        rhosw = stateFluid[:,0,0]*amutokg*invcm3toinvm3
        rhoio = stateFluid[:,0,1]*amutokg*invcm3toinvm3
        vdiff = (stateFluid[:,1:4,1]-stateFluid[:,1:4,0])*1.e3
        vdiffmagnorm = [np.sqrt(np.dot(vdiff[i,0:3],vdiff[i,0:3]))/uc for i in range(0,len(vdiff))]
        vdiffmagnorm = np.array(vdiffmagnorm)
        friction[:,0] =np.minimum(rhosw,rhoio)*(vdiffmagnorm)**alpha/tau*vdiff[:,0]
        friction[:,1] =np.minimum(rhosw,rhoio)*(vdiffmagnorm)**alpha/tau*vdiff[:,1]
        friction[:,2] =np.minimum(rhosw,rhoio)*(vdiffmagnorm)**alpha/tau*vdiff[:,2]
    if (GM=='SwHIonoO'):
        rhosw = stateFluid[:,0,0]*amutokg*invcm3toinvm3
        rhoO = stateFluid[:,0,1]*amutokg*invcm3toinvm3
        vdiff = (stateFluid[:,1:4,1]-stateFluid[:,1:4,0])*1.e3
        vdiffmagnorm = [np.sqrt(np.dot(vdiff[i,0:3],vdiff[i,0:3]))/uc for i in range(0,len(vdiff))]
        vdiffmagnorm = np.array(vdiffmagnorm)
        friction[:,0] =np.minimum(rhosw,rhoO)*(vdiffmagnorm)**alpha/tau*vdiff[:,0]
        friction[:,1] =np.minimum(rhosw,rhoO)*(vdiffmagnorm)**alpha/tau*vdiff[:,1]
        friction[:,2] =np.minimum(rhosw,rhoO)*(vdiffmagnorm)**alpha/tau*vdiff[:,2]
    if (GM=='SwIonoO' or GM=='SwIonoO28amu'):
        rhosw = stateFluid[:,0,0]*amutokg*invcm3toinvm3
        rhoio = stateFluid[:,0,1]*amutokg*invcm3toinvm3
        rhoO = stateFluid[:,0,2]*amutokg*invcm3toinvm3
        vdiff = (stateFluid[:,1:4,1]-stateFluid[:,1:4,0])*1.e3
        vdiffmagnorm = [np.sqrt(np.dot(vdiff[i,0:3],vdiff[i,0:3]))/uc for i in range(0,len(vdiff))]
        vdiffmagnorm = np.array(vdiffmagnorm)
        vdiff2 = (stateFluid[:,1:4,2]-stateFluid[:,1:4,0])*1.e3
        vdiffmagnorm2 = [np.sqrt(np.dot(vdiff2[i,0:3],vdiff2[i,0:3]))/uc for i in range(0,len(vdiff2))]
        vdiffmagnorm2 = np.array(vdiffmagnorm2)
        friction[:,0] =np.minimum(rhosw,rhoio)*(vdiffmagnorm)**alpha/tau*vdiff[:,0]
        friction[:,0] =friction[:,0]+np.minimum(rhosw,rhoO)*(vdiffmagnorm2)**alpha/tau*vdiff2[:,0]
        friction[:,1] =np.minimum(rhosw,rhoio)*(vdiffmagnorm)**alpha/tau*vdiff[:,1]
        friction[:,1] =friction[:,1]+np.minimum(rhosw,rhoO)*(vdiffmagnorm2)**alpha/tau*vdiff2[:,1]
        friction[:,2] =np.minimum(rhosw,rhoio)*(vdiffmagnorm)**alpha/tau*vdiff[:,2]
        friction[:,2] =friction[:,2]+np.minimum(rhosw,rhoO)*(vdiffmagnorm2)**alpha/tau*vdiff2[:,2]
    print(friction[:,0].max())
    friction = friction/amutokg/1.e4
    p = connectivity-1
    row,column=connectivity.shape
    frictioncent = np.zeros((row,3))
    for d in range(0,3):
        for i in range(0,row):
            frictioncent[i,d] = 0.25*(friction[p[i,0],d]+ \
                                      friction[p[i,1],d]+ \
                                      friction[p[i,2],d]+ \
                                      friction[p[i,3],d])
    return(xyz,frictioncent)

def gyrationsw(k,GM,reversal,cut,mode):
    home=os.getcwd()+'/'
    run =GM+'/'+reversal+'/'
    datalocation = home+run
    hr,mn=timestring(k)
    time = 't00'+hr+mn
    if not os.path.exists(datalocation):
        sys.exit('path does not exist')
    filelist=[]
    filelist.extend(glob.glob(datalocation+cut+'*'+time+'*'))
    filepath = filelist[0]
    filename = PlotProperties(filepath,-40.0,20.0,-20.0,20.0,True,cut,mode)
    data,connectivity,timetick = filename.get_data()
    #xyz,bfield,J,stateFluid,nFluid=filename.data_readable(data)
    xyz,bfield,J,stateFluid,nFluid=filename.data_readableSI(data)
    qe = 1.602176621e-19
    amutokg = 1.660539e-27
    kgtoamu = 1./amutokg
    invcm3toinvm3 = 1.e6
    
    uplus = np.zeros((len(xyz),3)) #initialize to size of gyration
    udiff = np.zeros((len(xyz),3)) #initialize to size of gyration
    #bfield = bfield*1.e-9
    p = connectivity-1
    row,column=connectivity.shape
    
    '''
    if (GM=='SwIono'):
        swi = stateFluid[:,0,0]*invcm3toinvm3
        iono= stateFluid[:,0,1]*invcm3toinvm3
        ntot = swi+iono
        for d in range(0,3): #d+1 -> goes from 1 to 3 because of velocity
            usw,uio = stateFluid[:,d+1,0]*1.e3,stateFluid[:,d+1,1]*1.e3
            uplus[:,d] = np.multiply(swi,usw)+np.multiply(iono,uio)
            uplus[:,d] = np.divide(uplus[:,d],ntot)
            udiff[:,d] = usw-uplus[:,d]        
    if (GM=='SwHIonoO'):
        swi = stateFluid[:,0,0]*invcm3toinvm3
        rhoO = stateFluid[:,0,1]*invcm3toinvm3
        nO= rhoO/16.
        ntot = swi+nO
        for d in range(0,3): #d+1 -> goes from 1 to 3 because of velocity
            usw,uO = stateFluid[:,d+1,0]*1.e3,stateFluid[:,d+1,1]*1.e3
            uplus[:,d] = np.multiply(swi,usw)+np.multiply(nO,uO)
            uplus[:,d] = np.divide(uplus[:,d],ntot)
            udiff[:,d] = usw-uplus[:,d]
    if (GM=='SwIonoO' or GM =='SwIonoO28amu'):
        swi = stateFluid[:,0,0]*invcm3toinvm3
        iono= stateFluid[:,0,1]*invcm3toinvm3
        rhoO=stateFluid[:,0,2]*invcm3toinvm3
        nO= rhoO/16.
        ntot = swi+iono+nO
        for d in range(0,3):
            usw,uO = stateFluid[:,d+1,0]*1.e3,stateFluid[:,d+1,2]*1.e3
            uio = stateFluid[:,d+1,1]*1.e3
            uplus[:,d] = np.multiply(swi,usw)+\
                         np.multiply(nO,uO) +\
                         np.multiply(iono,uio)
            uplus[:,d] = np.divide(uplus[:,d],ntot)
            udiff[:,d] = usw-uplus[:,d]
    '''
    
    if (GM=='SwIono'):
        swi = stateFluid[:,0,0]*kgtoamu
        iono= stateFluid[:,0,1]*kgtoamu
        ntot = swi+iono
        for d in range(0,3): #d+1 -> goes from 1 to 3 because of velocity
            usw,uio = stateFluid[:,d+1,0],stateFluid[:,d+1,1]
            uplus[:,d] = np.multiply(swi,usw)+np.multiply(iono,uio)
            uplus[:,d] = np.divide(uplus[:,d],ntot)
            udiff[:,d] = usw-uplus[:,d]        
    if (GM=='SwHIonoO'):
        swi = stateFluid[:,0,0]*kgtoamu
        rhoO = stateFluid[:,0,1]*kgtoamu
        nO= rhoO/16.
        ntot = swi+nO
        for d in range(0,3): #d+1 -> goes from 1 to 3 because of velocity
            usw,uO = stateFluid[:,d+1,0],stateFluid[:,d+1,1]
            uplus[:,d] = np.multiply(swi,usw)+np.multiply(nO,uO)
            uplus[:,d] = np.divide(uplus[:,d],ntot)
            udiff[:,d] = usw-uplus[:,d]
    if (GM=='SwIonoO' or GM =='SwIonoO28amu'):
        swi = stateFluid[:,0,0]*kgtoamu
        iono= stateFluid[:,0,1]*kgtoamu
        rhoO=stateFluid[:,0,2]*kgtoamu
        nO= rhoO/16.
        ntot = swi+iono+nO
        for d in range(0,3):
            usw,uO = stateFluid[:,d+1,0],stateFluid[:,d+1,2]
            uio = stateFluid[:,d+1,1]
            uplus[:,d] = np.multiply(swi,usw)+\
                         np.multiply(nO,uO) +\
                         np.multiply(iono,uio)
            uplus[:,d] = np.divide(uplus[:,d],ntot)
            udiff[:,d] = usw-uplus[:,d]

    '''
    nDensity = np.zeros((len(xyz),nFluid))
    for i in range(0,nFluid):
        nDensity[:,i] = stateFluid[:,0,i]*kgtoamu
    if (GM == 'SwHIonoO'):
        nDensity[:,1] = nDensity[:,1]/16.  
    if (GM == 'SwIonoO' or GM == 'SwIonoO28amu'):
        nDensity[:,2] = nDensity[:,2]/16.
    ntot = np.sum(nDensity,1)
    
    for d in range(0,3):
        for i in range(0,nFluid):
            uplus[:,d] = uplus[:,d]+nDensity[:,i]*stateFluid[:,d+1,i]
        uplus[:,d] = np.divide(uplus[:,d],ntot)
    udiff = np.zeros((len(xyz),3,nFluid))
    '''
    udiffb = np.zeros((len(xyz),3))
    udiffb = np.cross(udiff,bfield)
    gyration = np.zeros((len(xyz),3))
    gyrationcent = np.zeros((row,3))
    for d in range(0,3):
        gyration[:,d] = qe*swi[:]*udiffb[:,d]
    gyration = gyration/amutokg/1.e4
  
    for i in range(0,row):
        gyrationcent[i,:] = 0.25*(gyration[p[i,0],:]+ \
                                      gyration[p[i,1],:]+ \
                                      gyration[p[i,2],:]+ \
                                      gyration[p[i,3],:])
    return(xyz,gyrationcent)

def bodyforcesw(k,GM,reversal,cut,mode):
    home=os.getcwd()+'/'
    run =GM+'/'+reversal+'/'
    datalocation = home+run
    hr,mn=timestring(k)
    time = 't00'+hr+mn
    if not os.path.exists(datalocation):
        sys.exit('path does not exist')
    filelist=[]
    filelist.extend(glob.glob(datalocation+cut+'*'+time+'*'))
    filepath = filelist[0]
    filename = PlotProperties(filepath,-40.0,20.0,-20.0,20.0,True,cut,mode)
    data,connectivity,timetick = filename.get_data()
    xyz,bfield,J,stateFluid,nFluid=filename.data_readable(data)
    qe = 1.602176621e-19
    amutokg = 1.660539e-27
    invcm3toinvm3 = 1.e6
    nPatoPa=1.e-9
    muAtoA=1.e-6
    J = J*muAtoA
    bfield = bfield*1.e-9
    jxb = np.zeros((len(xyz),3))
    jxb[:,0] = J[:,1]*bfield[:,2]-J[:,2]*bfield[:,1]
    jxb[:,1] = J[:,2]*bfield[:,0]-J[:,0]*bfield[:,2]
    jxb[:,2] = J[:,0]*bfield[:,1]-J[:,1]*bfield[:,0]
    row,column=connectivity.shape
    jxbcenterx = np.zeros(row)
    jxbcenterz = np.zeros(row)
    xcenter = np.zeros(row)
    zcenter = np.zeros(row)
    nratcenter = np.zeros(row)
    dpedx = np.zeros(row)
    dpedz = np.zeros(row)
    bodyx = np.zeros(row)
    bodyz = np.zeros(row)
    Re =6.3781e6
    p = connectivity-1 #recalibrate array index by 1 to match python
    if (GM=='SwIono'):
        swi = stateFluid[:,0,0]*invcm3toinvm3
        iono= stateFluid[:,0,1]*invcm3toinvm3
        psw = stateFluid[:,4,0]*nPatoPa
        pio = stateFluid[:,4,1]*nPatoPa
        ntot = swi+iono
        nrat = np.divide(swi,ntot)
        ptot = psw+pio
    if (GM=='SwHIonoO'):
        swi = stateFluid[:,0,0]*invcm3toinvm3
        rhoO = stateFluid[:,0,1]*invcm3toinvm3
        psw = stateFluid[:,4,0]*nPatoPa
        pO = stateFluid[:,4,1]*nPatoPa
        nO= rhoO/16.
        ntot = swi+nO
        nrat = np.divide(swi,ntot)
        ptot = psw+pO
    if (GM=='SwIonoO' or GM =='SwIonoO28amu'):
        swi = stateFluid[:,0,0]*invcm3toinvm3
        iono= stateFluid[:,0,1]*invcm3toinvm3
        rhoO=stateFluid[:,0,2]*invcm3toinvm3
        psw = stateFluid[:,4,0]*nPatoPa
        pio = stateFluid[:,4,1]*nPatoPa
        pO = stateFluid[:,4,2]*nPatoPa
        nO= rhoO/16.
        ntot = swi+iono+nO
        nrat = np.divide(swi,ntot)
        ptot = psw+pio+pO
    pe = 0.2*ptot #electron pressure
    for i in range(0,row):
        jxbcenterbottom = (jxb[p[i,1],0]+jxb[p[i,0],0])/2.
        jxbcenterup =  (jxb[p[i,2],0]+jxb[p[i,3],0])/2.
        jxbcenterx[i] = (jxbcenterbottom+jxbcenterup)/2.
        jxbcenterbottom = (jxb[p[i,1],2]+jxb[p[i,0],2])/2.
        jxbcenterup =  (jxb[p[i,2],2]+jxb[p[i,3],2])/2.
        jxbcenterz[i] = (jxbcenterbottom+jxbcenterup)/2.
        xcenter[i] = 0.5*(xyz[p[i,1],0]+xyz[p[i,0],0])
        zcenter[i] = 0.5*(xyz[p[i,1],2]+xyz[p[i,2],2])
        nratcenter[i] = 0.25*(nrat[p[i,0]]+nrat[p[i,1]]+\
                              nrat[p[i,2]]+nrat[p[i,3]])
        ''' electron pressure gradient calculation
        calculate grad p in x direction by computing the 
        slope of the bottom row and top row of the cell and average them
        calculate grad p in z direction by computing the slope
        of the left side and right side of the cell and average them
        '''
        dpedx[i] = 0.5*(pe[p[i,1]]-pe[p[i,0]]+pe[p[i,2]]-pe[p[i,3]])/\
               ((xyz[p[i,1],0]-xyz[p[i,0],0])*Re)
        dpedz[i] = 0.5*(pe[p[i,2]]-pe[p[i,1]]+pe[p[i,3]]-pe[p[i,0]])/\
               ((xyz[p[i,2],2]-xyz[p[i,1],2])*Re)
        bodyx[i] = nratcenter[i]*(jxbcenterx[i]-dpedx[i])
        bodyz[i] = nratcenter[i]*(jxbcenterz[i]-dpedz[i])
    bodyx = bodyx/amutokg/1.e4
    bodyz = bodyz/amutokg/1.e4
    return(xcenter,zcenter,bodyx,bodyz)

def plotlastclosedcomponent(ax,k,GM,reversal,cut):
    '''
    plotlastclosed is a function that draws the last closed field line in the meridional plane
    
    '''
    home=os.getcwd()+'/'
    Title=['Magnetopause']
    quantity = 'lastclosed'
    root=home+'Data/Geopause/'+GM+'/'+reversal+'/'
    file1 = []
    hr,mn=timestring(k)
    time = 't'+hr+mn
    file1.extend(glob.glob(root+quantity+cut+'*'+time+'*'))
    filepath = file1[0]
    with open(filepath,'r') as infile:
        xnight1,ynight1,xday1,yday1=readlastclosed(infile)
    l, =ax.plot(xnight1,ynight1,color='cyan',label=Title[0])
    ax.plot(xday1,yday1,color='cyan')
    return(l)

def efield(k,GM,reversal,cut,mode):
    home=os.getcwd()+'/'
    run =GM+'/'+reversal+'/'
    datalocation = home+run
    hr,mn=timestring(k)
    time = 't00'+hr+mn
    if not os.path.exists(datalocation):
        sys.exit('path does not exist')
    filelist=[]
    filelist.extend(glob.glob(datalocation+cut+'*'+time+'*'))
    filepath = filelist[0]
    filename = PlotProperties(filepath,-40.0,20.0,-20.0,20.0,True,cut,mode)
    data,connectivity,timetick = filename.get_data()
    xyz,bfield,J,stateFluid,nFluid=filename.data_readable(data)
    qe = 1.602176621e-19
    amutokg = 1.660539e-27
    invcm3toinvm3 = 1.e6
    gyration = np.zeros((len(xyz),3))
    uplus = gyration #initialize to size of gyration
    udiff = gyration #initialize to size of gyration
    bfield = bfield*1.e-9
    p = connectivity-1
    row,column=connectivity.shape
    efieldcent = np.zeros((row,3))
    xcenter = np.zeros(row)
    zcenter = np.zeros(row)
    if (GM=='SwIono'):
        swi = stateFluid[:,0,0]*invcm3toinvm3
        iono= stateFluid[:,0,1]*invcm3toinvm3
        ntot = swi+iono
        for d in range(0,3): #d+1 -> goes from 1 to 3 because of velocity
            usw,uio = stateFluid[:,d+1,0]*1.e3,stateFluid[:,d+1,1]*1.e3
            uplus[:,d] = np.multiply(swi,usw)+np.multiply(iono,uio)
            uplus[:,d] = np.divide(uplus[:,d],ntot)
    if (GM=='SwHIonoO'):
        swi = stateFluid[:,0,0]*invcm3toinvm3
        rhoO = stateFluid[:,0,1]*invcm3toinvm3
        nO= rhoO/16.
        ntot = swi+nO
        for d in range(0,3): #d+1 -> goes from 1 to 3 because of velocity
            usw,uO = stateFluid[:,d+1,0]*1.e3,stateFluid[:,d+1,1]*1.e3
            uplus[:,d] = np.multiply(swi,usw)+np.multiply(nO,uO)
            uplus[:,d] = np.divide(uplus[:,d],ntot)
    if (GM=='SwIonoO' or GM =='SwIonoO28amu'):
        swi = stateFluid[:,0,0]*invcm3toinvm3
        iono= stateFluid[:,0,1]*invcm3toinvm3
        rhoO=stateFluid[:,0,2]*invcm3toinvm3
        nO= rhoO/16.
        ntot = swi+iono+nO
        for d in range(0,3):
            usw,uO = stateFluid[:,d+1,0]*1.e3,stateFluid[:,d+1,2]*1.e3
            uio = stateFluid[:,d+1,1]*1.e3
            uplus[:,d] = np.multiply(swi,usw)+\
                         np.multiply(nO,uO) +\
                         np.multiply(iono,uio)
            uplus[:,d] = np.divide(uplus[:,d],ntot)
    efield = -np.cross(uplus,bfield)
    for i in range(0,row):
        for j in range(0,3):
            efieldcent[i,j] = 0.25*(efield[p[i,0],j]+\
                                    efield[p[i,1],j]+\
                                    efield[p[i,2],j]+\
                                    efield[p[i,3],j])
        xcenter[i] = 0.5*(xyz[p[i,1],0]+xyz[p[i,0],0])
        zcenter[i] = 0.5*(xyz[p[i,1],2]+xyz[p[i,2],2])

    return(xcenter,zcenter,efieldcent)
    

def returnstate(k,GM,reversal,cut,mode):
    home=os.getcwd()+'/'
    run =GM+'/'+reversal+'/'
    datalocation = home+run
    hr,mn=timestring(k)
    time = 't00'+hr+mn
    if not os.path.exists(datalocation):
        sys.exit('path does not exist')
    filelist=[]
    filelist.extend(glob.glob(datalocation+cut+'*'+time+'*'))
    filepath = filelist[0]
    filename = PlotProperties(filepath,-40.0,20.0,-20.0,20.0,True,cut,mode)
    data,connectivity,timetick = filename.get_data()
    xyz,bfield,J,stateFluid,nFluid=filename.data_readable(data)    
    return(xyz,stateFluid,nFluid,bfield,J)


def clbstatelabels(GM,ns):

    if (GM=='SwIono'):
        if (ns==0):
            sub = 'sw,H}'
            species = 'Sw'
        if (ns==1):
            sub = 'io,H}'
            species = 'Io'
    if (GM=='SwHIonoO'):
        if (ns==0):
            sub = 'sw,H}'
            species = 'Sw'
        if (ns==1):
            sub = 'O}'
            species = 'O'
    if (GM=='SwIonoO'or GM=='SwIonoO28amu'):
        if (ns==0):
            sub = 'sw,H}'
            species = 'Sw'
        if (ns==1):
            sub = 'io,H}'
            species = 'Io'
        if (ns==2):
            sub = 'O}'
            species = 'O'
    clblabels = np.array([[r'$\rho_{'+sub+'$ [amu/cm$^3$]',
                           r'u$_{x,'+sub+'$ [km/s]',
                           r'($\rho u_x$)'+'$_{'+sub+'$'+' [amu/cm$^2$/sec]'],
                          [r'$p_{'+sub+'$ [nPa]',
                           r'u$_{z,'+sub+'$ [km/s]',
                           r'$(\rho u_z)$'+'$_{'+sub+'$'+' [amu/cm$^2$/sec]']])
    return (clblabels,species)

def clbtemplabels(GM):
    if (GM=='SwIono'):
        templabels = np.array([r'$T_{sw}$ [eV]',r'$T_{io}$ [eV]'])
    if (GM=='SwHIonoO'):
        templabels = np.array([r'$T_{sw}$ [eV]',r'$T_{O}$ [eV]'])
    if (GM=='SwIonoO'or GM=='SwIonoO28amu'):
        templabels = np.array([r'$T_{sw}$ [eV]',r'$T_{io}$ [eV]',r'$T_{O}$ [eV]'])
    return (templabels)


def plotstate(o):
    GM = ['SwIono','SwHIonoO','SwIonoO','SwIonoO28amu']
    Geopause = ['Density Geopause','Mass Density Geopause','Pressure Geopause']
    Title=['Sw H + Iono H','Sw H + Iono O','Sw H + Iono H + O','Sw H + Iono H + O at 28 amu/cc']
    IMF = ['South-to-North IMF','North-to-South IMF']
    reversal = ['runFlipNorth','runFlipSouth']
    quantity = ['geopause','rho','pressure']
    mode = '2d'
    minimum = 720
    plotnumber = 1+12*60
    colorsquant= ['sienna','green','k']
    cut = 'y'
    p = ordinate(cut)
    N = 1000

    if(p==1):
        ylabel = r'Y [$R_E$]'
    if(p==2):
        ylabel = r'Z [$R_E$]'
    rc('font',family='serif')
    rc('text', usetex=True)
    colormaps = np.array([['viridis','seismic','seismic'],['plasma','seismic','seismic']])
    for j in range(0,len(reversal)):
        for k in range(minimum,plotnumber):
             xyz,stateFluid,nFluid,bfield,J=returnstate(k,GM[o],reversal[j],cut,mode)
             xi,yi = gridpoints(xyz[:,0],N),gridpoints(xyz[:,p],N)
             for ns in range(0,nFluid):
                 fig,ax=plt.subplots(2,3,figsize=(15,10))
                 clblabels,species = clbstatelabels(GM[o],ns)
                 rhoi = twodinterpolate(xyz[:,0],xyz[:,p],xi,yi,stateFluid[:,0,ns])
                 vxi = twodinterpolate(xyz[:,0],xyz[:,p],xi,yi,stateFluid[:,1,ns])
                 voi = twodinterpolate(xyz[:,0],xyz[:,p],xi,yi,stateFluid[:,1+p,ns])
                 pti = twodinterpolate(xyz[:,0],xyz[:,p],xi,yi,stateFluid[:,4,ns])
                 rhovx,rhovoi = rhoi*vxi*1.e5, rhoi*voi*1.e5
                 z = np.array([[rhoi,vxi,rhovx],[pti,voi,rhovoi]])
                 vmins = np.array([[1.e-2,-600,-1.e9],[1.e-3,-200,-5.e8]])
                 vmaxs = np.array([[1.e3,600,1.e9],[1.e1,200,5.e8]])
                 for column in range(0,3):
                     for row in range(0,2):
                         if (column==0):
                             cb=ax[row,column].pcolormesh(xi,yi,z[row,column],
    #                                                      vmin = vmins[row,column],
    #                                                      vmax = vmaxs[row,column],
                                                          norm = colors.LogNorm(vmin=vmins[row,column],vmax = vmaxs[row,column]),
                                                      cmap=colormaps[row,column])
                         else:
                             cb=ax[row,column].pcolormesh(xi,yi,z[row,column],
                                                        vmin = vmins[row,column],
                                                         vmax = vmaxs[row,column],
                                                          #norm = colors.SymLogNorm(linthresh = 10, linscale = 1,vmin = vmins[row,column],vmax = vmaxs[row,column]),
                                                      cmap=colormaps[row,column])
                         drawgeopause(ax[row,column],k,GM[o],reversal[j],cut)
                         plotlastclosedcomponent(ax[row,column],k,GM[o],reversal[j],cut)
                         clb=fig.colorbar(cb, ax=ax[row,column])
                         clb.set_label(clblabels[row,column])
                         draw_earth(ax[row,column])
                         ax[row,column].set_xlim(-80,20)
                         ax[row,column].set_ylim(-50,50)
                         if (row ==1):
                             ax[row,column].set_xlabel(r'X [$R_E$]')
                         if(column == 0):
                             ax[row,column].set_ylabel(ylabel)
                 Suptitle = Title[o]+' '+IMF[j]
                 fig.suptitle(Suptitle,x=0.5,y=0.92)
                 home=os.getcwd()+'/'
                 plotlocation = home+'State/'+GM[o]+'/'+reversal[j]+'/'
                 if not os.path.exists(plotlocation):
                     os.makedirs(plotlocation)
                     print('making a plot directory')
                 plotpath = plotlocation+'state'+species+'.png'
                 fig.savefig(plotpath,dpi=300)
                 plt.close(fig)
    return


def plotfriction(o):
    GM = ['SwIono','SwHIonoO','SwIonoO','SwIonoO28amu']
    Geopause = ['Density Geopause','Mass Density Geopause','Pressure Geopause']
    Title=['Sw H + Iono H','Sw H + Iono O','Sw H + Iono H + O','Sw H + Iono H + O at 28 amu/cc']
    reversal = ['runFlipNorth','runFlipSouth']
    quantity = ['geopause','rho','pressure']
    mode = '2d'
    minimum = 0
    plotnumber = 1+12*60
    colorsquant= ['sienna','green','k']
    cut = 'y'
    N = 500
    rc('font',family='serif')
    rc('text', usetex=True)
    Title = [r'Solar Wind Friction in the X-direction', 
             r'Solar Wind Friction in the Y-direction', 
             r'Solar Wind Friction in the Z-direction']
    
    for k in range(minimum,plotnumber):
        fig,ax=plt.subplots(3,2,figsize=(10,10))
        plt.subplots_adjust(hspace=0.1)
        p = ordinate(cut)
        if(p==1):
            ylabel = r'Y [$R_E$]'
        if(p==2):
            ylabel = r'Z [$R_E$]'
        for j in range(0,2):
            xyz,friction=artificialfrictionsw(k,GM[o],reversal[j],cut,mode)
            xi,yi = gridpoints(xyz[:,0],N),gridpoints(xyz[:,p],N)
            for i in range(0,3):
                frictioni=twodinterpolate(xyz[:,0],xyz[:,p],xi,yi,friction[:,i])
                cs = ax[i,j].pcolormesh(xi,yi,frictioni,vmin=-15000,vmax=15000,cmap='seismic')
                drawgeopause(ax[i,j],k,GM[o],reversal[j],cut)
                ax[i,j].set_xlim(-80,20)
                ax[i,j].set_ylim(-20,20)
                draw_earth(ax[i,j])
                if (i==2):
                    ax[i,j].set_xlabel(r'X [$R_E$]')
                if (j==0):
                    ax[i,j].set_ylabel(ylabel)
                ax[i,j].set_title(Title[i])
                ax[i,j].set_aspect('equal')
                ax[i,j].grid()
        cbar_ax=fig.add_axes([0.15,0.06,0.7,0.01])
        clb=fig.colorbar(cs,cax=cbar_ax,aspect=20,orientation ='horizontal')
        clb.set_label('[amu/cm$^2$/sec$^2$]')
        hr,mn=timestring(k)
        home=os.getcwd()+'/'
        fig.text(0.47,0.08,'T='+hr+':'+mn)
        plotlocation = home+'Friction/'+GM[o]+'/'
        if not os.path.exists(plotlocation):
            os.makedirs(plotlocation)
            print('making a plot directory')
        plotpath = plotlocation+'meridionalfrictiont'+hr+mn+'.png'
        print(hr,mn,GM[o])
        fig.savefig(plotpath,dpi=300)
        #plt.show()
        plt.close(fig)

def plotgyration(o):
    GM = ['SwIono','SwHIonoO','SwIonoO','SwIonoO28amu']
    Geopause = ['Density Geopause','Mass Density Geopause','Pressure Geopause']
    Title=['Sw H + Iono H','Sw H + Iono O','Sw H + Iono H + O','Sw H + Iono H + O at 28 amu/cc']
    reversal = ['runFlipNorth','runFlipSouth']
    Legends = ['Number Density Geopause','Mass Density Geopause','Pressure Geopause','Magnetopause']
    IMF = [' for South-to-North IMF',' for North-to-South IMF']
    quantity = ['geopause','rho','pressure']
    mode = '2d'
    minimum = 12*60#0
    plotnumber = 1+12*60
    colorsquant= ['sienna','green','k']
    cmaps = ['seismic','seismic','bwr']
    loffset = [0.25,0.25,0.25,0.25]
    cut = 'y'
    N = 1000
    rc('font',family='serif')
    rc('text', usetex=True)
    Momentum = np.array([['Gyration in the x-direction','Body force in the x-direction','Friction in the x-direction'],['Gyration in the z-direction','Body force in the z-direction','Friction in the z-direction']])
    linthreshold = np.array([[1.e4,1.e4,1e2],[1e4,1e4,1e2]])
    vmins=np.array([[-1.e7,-1.e7,-1.e5],[-1.e7,-1.e7,-1.e5]])
    vmaxs=np.array([[1.e7,1.e7,1.e5],[1e7,1.e7,1.e5]])
    for a in range(0,2):
        for k in range(minimum,plotnumber):
            fig,ax=plt.subplots(2,3,figsize=(15,10))
            p = ordinate(cut)
            if(p==1):
                ylabel = r'Y [$R_E$]'
            if(p==2):
                ylabel = r'Z [$R_E$]'
            xyz,friction=artificialfrictionsw(k,GM[o],reversal[a],cut,mode)
            xyz,gyration =gyrationsw(k,GM[o],reversal[a],cut,mode)
            xc,zc,bodyx,bodyz=bodyforcesw(k,GM[o],reversal[a],cut,mode)
            xci,zci = gridpoints(xc,N),gridpoints(zc,N)
            gyrationxi=twodinterpolate(xc,zc,xci,zci,gyration[:,0])
            gyrationzi=twodinterpolate(xc,zc,xci,zci,gyration[:,2])
            bodyxi = twodinterpolate(xc,zc,xci,zci,bodyx)
            bodyzi = twodinterpolate(xc,zc,xci,zci,bodyz)
            frictionxi =twodinterpolate(xc,zc,xci,zci,friction[:,0])
            frictionzi = twodinterpolate(xc,zc,xci,zci,friction[:,2])
            xx = np.array([[xci,xci,xci],[xci,xci,xci]])
            yy = np.array([[zci,zci,zci],[zci,zci,zci]])
            zz = np.array([[gyrationxi,bodyxi,frictionxi],[gyrationzi,bodyzi,frictionzi]])

            for i in range(0,2):
                for j in range(0,3):
                    pcm=ax[i,j].pcolormesh(xx[i,j],yy[i,j],zz[i,j],norm=colors.SymLogNorm(linthresh=linthreshold[i,j],linscale=1,vmin=vmins[i,j],vmax=vmaxs[i,j]),cmap=cmaps[j])
                    #pcm=ax[i,j].pcolormesh(xx[i,j],yy[i,j],zz[i,j],vmin=vmins[i,j],vmax=vmaxs[i,j],cmap=cmaps[j])
                    l=drawgeopause(ax[i,j],k,GM[o],reversal[a],cut)
                    lc=plotlastclosedcomponent(ax[i,j],k,GM[o],reversal[a],cut)
                    ax[i,j].set_xlim(-80,20)
                    ax[i,j].set_ylim(-50,50)
                    clb=fig.colorbar(pcm, ax=ax[i,j])
                    clb.set_label('[amu/cm$^2$/sec$^2$]')
                    draw_earth(ax[i,j])
                    if(j==0):
                        ax[i,j].set_ylabel(r'Z [$R_E$]')
                    if(i==1):
                        ax[i,j].set_xlabel(r'X [$R_E$]')
                    ax[i,j].set_title(Momentum[i,j])
            Suptitle = Title[o]
            hr,mn=timestring(k)
            plt.suptitle(Suptitle+IMF[a]+' at T='+hr+':'+mn) 
            print(hr,mn,GM[o],reversal[a])
            home=os.getcwd()+'/'
            plt.figlegend((l[0],l[1],l[2],lc),Legends,(loffset[o],0.01),markerscale=3,ncol=4,fontsize=10,frameon=False)
            plotlocation = home+'Momentumwithlastclosed/'+GM[o]+'/'+reversal[a]+'/'
            if not os.path.exists(plotlocation):
                os.makedirs(plotlocation)
                print('making a plot directory')
            plotpath = plotlocation+'meridionalmomentumt'+hr+mn+'.png'
            fig.savefig(plotpath,dpi=300)
            plt.close(fig)
    return


def plotratio(o):
    GM = ['SwIono','SwHIonoO','SwIonoO','SwIonoO28amu']
    Geopause = ['Density Geopause','Mass Density Geopause','Pressure Geopause']
    Title=['Sw H + Iono H','Sw H + Iono O','Sw H + Iono H + O','Sw H + Iono H + O at 28 amu/cc']
    Legends = ['Number Density Geopause','Mass Density Geopause','Pressure Geopause','Magnetopause']
    reversal = ['runFlipNorth','runFlipSouth']
    loffset = [0.25,0.25,0.25,0.25]
    IMF = [' for South-to-North IMF',' for North-to-South IMF']
    quantity = ['geopause','rho','pressure']
    mode = '2d'
    minimum = 12*60#0
    plotnumber = 1+12*60
    colorsquant= ['sienna','green','k']
    cmaps = ['coolwarm','bwr','seismic']
    cut = 'y'
    N = 1000 
    rc('font',family='serif')
    rc('text', usetex=True)
    Momentum = np.array([['Gyration in the x-direction','Body force in the x-direction','Friction in the x-direction'],['Gyration in the z-direction','Body force in the z-direction','Friction in the z-direction']])
    linthreshold = np.array([[1.e5,1.e5,10],[1e-13,1e5,10]])
    vmins=np.array([[-1.e8,-1.e8,-1.e5],[-1.e-8,-1.e8,-1.e5]])
    vmaxs=np.array([[1.e8,1.e8,1.e5],[1e-8,1.e8,1.e5]])
    for a in range(0,2):
        for k in range(minimum,plotnumber):
            fig,ax=plt.subplots(2,3,figsize=(15,10))
            p = ordinate(cut)
            if(p==1):
                ylabel = r'Y [$R_E$]'
            if(p==2):
                ylabel = r'Z [$R_E$]'
            xyz,friction=artificialfrictionsw(k,GM[o],reversal[a],cut,mode)
            xyz,gyration =gyrationsw(k,GM[o],reversal[a],cut,mode)
            xc,zc,bodyx,bodyz=bodyforcesw(k,GM[o],reversal[a],cut,mode)
            xi,yi = gridpoints(xyz[:,0],N),gridpoints(xyz[:,p],N)
            xci,zci = gridpoints(xc,N),gridpoints(zc,N)
            gyrationxi=np.absolute(twodinterpolate(xc,zc,xci,zci,gyration[:,0]))
            gyrationzi=np.absolute(twodinterpolate(xc,zc,xci,zci,gyration[:,2]))
            bodyxi = np.absolute(twodinterpolate(xc,zc,xci,zci,bodyx))
            bodyzi = np.absolute(twodinterpolate(xc,zc,xci,zci,bodyz))
            frictionxi =np.absolute(twodinterpolate(xc,zc,xci,zci,friction[:,0]))
            frictionzi = np.absolute(twodinterpolate(xc,zc,xci,zci,friction[:,2]))
            momx = gyrationxi+bodyxi+frictionxi
            momz = gyrationzi+bodyzi+frictionzi
            gyrfracxi,gyrfraczi = np.divide(gyrationxi,momx),np.divide(gyrationzi,momz)
            bodyfracxi,bodyfraczi = np.divide(bodyxi,momx),np.divide(bodyzi,momz)
            fricfracxi,fricfraczi = np.divide(frictionxi,momx),np.divide(frictionzi,momz)
            xx = np.array([[xci,xci,xci],[xci,xci,xci]])
            yy = np.array([[zci,zci,zci],[zci,zci,zci]])
            zz = np.array([[gyrfracxi,bodyfracxi,fricfracxi],[gyrfraczi,bodyfraczi,fricfraczi]])

            for i in range(0,2):
                for j in range(0,3):
                    pcm=ax[i,j].pcolormesh(xx[i,j],yy[i,j],zz[i,j],vmin=0,vmax=1,cmap = 'magma')
                    l=drawgeopause(ax[i,j],k,GM[o],reversal[a],cut)
                    lc=plotlastclosedcomponent(ax[i,j],k,GM[o],reversal[a],cut)
                    ax[i,j].set_xlim(-80,20)
                    ax[i,j].set_ylim(-50,50)
                    clb=fig.colorbar(pcm, ax=ax[i,j])
                    clb.set_label('Abundance')
                    draw_earth(ax[i,j])
                    if(j==0):
                        ax[i,j].set_ylabel(r'Z [$R_E$]')
                    if(i==1):
                        ax[i,j].set_xlabel(r'X [$R_E$]')
                    ax[i,j].set_title(Momentum[i,j])
            Suptitle = Title[o]
            hr,mn=timestring(k)
            plt.figlegend((l[0],l[1],l[2],lc),Legends,(loffset[o],0.01),markerscale=3,ncol=4,fontsize=10,frameon=False)
            plt.suptitle(Suptitle+IMF[a]+' at T='+hr+':'+mn) 
            print(hr,mn,GM[o])
            home=os.getcwd()+'/'
            plotlocation = home+'MomentumRatio/'+GM[o]+'/'+reversal[a]+'/'
            if not os.path.exists(plotlocation):
                os.makedirs(plotlocation)
                print('making a plot directory')
            
            plotpath = plotlocation+'meridionalratiot'+hr+mn+'.png'
            fig.savefig(plotpath,dpi=300)
            plt.close(fig)
    return

def massfluid(GM,ns):
    if (GM == 'SwIono'):
        if (ns == 0):
            mass = 1.
        if (ns == 1):
            mass = 1.
    if (GM == 'SwHIonoO'):
        if (ns == 0):
            mass = 1.
        if (ns == 1):
            mass = 16.
    if (GM == 'SwIonoO' or GM == 'SwIonoO28amu'):
        if (ns == 0):
            mass = 1.
        if (ns == 1):
            mass = 1.
        if (ns == 2):
            mass = 16.
    return(mass)
    
def plottemperature(o):
    GM = ['SwIono','SwHIonoO','SwIonoO','SwIonoO28amu']
    numberfluid = [2,2,3,3]
    loffset = [0.12,0.12,0.25,0.25]
    Geopause = ['Density Geopause','Mass Density Geopause','Pressure Geopause']
    Legends = ['Number Density Geopause','Mass Density Geopause','Pressure Geopause','Magnetopause']
    Title=['Sw H + Iono H','Sw H + Iono O','Sw H + Iono H + O','Sw H + Iono H + O at 28 amu/cc']
    IMF = ['South-to-North IMF','North-to-South IMF']
    reversal = ['runFlipNorth','runFlipSouth']
    quantity = ['geopause','rho','pressure']
    mode = '2d'
    minimum = 720
    plotnumber = 1+12*60
    colorsquant= ['sienna','green','k']
    cut = 'y'
    p = ordinate(cut)
    N = 1000
    mpro = 1.6726219e-27
    invcm3tom3 = 1.e6
    nPatoPa = 1.e-9
    kb = 1.38064852e-23
    if(p==1):
        ylabel = r'Y [$R_E$]'
    if(p==2):
        ylabel = r'Z [$R_E$]'
    rc('font',family='serif')
    rc('text', usetex=True)
    clblabels = clbtemplabels(GM[o])
    colormaps = np.array([['viridis','seismic','seismic'],['plasma','seismic','seismic']])
    lengthsize = 5*numberfluid[o]
    vmins = np.array([[1.,1.,1.],[1.,1.,1.]])
    vmaxs = np.array([[1.e6,1.e6,1.e6],[1.e6,1.e6,1.e6]])
    for k in range(minimum,plotnumber):
        fig,ax=plt.subplots(len(reversal),numberfluid[o],figsize=(lengthsize,10))
        for j in range(0,len(reversal)):
             xyz,stateFluid,nFluid,bfield,J=returnstate(k,GM[o],reversal[j],cut,mode)
             xi,yi = gridpoints(xyz[:,0],N),gridpoints(xyz[:,p],N)
             for ns in range(0,nFluid):
                 ms = massfluid(GM[o],ns)
                 ni = twodinterpolate(xyz[:,0],xyz[:,p],xi,yi,stateFluid[:,0,ns])/ms*invcm3tom3
                 pti = twodinterpolate(xyz[:,0],xyz[:,p],xi,yi,stateFluid[:,4,ns])*nPatoPa
                 temper = np.divide(pti,ni)/kb/11600.
                 cb=ax[j,ns].pcolormesh(xi,yi,temper,
                                        norm = colors.LogNorm(
                                            vmin=vmins[j,ns],
                                            vmax = vmaxs[j,ns]),cmap = 'magma')
                 l=drawgeopause(ax[j,ns],k,GM[o],reversal[j],cut)
                 lc=plotlastclosedcomponent(ax[j,ns],k,GM[o],reversal[j],cut)
                 draw_earth(ax[j,ns])
                 ax[j,ns].set_xlim(-80,20)
                 ax[j,ns].set_ylim(-50,50)
                 if (j ==1):
                     ax[j,ns].set_xlabel(r'X [$R_E$]')
                 if(ns == 0):
                     ax[j,ns].set_ylabel(ylabel)
                 clb=fig.colorbar(cb, ax=ax[j,ns])
        clb.set_label(clblabels[ns])
        plt.figlegend((l[0],l[1],l[2],lc),Legends,(loffset[o],0.01),ncol=4,fontsize=10,frameon=False)
        Suptitle = Title[o]
        fig.suptitle(Suptitle,x=0.5,y=0.92)
        home=os.getcwd()+'/'
        plotlocation = home+'Temperature/'+GM[o]+'/'
        if not os.path.exists(plotlocation):
            os.makedirs(plotlocation)
            print('making a plot directory')
        plotpath = plotlocation+'temperature'+'.png'
        fig.savefig(plotpath,dpi=300)
        plt.close(fig)
    return

def plotefield(o):
    GM = ['SwIono','SwHIonoO','SwIonoO','SwIonoO28amu']
    Geopause = ['Density Geopause','Mass Density Geopause','Pressure Geopause']
    Title=['Sw H + Iono H','Sw H + Iono O','Sw H + Iono H + O','Sw H + Iono H + O at 28 amu/cc']
    IMF = ['South-to-North IMF','North-to-South IMF']
    reversal = ['runFlipNorth','runFlipSouth']
    Legends = ['Number Density Geopause','Mass Density Geopause','Pressure Geopause','Magnetopause']
    quantity = ['geopause','rho','pressure']
    mode = '2d'
    minimum = 720
    plotnumber = 1+12*60
    colorsquant= ['sienna','green','k']
    cut = 'y'
    p = ordinate(cut)
    N = 1000
    mpro = 1.6726219e-27
    invcm3tom3 = 1.e6
    nPatoPa = 1.e-9
    kb = 1.38064852e-23
    if(p==1):
        ylabel = r'Y [$R_E$]'
    if(p==2):
        ylabel = r'Z [$R_E$]'
    rc('font',family='serif')
    rc('text', usetex=True)
    clblabels = np.array([r'$E_x$ [mV/m]',r'$E_y$ [mV/m]',r'$E_z$ [mV/m]'])
    #lengthsize = 5*numberfluid[o]
    #vmins = np.array([[-10,-100.,-10.],[-10.,-100.,-10.]])
    #vmaxs = np.array([[10,100.,10.],[10.,100.,10.]])
    vmins = np.array([[-5,-20.,-10.],[-5,-20.,-10.]])
    vmaxs= np.array([[5.,20.,10.],[5.,20.,10.]])
    for k in range(minimum,plotnumber):
        fig,ax=plt.subplots(len(reversal),3,figsize=(15,10))
        for j in range(0,len(reversal)):
             xc,zc,Ef=efield(k,GM[o],reversal[j],cut,mode)
             xi,yi = gridpoints(xc,N),gridpoints(zc,N)
             for ns in range(0,3):
                 Efi = twodinterpolate(xc,zc,xi,yi,Ef[:,ns])*1.e3 #in mV/m
                 cb=ax[j,ns].pcolormesh(xi,yi,Efi,vmin=vmins[j,ns],
                                         vmax=vmaxs[j,ns],
                                        #norm = colors.SymLogNorm(linthresh = 1.e-2, linscale = 0.5,vmin = vmins[j,ns],vmax = vmaxs[j,ns]),
                                        cmap = 'seismic')
                 draw_earth(ax[j,ns])
                 l=drawgeopause(ax[j,ns],k,GM[o],reversal[j],cut)
                 lc=plotlastclosedcomponent(ax[j,ns],k,GM[o],reversal[j],cut)
                 ax[j,ns].set_xlim(-80,20)
                 ax[j,ns].set_ylim(-50,50)
                 if (j ==1):
                     ax[j,ns].set_xlabel(r'X [$R_E$]')
                 if(ns == 0):
                     ax[j,ns].set_ylabel(ylabel)
                 clb=fig.colorbar(cb, ax=ax[j,ns])
                 clb.set_label(clblabels[ns])
        plt.figlegend((l[0],l[1],l[2],lc),Legends,(0.25,0.01),ncol=4,fontsize=10,frameon=False)
        Suptitle = Title[o]
        fig.suptitle(Suptitle,x=0.5,y=0.92)
        home=os.getcwd()+'/'
        plotlocation = home+'Efield/'+GM[o]+'/'
        if not os.path.exists(plotlocation):
            os.makedirs(plotlocation)
            print('making a plot directory')
        plotpath = plotlocation+'efield'+'.png'
        fig.savefig(plotpath,dpi=300)
        plt.close(fig)
    return


def plotJy(core):
    o = core%4  #divide by 4 for number of ionospheric compositions
    j = core/4  #even means northward IMF, odd means southward IMF, division by 4 (as an integer) for number of IMF orientations
    GM = ['SwIono','SwHIonoO','SwIonoO','SwIonoO28amu']
    Geopause = ['Density Geopause','Mass Density Geopause','Pressure Geopause']
    Title=['Sw H + Iono H','Sw H + Iono O','Sw H + Iono H + O','Sw H + Iono H + O at 28 amu/cc']
    IMF = ['South-to-North IMF','North-to-South IMF']
    reversal = ['runFlipNorth','runFlipSouth']
    quantity = ['geopause','rho','pressure']
    Legends = ['Number Density Geopause','Mass Density Geopause','Pressure Geopause']
    mode = '2d'
    minimum = 0
    plotnumber = 1+12*60
    colorsquant= ['sienna','green','k']
    cut = 'y'
    p = ordinate(cut)
    N = 1000

    if(p==1):
        ylabel = r'Y [$R_E$]'
    if(p==2):
        ylabel = r'Z [$R_E$]'
    rc('font',family='serif')
    rc('text', usetex=True)
    for k in range(minimum,plotnumber):
        fig,ax=plt.subplots(1,1)
        hr,mn=timestring(k)
        time = 't'+hr+mn
        xyz,stateFluid,nFluid,bfield,Jdens=returnstate(k,GM[o],reversal[j],cut,mode)
        xi,yi = gridpoints(xyz[:,0],N),gridpoints(xyz[:,p],N)
        Jyi = twodinterpolate(xyz[:,0],xyz[:,p],xi,yi,Jdens[:,1])
        Jyi = 1.e3*Jyi # in nA/m^2
        cb = ax.pcolormesh(xi,yi,Jyi,vmin = -10,vmax=10,cmap='seismic')
        clb=fig.colorbar(cb, ax=ax)
        clb.set_label(r'$J_y$ [nA/m$^2$]')
        l = drawgeopause(ax,k,GM[o],reversal[j],cut)
        draw_earth(ax)
        ax.set_xlim(-80,20)
        ax.set_ylim(-50,50)
        ax.set_xlabel(r'X [$R_E$]')
        ax.set_ylabel(ylabel)
        Suptitle = Title[o]+' at '+'T='+hr+':'+mn+' for '+IMF[j]
        fig.suptitle(Suptitle,x=0.5,y=0.92)
        home=os.getcwd()+'/'
        plotlocation = home+'Jy/'+GM[o]+'/'+reversal[j]+'/'
        plt.figlegend((l[0],l[1],l[2]),Legends,(0.03,-0.004),ncol=3,fontsize=8,frameon=False)
        if not os.path.exists(plotlocation):
            os.makedirs(plotlocation)
            print('making a plot directory')
        plotpath = plotlocation+'Jy'+time+'.png'
        fig.savefig(plotpath,dpi=300)
        plt.close(fig)
        print(hr,mn,GM[o],reversal[j],o,j)
    return




'''
Plots Jy
'''


if __name__ == "__main__":
    multiprocessing.freeze_support()
    some_list = range(0,8)
    num_proc = 8
    p = multiprocessing.Pool(num_proc)
    p.map(plotJy, some_list)


#plotefield(3)
#plottemperature(3)

'''
if __name__ == "__main__":
    multiprocessing.freeze_support()
    some_list = range(0,4)
    num_proc = 4
    p = multiprocessing.Pool(num_proc)
    p.map(plotgyration, some_list)
'''


#plotratio(2)
#plotlinegeopauseall()
#plotlinegeopausecombined()
#plotpanelgeopause()
#plotpanelall()
#plotpanelmagnetopause()
'''
if __name__ == "__main__":
    plotgyration(2)

'''
