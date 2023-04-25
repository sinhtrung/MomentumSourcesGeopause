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
from matplotlib import ticker
import cProfile
import multiprocessing
from classes import *
import plotting as myplot



def momentum(k,GM,reversal,cut,mode):
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
    xyz,bfield,J,stateFluid,nFluid=filename.data_readableSI(data) #everything is in SI Units
    lines = len(xyz)
    jxb = np.zeros((len(xyz),3))
    gyration,udiff = np.zeros((len(xyz),3,nFluid)),np.zeros((len(xyz),3,nFluid))
    jxb = np.cross(J,bfield)
    row,column=connectivity.shape
    # initialize arrays
    jxbcenterx,jxbcenterz = np.zeros(row),np.zeros(row)
    xcenter,zcenter = np.zeros(row),np.zeros(row)
    nratcenter = np.zeros((row,nFluid))
    dpedx,dpedz = np.zeros(row),np.zeros(row)    
    nrat = np.zeros((len(xyz),nFluid))
    bodyx,bodyz = np.zeros((row,nFluid)),np.zeros((row,nFluid))
    friction = np.zeros((len(xyz),3,nFluid))
    gyrationcent = np.zeros((row,3,nFluid))
    frictioncent = np.zeros((row,3,nFluid))
    vdiffmagnorm = np.zeros((len(xyz),nFluid,nFluid))
    vdiff = np.zeros((len(xyz),3,nFluid,nFluid))
    uplus = np.zeros((len(xyz),3))
    qe = 1.602176621e-19
    amutokg = 1.660539e-27
    kgtoamu = 1./amutokg
    Re = 6.3781e6
    tau = 1000.     #cutoff time scale
    uc = 100.*1.e3  #100 km/s to m/s cutoff velocity
    alpha = 2.
    p = connectivity-1 #recalibrate array index by 1 to match python from Fortran
    nDensity = np.zeros((len(xyz),nFluid))
    pressure =  np.zeros((len(xyz),nFluid))
    #Assign to number density array
    for i in range(0,nFluid):
        nDensity[:,i] = stateFluid[:,0,i]*kgtoamu #convert mass density in SI units to amu/m^3
        pressure[:,i] = stateFluid[:,4,i]
       
    if (GM == 'SwHIonoO'): #change mass density to number density for oxygen by dividing by the appropriate mass unit
        nDensity[:,1] = nDensity[:,1]/16.  
    if (GM == 'SwIonoO' or GM == 'SwIonoO28amu'):
        nDensity[:,2] = nDensity[:,2]/16.

    ntot = np.sum(nDensity,1)
    ptot = np.sum(pressure,1)

    #species velocity
    uspecies = np.zeros((len(xyz),3,nFluid))
    ###### calculate the gyration term
    # first calculate u+ (we are ignoring charge since we have a singly ionized plasma)
    uspecies = stateFluid[:,1:4,:]
    # calculate the numerator by calculating n_i*u_i for each direction
    for d in range(0,3):
        for i in range(0,nFluid):
            uplus[:,d] = uplus[:,d]+nDensity[:,i]*uspecies[:,d,i]

    # calculate u+ by dividing by the total number density
    for d in range(0,3):
        uplus[:,d] = np.divide(uplus[:,d],ntot[:])

    # there is no charge given that the ion and electron species have charge e

    # Calculate the velocity in the charge averaged frame
    for d in range(0,3):
        for i in range(0,nFluid):
            udiff[:,d,i] = uspecies[:,d,i]-uplus[:,d]
    ucrossb = np.zeros((len(xyz),3,nFluid))
    # Calculate the u x B term
    for i in range(0,nFluid):
        ucrossb[:,:,i] = np.cross(udiff[:,:,i],bfield)

    # Multiply by the charge and density of the ion to obtain the gyration term
    for d in range(0,3):
        for i in range(0,nFluid):
            gyration[:,d,i] = qe*nDensity[:,i]*ucrossb[:,d,i]
    ##### End of calculation of the gyration term
            
  
    ##### calculate friction
    # calculate difference between ion velocities, and normalized velocities
    for i in range(0,nFluid):
        for j in range(0,nFluid):
            vdiff[:,0:3,i,j] = stateFluid[:,1:4,j]-stateFluid[:,1:4,i]
            vdiffmagnorm[:,i,j] = [(np.sqrt(np.dot(vdiff[k,0:3,i,j],vdiff[k,0:3,i,j]))/uc)**alpha/tau for k in range(0,len(xyz))]

    #statefluid is used because it is the mass density
    for d in range(0,3):
        for i in range(0,nFluid):
            for j in range(0,nFluid):
                friction[:,d,i] = friction[:,d,i]+np.minimum(stateFluid[:,0,j],stateFluid[:,0,i])*vdiff[:,d,i,j]*vdiffmagnorm[:,i,j]

    
    #### End of friction calculation

    ##### This is the start of the calculation for the body forces
    # Calculate the fraction of each ion number density to total density
    for i in range(0,nFluid):
        nrat[:,i]=np.divide(nDensity[:,i],ntot)
    pe = 0.2*ptot #electron pressure
    

    
    for i in range(0,row):
        jxbcenterbottom = (jxb[p[i,1],0]+jxb[p[i,0],0])/2.
        jxbcenterup =  (jxb[p[i,2],0]+jxb[p[i,3],0])/2.
        jxbcenterx[i] = (jxbcenterbottom+jxbcenterup)/2.
        jxbcenterbottom = (jxb[p[i,1],2]+jxb[p[i,0],2])/2.
        jxbcenterup =  (jxb[p[i,2],2]+jxb[p[i,3],2])/2.
        jxbcenterz[i] = (jxbcenterbottom+jxbcenterup)/2.
        xcenter[i] = 0.5*(xyz[p[i,1],0]+xyz[p[i,0],0])/Re     #center values of x and z in each cell
        zcenter[i] = 0.5*(xyz[p[i,1],2]+xyz[p[i,2],2])/Re
        
        ''' electron pressure gradient calculation
        calculate grad p in x direction by computing the         slope of the bottom row and top row of the cell and average them
        calculate grad p in z direction by computing the slope
        of the left side and right side of the cell and average them
        '''
        dpedx[i] = 0.5*(pe[p[i,1]]-pe[p[i,0]]+pe[p[i,2]]-pe[p[i,3]])/\
               ((xyz[p[i,1],0]-xyz[p[i,0],0]))
        dpedz[i] = 0.5*(pe[p[i,2]]-pe[p[i,1]]+pe[p[i,3]]-pe[p[i,0]])/\
               ((xyz[p[i,2],2]-xyz[p[i,1],2]))
        for j in range(0,nFluid):
            nratcenter = 0.25*(nrat[p[i,0],j]+nrat[p[i,1],j]+\
                              nrat[p[i,2],j]+nrat[p[i,3],j])
            bodyx[i,j] = nratcenter*(jxbcenterx[i]-dpedx[i])
            bodyz[i,j] = nratcenter*(jxbcenterz[i]-dpedz[i])
    #### End of body force term calculations

    # Calculate center values of gyration and friction terms
    for i in range(0,row):
        for j in range(0,nFluid):
            gyrationcent[i,:,j] = 0.25*(gyration[p[i,0],:,j]+ \
                                      gyration[p[i,1],:,j]+ \
                                      gyration[p[i,2],:,j]+ \
                                      gyration[p[i,3],:,j])
            frictioncent[i,:,j] = 0.25*(friction[p[i,0],:,j]+ \
                                      friction[p[i,1],:,j]+ \
                                      friction[p[i,2],:,j]+ \
                                      friction[p[i,3],:,j])
    
    bodyx = bodyx/amutokg/1.e4 #convert to N/m^3 to amu/cm^2/s^2
    bodyz = bodyz/amutokg/1.e4
    gyrationcent = gyrationcent/amutokg/1.e4
    frictioncent = frictioncent/amutokg/1.e4
    return(xcenter,zcenter,bodyx,bodyz,gyrationcent,frictioncent,nFluid)

def specieslabelpaper(GM):
    if (GM == 'SwIono'):
        Title = [r'Sw H$^+$ ',r'Iono H$^+$ ']
        Filename = ['SwH','IonoH']
    if (GM == 'SwHIonoO'):
        Title = [r'Sw H$^+$ ',r'Iono O$^+$ ']
        Filename = ['SwH','O']
    if (GM == 'SwIonoO' or GM =='SwIonoO28amu'):
        Title = [r'Sw H$^+$ ', r'Iono H$^+$ ',r'Iono O$^+$ ']
        Filename = ['SwH','IonoH','O']
    return(Title,Filename)

def plotmomentum(o): #perform 2d interpolation to calculate the force density
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
    a = 0
    colorsquant= ['sienna','green','k']
    cmaps = ['seismic','seismic','bwr']
    loffset = [0.25,0.25,0.25,0.25]
    cut = 'y'
    N = 1000
    rc('font',family='serif')
    rc('text', usetex=True)
    Momentum = np.array([['gyration in the x-direction','bulk terms in the x-direction','friction in the x-direction'],['gyration in the z-direction','bulk terms in the z-direction','friction in the z-direction']])
    linthreshold = np.array([[1.e4,1.e4,1e2],[1e4,1e4,1e2]])
    vmins=np.array([[-1.e7,-1.e7,-1.e5],[-1.e7,-1.e7,-1.e5]])
    vmaxs=np.array([[1.e7,1.e7,1.e5],[1e7,1.e7,1.e5]])
    letter = np.array([['(a)','(b)','(c)'],['(d)','(e)','(f)']])
    for a in range(0,2):
        for k in range(minimum,plotnumber):
           
            p = ordinate(cut)
            if(p==1):
                ylabel = r'Y [$R_E$]'
            if(p==2):
                ylabel = r'Z [$R_E$]'    
            xc,zc,bodyxx,bodyzz,gyrationcent,frictioncent,nFluid=momentum(k,GM[o],reversal[a],cut,mode)
            #xold,zold,bodyx,bodyz = bodyforcespecies(k,GM[o],reversal[a],cut,mode)
            xci,zci = gridpoints(xc,N),gridpoints(zc,N)
            Titleadd,Extension = specieslabelpaper(GM[o])
            for m in range(0,nFluid):
                fig,ax=plt.subplots(2,3,figsize=(15,10))
                gyrationxi =twodinterpolate(xc,zc,xci,zci,gyrationcent[:,0,m])
                gyrationzi = twodinterpolate(xc,zc,xci,zci,gyrationcent[:,2,m])
                bodyxi = twodinterpolate(xc,zc,xci,zci,bodyxx[:,m])
                bodyzi = twodinterpolate(xc,zc,xci,zci,bodyzz[:,m])
                frictionxi =twodinterpolate(xc,zc,xci,zci,frictioncent[:,0,m])
                frictionzi = twodinterpolate(xc,zc,xci,zci,frictioncent[:,2,m])
                xx = np.array([[xci,xci,xci],[xci,xci,xci]])
                yy = np.array([[zci,zci,zci],[zci,zci,zci]])
                zz = np.array([[gyrationxi,bodyxi,frictionxi],[gyrationzi,bodyzi,frictionzi]])
                for i in range(0,2):
                    for j in range(0,3):
                        pcm=ax[i,j].pcolormesh(xx[i,j],yy[i,j],zz[i,j],shading='nearest',norm=colors.SymLogNorm(linthresh=linthreshold[i,j],linscale=1,vmin=vmins[i,j],vmax=vmaxs[i,j],base=10),cmap=cmaps[j])
                        l=myplot.drawgeopause(ax[i,j],k,GM[o],reversal[a],cut)
                        lc=myplot.plotlastclosedcomponent(ax[i,j],k,GM[o],reversal[a],cut)
                        ax[i,j].set_xlim(-80,20)
                        ax[i,j].set_ylim(-50,50)
                        clb=fig.colorbar(pcm, ax=ax[i,j])
                        clb.set_label('[amu/cm$^2$/sec$^2$]')
                        draw_earth(ax[i,j])
                        if(j==0):
                            ax[i,j].set_ylabel(r'Z [$R_E$]')
                        if(i==1):
                            ax[i,j].set_xlabel(r'X [$R_E$]')
                        ax[i,j].set_title(letter[i,j]+' '+Titleadd[m]+Momentum[i,j])
                        #ax[i,j].text(0.9,0.5,letter[i,j],transform=ax[i,j].transAxes,fontsize=10,color='black')
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
                plotpath = plotlocation+'meridionalmomentum'+Extension[m]+'t'+hr+mn+'.png'
                fig.savefig(plotpath,dpi=300)
                plt.close(fig)
                print('plot saved to:',plotpath)
            
    return

def plotratiorelative(o): #plot the relative ratio of the magnitude of the force densities to each other
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
    a = 0
    colorsquant= ['sienna','green','k']
    cmaps = ['seismic','seismic','bwr']
    loffset = [0.25,0.25,0.25,0.25]
    cut = 'y'
    N = 1000
    rc('font',family='serif')
    rc('text', usetex=True)
    Momentum = np.array([['gyration/body force in x-direction','friction/body force in the x-direction','friction/gyration in the x-direction'],['gyration/body force in the z-direction','friction/body force in the z-direction','gyration/friction in the z-direction']])
    linthreshold = np.array([[1.e4,1.e4,1e2],[1e4,1e4,1e2]])

    for a in range(0,2):
        for k in range(minimum,plotnumber):
           
            p = ordinate(cut)
            if(p==1):
                ylabel = r'Y [$R_E$]'
            if(p==2):
                ylabel = r'Z [$R_E$]'    
            xc,zc,bodyxx,bodyzz,gyrationcent,frictioncent,nFluid=momentum(k,GM[o],reversal[a],cut,mode)
            #xold,zold,bodyx,bodyz = bodyforcespecies(k,GM[o],reversal[a],cut,mode)
            xci,zci = gridpoints(xc,N),gridpoints(zc,N)
            Titleadd,Extension = specieslabelpaper(GM[o])
            for m in range(0,nFluid):
                fig,ax=plt.subplots(2,3,figsize=(15,10))
                gyrationxi =np.absolute(twodinterpolate(xc,zc,xci,zci,gyrationcent[:,0,m]))
                gyrationzi = np.absolute(twodinterpolate(xc,zc,xci,zci,gyrationcent[:,2,m]))
                bodyxi = np.absolute(twodinterpolate(xc,zc,xci,zci,bodyxx[:,m]))
                bodyzi = np.absolute(twodinterpolate(xc,zc,xci,zci,bodyzz[:,m]))
                frictionxi =np.absolute(twodinterpolate(xc,zc,xci,zci,frictioncent[:,0,m]))
                frictionzi = np.absolute(twodinterpolate(xc,zc,xci,zci,frictioncent[:,2,m]))
                gyrtobfxi,gyrtobfzi = np.divide(gyrationxi,bodyxi),np.divide(gyrationzi,bodyzi)
                frictobfxi,frictobfzi = np.divide(frictionxi,bodyxi),np.divide(frictionzi,bodyzi)
                frictogyrxi,frictogyrzi = np.divide(frictionxi,gyrationxi),np.divide(frictionzi,gyrationzi)
                xx = np.array([[xci,xci,xci],[xci,xci,xci]])
                yy = np.array([[zci,zci,zci],[zci,zci,zci]])
                zz = np.array([[gyrtobfxi,frictobfxi,frictogyrxi],[gyrtobfzi,frictobfzi,frictogyrzi]])
                for i in range(0,2):
                    for j in range(0,3):
                        #pcm=ax[i,j].pcolormesh(xx[i,j],yy[i,j],zz[i,j],shading='nearest',vmin=0,vmax=10,cmap='magma')
                        pcm=ax[i,j].pcolormesh(xx[i,j],yy[i,j],zz[i,j],shading='nearest',norm=colors.LogNorm(vmin=1.e-2,vmax=1.e2),cmap='magma')
                        l=myplot.drawgeopause(ax[i,j],k,GM[o],reversal[a],cut)
                        lc=myplot.plotlastclosedcomponent(ax[i,j],k,GM[o],reversal[a],cut)
                        ax[i,j].set_xlim(-80,20)
                        ax[i,j].set_ylim(-50,50)
                        clb=fig.colorbar(pcm, ax=ax[i,j])
                        clb.set_label('Fractional Force Density')
                        draw_earth(ax[i,j])
                        if(j==0):
                            ax[i,j].set_ylabel(r'Z [$R_E$]')
                        if(i==1):
                            ax[i,j].set_xlabel(r'X [$R_E$]')
                        ax[i,j].set_title(Titleadd[m]+Momentum[i,j])
                Suptitle = Title[o]
                hr,mn=timestring(k)
                plt.suptitle(Suptitle+IMF[a]+' at T='+hr+':'+mn) 
                print(hr,mn,GM[o],reversal[a])
                home=os.getcwd()+'/'
                plt.figlegend((l[0],l[1],l[2],lc),Legends,(loffset[o],0.01),markerscale=3,ncol=4,fontsize=10,frameon=False)
                plotlocation = home+'MomentumRatio/'+GM[o]+'/'+reversal[a]+'/'
                if not os.path.exists(plotlocation):
                    os.makedirs(plotlocation)
                    print('making a plot directory')
                plotpath = plotlocation+'meridionalmomentumratiorelative'+Extension[m]+'t'+hr+mn+'.png'
                fig.savefig(plotpath,dpi=300)
                plt.close(fig)
                print('plot saved to:',plotpath)
    return

def plotforceratio(o): #plot the relative ratio of the magnitude of the force densities to the magnitude of the total force density
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
    a = 0
    colorsquant= ['sienna','green','k']
    cmaps = ['seismic','seismic','bwr']
    loffset = [0.25,0.25,0.25,0.25]
    cut = 'y'
    N = 1000
    rc('font',family='serif')
    rc('text', usetex=True)
    Momentum = np.array([['gyration/total force in x-direction','body force/total force in the x-direction','friction/total force in the x-direction'],['gyration/total force in the z-direction','body force/total force in the z-direction','friction/total force in the z-direction']])

    for a in range(0,2):
        for k in range(minimum,plotnumber):
           
            p = ordinate(cut)
            if(p==1):
                ylabel = r'Y [$R_E$]'
            if(p==2):
                ylabel = r'Z [$R_E$]'    
            xc,zc,bodyxx,bodyzz,gyrationcent,frictioncent,nFluid=momentum(k,GM[o],reversal[a],cut,mode)
            #xold,zold,bodyx,bodyz = bodyforcespecies(k,GM[o],reversal[a],cut,mode)
            xci,zci = gridpoints(xc,N),gridpoints(zc,N)
            Titleadd,Extension = specieslabelpaper(GM[o])
            for m in range(0,nFluid):
                fig,ax=plt.subplots(2,3,figsize=(15,10))
                gyrationxi =twodinterpolate(xc,zc,xci,zci,gyrationcent[:,0,m])
                gyrationzi = twodinterpolate(xc,zc,xci,zci,gyrationcent[:,2,m])
                bodyxi = twodinterpolate(xc,zc,xci,zci,bodyxx[:,m])
                bodyzi = twodinterpolate(xc,zc,xci,zci,bodyzz[:,m])
                frictionxi =twodinterpolate(xc,zc,xci,zci,frictioncent[:,0,m])
                frictionzi = twodinterpolate(xc,zc,xci,zci,frictioncent[:,2,m])
                totalfxi = gyrationxi+bodyxi+frictionxi
                totalfzi = gyrationzi+bodyzi+frictionzi
                
                gyrtototfxi,gyrtototfzi = np.absolute(np.divide(gyrationxi,totalfxi)),np.absolute(np.divide(gyrationzi,totalfzi))
                bodytototfxi,bodytototfzi = np.absolute(np.divide(bodyxi,totalfxi)),np.absolute(np.divide(bodyzi,totalfzi))
                frictototfxi,frictototfzi = np.absolute(np.divide(frictionxi,totalfxi)),np.absolute(np.divide(frictionzi,totalfzi))
                xx = np.array([[xci,xci,xci],[xci,xci,xci]])
                yy = np.array([[zci,zci,zci],[zci,zci,zci]])
                zz = np.array([[gyrtototfxi,bodytototfxi,frictototfxi],[gyrtototfzi,bodytototfzi,frictototfzi]])
                for i in range(0,2):
                    for j in range(0,3):
                        pcm=ax[i,j].pcolormesh(xx[i,j],yy[i,j],zz[i,j],shading='nearest',vmin=0,vmax=2,cmap='plasma')
                        l=myplot.drawgeopause(ax[i,j],k,GM[o],reversal[a],cut)
                        lc=myplot.plotlastclosedcomponent(ax[i,j],k,GM[o],reversal[a],cut)
                        ax[i,j].set_xlim(-80,20)
                        ax[i,j].set_ylim(-50,50)
                        clb=fig.colorbar(pcm, ax=ax[i,j])
                        #clb.set_label('Fractional Force Density')
                        draw_earth(ax[i,j])
                        if(j==0):
                            ax[i,j].set_ylabel(r'Z [$R_E$]')
                        if(i==1):
                            ax[i,j].set_xlabel(r'X [$R_E$]')
                        ax[i,j].set_title(Titleadd[m]+Momentum[i,j])
                Suptitle = Title[o]
                hr,mn=timestring(k)
                plt.suptitle(Suptitle+IMF[a]+' at T='+hr+':'+mn) 
                print(hr,mn,GM[o],reversal[a])
                home=os.getcwd()+'/'
                plt.figlegend((l[0],l[1],l[2],lc),Legends,(loffset[o],0.01),markerscale=3,ncol=4,fontsize=10,frameon=False)
                plotlocation = home+'MomentumRatio/'+GM[o]+'/'+reversal[a]+'/'
                if not os.path.exists(plotlocation):
                    os.makedirs(plotlocation)
                    print('making a plot directory')
                plotpath = plotlocation+'meridionalmomentumratioabsolute'+Extension[m]+'t'+hr+mn+'.png'
                fig.savefig(plotpath,dpi=300)
                plt.close(fig)
                print('plot saved to:',plotpath)
    return



def plotratiototal(o): #plot the ratio of the magnitude of the force densities relative to the sum of the magnitude of the force densities
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
    a = 0
    colorsquant= ['sienna','green','k']
    cmaps = ['seismic','seismic','bwr']
    loffset = [0.25,0.25,0.25,0.25]
    cut = 'y'
    N = 1000
    rc('font',family='serif')
    rc('text', usetex=True)
    Momentum = np.array([['gyration in the x-direction','bulk terms in the x-direction','friction in the x-direction'],['gyration in the z-direction','bulk terms in the z-direction','friction in the z-direction']])
    linthreshold = np.array([[1.e4,1.e4,1e2],[1e4,1e4,1e2]])
    vmins=np.array([[-1.e7,-1.e7,-1.e5],[-1.e7,-1.e7,-1.e5]])
    vmaxs=np.array([[1.e7,1.e7,1.e5],[1e7,1.e7,1.e5]])
    letter = np.array([['(a)','(b)','(c)'],
              ['(d)','(e)','(f)']
    ])
    for a in range(0,2):
        for k in range(minimum,plotnumber):
           
            p = ordinate(cut)
            if(p==1):
                ylabel = r'Y [$R_E$]'
            if(p==2):
                ylabel = r'Z [$R_E$]'    
            xc,zc,bodyxx,bodyzz,gyrationcent,frictioncent,nFluid=momentum(k,GM[o],reversal[a],cut,mode)
            #xold,zold,bodyx,bodyz = bodyforcespecies(k,GM[o],reversal[a],cut,mode)
            xci,zci = gridpoints(xc,N),gridpoints(zc,N)
            Titleadd,Extension = specieslabelpaper(GM[o])
            for m in range(0,nFluid):
                fig,ax=plt.subplots(2,3,figsize=(15,10))
                gyrationxi =np.absolute(twodinterpolate(xc,zc,xci,zci,gyrationcent[:,0,m]))
                gyrationzi = np.absolute(twodinterpolate(xc,zc,xci,zci,gyrationcent[:,2,m]))
                bodyxi = np.absolute(twodinterpolate(xc,zc,xci,zci,bodyxx[:,m]))
                bodyzi = np.absolute(twodinterpolate(xc,zc,xci,zci,bodyzz[:,m]))
                frictionxi =np.absolute(twodinterpolate(xc,zc,xci,zci,frictioncent[:,0,m]))
                frictionzi = np.absolute(twodinterpolate(xc,zc,xci,zci,frictioncent[:,2,m]))
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
                        pcm=ax[i,j].pcolormesh(xx[i,j],yy[i,j],zz[i,j],shading='nearest',vmin=0,vmax=1,cmap='plasma')
                        #pcm=ax[i,j].pcolormesh(xx[i,j],yy[i,j],zz[i,j],vmin=vmins[i,j],vmax=vmaxs[i,j],cmap=cmaps[j])
                        l=myplot.drawgeopause(ax[i,j],k,GM[o],reversal[a],cut)
                        lc=myplot.plotlastclosedcomponent(ax[i,j],k,GM[o],reversal[a],cut)
                        ax[i,j].set_xlim(-80,20)
                        ax[i,j].set_ylim(-50,50)
                        clb=fig.colorbar(pcm, ax=ax[i,j])
                        clb.set_label('Absolute Abundance')
                        draw_earth(ax[i,j])
                        if(j==0):
                            ax[i,j].set_ylabel(r'Z [$R_E$]')
                        if(i==1):
                            ax[i,j].set_xlabel(r'X [$R_E$]')
                        ax[i,j].set_title(letter[i,j]+' '+Titleadd[m]+Momentum[i,j])
                Suptitle = Title[o]
                hr,mn=timestring(k)
                plt.suptitle(Suptitle+IMF[a]+' at T='+hr+':'+mn) 
                print(hr,mn,GM[o],reversal[a])
                home=os.getcwd()+'/'
                plt.figlegend((l[0],l[1],l[2],lc),Legends,(loffset[o],0.01),markerscale=3,ncol=4,fontsize=10,frameon=False)
                plotlocation = home+'MomentumRatio/'+GM[o]+'/'+reversal[a]+'/'
                if not os.path.exists(plotlocation):
                    os.makedirs(plotlocation)
                    print('making a plot directory')
                plotpath = plotlocation+'meridionalmomentumratio'+Extension[m]+'t'+hr+mn+'.png'
                fig.savefig(plotpath,dpi=300)
                plt.close(fig)
                print('plot saved to:',plotpath)
    return


def plotratioscatter(o):
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
    a = 0
    colorsquant= ['sienna','green','k']
    cmaps = ['seismic','seismic','bwr']
    loffset = [0.25,0.25,0.25,0.25]
    cut = 'y'
    N = 1000
    rc('font',family='serif')
    rc('text', usetex=True)
    Momentum = np.array([['gyration in the x-direction','body force fraction in the x-direction','friction in the x-direction'],['gyration in the z-direction','body force in the z-direction','friction in the z-direction']])
    linthreshold = np.array([[1.e4,1.e4,1e2],[1e4,1e4,1e2]])
    vmins=np.array([[-1.e7,-1.e7,-1.e5],[-1.e7,-1.e7,-1.e5]])
    vmaxs=np.array([[1.e7,1.e7,1.e5],[1e7,1.e7,1.e5]])
    for a in range(0,2):
        for k in range(minimum,plotnumber):
            p = ordinate(cut)
            if(p==1):
                ylabel = r'Y [$R_E$]'
            if(p==2):
                ylabel = r'Z [$R_E$]'    
            xc,zc,bodyxx,bodyzz,gyrationcent,frictioncent,nFluid=momentum(k,GM[o],reversal[a],cut,mode)
            Titleadd,Extension = specieslabelpaper(GM[o])
            for m in range(0,nFluid):
                fig,ax=plt.subplots(2,3,figsize=(15,10))
                momx = gyrationcent[:,0,m]+bodyxx[:,m]+frictioncent[:,0,m]
                momz = gyrationcent[:,2,m]+bodyzz[:,m]+frictioncent[:,2,m]
                gyrfracx,gyrfracz = np.divide(gyrationcent[:,0,m],momx), np.divide(gyrationcent[:,2,m],momz)
                bodyfracx,bodyfracz = np.divide(bodyxx[:,m],momx), np.divide(bodyzz[:,m],momz)
                fricfracx,fricfracz = np.divide(frictioncent[:,0,m],momx), np.divide(frictioncent[:,2,m],momz)
                xx = np.array([[xc,xc,xc],[xc,xc,xc]])
                yy = np.array([[zc,zc,zc],[zc,zc,zc]])
                zz = np.array([[gyrfracx,bodyfracx,fricfracx],[gyrfracz,bodyfracz,fricfracz]])
                for i in range(0,2):
                    for j in range(0,3):
                        pcm=ax[i,j].scatter(xx[i,j],yy[i,j],c = zz[i,j],s=1,vmin=0,vmax=1,cmap='magma')
                        l=myplot.drawgeopause(ax[i,j],k,GM[o],reversal[a],cut)
                        lc=myplot.plotlastclosedcomponent(ax[i,j],k,GM[o],reversal[a],cut)
                        ax[i,j].set_xlim(-80,20)
                        ax[i,j].set_ylim(-50,50)
                        clb=fig.colorbar(pcm, ax=ax[i,j])
                        clb.set_label('Absolute Abundance')
                        draw_earth(ax[i,j])
                        if(j==0):
                            ax[i,j].set_ylabel(r'Z [$R_E$]')
                        if(i==1):
                            ax[i,j].set_xlabel(r'X [$R_E$]')
                        ax[i,j].set_title(Titleadd[m]+Momentum[i,j])
                Suptitle = Title[o]
                hr,mn=timestring(k)
                plt.suptitle(Suptitle+IMF[a]+' at T='+hr+':'+mn) 
                print(hr,mn,GM[o],reversal[a])
                home=os.getcwd()+'/'
                plt.figlegend((l[0],l[1],l[2],lc),Legends,(loffset[o],0.01),markerscale=3,ncol=4,fontsize=10,frameon=False)
                plotlocation = home+'MomentumRatio/'+GM[o]+'/'+reversal[a]+'/'
                if not os.path.exists(plotlocation):
                    os.makedirs(plotlocation)
                    print('making a plot directory')
                plotpath = plotlocation+'meridionalmomentumratio'+Extension[m]+'t'+hr+mn+'.png'
                fig.savefig(plotpath,dpi=300)
                plt.close(fig)
            
    return
plotratiototal(1)
plotratiototal(0)


plotmomentum(0)
plotmomentum(1)
