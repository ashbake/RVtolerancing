
import pyzdde.zdde as pyz
import sys

from astropy.io import fits
import matplotlib.pylab as plt
import numpy as np
from matplotlib.patches import Rectangle 
plt.ion()


if len(sys.argv) > 1:
    filename = sys.argv[1]
else:
    filename = 'psf_locs.txt'



def load_orders_wl():
    """
    load orders and wavelength
    """
    f = np.loadtxt('./orders_20220608C_HISPEC_SPECTRO_HK_pyechelle.csv',delimiter=',',skiprows=1).T
    orders, minwl,maxwl = f[0],f[1],f[2]
    order_sub = [2,15,37]
    isort = np.argsort(orders[order_sub])
    return orders[order_sub][isort],minwl[order_sub][isort]+10,maxwl[order_sub][isort]-10



def load_tolerance_data(tol_filename='./tolerance_data_file.txt'):
    """
    tolerance data is in following form:
    #surface_name, surface_number, axis, nominal, min, max
    # 0 1 2 3 4 5: X Y Z TX TY TZ

    load this to make appropriate zemax mods
    """
    f = np.loadtxt(tol_filename,delimiter=',',dtype='str').T
    surf_name = f[0]
    surf_num  = f[1].astype('int')
    axis      = f[2].astype('int')
    nom_val   = f[3].astype('float')
    val       = f[4].astype('float')

    return surf_name, surf_num, axis,nom_val,val/1e6


def getpsf_orders_fibers(ln,orders,minwls,maxwls,NWAVE=15,fiber_locs=[ -0.065, 0.065]):
    """
    loop through orders and fibers getting the psf

    input: fiber_locs
        locations of fiber positions to get PSF for
    """
    NFIBER =  len(fiber_locs)
    NORDER =  len(orders)
    psf_cents = np.zeros((NORDER,NFIBER,NWAVE,3))

    fiboffsetsurf  = 1   # surface of fiber offset
    fiboffsetparam = 1   # parameter for fiber offset position x decenter
    echelle_surf   = 28

    # step through offset
    for io, o in enumerate(list(orders)):
        wl = np.linspace(minwls[io]/1000,maxwls[io]/1000, NWAVE)
        ln.zSetSurfaceParameter(echelle_surf, 2, int(o))
        print('Trace order', o)
        for ifo, fiber_offset in enumerate(fiber_locs):
            ln.zSetSurfaceParameter(fiboffsetsurf, fiboffsetparam, fiber_offset) # set fiber offset
            for iw, w in enumerate(wl):
                ln.zSetWave(1, w, 1.)
                ln.zPushLens()
                #print('WAVELENGTH=' + str(w))
                try: 
                    #psfparams,test = ln.zGetPSF(which='huygens',timeout=1000)
                    #cehx = ln.zOperandValue('CEHX',1,1,0,2,2,0) # 2,2, refers to 64x64 grid size for pupil and image sampling. could increase these
                    #cehy = ln.zOperandValue('CEHY',1,1,0,2,2,0) # found that cehx/y doesn make a huge difference
                    #print('CEHX' + str(cehx))
                    cehx = ln.zOperandValue('CENX',0,1,1,0,15) # 2,2, refers to 64x64 grid size for pupil and image sampling. could increase these
                    cehy = ln.zOperandValue('CENY',0,1,1,0,15)
                    # psf_cents[io,ifo,iw] = [w,psfparams.centerCoordX, psfparams.centerCoordY] # wrong! gives center of box not the centroid
                    psf_cents[io,ifo,iw] = [w,cehx,cehy] # should be right
                except TypeError:
                    print('%s %s %s' %(o,fiber_offset, w))#psf_cents[io,ifo,iw] = [w,0,0]
                    
    return psf_cents


def save_psf_cents(savename,psf_cents,NWAVE,NFIBER,NORDER):
    #savetofits        
    hdr = fits.Header()
    hdr['NAXIS1'] = (3, 'wave, X,Y')
    hdr['NAXIS2'] = (NWAVE,'number of wavelengths')
    hdr['NAXIS3'] = (NFIBER,'number of fibers')
    hdr['NAXIS4'] = (NORDER, 'number of orders')
    hdu = fits.PrimaryHDU(psf_cents,header=hdr)
    hdu.writeto(savename,overwrite=True)


def plot_psf_shifts(X1,Y1,Xshift,Yshift,RVerr,surname,value):
    """
    take psf shift data and turn into RV vector field map

    plot vectors as difference in shift of fiber 1 minus difference in shift of fiber 2 (difference: before/after perturbation)
    """ 
    ifib1, ifib2 = 0, 1# which fiber inds to take
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for io, o in enumerate(list(orders)):
        ax = plt.quiver(X1,Y1,Xshift,Yshift,RVerr,
                   width=0.003,
                   headwidth=4,
                   headlength=6)
        plt.plot(X1,Y1,'gray',alpha=0.3,zorder=-1)

    cbar = plt.colorbar()
    cbar.set_label('RV (m/s)')
    plt.title('HK %s, +/-%s' %(surname,value))
    plt.xlabel('Detector Y Position (mm)')
    plt.ylabel('Detector X Position (mm)')
    plt.savefig('testplot.png')
    # draw detector outline
    detector_halfwidth = 30.7
    #rect = Rectangle((-1*detector_halfwidth,-1*detector_halfwidth), 2*detector_halfwidth, 2*detector_halfwidth,linewidth=1,edgecolor='gray',facecolor='gray')
    #plt.axvline()

    #return X1,Y1,U,V,U*mm_to_ms

def get_psf_shifts(psf_cents_min,psf_cents_max):
    """
    """
    ifib1=0;ifib2=1
    mm_to_ms = 1000/0.015
 
    X1 = psf_cents_min[:,ifib1,:,1]
    X2 = psf_cents_min[:,ifib2,:,1]
    Y1 = psf_cents_min[:,ifib1,:,2]
    Y2 = psf_cents_min[:,ifib2,:,2]
    A1 = psf_cents_max[:,ifib1,:,1]
    A2 = psf_cents_max[:,ifib2,:,1]
    B1 = psf_cents_max[:,ifib1,:,2]
    B2 = psf_cents_max[:,ifib2,:,2]
    U = (X1 - X2) - (A1 - A2) # x vector, RV direction
    V = (Y1 - Y2) - (B1 - B2) # y vector
    dX = X1 - A1
    dY = Y1 - B1
    return X1,Y1,U*mm_to_ms,V*mm_to_ms,dX*mm_to_ms,dY*mm_to_ms #save x shift of one fiber to compare to Jason's result



def get_one_shift(surf_num, axis,  val, nom_val, orders, minwls,maxwls, surf_name=None, saveme=False, plotme=False):
    """
    """
    NWAVE = 3
    fiber_locs=[ -0.065, 0.065, 0.195]
    NFIBER=len(fiber_locs)
    if axis >=3: val *= 1e-3 * 180/np.pi # convert to deg if offset is an angle

    # shift to min value, get psfs
    if axis > 0: ln.zSetSurfaceParameter(surf_num, axis,  nom_val + -1*val)
    if axis==0:  ln.zSetThickness(surf_num, nom_val + -1*val)
    ln.zPushLens()
    psf_cents_min = getpsf_orders_fibers(ln,orders,minwls,maxwls,NWAVE=NWAVE,fiber_locs=fiber_locs)
    if saveme: save_psf_cents('./testsavepsf_min.fits',psf_cents_min,NWAVE,NFIBER,NORDER)

    # shift to max value, get psfs
    if axis > 0: ln.zSetSurfaceParameter(surf_num, axis,  nom_val)
    if axis==0:  ln.zSetThickness(surf_num, nom_val)
    ln.zPushLens()
    psf_cents_max = getpsf_orders_fibers(ln,orders,minwls,maxwls,NWAVE=NWAVE,fiber_locs=fiber_locs)
    if saveme: save_psf_cents('./testsavepsf_max.fits',psf_cents_max,NWAVE,NFIBER,NORDER)
    
    #return to nominal before editing next param
    ln.zSetSurfaceParameter(surf_num, axis,  nom_val)
    if axis==0:
        ln.zSetThickness(surf_num, nom_val)
    else:
        ln.zSetSurfaceParameter(surf_num, axis,  nom_val)
    ln.zPushLens()

    if plotme: plot_psf_shifts(orders, psf_cents_min,psf_cents_max,name=surf_name,axis=axis,value=val)

    X,Y,U,V,Xshift,Yshift= get_psf_shifts(psf_cents_min,psf_cents_max)

    return X,Y,U,V,Xshift,Yshift



if __name__=='__main__':
    fiber_locs=[ -0.065, 0.065, 0.195]
    NFIBER=len(fiber_locs)
    try:
        ln
    except NameError:
        ln = pyz.createLink()

    orders,minwls,maxwls = load_orders_wl()
    NORDER = len(orders)
    surf_name, surf_num, axis,nom_val,val = load_tolerance_data(tol_filename='./tol_data_file.csv') # just use maxval bc minval is just neg right now

    # make sure all are set to nominal values before doing individual tweaks
    for itolfix in np.arange(len(surf_name)):
        if axis[itolfix]==0:
            ln.zSetThickness(surf_num[itolfix], nom_val[itolfix])
        else:
            ln.zSetSurfaceParameter(surf_num[itolfix], axis[itolfix],  nom_val[itolfix])
    ln.zPushLens()

    # save for each tolerance dimension
    #f = open('tol_rv_rate.csv','w')
    #f.write('itol,avg_rate,std_rate\n')
    for itol in np.arange(50,51):#,len(surf_name)):
        if surf_num[itol]<0:
            continue
        # shape out: norders, nwavelengths ( take zero index order - already doing just one order to start)
        print(itol)
        X1,Y1,RVerr1,V1,Xshift1,Yshift1 = get_one_shift(surf_num[itol], axis[itol],  val[itol], nom_val[itol], orders, minwls, maxwls,surf_name=surf_name[itol])
        X2,Y2,RVerr2,V2,Xshift2,Yshift2 = get_one_shift(surf_num[itol], axis[itol],  5*val[itol], nom_val[itol], orders, minwls, maxwls,surf_name=surf_name[itol])
        RVerr_rate = 1e-6 * (RVerr2 - RVerr1) / (5*val[itol] - val[itol]) #1e6 converts back to nano units
        avg_rate   = np.mean(np.abs(RVerr_rate))
        std_rate   = np.std(np.abs(RVerr_rate))
        jason_comp = np.sqrt(np.sum(Xshift1**2 + Yshift1**2)/(2*len(Xshift1.flatten()))) # mimic zemax error fxn
        plot_psf_shifts(X1,Y1,Xshift1,Yshift1,RVerr1,surf_name[itol],val[itol])
        #f.write('%s,%s,%s,%s,%s\n' %(surf_name[itol],itol,avg_rate,std_rate,jason_comp))

    #f.close()


