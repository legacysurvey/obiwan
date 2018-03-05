"""Adapted from Dustin's legacyanalysis/brick-summary.py

Makes just the plots
"""

from __future__ import print_function
import os
import sys
import fitsio
import numpy as np
from glob import glob
from collections import Counter

import matplotlib
matplotlib.use('Agg')
matplotlib.rc('text') #, usetex=True)
matplotlib.rc('font', family='serif')
from matplotlib.cm import get_cmap
import pylab as plt

from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage.morphology import binary_dilation

from astrometry.util.fits import fits_table,merge_tables
from astrometry.util.plotutils import antigray
from astrometry.util.starutil_numpy import lbtoradec
from astrometry.libkd.spherematch import match_radec

from tractor.sfd import SFDMap

def colorbar_axes(parent, frac=0.12, pad=0.03, aspect=20):
	pb = parent.get_position(original=True).frozen()
	# new parent box, padding, child box
	(pbnew, padbox, cbox) = pb.splitx(1.0-(frac+pad), 1.0-frac)
	cbox = cbox.anchored('C', cbox)
	parent.set_position(pbnew)
	parent.set_anchor((1.0, 0.5))
	cax = parent.get_figure().add_axes(cbox)
	cax.set_aspect(aspect, anchor=((0.0, 0.5)), adjustable='box')
	parent.get_figure().sca(parent)
	return cax

# From http://scipy-cookbook.readthedocs.io/items/Matplotlib_ColormapTransformations.html
def cmap_discretize(cmap, N):
    """Return a discrete colormap from the continuous colormap cmap.

        cmap: colormap instance, eg. cm.jet.
        N: number of colors.

    Example
        x = resize(arange(100), (5,100))
        djet = cmap_discretize(cm.jet, 5)
        imshow(x, cmap=djet)
    """
    if type(cmap) == str:
        cmap = get_cmap(cmap)
    colors_i = np.concatenate((np.linspace(0, 1., N), (0.,0.,0.,0.)))
    colors_rgba = cmap(colors_i)
    indices = np.linspace(0, 1., N+1)
    cdict = {}
    for ki,key in enumerate(('red','green','blue')):
        cdict[key] = [ (indices[i], colors_rgba[i-1,ki], colors_rgba[i,ki]) for i in range(N+1) ]
    # Return colormap object.
    return matplotlib.colors.LinearSegmentedColormap(cmap.name + "_%d"%N, cdict, 1024)



def plots(heatmap_table_fn, outdir='./',
          ramin=None,ramax=None,decmin=None,decmax=None,
          raticks=[]):
    T = fits_table(heatmap_table_fn)

    B = fits_table(os.path.join(os.environ['LEGACY_SURVEY_DIR'],
                                'survey-bricks.fits.gz'))
    print('Looking up brick bounds')
    ibrick = dict([(n,i) for i,n in enumerate(B.brickname)])
    bi = np.array([ibrick[n] for n in T.brickname])
    T.ra1 = B.ra1[bi]
    T.ra2 = B.ra2[bi]
    T.dec1 = B.dec1[bi]
    T.dec2 = B.dec2[bi]
    assert(np.all(T.ra2 > T.ra1))
    T.area = ((T.ra2 - T.ra1) * (T.dec2 - T.dec1) *
              np.cos(np.deg2rad((T.dec1 + T.dec2) / 2.)))
    del B
    del bi
    del ibrick

    print('Total sources:', sum(T.nobjs))
    print('Approx area:', len(T)/16., 'sq deg')
    print('Area:', np.sum(T.area))
    print('g,r,z coverage:', sum((T.nexp_g > 0) * (T.nexp_r > 0) * (T.nexp_z > 0)) / 16.)

    decam = True
    # vs MzLS+BASS
    #release = 'MzLS+BASS DR4'
    #release = 'DECaLS DR5'
    release = 'DECaLS DR3'
    
    if decam:
        # DECam
        #ax = [360, 0, -21, 36]
        ax = [300, -60, -21, 36]
        if ramin:
            ax[1]= ramin
        if ramax:
            ax[0]= ramax
        if decmin:
            ax[2]= decmin
        if decmax:
            ax[3]= decmax

        def map_ra(r):
                return r + (-360 * (r > 300))


        def map_ra(r):
                return r

    udec = np.unique(T.dec)
    print('Number of unique Dec values:', len(udec))
    print('Number of unique Dec values in range', ax[2],ax[3],':',
          np.sum((udec >= ax[2]) * (udec <= ax[3])))
        
    def radec_plot():
        plt.axis(ax)
        plt.xlabel('RA (deg)')
        if decam:
            # plt.xticks(np.arange(0, 361, 45))
            #tt = np.arange(0, 361, 60)
            #plt.xticks(tt, map_ra(tt))
            if len(raticks) > 0:
                plt.xticks(raticks)
            else:
                raise ValueError
                plt.xticks([-60,0,60,120,180,240,300], [300,0,60,120,180,240,300])
        else:
            plt.xticks(np.arange(90, 311, 30))

        plt.ylabel('Dec (deg)')

        def plot_broken(rr, dd, *args, **kwargs):
            dr = np.abs(np.diff(rr))
            I = np.flatnonzero(dr > 90)
            #print('breaks:', rr[I])
            #print('breaks:', rr[I+1])
            if len(I) == 0:
                plt.plot(rr, dd, *args, **kwargs)
                return
            for lo,hi in zip(np.append([0], I+1), np.append(I+1, -1)):
                #print('Cut:', lo, ':', hi, '->', rr[lo], rr[hi-1])
                plt.plot(rr[lo:hi], dd[lo:hi], *args, **kwargs)

        # Galactic plane lines
        gl = np.arange(361)
        gb = np.zeros_like(gl)
        rr,dd = lbtoradec(gl, gb)
        plot_broken(map_ra(rr), dd, 'k-', alpha=0.5, lw=1)
        rr,dd = lbtoradec(gl, gb+10)
        plot_broken(map_ra(rr), dd, 'k-', alpha=0.25, lw=1)
        rr,dd = lbtoradec(gl, gb-10)
        plot_broken(map_ra(rr), dd, 'k-', alpha=0.25, lw=1)
        
    plt.figure(1, figsize=(8,5))
    plt.subplots_adjust(left=0.1, right=0.98, top=0.93)

    plt.figure(2, figsize=(8,4))
    #plt.subplots_adjust(left=0.06, right=0.98, top=0.98)
    plt.subplots_adjust(left=0.08, right=0.98, top=0.98)
    plt.figure(1)
    
    # Map of the tile centers we want to observe...
    if decam:
        O = fits_table(os.path.join(os.environ['obiwan_code'],
                                    'svn_decam/obstatus/decam-tiles_obstatus.fits'))
    else:
        O = fits_table(os.path.join(os.environ['obiwan_code'],
                                    'svn_mosaic/obstatus/mosaic-tiles_obstatus.fits'))
    O.cut(O.in_desi == 1)
    rr,dd = np.meshgrid(np.linspace(ax[1],ax[0], 700),
                        np.linspace(ax[2],ax[3], 200))
    I,J,d = match_radec(O.ra, O.dec, rr.ravel(), dd.ravel(), 1.)
    desimap = np.zeros(rr.shape, bool)
    desimap.flat[J] = True

    # Smoothed DESI boundary contours
    C = plt.contour(gaussian_filter(
        binary_dilation(desimap).astype(np.float32), 2),
        [0.5], extent=[ax[1],ax[0],ax[2],ax[3]])
    plt.clf()
    desi_map_boundaries = C.collections[0]
    def desi_map_outline():
        segs = desi_map_boundaries.get_segments()
        for seg in segs:
            plt.plot(seg[:,0], seg[:,1], 'b-')
    
    def desi_map():
        # Show the DESI tile map in the background.
        plt.imshow(desimap, origin='lower', interpolation='nearest',
                   extent=[ax[1],ax[0],ax[2],ax[3]], aspect='auto',
                   cmap=antigray, vmax=8)

    base_cmap = 'viridis'

    # Dust map -- B&W version
    nr,nd = 610,350
    plt.figure(2)
    plt.clf()
    dmap = np.zeros((nd,nr))
    rr = np.linspace(ax[0], ax[1], nr)
    dd = np.linspace(ax[2], ax[3], nd)
    rr = rr[:-1] + 0.5*(rr[1]-rr[0])
    dd = dd[:-1] + 0.5*(dd[1]-dd[0])
    rr,dd = np.meshgrid(rr,dd)
    I,J,d = match_radec(rr.ravel(), dd.ravel(),
                        O.ra, O.dec, 1.0, nearest=True)
    iy,ix = np.unravel_index(I, rr.shape)
    #dmap[iy,ix] = O.ebv_med[J]
    sfd = SFDMap()
    ebv = sfd.ebv(rr[iy,ix], dd[iy,ix])
    dmap[iy,ix] = ebv
    mx = np.percentile(dmap[dmap > 0], 98)
    plt.imshow(dmap, extent=[ax[0],ax[1],ax[2],ax[3]], interpolation='nearest', origin='lower',
                   aspect='auto', cmap='Greys', vmin=0, vmax=mx)
    #desi_map_outline()
    radec_plot()
    cax = colorbar_axes(plt.gca(), frac=0.12)        
    cbar = plt.colorbar(cax=cax)
    cbar.set_label('Extinction E(B-V)')
    plt.savefig(os.path.join(outdir,'ext-bw.png'))
    plt.clf()
    dmap = sfd.ebv(rr.ravel(), dd.ravel()).reshape(rr.shape)
    plt.imshow(dmap, extent=[ax[0],ax[1],ax[2],ax[3]],
               interpolation='nearest', origin='lower',
               aspect='auto', cmap='Greys', vmin=0, vmax=0.25)
    desi_map_outline()
    radec_plot()
    cax = colorbar_axes(plt.gca(), frac=0.12)        
    cbar = plt.colorbar(cax=cax)
    cbar.set_label('Extinction E(B-V)')
    plt.savefig(os.path.join(outdir,'ext-bw-2.png'))
    plt.figure(1)

    #sys.exit(0)
    
    plt.clf()
    depthlo,depthhi = 21.5, 25.5
    for band in 'grz':
        depth = T.get('galdepth_%s' % band)
        ha = dict(histtype='step',  bins=50, range=(depthlo,depthhi))
        ccmap = dict(g='g', r='r', z='m')
        plt.hist(depth[depth>0], label='%s band' % band,
                 color=ccmap[band], **ha)
    plt.xlim(depthlo, depthhi)
    plt.xlabel('Galaxy depth (median per brick) (mag)')
    plt.ylabel('Number of Bricks')
    plt.title(release)
    plt.savefig(os.path.join(outdir,'galdepths.png'))

    for band in 'grz':
        depth = T.get('galdepth_%s' % band)
        nexp = T.get('nexp_%s' % band)
        #lo,hi = 22.0-0.05, 24.2+0.05
        lo,hi = depthlo-0.05, depthhi+0.05
        nbins = 1 + int((depthhi - depthlo) / 0.1)
        ha = dict(histtype='step',  bins=nbins, range=(lo,hi))
        ccmap = dict(g='g', r='r', z='m')
        area = 0.25**2
        plt.clf()
        I = np.flatnonzero((depth > 0) * (nexp == 1))
        plt.hist(depth[I], label='%s band, 1 exposure' % band,
                 color=ccmap[band], lw=1,
                 weights=area * np.ones_like(depth[I]),
                 **ha)
        I = np.flatnonzero((depth > 0) * (nexp == 2))
        plt.hist(depth[I], label='%s band, 2 exposures' % band,
                 color=ccmap[band], lw=2, alpha=0.5,
                 weights=area * np.ones_like(depth[I]),
                 **ha)
        I = np.flatnonzero((depth > 0) * (nexp >= 3))
        plt.hist(depth[I], label='%s band, 3+ exposures' % band,
                 color=ccmap[band], lw=3, alpha=0.3,
                 weights=area * np.ones_like(depth[I]),
                 **ha)
        plt.title('%s: galaxy depths, %s band' % (release, band))
        plt.xlabel('5-sigma galaxy depth (mag)')
        plt.ylabel('Square degrees')
        plt.xlim(lo, hi)
        plt.xticks(np.arange(depthlo, depthhi+0.01, 0.2))
        plt.legend(loc='upper right')
        plt.savefig(os.path.join(outdir,'depth-hist-%s.png' % band))

    for band in 'grz':
        plt.clf()
        desi_map()
        N = T.get('nexp_%s' % band)
        I = np.flatnonzero(N > 0)
        #cm = matplotlib.cm.get_cmap('jet', 6)
        #cm = matplotlib.cm.get_cmap('winter', 5)

        mx = 10
        cm = cmap_discretize(base_cmap, mx)
        plt.scatter(map_ra(T.ra[I]), T.dec[I], c=N[I], s=3,
                    edgecolors='none',
                    vmin=0.5, vmax=mx + 0.5, cmap=cm)
        radec_plot()
        cax = colorbar_axes(plt.gca(), frac=0.08)
        plt.colorbar(cax=cax, ticks=range(mx+1))
        plt.title('%s: Number of exposures in %s' % (release, band))
        plt.savefig(os.path.join(outdir,'nexp-%s.png' % band))

        #cmap = cmap_discretize(base_cmap, 15)
        cmap = cmap_discretize(base_cmap, 10)
        plt.clf()
        desi_map()
        psf = T.get('psfsize_%s' % band)
        I = np.flatnonzero(psf > 0)
        plt.scatter(map_ra(T.ra[I]), T.dec[I], c=psf[I], s=3,
                    edgecolors='none', cmap=cmap,
                    vmin=0.5, vmax=2.5)
        #vmin=0, vmax=3.)
        radec_plot()
        plt.colorbar()
        plt.title('%s: PSF size, band %s' % (release, band))
        plt.savefig(os.path.join(outdir,'psfsize-%s.png' % band))

        plt.clf()
        desi_map()

        depth = T.get('galdepth_%s' % band) - T.get('ext_%s' % band)
        mn,mx = np.percentile(depth[depth > 0], [10,98])
        mn = np.floor(mn * 10) / 10.
        mx = np.ceil(mx * 10) / 10.
        cmap = cmap_discretize(base_cmap, 1+int((mx-mn+0.001)/0.1))
        I = (depth > 0)
        plt.scatter(map_ra(T.ra[I]), T.dec[I], c=depth[I], s=3,
                    edgecolors='none', vmin=mn-0.05, vmax=mx+0.05, cmap=cmap)
        radec_plot()
        plt.colorbar()
        plt.title('%s: galaxy depth, band %s, median per brick, extinction-corrected' % (release, band))
        radec_plot()
        plt.savefig(os.path.join(outdir,'galdepth-%s.png' % band))

        # B&W version
        plt.figure(2)
        plt.clf()
        mn,mx = np.percentile(depth[depth > 0], [2,98])
        print('Raw mn,mx', mn,mx)
        mn = np.floor((mn+0.05) * 10) / 10. - 0.05
        mx = np.ceil( (mx-0.05) * 10) / 10. + 0.05
        print('rounded mn,mx', mn,mx)
        nsteps = int((mx-mn+0.001)/0.1)
        print('discretizing into', nsteps, 'colormap bins')
        #nsteps = 1+int((mx-mn+0.001)/0.1)
        cmap = cmap_discretize(antigray, nsteps)
        nr,nd = 610,228
        dmap = np.zeros((nd,nr))
        rr = np.linspace(ax[0], ax[1], nr)
        dd = np.linspace(ax[2], ax[3], nd)
        rr = rr[:-1] + 0.5*(rr[1]-rr[0])
        dd = dd[:-1] + 0.5*(dd[1]-dd[0])
        rr,dd = np.meshgrid(rr,dd)
        I,J,d = match_radec(rr.ravel(), dd.ravel(),
                            T.ra, T.dec, 0.2, nearest=True)
        iy,ix = np.unravel_index(I, rr.shape)
        dmap[iy,ix] = depth[J]
        plt.imshow(dmap, extent=[ax[0],ax[1],ax[2],ax[3]], interpolation='nearest', origin='lower',
                   aspect='auto', cmap=cmap, vmin=mn, vmax=mx)
        desi_map_outline()
        radec_plot()
        cax = colorbar_axes(plt.gca(), frac=0.12)        
        cbar = plt.colorbar(cax=cax, ticks=np.arange(20, 26, 0.5)) #ticks=np.arange(np.floor(mn/5.)*5., 0.1+np.ceil(mx/5.)*5, 0.2))
        cbar.set_label('Depth (5-sigma, galaxy profile, AB mag)')
        plt.savefig(os.path.join(outdir,'galdepth-bw-%s.png' % band))
        plt.figure(1)
        
        plt.clf()
        desi_map()
        ext = T.get('ext_%s' % band)
        mn = 0.
        mx = 0.5
        cmap = 'hot'
        cmap = cmap_discretize(cmap, 10)
        #cmap = cmap_discretize(base_cmap, 1+int((mx-mn+0.001)/0.1))
        plt.scatter(map_ra(T.ra), T.dec, c=ext, s=3,
                    edgecolors='none', vmin=mn, vmax=mx, cmap=cmap)
        radec_plot()
        plt.colorbar()
        plt.title('%s: extinction, band %s' % (release, band))
        plt.savefig(os.path.join(outdir,'ext-%s.png' % band))


    T.ngal = T.nsimp + T.nrex + T.nexp + T.ndev + T.ncomp
        
    for col in ['nobjs', 'npsf', 'nsimp', 'nrex', 'nexp', 'ndev', 'ncomp', 'ngal']:
        if not col in T.get_columns():
            continue
        plt.clf()
        desi_map()
        N = T.get(col) / T.area
        mx = np.percentile(N, 99.5)
        plt.scatter(map_ra(T.ra), T.dec, c=N, s=3,
                    edgecolors='none', vmin=0, vmax=mx)
        radec_plot()
        cbar = plt.colorbar()
        cbar.set_label('Objects per square degree')
        tt = 'of type %s' % col[1:]
        if col == 'nobjs':
            tt = 'total'
        plt.title('%s: Number of objects %s' % (release, tt))
        plt.savefig(os.path.join(outdir,'nobjs-%s.png' % col[1:]))

        # B&W version
        plt.figure(2)
        plt.clf()
        # plt.scatter(map_ra(T.ra), T.dec, c=N, s=3,
        #             edgecolors='none', vmin=0, vmax=mx, cmap=antigray)
        # Approximate pixel size in PNG plot
        # This doesn't work correctly -- we've already binned to brick resolution, so get moire patterns
        # nobjs,xe,ye = np.histogram2d(map_ra(T.ra), T.dec, weights=T.get(col),
        #                              bins=(nr,nd), range=((ax[1],ax[0]),(ax[2],ax[3])))
        # nobjs = nobjs.T
        # area = np.diff(xe)[np.newaxis,:] * (np.diff(ye) * np.cos(np.deg2rad(ye[:-1])))[:,np.newaxis]
        # nobjs /= area
        # plt.imshow(nobjs, extent=[ax[1],ax[0],ax[2],ax[3]], interpolation='nearest', origin='lower',
        #           aspect='auto')
        #print('Computing neighbours for nobjs plot...')
        nr,nd = 610,228
        nobjs = np.zeros((nd,nr))
        rr = np.linspace(ax[0], ax[1], nr)
        dd = np.linspace(ax[2], ax[3], nd)
        rr = rr[:-1] + 0.5*(rr[1]-rr[0])
        dd = dd[:-1] + 0.5*(dd[1]-dd[0])
        rr,dd = np.meshgrid(rr,dd)
        I,J,d = match_radec(rr.ravel(), dd.ravel(),
                            T.ra, T.dec, 0.2, nearest=True)
        iy,ix = np.unravel_index(I, rr.shape)
        nobjs[iy,ix] = T.get(col)[J] / T.area[J]
        #print('done')

        #mx = 2. * np.median(nobjs[nobjs > 0])
        mx = np.percentile(N, 99)

        plt.imshow(nobjs, extent=[ax[0],ax[1],ax[2],ax[3]], interpolation='nearest', origin='lower',
                   aspect='auto', cmap='Greys', vmin=0, vmax=mx)
        desi_map_outline()
        radec_plot()
        #cax = colorbar_axes(plt.gca(), frac=0.08)
        cax = colorbar_axes(plt.gca(), frac=0.12)        
        cbar = plt.colorbar(cax=cax,
                            format=matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
        cbar.set_label('Objects per square degree')
        plt.savefig(os.path.join(outdir,'nobjs-bw-%s.png' % col[1:]))
        #plt.savefig(os.path.join(outdir,'nobjs-bw-%s.png' % col[1:]))
        plt.figure(1)
        
    Ntot = T.nobjs
    for col in ['npsf', 'nsimp', 'nrex', 'nexp', 'ndev', 'ncomp', 'ngal']:
        if not col in T.get_columns():
            continue
        plt.clf()
        desi_map()
        N = T.get(col) / (Ntot.astype(np.float32))
        N[Ntot == 0] = 0.
        print(col, 'max frac:', N.max())
        mx = np.percentile(N, 99.5)
        print('mx', mx)
        plt.scatter(map_ra(T.ra), T.dec, c=N, s=3,
                    edgecolors='none', vmin=0, vmax=mx)
        radec_plot()
        plt.colorbar()
        plt.title('%s: Fraction of objects of type %s' % (release, col[1:]))
        plt.savefig(os.path.join(outdir,'fobjs-%s.png' % col[1:]))

        # B&W version
        plt.figure(2)
        plt.clf()
        #plt.scatter(map_ra(T.ra), T.dec, c=N * 100., s=3,
        #            edgecolors='none', vmin=0, vmax=mx*100., cmap=antigray)

        fobjs = np.zeros((nd,nr))
        rr = np.linspace(ax[0], ax[1], nr)
        dd = np.linspace(ax[2], ax[3], nd)
        rr = rr[:-1] + 0.5*(rr[1]-rr[0])
        dd = dd[:-1] + 0.5*(dd[1]-dd[0])
        rr,dd = np.meshgrid(rr,dd)
        I,J,d = match_radec(rr.ravel(), dd.ravel(),
                            T.ra, T.dec, 0.2, nearest=True)
        iy,ix = np.unravel_index(I, rr.shape)
        fobjs[iy,ix] = N[J] * 100.

        #mx = 2. * np.median(fobjs[fobjs > 0])
        mx = np.percentile(N * 100., 99)

        plt.imshow(fobjs, extent=[ax[0],ax[1],ax[2],ax[3]], interpolation='nearest', origin='lower',
                   aspect='auto', cmap='Greys', vmin=0, vmax=mx)

        desi_map_outline()
        radec_plot()
        cax = colorbar_axes(plt.gca(), frac=0.12)        
        cbar = plt.colorbar(cax=cax,
                            format=matplotlib.ticker.FuncFormatter(lambda x, p: '%.2g' % x))
        cbar.set_label('Percentage of objects of type %s' % col[1:].upper())
        plt.savefig(os.path.join(outdir,'fobjs-bw-%s.png' % col[1:]))
        #plt.savefig(os.path.join(outdir,'fobjs-bw-%s.png' % col[1:]))
        plt.figure(1)
        
    return 0
       


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--heatmap_table', type=str)
    parser.add_argument('--outdir', default='./', type=str)
    opt = parser.parse_args()

    plots(heatmap_table_fn= opt.heatmap_table, 
          outdir=opt.outdir,
          ramin=10,ramax=300,decmin=-10,decmax=28,
          raticks=[10,120,240,300])

if __name__ == '__main__':
    main()
