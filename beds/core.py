#!/usr/bin/env python

# inst: university of bristol
# auth: jeison sosa
# mail: sosa.jeison@gmail.com / j.sosa@bristol.ac.uk

import os
import numpy as np
import pandas as pd
import gdalutils as gu
import xarray as xr
from glob import glob
from subprocess import call
from lfptools.buildmodel import write_bci
from lfptools.buildmodel import write_bdy
from lfptools.buildmodel import write_evap
from lfptools.buildmodel import write_par


def beds(widthtif, bnkfixtif, runoffcsv, date1, date2, bedtif, lisfloodfp):

    # Create a temp temporal work folder
    outfolder = os.path.dirname(bedtif) + '/beds-temp/'
    try:
        os.makedirs(outfolder + 'lfp/nc/')
    except FileExistsError:
        pass

    # Determine end of the simulation, how many days
    t = (pd.to_datetime(date2, format='%Y-%m-%d') -
         pd.to_datetime(date1, format='%Y-%m-%d')).days + 1

    # Create 1D DEM, synthetic
    demtif = outfolder + 'dem1d.tif'
    wdt = gu.get_data(widthtif)
    geo = gu.get_geo(widthtif)
    dem = np.where(wdt > 0, 10000, 0)
    gu.write_raster(dem, demtif, geo, 'Int16', 0)

    # Convert input files to ASCII
    widthasc = outfolder + 'width.asc'
    call(['gdal_translate',
          '-of', 'AAIGRID',
          widthtif, widthasc])

    demasc = outfolder + 'dem.asc'
    call(['gdal_translate',
          '-of', 'AAIGRID',
          demtif, demasc])

    bnkfixasc = outfolder + 'bnkfix.asc'
    call(['gdal_translate',
          '-of', 'AAIGRID',
          bnkfixtif, bnkfixasc])

    # Write LISFLOOD-FP files
    bcilfp = outfolder + 'lfp.bci'
    write_bci(bcilfp, runoffcsv)

    bdylfp = outfolder + 'lfp.bdy'
    write_bdy(bdylfp, runoffcsv, t)

    evaplfp = outfolder + 'lfp.evap'
    write_evap(evaplfp, t)

    parlfp = outfolder + 'lfp.par'
    write_par(parlfp=parlfp,
              bcilfp=bcilfp,
              bdylfp=bdylfp,
              evaplfp=evaplfp,
              gaugelfp='./none',
              stagelfp='./none',
              dembnktif=demasc,
              wdttif=widthasc,
              bedtif=bnkfixasc,
              t=t)

    # Run simulation
    call([lisfloodfp, '-v', 'lfp.par'], cwd=outfolder)

    # Write netCDFs for WATER DEPTHS
    myfiles = sorted(glob(outfolder + '/lfp/*.wd'))
    for myfile in myfiles:
        fname = outfolder + 'lfp/nc/' + os.path.basename(myfile) + '.nc'
        xr.open_rasterio(myfile).to_dataset(name='myvar').to_netcdf(fname)

    # Read netCDFs
    ds = xr.open_mfdataset(outfolder + 'lfp/nc/*.nc',
                           concat_dim='band',
                           autoclose=True,
                           parallel=False,
                           chunks={'band': 10})

    # Calculating mean
    method = ds.where(ds > 0, 0).myvar.mean('band')

    # Saving result in netCDF
    method.to_netcdf(outfolder + 'mean.nc')

    # Reading banks
    bnkfix = gu.get_data(bnkfixtif)
    bnkfix = np.where(bnkfix > 0, bnkfix, 0)

    # Calculating bed
    bed = bnkfix - method.compute().data

    # Write final raster
    gu.write_raster(bed, bedtif, geo, 'Float64', 0)
