#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri April 24 05:54:18 2025

This code is a workaround for a composite ATL23 made up of combinations of 
dot_icefree_avg_albm, dot_avg_inice_7in12, and the open water DOT, dot_avg_albm.

@author: suzanne
"""

import glob, re
import h5py
import numpy as np
import matplotlib.pyplot as plt
from pyproj import CRS
from pyproj import Transformer
from scipy.ndimage import gaussian_filter
import datetime as dt

def gaussian_filter_nan(A, sigma):
    # Mask of valid values
    M = np.isfinite(A).astype(float)

    # Replace NaNs with 0 for convolution
    A0 = np.nan_to_num(A, nan=0.0)

    # Convolve data and mask
    num = gaussian_filter(A0, sigma=sigma)
    den = gaussian_filter(M, sigma=sigma)

    # Normalize
    out = num / den
    out[den == 0] = np.nan  # where no valid neighbors
    
    return out

path23 = '/Volumes/ATL12_Data/atl23_rel002'
a23files = sorted(glob.glob(f'{path23}/ATL23_*.h5'))

ymstr = [re.search(r"_(\d{6})", afile).group(1) for afile in a23files]
monstr = [dt.datetime(int(ymstr[mm][:4]),int(ymstr[mm][4:]),15).strftime('%b%Y') for mm in np.arange(len(ymstr))]

# establish empty list
dot23np = []

for ia,afile in enumerate(a23files):
    print(afile)
    print(re.search(r"_(\d{6})", afile).group(1))
    
    with h5py.File(afile,'r') as hf:
        
        if ia==0:
            # get projection
            gridx = hf['north_polar']['ds_grid_x'][:]
            gridy = hf['north_polar']['ds_grid_y'][:]
            xg23,yg23 = np.meshgrid(gridx,gridy)
            # lat = f['grid_lat'][:]
            # lon = f['grid_lon'][:]
            # print(f['crs'].attrs['proj4text'])
            proj4str = hf['north_polar']['crs'].attrs['proj4text'].decode("utf-8")
            crs23 = CRS.from_proj4(proj4str)
            transformer = Transformer.from_crs(crs23,"EPSG:4326")
            lat23,lon23 = transformer.transform(xg23,yg23)  # center lon/lat
            landmasknp = hf['north_polar']['landmask'][:]
            landmasknp = np.where(landmasknp==0,np.nan,1)
            # fig,ax = plt.subplots(1,1,figsize=(12,12),subplot_kw=dict(projection=ccrs.NorthPolarStereo(-45)))
            # ax.coastlines(zorder=3)
            # ax.gridlines(crs=ccrs.PlateCarree(),xlocs=np.arange(-180,180,45),ylocs=np.arange(70,90,5),color='gray',linewidth=0.3)
            # ax.set_extent([-180,180,68,90],crs=ccrs.PlateCarree())
            # pcm = ax.pcolormesh(lon23,lat23,landmasknp,cmap='jet',transform=ccrs.PlateCarree())
            # cbar = fig.colorbar(pcm,shrink=0.8)
            # ax.set_title(f'ATL23 landmask')

        tmpd = hf['north_polar']['dot_avg_albm'][:]
        fv = hf['north_polar']['dot_avg_albm'].attrs['_FillValue']
        tmpd = np.where(tmpd == fv,np.nan,tmpd)
        
        tmpi = hf['north_polar']['dot_icefree_avg_albm'][:]
        fv = hf['north_polar']['dot_icefree_avg_albm'].attrs['_FillValue']
        tmpi = np.where(tmpi == fv,np.nan,tmpi)
                
        tmp12_22 = hf['north_polar']['dot_avginice_7in12_albm'][:]
        fv = hf['north_polar']['dot_avginice_7in12_albm'].attrs['_FillValue']
        tmp12_22 = np.where(tmp12_22 == fv,np.nan,22) # 22 is a placeholder
                
        # where dot_icefree_avg_albm is valid, use that
        dot = tmpi   
        
        # where not dot_icefree_avg_albm, yet dot_avginice_7in12_albm is valid, use placeholder (here = 22)
        dot = np.where(np.isnan(dot),tmp12_22,dot)  
        
        # where nans still exist, use dot_avg_albm
        dot = np.where(np.isnan(dot),tmpd,dot)  
        
        # set placeholder (22) to nan
        dot[dot==22] = np.nan
 
        dot23np.append(dot)
    
dot23np = np.array(dot23np)        
       
# smooth
dot23npsm = []
for im,dot in enumerate(dot23np):
    dot23npsm.append(gaussian_filter_nan(dot23np[im],5))
dot23npsm *= landmasknp
prettydates = [dt.datetime.strftime(dt.datetime.strptime(ymstr[tt],'%Y%m'),'%b %Y') for tt in np.arange(len(ymstr))]
dates = [dt.datetime(int(tt[:4]),int(tt[4:]),15) for tt in ymstr]

    

