#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script that provides basic tile-creation for a directory of ATL03 and ATL08
track files

Copyright 2021 Applied Research Laboratories, University of Texas at Austin

This package is free software; the copyright holder gives unlimited
permission to copy and/or distribute, with or without modification, as
long as this notice is preserved.

Authors:
    Max Daniller-Varghese
    
Created on Mon Jan 25 09:37:09 2021

"""
#Basics
import os
import sys
import ogr
from math import ceil
import pandas as pd
import numpy as np
import datetime
import time
import geopandas as gpd
import geohash
import subprocess
# Laz generator function
import pylas
#ARL functions
from icesatIO import selectwkt
from icesatReader import get_atl03_struct, read_atl03_geolocation, get_atl08_struct, append_atl03_geolocation, get_attribute_info

# Depreciation warning supressor. Not a lts
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

pd.options.mode.chained_assignment = None  # default='warn'

#%
def file_sizer(filelist, path):
    # Loop and add files to list. This sorts the files in a list by size
    pairs = []
    for file in filelist:

        # Use join to get full file path.
        location = path + file
    
        # Get size and add to list of tuples.
        size = os.path.getsize(location)
        pairs.append((size, file))
    
    # Sort list of tuples by the first element, size.
    pairs.sort(key=lambda s: s[0])
    pairs=np.array(pairs)
    orderedlist=list(pairs[:,1])
    return orderedlist
#
def coordrounder(d, tile_size):
    # rounds coordinates to the increment of the size of the tiles being generated
    # it is robust for both the northern and southern hemisphere
    d[d>0] = tile_size * np.floor(d[d>0]/tile_size)
    d[d<0] = tile_size * np.ceil(d[d<0]/tile_size)
    return d


def geohash_ph(df,tile_size):
    # encodes the rounded photon coordinates to a geohashed region on earth for fast sorting
    start_time = time.time()
    # df['geokey_rounded'] = df.apply(lambda x: geohash.encode(latitude=coordrounder(x.lat_ph,tile_size),longitude=coordrounder(x.lon_ph,tile_size), precision=4), axis=1)
    df['lonround']=coordrounder(df.lon_ph.copy(),tile_size)
    df['latround']=coordrounder(df.lat_ph.copy(),tile_size)
    df["geokey_rounded"] = df.apply(lambda x: geohash.encode(latitude=x.latround, longitude=x.lonround,precision=4),axis=1)
    print('geohash_ph took: ' +str(time.time() - start_time))
    return df
    
def geohasher(grid):
    # encodes the corners of the grid that maps the tiles into geohashed regions
    # It requires an extremenly small amount of wiggle room to avoid errors at the poles
    dither = 0.00000001
    grid['NE'] = grid.geometry.bounds.apply(lambda x: geohash.encode(longitude=x.maxx-dither, 
                               latitude=x.maxy-dither, precision=4), axis=1)
    grid['NW'] = grid.geometry.bounds.apply(lambda x: geohash.encode(longitude=x.minx+dither, 
                               latitude=x.maxy-dither, precision=4), axis=1)
    grid['SW'] = grid.geometry.bounds.apply(lambda x: geohash.encode(longitude=x.minx+dither, 
                               latitude=x.miny+dither, precision=4), axis=1)
    grid['SE'] = grid.geometry.bounds.apply(lambda x: geohash.encode(longitude=x.maxx-dither, 
                               latitude=x.miny+dither, precision=4), axis=1)
    
    return grid
       
def log_writer(output_directory,file_name, idx, time):
    #logs the file name, it's index in the ordered list, and the time of initiation
    global d, row_out
    d={'file_name':[file_name], 'index':[idx], 'time':[time]}
    row_out = pd.DataFrame(data=d)
    if os.path.exists(output_directory+'log.csv')== True:
        row_out.to_csv(output_directory+'log.csv',mode='a',index=False, header=False)
    else:
        # create file and save it
        row_out.to_csv(output_directory+'log.csv',mode='w', index=False, header=True)   
        

def csv_writer(df_cell,putfile):
    #writes the dataframe into a csv. This format is bloated, but is included for completeness and flexibility
    # global df_cell_out, outputfile
    outputfile = putfile[:-3] + "csv"
    df_cell_out = df_cell[["lon_ph", "lat_ph", "h_ph","classification","signal_conf_ph","hemi","zone","geoid", "beam", "segment_snowcover", "segment_landcover", "file_name", "ph_h", "solar_elevation"]]
    if os.path.exists(outputfile) == True:
        # append the data to the file
        df_cell_out.to_csv(outputfile,mode='a',index=False, header=False)
    else:
        # create file and save it
        df_cell_out.to_csv(outputfile,mode='w', index=False, header=True)


def yrdoy(df):
    # create year and day of year field based on file name. Format is yyddd.
    # For example January 1, 2021 would be '21001' 
    datecode = pd.to_datetime(df.YYMMDD, format="%y%m%d")
    datecode_converted = datecode.dt.strftime("%y%j") 
    datecode_converted = datecode_converted.astype(int)
    return datecode_converted  


        
def write_gridded_las_pylas(df,outfilepath): 
    # write the gridded file into LAZ format. If you're having any trouble at all,
    # I recommned uninstalling pylas and lazrs, and "pip install pylas[lazrs]"
    # It's very sensitive to pylas and lazrs being built together.
    global wkt, new_vlr, proj,hemi,zone, datecode_converted, old_file
    global las, zz, delta_time, YYMMDD, track_num
    
    # read in the old file if it exists and append data together
    if os.path.exists(outfilepath[:-4] + '.laz')==True:
        print('readin')
        old_file = pylas.read(outfilepath[:-4] + '.laz')
        xx = np.concatenate([np.array(old_file.x),np.array(df.lon_ph)])
        yy = np.concatenate([np.array(old_file.y),np.array(df.lat_ph)])
        zz = np.concatenate([np.array(old_file.z),np.array(df.h_ph)])
        cc = np.concatenate([np.array(old_file.classification),np.array(df.classification)])
        delta_time = np.concatenate([np.array(old_file.delta_time, dtype='double'),np.array(df.delta_time, dtype='double')])
        YYMMDD = np.concatenate([np.array(old_file.YYMMDD, dtype='int'),np.array(df.YYMMDD, dtype='int')])
        track_num = np.concatenate([np.array(old_file.track_num, dtype='int'),np.array(df.track_num, dtype='int')])

        print('concat')
        gg = np.concatenate([np.array(old_file.geoid),np.array(df.geoid)])
        bb = np.concatenate([np.array(old_file.beam_number),np.array(df.beam)]) 
        sf = np.concatenate([np.array(old_file.snow_flag),np.array(df.segment_snowcover)])
        lf = np.concatenate([np.array(old_file.landcover),np.array(df.segment_landcover)])
        mf = np.concatenate([np.array(old_file.msw_flag),np.array(df.msw_flag)])
        datecode_converted = np.concatenate([np.array(old_file.yrdoy),np.array(yrdoy(df), dtype='int')])
        zc = np.concatenate([np.array(old_file.ph_h),np.array(df.ph_h)])
        se = np.concatenate([np.array(old_file.solar_elevation),np.array(df.solar_elevation)])
        sc= np.concatenate([np.array(old_file.signal_conf),np.array(df.signal_conf_ph)])
        
    else:
        xx = np.array(df.lon_ph)
        yy = np.array(df.lat_ph)
        zz = np.array(df.h_ph)
        cc = np.array(df.classification)
        delta_time = np.array(df.delta_time, dtype='double')
        YYMMDD = np.array(df.YYMMDD, dtype='double')
        track_num = np.array(df.track_num, dtype='double')
        gg = np.array(df.geoid)
        bb = np.array(df.beam) 
        sf = np.array(df.segment_snowcover)
        lf = np.array(df.segment_landcover)
        mf = np.array(df.msw_flag)
        datecode_converted = np.array(yrdoy(df), dtype='int')
        zc = np.array(df.ph_h)
        se = np.array(df.solar_elevation)
        sc= np.array(df.signal_conf_ph)
        
    hemi = df.loc[df.hemi.first_valid_index()].hemi
    zone = df.loc[df.zone.first_valid_index()].zone    
    
    # Set projection and write .las file
    proj = 'wgs84'    
    wkt = selectwkt(proj,hemi,zone) #well-known-text
    
    #Create new VLR
    new_vlr = pylas.vlrs.VLR(user_id = "LASF_Projection",
                                record_id = 2112,
                               description = "OGC Coordinate System WKT")
    new_vlr.record_data = wkt
    
    #Create new Header
    las = pylas.create(point_format_id=6,file_version="1.4")
    
    las.vlrs.append(new_vlr)
    las.header.global_encoding.wkt = 1
    las.header.x_scale = 0.000001
    las.header.y_scale = 0.000001
    las.header.z_scale = 0.000001

    #Add extra data fields
    names=["geoid", "beam_number", "snow_flag", "landcover", "yrdoy", "ph_h","solar_elevation", "signal_conf", "delta_time", "YYMMDD" , "track_num", "msw_flag"]
    data_types= ["double","uint8","uint8","uint8","uint32","double","double", "uint8", "double", "uint32","uint32", "uint8"]
    descriptions = ["Geoid (m)", "Beam Descriptor", "Modis Snowcover Flag", "Landcover class", "Year and Day of Year", "Canopy Height (m)", "Solar Elevation", "Photon Signal Confidence", "Delta Time", "YYMMDD",  "ATL Track Number", "Multiple Scatter Warning Flag"]

    for idx, nm in enumerate(names):
        las.add_extra_dim(name=nm,
                type=data_types[idx],
                description=descriptions[idx])
    
    # Write data
    las.x = xx
    las.y = yy
    las.z = zz
    las.classification = cc
    las.geoid = gg
    las.beam_number = bb
    las.snow_flag = sf
    las.landcover = lf
    las.msw_flag = mf
    las.yrdoy = datecode_converted
    las.ph_h = zc
    las.solar_elevation = se
    las.signal_conf = sc
    las.delta_time = delta_time
    las.YYMMDD = YYMMDD
    las.track_num = track_num

    las.write(outfilepath[:-4] + '.laz',do_compress=True)

    return las    


def grid_file_namer(grid, releasenum, cyclenum):
    #Assign a file name for each tile overlapping the grid
    global cell_name
    st=time.time()
    minx = grid.geometry.bounds.minx
    miny = grid.geometry.bounds.miny
    maxx = grid.geometry.bounds.maxx
    maxy = grid.geometry.bounds.maxy
    
    cell_name = []
    for ii in range(len(grid)):
        cell_name.append('ATL_tile_'+'r'+str(releasenum) +'_cycle' +str(cyclenum) + '_' + str(round(maxy[ii]))+'N'+str(round(maxx[ii]))+'E_'+str(round(miny[ii]))+'N'+str(round(minx[ii]))+'E_'+str(round(maxy[ii]-miny[ii])) + '.laz')
    
    grid = grid.assign(cell_name=cell_name)
    print('grid takes ' + str(time.time()-st))
    return grid


def track_to_grid(df, output_directory, shapefilenamepath):
    #Determine the overlap of the ATL track on the grid
    global grid, gdf, union, df_cell, grid_file_name

    proj = 'wgs84'
    releasenum = df.loc[df.releasenum.first_valid_index()].releasenum
    cyclenum = df.loc[df.cyclenum.first_valid_index()].cyclenum
    grid = gpd.read_file(shapefilenamepath)
    grid = grid_file_namer(grid, releasenum, cyclenum)
    grid.crs = {'init' :'epsg:4326'}
    grid=geohasher(grid)

    ##%% Union and identification of grid
    start_time = time.time()
    tile_size = grid.geometry.bounds.maxx[0] - grid.geometry.bounds.minx[0]
    union = geohash_ph(df,tile_size)
    print("postoverlay took " + str(time.time() - start_time))
    start_time = time.time()

    for ii in np.unique(union.geokey_rounded): # Cycle through each grid cell in the dataframe
        subgrid = grid[grid.SW == ii]
        if subgrid.empty == False:
            grid_file_name = grid[grid.SW == ii].cell_name.to_string(index=False).strip()
            df_cell = union.loc[union.geokey_rounded == ii]
            print("postunion " + grid_file_name)
            # csv_writer(df_cell, output_directory + grid_file_name)
            write_gridded_las_pylas(df_cell,output_directory + grid_file_name)
        print("postunion took " + str(time.time() - start_time))

    
def append_atl03_atl08features(heights_i, segment_landcover):
    #Add atl08 features to atl03
    # global expanded_df
    st=time.time()
    
    segment_landcover['segment_id'] = segment_landcover.apply(lambda x:  np.arange(x.segment_id_beg, x.segment_id_end+1), axis=1)
    expanded_df = segment_landcover.explode('segment_id')
    heights_expanded = pd.merge(heights_i, expanded_df, how='outer', left_on='segment_id', right_on='segment_id')
    print('join took' +str(time.time()-st))
    return heights_expanded    


def read_track(dir_03, dir_08, shapefilenamepath, output_directory):
    #Read the track in, add atl08 features, and process it into a LAZ tile
    global heights_i, heights, atl03, atl08, canopy08, list_of_atl03_files, info, file_name
    from icesatReader import get_atl03_struct, read_atl03_geolocation, get_atl08_struct, append_atl03_geolocation, get_attribute_info, get_atl08canopy_struct, read_atl08_canopy_h, geoid_corrector
        
    list_in_atl03_path = os.listdir(dir_03)
    list_of_atl03_files = [f for f in list_in_atl03_path if f.startswith('ATL03') & f.endswith('h5')]    
    orderedlist_03 = file_sizer(list_of_atl03_files, dir_03) # order by size
    
    list_in_atl08_path = os.listdir(dir_08)
    list_of_atl08_files = [f for f in list_in_atl08_path if f.startswith('ATL08') & f.endswith('h5')]
    
    for idx, file_name in enumerate(orderedlist_03, start=0):
        log_writer(output_directory,file_name, idx, time.time())
        atl03filepath = dir_03 + file_name
        atl08filepath = dir_08 + str.replace(file_name, 'ATL03', 'ATL08')
        gt_all = ['gt1r', 'gt1l', 'gt2r', 'gt2l', 'gt3r', 'gt3l']

        for gt in gt_all:
            info = get_attribute_info(atl03filepath,gt)
            try:
                atl03 = get_atl03_struct(atl03filepath, gt, atl08filepath)                    
                geolocation = read_atl03_geolocation(atl03filepath, gt)
                    
                atl08 = get_atl08_struct(atl08filepath, gt, atl03)
                canopy08 = read_atl08_canopy_h(atl08filepath,gt)
                canopy08.columns = ['classed_pc_flag', 'classed_pc_indx', 'd_flag', 'delt_time', 'ph_h','ph_segment_id']
                segment_landcover = atl08.df[['segment_landcover','segment_snowcover','segment_id_beg','segment_id_end', 'dem_h','solar_elevation', 'msw_flag']]
                heights_i = atl03.df
                heights_i = append_atl03_geolocation(heights_i, geolocation, 
                                                    fields = ['segment_id'])
                heights_i = heights_i[heights_i.classification>0]
                heights_i = append_atl03_atl08features(heights_i, segment_landcover) 
                heights_i[heights_i.geoid>10e10] = np.nan
                
                heights_i = heights_i.assign(beam=info['atlas_spot_number'], hemi=atl03.hemi, 
                                             zone=atl03.zone, epsg=atl03.epsg, 
                                             file_name=atl03.atlFileName, 
                                             releasenum=atl03.releaseNum,
                                             cyclenum = atl03.unknown[0:2],
                                             YYMMDD= atl03.atlFileName.split('_')[1][2:8],
                track_num = atl03.atlFileName.split('_')[2][0:4])
                
                heights_i = heights_i.join(canopy08,how='inner',lsuffix='delta_time',rsuffix='delt_time')
                nabool = pd.isna(heights_i.classification)
                heights_i = heights_i[~nabool]                
                heights_i[heights_i.ph_h>10e10] = np.nan
                
                #TODO Chunker for extremely large input files. Currently breaks in loop
                # from more_itertools import chunked
                # CHUNK_SIZE = 100000 # variable for available computer memory?
                # index_chunks = chunked(heights_i.index, CHUNK_SIZE)
                # for ii in index_chunks:
                #     track_to_grid(heights_i.iloc[ii], output_directory, shapefilenamepath)
                track_to_grid(heights_i, output_directory, shapefilenamepath)
                print('Success on ' + file_name + ' beam ' + gt+ 'index: ' +str(idx))
                
            except:
                print('Failure on ' + file_name + ' beam ' + gt + 'index: ' +str(idx))
                
    return

def grid_maker(outputGridfn,xmin,xmax,ymin,ymax,gridHeight,gridWidth):
    #Make the grid that defines the tile edges
    from osgeo import osr, gdal, ogr
    # convert sys.argv to float
    xmin = float(xmin)
    xmax = float(xmax)
    ymin = float(ymin)
    ymax = float(ymax)
    gridWidth = float(gridWidth)
    gridHeight = float(gridHeight)

    # get rows
    rows = ceil((ymax-ymin)/gridHeight)
    # get columns
    cols = ceil((xmax-xmin)/gridWidth)

    # start grid cell envelope
    ringXleftOrigin = xmin
    ringXrightOrigin = xmin + gridWidth
    ringYtopOrigin = ymax
    ringYbottomOrigin = ymax-gridHeight

    # create output file
    outDriver = ogr.GetDriverByName('ESRI Shapefile')
    if os.path.exists(outputGridfn):
        os.remove(outputGridfn)
    outDataSource = outDriver.CreateDataSource(outputGridfn)
    outLayer = outDataSource.CreateLayer(outputGridfn,geom_type=ogr.wkbPolygon )
    featureDefn = outLayer.GetLayerDefn()

    # create grid cells
    countcols = 0
    while countcols < cols:
        countcols += 1

        # reset envelope for rows
        ringYtop = ringYtopOrigin
        ringYbottom =ringYbottomOrigin
        countrows = 0

        while countrows < rows:
            countrows += 1
            ring = ogr.Geometry(ogr.wkbLinearRing)
            ring.AddPoint(ringXleftOrigin, ringYtop)
            ring.AddPoint(ringXrightOrigin, ringYtop)
            ring.AddPoint(ringXrightOrigin, ringYbottom)
            ring.AddPoint(ringXleftOrigin, ringYbottom)
            ring.AddPoint(ringXleftOrigin, ringYtop)
            poly = ogr.Geometry(ogr.wkbPolygon)
            poly.AddGeometry(ring)

            # add new geom to layer
            outFeature = ogr.Feature(featureDefn)
            outFeature.SetGeometry(poly)
            outLayer.CreateFeature(outFeature)
            outFeature = None

            # new envelope for next poly
            ringYtop = ringYtop - gridHeight
            ringYbottom = ringYbottom - gridHeight

        # new envelope for next poly
        ringXleftOrigin = ringXleftOrigin + gridWidth
        ringXrightOrigin = ringXrightOrigin + gridWidth

    # Save and close DataSources
    outDataSource = None

def main(dir_03, dir_08, dir_out, shapefilenamepath, xmin, xmax, ymin, ymax, gridHeight, gridWidth):
    # Make a gridded shapefile based on user specified 
    print("working...")
    grid = grid_maker(shapefilenamepath,xmin,xmax,ymin,ymax,gridHeight,gridWidth)
    print("made the grid")
    read_track(dir_03, dir_08, shapefilenamepath, dir_out)
    print("read the track")
        
    print("finished!")
#%

if __name__ == "__main__":
        
    inpFile = 'example_local.inp'
    with open(inpFile) as inputFile:
        exec(inputFile.read(), globals())
        main(dir_03, dir_08, dir_out, shapefilenamepath, xmin, xmax, ymin, ymax, gridHeight, gridWidth)
