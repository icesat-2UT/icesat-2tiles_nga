# -*- coding: utf-8 -*-
"""
Script that contains utliity functions for ICESat-2 Mega Tiles

Copyright 2021 Applied Research Laboratories, University of Texas at Austin

This package is free software; the copyright holder gives unlimited
permission to copy and/or distribute, with or without modification, as
long as this notice is preserved.

Authors:
    Mike Alonzo
    
Date: June 2021
"""

# Import Python modules
import numpy as np
import sys
import pyproj as proj
import os
import glob
import pandas as pd
import pylas # pip install lazrs for .laz files
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as mcol
import warnings
import time
import simplekml
from osgeo import gdal, ogr

#from mpl_toolkits.basemap import Basemap

verNum = 'v1.0'

EPSG_ARCTIC = '3413'
EPSG_ANTARCTIC = '3976'
EPSG_ERR = '-1'
DEG_ARCTIC = 84.0
DEG_ANTARCTIC = -80.0

#EPSG_CEA = '6933' # Lambert Cylindrical Equal Area
#EPSG_LAEA_N = '6931' # Lambert Azimuthal Equal Area
#EPSG_LAEA_S = '6932'
#GT_NUMS = ['gt1r', 'gt2r', 'gt3r', 'gt1l', 'gt2l', 'gt3l']


# Declare global variables
class globalVars(): 
    
    # Set method definitions and usages
    homeDef = 'To return to this home window'
    exitDef = "To exit this program (can also type 'quit' or 'stop')"
    helpDef = 'To show help info for any method below'
    helpUse = 'Usage: help methodName'
    helpExp = 'Example: help process'
    listDef = 'To list files with a specified extension'
    listUse = 'Usage: list extName'
    listExp = 'Example: list inp'
    loadDef = 'To load previously processed Mega Tile data from .gzip file'
    loadUse = 'Usage: load processedFileName.gzip'
    loadExp = 'Example: load file1.gzip'
    processDef = 'To process ICESat-2 Mega Tile .laz files from .inp file'
    processUse = 'Usage: process inputFileName.inp'
    processExp = 'Example: process test.inp'
    filterDef = 'To filter processed ICESat-2 Mega Tile data'
    filterUse = 'Usage: filter colNameToFilter filterOperation filterValue'
    filterExp = 'Example: filter classification == 1 OR filter z_hae > 20'
    computeDef = 'To compute stats on processed ICESat-2 Mega Tile data'
    computeUse = 'Usage: compute -g groupByColumn -c columnToOperateOn -o operationToPerform'
    computeExp = "Example: compute -g classification -c ['z_hae'] -o mean"
    plotDef = 'To plot processed ICESat-2 Mega Tile data'
    plotUse = 'Usage: plot plotTypeName'
    plotExp = 'Example: plot xy'
    exportDef = 'To export data to .gzip, .csv, .laz, .kml, or .tif files'
    exportUse = 'Usage: export fileName.extName'
    exportExp = 'Example: export fileName.csv OR export fileName.laz latlon'

# endClass


# Class to define global color scheme
class globalColors():
    
    # Set colors
    brown = np.array([[ 0.6275, 0.3216, 0.1765]])
    darkGreen = np.array([[0, 0.5, 0]])
    liteGreen = np.array([[0, 1.0, 0]])
    
# endClass
    

# Object for gridMetricNew function
class GridStruct:
    def __init__(self, x, y, grid, time):
        self.x = x
        self.y = y
        self.grid = grid
        self.t = time
    # endDef
# endClass
    
    
# Function to get lat filter
def getLatFilter(inpLat):
    
    if(inpLat>=0):
        outFilter = str(inpLat) + 'N'
    else:
        outFilter = str(inpLat) + 'S'
    # endIf
    
    return outFilter

# endDef
    
# Function to get lon filter
def getLonFilter(inpLon):
    
    if(inpLon>=0):
        outFilter = str(inpLon) + 'E'
    else:
        outFilter = str(inpLon) + 'W'
    # endIf
    
    return outFilter

# endDef    
    

# Find UTM zone for individual lon/lat.
def find_utm_zone_point(lon, lat, module=None, m=None):

    """
    Input:
    longitude, latitude, {opt} module, {opt} m


    Output:
    epsg_code string


    Desc:
        longitude, latitude - in degrees
            No checks on lon/lat bounds

        module - options: [None, 'mgrs', 'pygeodesy']
            optional
            If nothing is given, it defaults to None
            If a wrong module is given, it defaults to None condition
            If mgrs or pygeodesy are given, it uses the 
            respective module

        m - mgrs object, i.e. m = mgrs.MGRS()
            optional input to allow user to predefine
            MGRS object since it doesn't have to be
            defined for every loop

            If it defined as anything other than
            mgrs.MGRS(), code will throw an error


    Latitude Bounds
    UTM exists between -80 and 84 deg latitude. UPS or
    NSIDC Sea Ice Polar Stereographic exists either
    south of -80 deg or north of 84 deg. The exact
    latitude bounds differ by method, as shown below
    in interval set notation:

    if module == None:
        antarctic:  [-90, -80]
        UTM:        (-80, 84)
        arctic:     [84,90]
        Does not take Norway/Svalbard UTM zone changes
        into account

    if module == 'mgrs' or module == 'pygeodesy':
        antarctic:  [-90, -80)
        UTM:        [-80, 84)
        arctic:     [84,90]
        Includes all one-offs in UTM zones
        

    Precision
        None has a precision of round-off
        mgrs has a precision of 1m
        pygeodesy has a precision of round-off

    """


    if not (module == 'mgrs' or module == 'pygeodesy'): # None condition

        utm_band = str(int((np.floor((lon + 180) / 6 ) % 60) + 1))
        if len(utm_band) == 1:
            utm_band = '0'+utm_band
            # print(utm_band)
        if lat >= 0 and lat < 84:
            epsg_code = '326' + utm_band
        elif lat >= 84:
            # epsg_code = '3411'
            epsg_code = EPSG_ARCTIC
        elif lat <= -80:
            # epsg_code = '3412'
            epsg_code = EPSG_ANTARCTIC
        else:
            epsg_code = '327' + utm_band
        return epsg_code


    elif module == 'mgrs':
        import mgrs
        if m == None:
            m = mgrs.MGRS()

        # UTM defined in [-80, 84)
        if DEG_ANTARCTIC <= lat < DEG_ARCTIC:
            # [-80,84)
            mgrs_code = m.toMGRS(lat, lon, MGRSPrecision=5) # 1m precision
            UTM_code = m.MGRSToUTM(mgrs_code)
            utm_band = str(UTM_code[0]).zfill(2)
            hemi = UTM_code[1].decode() # 'N' or 'S'

            if not (hemi == 'N' or hemi == 'S'):
                # debug check since code is new
                print('warning ({}): hemi={} is something other than N or S'.format(module, hemi))

            epsg_code = '326' + utm_band # N
            if hemi == 'S':
                epsg_code = '327' + utm_band

        elif lat < DEG_ANTARCTIC:
            # [-90, -80), polar antarctic
            epsg_code = EPSG_ANTARCTIC

        elif lat >= DEG_ARCTIC:
            # [84, 90], polar arctic
            epsg_code = EPSG_ARCTIC

        return epsg_code


    elif module == 'pygeodesy':
        import pygeodesy as pg

        if DEG_ANTARCTIC <= lat < DEG_ARCTIC:
             # [-80,84)
            z = pg.utm.utmZoneBand5(lat, lon)
            utm_band = str(z.zone).zfill(2)
            hemi = z.hemipole # 'N' or 'S'

            if not (hemi == 'N' or hemi == 'S'):
                # debug check since code is new
                print('warning ({}): hemi={} is something other than N or S'.format(module, hemi))

            epsg_code = '326' + utm_band # N
            if hemi == 'S':
                epsg_code = '327' + utm_band


        elif lat < DEG_ANTARCTIC:
            # [-90, -80), polar antarctic
            epsg_code = EPSG_ANTARCTIC

        elif lat >= DEG_ARCTIC:
            # [84, 90], polar arctic
            epsg_code = EPSG_ARCTIC

        return epsg_code

# endDef


# Find UTM zone for numpy array of lon/lat.
def find_utm_zone_arr(lon, lat, mode=True, module=None, m=None):

    """
    Input:
    longitude, latitude, {opt} mode, {opt} module, {opt} m
    

    Output:
    if mode:
        return a single epsg_code
    else:
        return a list of epsg_codes relative to lon/lat


    Desc:
        longitude, latitude - in degrees
            No checks on lon/lat bounds
            These must be the same length

        mode - whether to return the mode of zones or not
            default is True

        module - options: [None, 'mgrs', 'pygeodesy']
            optional
            If nothing is given, it defaults to None
            If a wrong module is given, it defaults to None condition
            If mgrs or pygeodesy are given, it uses the 
            respective module

        m - mgrs object, i.e. m = mgrs.MGRS()
            optional input to allow user to predefine
            MGRS object since it doesn't have to be
            defined for every loop

            If it defined as anything other than
            mgrs.MGRS(), code will throw an error

    See find_utm_zone_point() for additional details.

    """
    default_mod_bool = not (module == 'mgrs' or module == 'pygeodesy')

    if type(lon) != np.ndarray:
        lon = np.array(lon)
    if type(lat) != np.ndarray:
        lat = np.array(lat)

    if lon.size == 1:
        lon = np.array([lon])
    if lat.size == 1:
        lat = np.array([lat])

    if len(lon) != len(lat):
        print('error: len(lon) != len(lat)')
        return EPSG_ERR
        
    if default_mod_bool and mode:
        # OG way to find zone via the mode of all points

        arctic = len(lat[lat >= 84.0])
        nhem = len(lat[(lat < 84.0) & (lat >= 0)])
        shem = len(lat[(lat > -80.0) & (lat < 0)])
        antarctic = len(lat[lat < -80.0])
        if arctic > nhem: #Arctic case
            epsg_code = EPSG_ARCTIC # '3413'
        elif antarctic > shem: #Antarctic case
            epsg_code = EPSG_ANTARCTIC #'3976'
        else: #If not Arctic or Antarctic it is within UTM Zone
            tz = np.floor((((lon + 180) / 6 ) % 60) + 1).astype(int)
            zone = np.unique(tz)
            if len(zone) == 1:
                zone = zone[0]
            elif len(zone) == 2:
                z1 = len(tz[tz == zone[0]])
                z2 = len(tz[tz == zone[1]])
                if z1 > z2:
                    zone = zone[0]
                else:
                    zone = zone[1]
            elif len(zone) == 3:
                zone = zone[1]
            elif len(zone) > 3:
                from scipy.stats import mode as mode_func
                zone = mode_func(zone)[0][0]
                print("Warning: Input ground track present in more than 3 UTM zones. \
                         \nRecommend manually selecting GCS.")
            else:
                # len(zone) == 0
                print("Warning: zone == [], lon/lat may not have values")
                sys.exit()

            if nhem >= shem:
                if len(str(zone)) == 1:
                    zone = "0" + str(zone)
                epsg_code = '326' + str(zone)
            else:
                if len(str(zone)) == 1:
                    zone = "0" + str(zone)
                epsg_code = '327' + str(zone)
        return epsg_code

    else:
        # OG method and mode or
        # OG method and not mode or
        # different method and mode or
        # different method and not mode

        epsg_code_arr = []
        for i in range(len(lon)):
            out = find_utm_zone_point(lon[i], lat[i], module=module, m=m)
            epsg_code_arr.append(out)

        if mode:
            from scipy.stats import mode as mode_func
            epsg_code = mode_func(epsg_code_arr)[0][0]
            return epsg_code

        else:
            return epsg_code_arr

# endDef
            

# Transform GCS/PCS based on EPSG and x/y. 
def transform(epsg_in, epsg_out, x, y, use_old_version=False):
    try:
        # Using pyproj version 2 and above
        # https://pyproj4.github.io/pyproj/stable/gotchas.html#upgrading-to-pyproj-2-from-pyproj-1
        transformer = proj.Transformer.from_crs(epsg_in, epsg_out, always_xy=True)
        xx, yy = transformer.transform(x, y)
    except:
        print('PYPROJ failed, attemping with GDAL...')
    # endTry
    
    return xx,yy
    # endIf
# endDef
    
    
# Transform from lon/lat to given EPSG. 
def wgs84_to_epsg_transform(epsg_out, lon, lat):
    epsg_in = 'epsg:4326'
    epsg_out = ('epsg:{0}'.format(str(epsg_out)))
    xx, yy = transform(epsg_in, epsg_out, lon, lat)
    return xx,yy

# endDef
    

def identifyEPSG(hemi,zone):
    if hemi == 'N':
        outstring = 'epsg:326'
    elif hemi == "S":
        outstring = 'epsg:327'
    else:
        outstring = ''
        print('error')
    
    outstring = outstring + (str(zone).zfill(2))
    
    return outstring
# endDef
    

def identify_hemi_zone(epsg):
    epsg = str(epsg)
    epsg = epsg.split(':')[-1]
    
    if (epsg[0:3] == '326'):
        hemi = 'N'
        zone = epsg[3:5]
    elif (epsg[0:3] == '327'):
        hemi = 'N'
        zone = epsg[3:5]
    elif(epsg =='3413'):
        zone = epsg 
        hemi = 'arctic'
    elif(epsg =='3976'):
        zone = epsg 
        hemi = 'arctic'
    else:
        print('Could not read EPSG for hemi/zone')
        zone = epsg 
        hemi = ''
    return hemi, zone

# endDef
    

# Calls functions to find EPSG code and perform GCS transform automatically.
def wgs84_to_utm_find_and_transform(lon, lat):
    epsg_out = find_utm_zone_arr(lon, lat)
    xx, yy = wgs84_to_epsg_transform(epsg_out, lon, lat)
    return xx,yy,epsg_out

# endDef
    

# Inputs: lon, lat, (UTM zone), (UTM hemisphere)
def getLatLon2UTM(*args):

    # Set EPSG code for lat/lon coords
    epsg_in = 'epsg:4326'
    
    # Get lats/lons
    lon = args[0]
    lat = args[1]
    
    # Call function based on number of input args
    if(len(args) > 2):
        
        # Get zone/hemi
        zone = args[2]
        hemi = args[3]
        
        # Get EPSG out code for UTM coords
        if(hemi=='N'):
            if len(zone) == 1:
                zone = '0' + zone
            epsg_out = 'epsg:326' + zone
        else:
            if len(zone) == 1:
                zone = '0' + zone 
            epsg_out = 'epsg:327' + zone
        # endif
        
        # Call transform function
        if len(lon) > 0 and len(lat) > 0:
            # print(epsg_in)
            # print(epsg_out)
            xx, yy = transform(epsg_in, epsg_out, lon, lat)
        else:
            xx, yy = np.array([]), np.array([])

    else:
        
        # Get UTM coords
        xx, yy, epsg_out = wgs84_to_utm_find_and_transform(lon, lat)
        
        if(epsg_out=='3413'):
            
            zone = epsg_out 
            hemi = 'arctic'
            
        elif(epsg_out=='3976'):
            
            zone = epsg_out 
            hemi = 'antarctic'
        
        else:
        
            # Store zone
            zone = epsg_out[3:]
            
            # Store hemisphere
            if(epsg_out[0:3]=='326'):   
                hemi = 'N'  
            else:
                hemi = 'S'
            # endIf
            
        # endIf
        
    # endIf
    
    # Return output
    return xx, yy, zone, hemi

# endDef
    

# Function to convert UTM to lat/lon
def getUTM2LatLon(x,y,zone,hemi):
    
    # Set EPSG code for lat/lon coords
    epsg_out = 'epsg:4326'
    
    # Get EPSG code for UTM coords
    if(hemi=='N'):
        if len(zone) == 1:
            zone = "0" + zone
        epsg_in = 'epsg:326' + zone
    elif(hemi=='S'):
        if len(zone) == 1:
            zone = "0" + zone
        epsg_in = 'epsg:327' + zone
    elif(hemi=='arctic'):
        epsg_in = 'epsg:3413'
    elif(hemi=='antarctic'):
        epsg_in = 'epsg:3976'
    # endif
    
    # Call transform function
    lon, lat = transform(epsg_in, epsg_out, x, y)
        
    return lat, lon

# endDef
    

# Function to process mega tiles
def process_mega_tiles(inpFile):
    
    # Initialize inputs
    df_all = []
    verbose = True
    
    # inpFile = 'test_1.inp'
    with open(inpFile) as inputFile:
        
        try:        
            exec(inputFile.read(), globals())
        except:
            print('ERROR: Could not read .inp file, check contents.\n')
        # endTry
    # endWith
        
    # Print user input 
    if(verbose):
        print()
        print('USER INPUT:')
        print('------------')
        print('Min Lat: %0.0f' %latMin)
        print('Max Lat: %0.0f' %latMax)
        print('Min Lon: %0.0f' %lonMin)
        print('Max Lon: %0.0f' %lonMax)
        print()
        print('Cycle(s): %s' %cycle)
        print('-----------------')
        print()
    # endIf
    
    # Arrays with lat/lon bounds
    lats = np.arange(latMin, latMax+1, 1)
    lons = np.arange(lonMin, lonMax+1, 1)
    
    # Filter file paths
    inpPaths = np.array([])
    inpPaths_curr = []
    for i in range(0,len(cycle)):
        
        # Get current cycle number
        cycle_curr = cycle[i]
        
        # Get LAZ files that meet lat/lon bounds from user
        inpFilesLatLon = getFileNamesInLatLon(inpDir, latMin, latMax, lonMin, lonMax)
        
        # Get filtered file paths with matching cycle number
        if(cycle_curr.lower()=='all'):
            inpPaths_curr = inpFilesLatLon
        else:
            cycleSearch = 'cycle' + cycle_curr
            for match in inpFilesLatLon:
                if(cycleSearch in match):
                    inpPaths_curr.append(match)
                # endIf
            # endFor
        # endIf
        
        # Concatenate filtered file paths
        inpPaths = np.concatenate([inpPaths, inpPaths_curr])
        
    # endFor
        
    # Sort filtered file paths
    inpPaths = np.sort(inpPaths)
        
    # Initialize parameters
    numFiles = len(inpPaths)
    df_all = pd.DataFrame()
    primary_zone = []
    zoneCount = 2
        
    if(numFiles>0):
        
        # Load each mega tile (.laz) file and store as a Pandas dataframe
        for i in range(0,numFiles):
            
            # Get path to current file
            inpPath = inpPaths[i]
        
            # Read .laz file
            with pylas.open(inpPath) as lasFile:
                
                # Get file info
                fileType = (inpPath.split('.')[-1]).upper()
                inpFile = os.path.basename(inpPath)
                
                # Read las/laz file
                if(verbose):
                    print('Reading %s file #%d of %d...' %(fileType, i+1, numFiles))
                    time.sleep(0.1)
                # endIf
                  
                # Read las/laz file
                las = lasFile.read()
                
                # Get meta-data
                numPoints = las.header.point_count
                verMajor = las.header.version_major
                verMinor = las.header.version_minor
                
                # Print to screen
                if(verbose):
                    print('File Name: %s' %inpFile)
                    print('%s Version: %s.%s' %(fileType, verMajor, verMinor))
                    print('Number of Points: %s' %numPoints)
                # endIf
                
            # endWith
                
            # Create data frame
            df = pd.DataFrame()
            
            # Convert lat/lon to UTM
            if(i==0):
                x, y, zone, hemi = getLatLon2UTM(las.x, las.y)
                primary_zone = zone
                primary_hemi = hemi
            else:
                x, y, zone, hemi = getLatLon2UTM(las.x, las.y)
                if(zone!=primary_zone or hemi!=primary_hemi):
                    print('WARNING: Data spans %d UTM zones. Forcing to use UTM Zone = %s %s.' %(zoneCount,primary_zone,primary_hemi))
                    x, y, zone, hemi = getLatLon2UTM(las.x, las.y, primary_zone, primary_hemi)
                    zoneCount += 1
                # endIf
            # endIf
            
            if(verbose):
                if(primary_hemi=='N' or primary_hemi=='S'):
                    print('UTM Zone/Hemi: %s %s' %(primary_zone, primary_hemi))
                else:
                    print('Polar Stereo EPSG/Region: %s (%s)' %(primary_zone, primary_hemi))
                # endIf
                print()
            # endIf
            
            # Store fields into dataframe 
            df['lon'] = las.x
            df['lat'] = las.y
            df['x'] = x
            df['y'] = y
            df['zone'] = zone
            df['hemi'] = hemi
            df['z_hae'] = las.z
            df['z_msl'] = las.z - las.geoid
            df['intensity'] = las.intensity
            df['return_number'] = las.return_number
            df['number_of_returns'] = las.number_of_returns
            df['synthetic'] = las.synthetic
            df['key_point'] = las.key_point
            df['withheld'] = las.withheld
            df['overlap'] = las.overlap
            df['scanner_channel'] = las.scanner_channel
            df['scan_direction_flag'] = las.scan_direction_flag
            df['edge_of_flight_line'] = las.edge_of_flight_line
            df['classification'] = las.classification
            df['user_data'] = las.user_data
            df['scan_angle_rank'] = las.scan_angle_rank
            df['point_source_id'] = las.point_source_id
            df['gps_time'] = las.gps_time
            df['geoid'] = las.geoid
            df['beam_number'] = las.beam_number
            df['snow_flag'] = las.snow_flag
            df['landcover'] = las.landcover
            # df['landcover'] = las.land_flag
            df['yrdoy'] = las.yrdoy            
            df['yymmdd'] = las.YYMMDD
            df['track_num'] = las.track_num 
            df['rel_canopy_height'] = las.ph_h
            df['solar_elev'] = las.solar_elevation
            df['signal_conf'] = las.signal_conf
            
            # Concatenate data frames together
            df_all = pd.concat([df_all, df])
            
        # endFor
        
        # Reset indices
        df_all.reset_index(drop=True, inplace=True)
        
        # Print default results
        print('Preparing default stats...')
        print_standard_results(df_all)
            
        print()
        print('Processing complete.\n')
    
    else:
        
        print('ERROR: No LAZ files meet your lat/lon/cycle search criteria.')
        print('       Please adjust your input file.\n')
        
    # endIf
    
    return df_all

    #if(plotData):
    #        
    #    # Plot Lat/Lon bounds on map
    #    plt.figure()
    #    roundNearest = 10
    #    offset = 10
    #    llcrnrlon = roundNearest*(np.round((lonMin-offset)/roundNearest))
    #    urcrnrlon = roundNearest*(np.round((lonMax+offset)/roundNearest))
    #    llcrnrlat = roundNearest*(np.round((latMin-offset)/roundNearest))
    #    urcrnrlat = roundNearest*(np.round((latMax+offset)/roundNearest))
    #    mymap = Basemap(projection='merc',
    #                     resolution = 'l', area_thresh = 1000.0,
    #                     llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat, 
    #                     urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat)
    #    mymap.drawcoastlines()
    #    mymap.drawcountries()
    #    mymap.fillcontinents(color = 'white', alpha = 0.3)
    #    mymap.shadedrelief()
    #    parallels = np.arange(0., 91, 5.)
    #    mymap.drawparallels(parallels,labels=[False,True,True,False])
    #    meridians = np.arange(0., 361., 5.)
    #    mymap.drawmeridians(meridians,labels=[True,False,False,True])
    #    lonBounds = np.array([lonMin, lonMin, lonMax, lonMax, lonMin])
    #    latBounds = np.array([latMin, latMax, latMax, latMin, latMin])
    #    map_xBounds, map_yBounds = mymap(lonBounds, latBounds)
    #    map_xPts, map_yPts = mymap(lonPts, latPts)
    #    mymap.pcolormesh(map_xPts, map_yPts, rasterData.grid, cmap=cmap, zorder=10)
    #    mymap.plot(map_xBounds, map_yBounds, color='r', linewidth='3', zorder=11)
    #    cbar = mymap.colorbar(pad=0.5)
    #    cbar.set_label(rasterOp + ' ' + rasterParam + ' (m)') 

# endDef
    
   
# Function to show main menu
def showHeader():
    
    # Clear terminal
    os.system('cls' if os.name == 'nt' else 'clear')
    
    print()
    print('===========================================================================')
    print('|                                                                         |')
    print('|   ICESAT-2 MEGA TILER (%s)                                            |' %verNum)
    print('|                                                                         |')
    print('|   Copyright 2021                                                        |')
    print('|   The Applied Research Laboratories, The University of Texas at Austin  |')
    print('|                                                                         |')
    print('|   This package is free software; the copyright holder gives unlimited   |')
    print('|   permission to copy and/or distribute, with or without modification,   |')
    print('|   as long as this notice is preserved.                                  |')
    print('|                                                                         |')
    print('===========================================================================')
    
# endDef
    
 
# Function to show Help help
def showHelpHelp():
    
    print()
    print('help --> ' + globalVars.helpDef)
    print('         ' + globalVars.helpUse)
    print('         ' + globalVars.helpExp)
    print()
    
# endDef   
 
    
# Function to show List help
def showListHelp():
    
    print()
    print('list --> ' + globalVars.listDef)
    print('         ' + globalVars.listUse)
    print('         ' + globalVars.listExp)
    print()
    
# endDef
    
    
# Function to show Load help
def showLoadHelp():
    
    print()
    print('load --> ' + globalVars.loadDef)
    print('         ' + globalVars.loadUse)
    print('         ' + globalVars.loadExp)
    print()
    
# endDef
    
    
# Function to show Process help
def showProcessHelp():
    
    print()
    print('process --> ' + globalVars.processDef)
    print('            ' + globalVars.processUse)
    print('            ' + globalVars.processExp)
    print()
    
# endDef
    

# Function to show Filter help
def showFilterHelp():
    
    print()
    print('filter --> ' + globalVars.filterDef)
    print('           ' + globalVars.filterUse)
    print('           ' + globalVars.filterExp)
    print('           Options:')
    print('             colNameToFilter -> column name to filter, ex: classification (for column names use: filter list)')
    print('             filterOperation--> filter operation to use, ex: ==, >, <, >=, <=')
    print('             filterValue -----> value to filter on')
    print('           Methods:')
    print('             list --> list dataframe columns after processing, ex: filter list')
    print()
    
# endDef

    
# Function to show Compute help
def showComputeHelp():
    
    print()
    print('compute --> ' + globalVars.computeDef)
    print('            ' + globalVars.computeUse)
    print('            ' + globalVars.computeExp)
    print('            Options:')
    print('              -g --> column to group output by (optional), ex: -g classification or -g beam_number')
    print("              -c --> column(s) to operate on, ex: -c ['z_hae'] or -c ['lat','lon','z_msl']")
    print('              -o --> operation to perform on data, ex: -o min, -o max, or -o mean')
    print('            Methods:')
    print('              list --> list dataframe columns after processing, ex: compute list')
    print()
    
# endDef
  
    
# Function to show Plot help
def showPlotHelp():
    
    print()
    print('plot --> ' + globalVars.plotDef)
    print('         ' + globalVars.plotUse)
    print('         ' + globalVars.plotExp)
    print('         Methods:')
    print('           xy -------> plot xy (UTM or Polar Stereo) data, ex: plot xy OR plot xy class')
    print('           Options:')
    print('             alt ----> color xy data by altitude (default)')
    print('             class --> color xy data by classification')
    print('')
    print('           latlon ---> plot lat/lon data, ex: plot latlon OR plot latlon class')
    print('           Options:')
    print('             alt ----> color lat/lon data by altitude (default)')
    print('             class --> color lat/lon data by classification')
    print('')
    print('           3d -------> plot 3d point cloud, ex: plot 3d OR plot 3d class')
    print('           Options:')
    print('             alt ----> color 3d data by altitude (default)')
    print('             class --> color 3d data by classification')
    print('')
    print('           hist -----> plot histogram of classifications, ex: plot hist')
    print('           Options:')
    print('             all ----> plot histogram of all fields')
    print()
    
# endDef
    
    
# Function to show Export help
def showExportHelp():
    
    print()
    print('export --> ' + globalVars.exportDef)
    print('           ' + globalVars.exportUse)
    print('           ' + globalVars.exportExp)
    print('           Methods:')
    print('             gzip --> export to .gzip file')
    print('             csv ---> export to .csv file')
    print('             laz ---> export to .laz 1.4 file')
    print('             Options:')
    print('               xy -----> use xy (UTM or Polar Stereo) coordinates in .laz file (default)')
    print('               latlon -> use Lat/Lon coordinates in .laz file, ex: export fileName.laz latlon')
    print('             kml ---> export to .kml file')
    print('             tif ---> export to .tif file')
    print('             Options:')
    print('               xy -----> use xy (UTM or Polar Stereo) coordinates in .tif file (default)')
    print('               latlon -> use Lat/Lon coordinates in .tif file, ex: export fileName.tif latlon')
    print()
    
# endDef
    
    
# Function to show main menu
def showMainMenu():
    
    # Print method definitions and usages
    print()
    print('COMMANDS:')
    print('---------')
    print('home -----> ' + globalVars.homeDef)
    print('exit -----> ' + globalVars.exitDef)
    print('help -----> ' + globalVars.helpDef)
    print('            ' + globalVars.helpUse)
    print('list -----> ' + globalVars.listDef)
    print('            ' + globalVars.listUse)
    print('load -----> ' + globalVars.loadDef)
    print('            ' + globalVars.loadUse)
    print('process --> ' + globalVars.processDef)
    print('            ' + globalVars.processUse)
    print('filter ---> ' + globalVars.filterDef)
    print('            ' + globalVars.filterUse)
    print('compute --> ' + globalVars.computeDef)
    print('            ' + globalVars.computeUse)
    print('plot -----> ' + globalVars.plotDef)
    print('            ' + globalVars.plotUse)
    print('export ---> ' + globalVars.exportDef)
    print('            ' + globalVars.exportUse)
    print()

# endDef
        
        
# Function to filter processed data
def computeData(df_all=[], groupBy='None', fieldsToOperateOn=['z_hae'], operation='mean'):
        
    fieldsToOperateOn = eval(fieldsToOperateOn)
    
    if(groupBy!='None'):
        evalStr = 'df_all.groupby(groupBy)[fieldsToOperateOn].' + operation + '()'
    elif(groupBy=='None'): 
        evalStr = 'df_all[fieldsToOperateOn].' + operation + '()'
    else:
        print('ERROR: Incorrect inputs.')
        return -1
    # endIf
    
    # Aggregate results for Column To Operate On by Column Name to Group By
    df_filt = eval(evalStr)
    
    # Print results
    print()
    print(df_filt)
    print()    
    
    return df_filt

# endDef
    
    
# Function to plot utm data
def plot_utm(df_all, rasterData, colorByAlt=True):
    
    # Ignore warnings
    warnings.filterwarnings('ignore')
    
    print()
    print('Loading figure...\n')
    
    # Raster UTM data
    rasterRes = 1000
    # if(len(rasterData.x)==0):
    if(colorByAlt):
        rasterParam = 'z_hae' 
        rasterOp = 'mean'
        rasterData = getRaster(df_all['x'], df_all['y'], df_all[rasterParam], resolution=rasterRes, method=rasterOp, fillValue = np.nan, time = [], xAllArray = [], yAllArray = [], origin=None)
    else:
        rasterParam = 'classification' 
        rasterOp = 'median'
        rasterData = getRaster(df_all['x'], df_all['y'], df_all[rasterParam], resolution=rasterRes, method=rasterOp, fillValue = np.nan, time = [], xAllArray = [], yAllArray = [], origin=None)
        rasterData.grid = np.ceil(rasterData.grid)
    # endIf
    # endIf
    
    # Get zone/hemi
    zone = df_all['zone'].iloc[0]
    hemi = df_all['hemi'].iloc[0]
    
    # Get x,y labels
    if(hemi=='arctic' or hemi=='antarctic'):
        xlabel = 'Polar Stereo X (m)'
        ylabel = 'Polar Stereo Y (m)'
    else:
        xlabel = 'UTME (m)'
        ylabel = 'UTMN (m)'
    # endIf
    
    # Plot UTME, UTMN, Z
    plt.figure()
    plt.ion()
    if(colorByAlt):
        cmap = plt.get_cmap('jet')
    else:
        cmap_np = np.concatenate([globalColors.brown, globalColors.darkGreen, globalColors.liteGreen])
        cmap_list = cmap_np.tolist()
        cmap = LinearSegmentedColormap.from_list('BrDgLg', cmap_list)
    # endIf
    fig = plt.pcolormesh(rasterData.x, rasterData.y, rasterData.grid, cmap=cmap)
    plt.grid(b = True, which = 'major', axis = 'both')
    plt.ticklabel_format(style='plain')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if(hemi=='arctic' or hemi=='antarctic'):
        titleEnd = 'Polar Stereo: ' + hemi
    else:
        titleEnd = 'UTM Zone: ' + zone + hemi
    # endIf
    plt.title('Raster ' + rasterOp + ' ' + rasterParam + ' (Resolution: ' + str(rasterRes) + ' m, ' + titleEnd + ')')
    if(colorByAlt):
        cbar = plt.colorbar()
    else:
        ticks = np.array([1,2,3])
        bounds = np.array([0.5,1.5,2.5,3.5])
        norm = mcol.BoundaryNorm(bounds,cmap.N)
        cbar = plt.colorbar(fig, cmap=cmap, ticks=ticks, boundaries=bounds, norm=norm)
        cbar.set_ticklabels(['Ground','Canopy','Canopy Top'])
    # endIf
    cbar.set_label(rasterOp + ' ' + rasterParam + ' (m)')
    plt.axis('equal')
    plt.tight_layout()
    plt.show()
        
    return rasterData
    
# endDef
    
    
# Function to plot lat/lon data
def plot_ll(df_all, rasterData, colorByAlt=True):
    
    # Ignore warnings
    warnings.filterwarnings('ignore')
    
    print()
    print('Loading figure...\n')
    
    # Raster UTM data
    rasterRes = 1000
    # if(len(rasterData.x)==0):
    if(colorByAlt):
        rasterParam = 'z_hae' 
        rasterOp = 'mean'
        rasterData = getRaster(df_all['x'], df_all['y'], df_all[rasterParam], resolution=rasterRes, method=rasterOp, fillValue = np.nan, time = [], xAllArray = [], yAllArray = [], origin=None)

    else:
        rasterParam = 'classification' 
        rasterOp = 'median'
        rasterData = getRaster(df_all['x'], df_all['y'], df_all[rasterParam], resolution=rasterRes, method=rasterOp, fillValue = np.nan, time = [], xAllArray = [], yAllArray = [], origin=None)
        rasterData.grid = np.ceil(rasterData.grid)
    # endIf
    # endIf
    
    # Get zone/hemi
    zone = df_all['zone'].iloc[0]
    hemi = df_all['hemi'].iloc[0]
    
    # Convert to lat/lon
    latPts, lonPts = getUTM2LatLon(rasterData.x, rasterData.y, zone, hemi)
    
    # Plot lat/lon, Z
    plt.figure()
    plt.ion()
    if(colorByAlt):
        cmap = plt.get_cmap('jet')
    else:
        cmap_np = np.concatenate([globalColors.brown, globalColors.darkGreen, globalColors.liteGreen])
        cmap_list = cmap_np.tolist()
        cmap = LinearSegmentedColormap.from_list('BrDgLg', cmap_list)
    # endIf
    fig = plt.pcolormesh(lonPts, latPts, rasterData.grid, cmap=cmap)
    plt.grid(b = True, which = 'major', axis = 'both')
    plt.ticklabel_format(style='plain')
    plt.xlabel('Lon (deg)')
    plt.ylabel('Lat (deg)')
    plt.title('Raster ' + rasterOp + ' ' + rasterParam + ' (Resolution: ' + str(rasterRes) + ' m)')
    if(colorByAlt):
        cbar = plt.colorbar()
    else:
        ticks = np.array([1,2,3])
        bounds = np.array([0.5,1.5,2.5,3.5])
        norm = mcol.BoundaryNorm(bounds,cmap.N)
        cbar = plt.colorbar(fig, cmap=cmap, ticks=ticks, boundaries=bounds, norm=norm)
        cbar.set_ticklabels(['Ground','Low Veg','High Veg'])
    # endIf
    cbar.set_label(rasterOp + ' ' + rasterParam + ' (m)')
    plt.axis('equal')
    plt.tight_layout()
    plt.show()
        
    return rasterData
    
# endDef
    
    
def plot_hist(df_all, plotAll=False):
    
    # Ignore warnings
    warnings.filterwarnings('ignore')
    
    print()
    print('Loading figure...\n')
    
    if(plotAll):
        
        df_all.hist()
        plt.tight_layout()
        plt.show()
        
    else:
        
        f, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)
        plt.ion()
        zVals = df_all['z_hae'][df_all['classification']==1].to_numpy()
        ax1.hist(zVals, density=True, color=globalColors.brown)
        ax1.grid(b = True, which = 'major', axis = 'both')
        ax1.set_ylabel('% Counts')
        ax1.set_title('Class 1 (Ground)')
        zVals = df_all['z_hae'][df_all['classification']==2].to_numpy()
        ax2.hist(zVals, density=True, color=globalColors.darkGreen)
        ax2.grid(b = True, which = 'major', axis = 'both')
        ax2.set_ylabel('% Counts')
        ax2.set_title('Class 2 (Low Veg)')
        zVals = df_all['z_hae'][df_all['classification']==3].to_numpy()
        ax3.hist(zVals, density=True, color=globalColors.liteGreen)
        ax3.grid(b = True, which = 'major', axis = 'both')
        ax3.set_ylabel('% Counts')
        ax3.set_title('Class 3 (High Veg)')
        ax3.set_xlabel('Z (m)')
        plt.tight_layout()
        plt.show()
        
    # endIf
        
# endDef
    

# Function to plot 3d data with pptk
def plot_3d(df_all, rasterData, pptkViewers, colorByAlt=True):
    
    # Import pptk viewer
    import pptk
    
    # Print to screen
    print()
    print('Loading figure...\n')
    
    # Raster UTM data
    rasterParam = 'z_hae' 
    rasterRes = 1000
    rasterOp = 'mean'
    # if(len(rasterData.x)==0):
    rasterData = getRaster(df_all['x'], df_all['y'], df_all[rasterParam], resolution=rasterRes, method=rasterOp, fillValue = np.nan, time = [], xAllArray = [], yAllArray = [], origin=None)
    # endIf
    
    # Get xyz points
    xyzPoints = df_all[['x','y','z_hae']].to_numpy()
    
    # Get ground/veg classifications
    ground_inds = df_all['classification'].to_numpy()==1
    loVeg_inds = df_all['classification'].to_numpy()==2
    hiVeg_inds = df_all['classification'].to_numpy()==3

    # Assign colors to ground/veg classes
    colors = np.empty((np.shape(xyzPoints)), dtype='float')
    colors[ground_inds,:] = globalColors.brown
    colors[loVeg_inds,:] = globalColors.darkGreen
    colors[hiVeg_inds,:] = globalColors.liteGreen
    
    # View point cloud with pptk (color by altitude or classification)
    if(colorByAlt):
        v = pptk.viewer(xyzPoints, xyzPoints[:,2])
        scaleMin = np.nanmin(np.ravel(rasterData.grid))
        scaleMax = np.nanmax(np.ravel(rasterData.grid))
        v.color_map('jet', scale=[scaleMin, scaleMax])
    else:
        v = pptk.viewer(xyzPoints, colors)
        v.set(point_size=0.1)
    # endIf
    v.set(point_size=0.1, phi=-np.pi/2, theta=np.pi/2)
        
    # Store pptk viewers
    pptkViewer = v
    pptkViewers.append(pptkViewer)
    
    # return rasterData, pptkViewer
    return pptkViewers
        
# endDef
    

# Function to handle User Input
def handle_user_input(userInp):
    
    # Split user input 
    userInpSplitRaw = userInp.split()
    
    # Make each user input lower case
    userInpSplit = [element.lower() for element in userInpSplitRaw]

    # Get user input key (first keyword)
    if(len(userInpSplit)==0):
        userInpKey = 'Empty'
    else:
        userInpKey = userInpSplit[0]
    # endIf
    
    return userInpSplit, userInpKey

# endDef
    

# Function to handle Home option
def handle_home_option():
    
    # Show header message
    showHeader()
    
    # Show main home menu
    showMainMenu()

# endDef
        
# Function to Exit option
def handle_exit_option(pptkViewers):
    
    # Clear terminal
    os.system('cls' if os.name == 'nt' else 'clear')
    
    # Close pptkViewer if open
    if(len(pptkViewers)>0):
        for i in range(0,len(pptkViewers)):
            pptkViewer = pptkViewers[i]
            pptkViewer.close()
        # endFor
    # endIf
    
    # Print exit message
    print()
    print('Exiting ICESat-2 Mega Tiler, bye!\n')
    
# endDef
    
    
# Function to handle Help option
def handle_help_option(userInpSplit):
    
    if(len(userInpSplit)==2):
        
        # Get user option
        userOption = userInpSplit[-1]
        
        # Get help info for user option
        if(userOption=='list'):
            
            showListHelp()
            
        elif(userOption=='load'):
            
            showLoadHelp()
            
        elif(userOption=='process'):
            
            showProcessHelp()
            
        elif(userOption=='filter'):
            
            showFilterHelp()
            
        elif(userOption=='compute'):
            
            showComputeHelp()
            
        elif(userOption=='plot'):
            
            showPlotHelp()
            
        elif(userOption=='export'):
            
            showExportHelp()
            
        else:
            
            print('ERROR: Unknown command.\n')
        
    else:
        
        print('ERROR: Incorrect number of inputs.\n')
        
    # endIf
    
# endDef


# Function to handle List option
def handle_list_option(userInpSplit):
        
    if(len(userInpSplit)==2):
        
        # Get user extension
        ext = userInpSplit[-1]
        
        # Get files in current directory with extension
        extFiles = glob.glob(os.path.normpath(os.curdir + '/*.' + ext))
            
        # Show files with extension
        if(len(extFiles)>0):
            print()
            print('%s Files:' %ext)
            print('------------')
            for i in range(0,len(extFiles)):
                print('%d) %s' %(i+1,extFiles[i]))
            # endFor
            print()
        else:
            print('WARNING: No files in current directory with .%s extension.\n' %ext)
        # endIf
        
    else: 
        print('ERROR: Incorrect number of inputs.\n')
    # endIf

# endDef
    
    
# Function to handle Load option
def handle_load_option(userInpSplit, df_all):
        
    if(len(userInpSplit)==2):
    
        # Get input file to read            
        outName = userInpSplit[-1]
        
        # Get extension
        ext = outName.split('.')[-1].lower()
        
        # Read input file
        if(ext=='gzip'):
            
            # Read gzip file
            print()
            print('Loading %s...' %outName)
            df_all = pd.read_parquet(outName)
            print('Complete.\n')
            
            # Print default results
            print('Preparing default stats...')
            print_standard_results(df_all)
            
        else:
            
            print('ERROR: Incorrect file extension.')
            print('File must be a .gzip file created by this tool.\n')
        
        #endIf
    
    else:
        
        print('ERROR: Incorrect number of inputs.\n')
        
    # endIf
    
    return df_all
        
# endDef
    

# Function to handle Process option
def handle_process_option(userInpSplit, df_all):
        
    if(len(userInpSplit)==2):
    
        # Get input file to process
        inpFile = userInpSplit[1]
        
        # Process input file
        try:
            df_all = process_mega_tiles(inpFile)
        except:
            print('ERROR: Could not process file, please check inputs.\n')
        #endTry
    
    else:
        
        print('ERROR: Incorrect number of inputs.\n')
        
    # endIf
    
    return df_all
        
# endDef
    

# Function to handle Filter option
def handle_filter_option(userInpSplit, df_all):
        
    if('list' in userInpSplit):
            
        # Show columns in dataframe if dataframe exists
        if(len(df_all)>0):
            
            colNames = df_all.columns.to_list()
            
            print()
            print('Column names:')
            print('-------------')
            for i in range(0,len(colNames)):
                print('%d) %s' %(i+1,colNames[i]))
            # endFor
            print()
        
        else:
                
            print('ERROR: Must process data first to view columns.\n')
            
        # endIf
        
    else:
        
        if(len(userInpSplit)==4):
                        
            columnName = userInpSplit[1]
            filterOperation = userInpSplit[2]
            filterValue = userInpSplit[3]
            
            # Filter processed data
            print()
            print('Filtering data...')
            try:
                df_all = filterData(df_all, columnName, filterOperation, filterValue)
            except:
                print('ERROR: Could not filter data, please check inputs.')
            # endTry
    
        else:
            
            print('ERROR: Incorrect number of inputs.\n')
            
        # endIf
    # endIf
    
    return df_all
            
# endDef
        

# Function to handle Compute option
def handle_compute_option(userInpSplit, df_all):
        
    if(any('list' in s for s in userInpSplit)):
            
        # Show columns in dataframe if dataframe exists
        if(len(df_all)>0):
            
            colNames = df_all.columns.to_list()
            
            print()
            print('Column names:')
            print('-------------')
            for i in range(0,len(colNames)):
                print('%d) %s' %(i+1,colNames[i]))
            # endFor
            print()
        
        else:
                
            print('ERROR: Must process data first to view columns.\n')
            
        # endIf
        
    else:
        
        if(len(userInpSplit)>1):
        
            # Initialize inputs
            groupBy = 'None'
            fieldsToOperateOn = 'None'
            operation = 'None'
            
            print()
            
            # Input: groupBy
            if(any('-g' in s for s in userInpSplit)):
                
                index = userInpSplit.index('-g')
                
                groupBy = userInpSplit[index+1]
                print('Group By: %s ' %groupBy)
                
            # endIf
            
            # Input: fieldsToOperateOn
            if(any('-c' in s for s in userInpSplit)):
                
                index = userInpSplit.index('-c')
                
                fieldsToOperateOn = userInpSplit[index+1]
                if(len(eval(fieldsToOperateOn))==1):
                    colText = 'Column'
                else:
                    colText = 'Columns'
                # endIf
                
                print('%s to Operate on: %s' %(colText, fieldsToOperateOn))
                
            # endIf
            
            # Input: operation
            if(any('-o' in s for s in userInpSplit)):
                
                index = userInpSplit.index('-o')
                
                operation = userInpSplit[index+1]
                print('Operation: %s' %operation)
                
            # endIf
        
            # Filter processed data
            if(fieldsToOperateOn!='None' and operation!='None'):
                print()
                print('Computed Output:')
                print('----------------')
                try:
                    _ = computeData(df_all=df_all, 
                                    groupBy=groupBy, 
                                    fieldsToOperateOn=fieldsToOperateOn, 
                                    operation=operation)
                except:
                    print('ERROR: Could not compute data, please check inputs.')
                # endTry
            # endIf
    
        else:
            
            print('ERROR: Incorrect number of inputs.\n')
            
        # endIf
    # endIf
            
# endDef
    

# Function to handle Plot option
def handle_plot_option(userInpSplit, df_all, rasterData, pptkViewers):
                
    if(len(df_all)>0):
        
        if('xy' in userInpSplit): # Input: utm
            
            if(any('class' in s for s in userInpSplit)):
                
                # Plot lat/lon data by classification                
                rasterData = plot_utm(df_all, rasterData, colorByAlt=False)
                
            else:
                
                # Plot lat/lon data by classification                
                rasterData = plot_utm(df_all, rasterData, colorByAlt=True)
                
            # endIf
            
        elif('latlon' in userInpSplit): # Input: lat/lon
            
            if(any('class' in s for s in userInpSplit)):
                
                # Plot lat/lon data by classification                
                rasterData = plot_ll(df_all, rasterData, colorByAlt=False)
                
            else:
                
                # Plot lat/lon data by classification                
                rasterData = plot_ll(df_all, rasterData, colorByAlt=True)
                
            # endIf
            
        # endIf
        
        elif('hist' in userInpSplit): # Input: histogram
            
            if(any('all' in s for s in userInpSplit)):
                
                # Plot histogram for all fields                
                plot_hist(df_all, plotAll=True)
                
            else:
                
                # Plot histogram                
                plot_hist(df_all)
            
            # endIf
            
        # endIf
        
        elif('3d' in userInpSplit): # Input: 3d point cloud
            
            if(any('class' in s for s in userInpSplit)):
                
                # Plot 3d point cloud by classification                
                pptkViewers = plot_3d(df_all, rasterData, pptkViewers, colorByAlt=False)
                
            else:
                
                # Plot 3d point cloud by altitude              
                pptkViewers = plot_3d(df_all, rasterData, pptkViewers, colorByAlt=True)
                
            # endIf
        
        else:
            
            print('ERROR: Unknown command (only utm, latlon, hist, or 3d). Use help plot for more info.\n')
        
        # endIf

    else:
        
        print('ERROR: No processed data to plot.\n')
        
    # endIf
    
    return rasterData, pptkViewers
        
# endDef
    

# Function to handle Export option
def handle_export_option(userInpSplit, df_all, rasterData):
        
    if(len(userInpSplit)>=2):
        
        if(len(df_all)>0):
        
            outName = userInpSplit[1]
            ext = outName.split('.')[-1].lower()
            
            # Determine export file type
            if(ext=='gzip'):
                
                # Write .gzip file
                print()
                print('Exporting data to %s...' %outName)
                df_all.to_parquet(outName, compression='gzip')
                print('Complete.\n')
                
            elif(ext=='csv'):
                
                # Write to .csv file
                print()
                print('Exporting data to %s...' %outName)
                df_all.to_csv(outName)
                print('Complete.\n')
    
            elif(ext=='laz'):
                
                try:
                    # Write to .laz file (v 1.4)
                    print()
                    if('latlon' in userInpSplit):
                        print('Exporting data to %s 1.4 file format (Lat/Lon coords)...' %outName)
                        writeLAZ(df_all, outName, utm=False)
                    else:
                        if(df_all.hemi[0]=='arctic' or df_all.hemi[0]=='antarctic'):
                            print('Exporting data to %s 1.4 file format (Polar Stereo coords)...' %outName)
                        else:
                            print('Exporting data to %s 1.4 file format (UTM coords)...' %outName)
                        # endIf
                        writeLAZ(df_all, outName, utm=True)
                    # endIf
                    print('Complete.\n')
                    
                except:
                    print('ERROR: Could not export to .laz file, please check inputs.')
                # endTry
                
            elif(ext=='kml'):
                
                try:
                    # Write to .kml file
                    print()
                    print('Exporting data to %s...' %outName)
                    writeKML(df_all, outName)
                    print('Complete.\n')
                except:
                    print('ERROR: Could not export to .kml file, please check inputs.')
                # endTry
                
            elif(ext=='tif'):
                try:
                    # Write to .tif file
                    print()
                    if('latlon' in userInpSplit):
                        print('Exporting data to %s (Lat/Lon coords)...' %outName)
                        writeTIF(df_all, outName, utm=False, colorByAlt=True)
                    else:
                        if(df_all.hemi[0]=='arctic' or df_all.hemi[0]=='antarctic'):
                            print('Exporting data to %s (Polar Stereo coords)...' %outName)
                        else:
                            print('Exporting data to %s (UTM coords)...' %outName)
                        # endIf
                        writeTIF(df_all, outName, utm=True, colorByAlt=True)
                    # endIf   
                    print('Complete.\n')
                except:
                    print('ERROR: Could not export to .tif file, please check inputs.')
                # endTry
                
            else:
                
                print('ERROR: Incorrect file extension.')
                print('Must be either: .gzip, .csv, .laz, .kml, or .tif.\n')
                
            # endIf
        
        else:
            
            print('ERROR: No data to export.\n')
            
        # endIf
    
    else:
        
        print('ERROR: Incorrect number of inputs.\n')
        
    # endIf
        
# endDef
        
        
# Function to grid point cloud data
def getRaster(x, y, z, resolution, method, fillValue = -999, time = [], xAllArray = [], yAllArray = [],
                origin=None):
    
    # USER INPUTS
    # ---------------------------
    # x = input x array of values
    # y = input y array of values
    # z = input z array of values
    # resolution = resolution of grid cells (N = N x N, [M,N] = M x N)
    # method = operation to perform in each grid cell
    #   - min
    #   - max
    #   - mean (default)
    #   - median
    #   - range
    #   - std (standard deviation)
    #   - numel (number of elements)
    # fillValue = value to fill in empty grid cells
    # time = secondary array (like z) to perform operation on in each grid cell
    # xAllArray = force output X grid cells (use np.arange(start, stop + 1, step))
    # yAllArray = force output Y grid cells (use np.arange(start, stop + 1, step))
    #
    # OUTPUTS
    # -------------
    # output.x = x (2D raster)
    # output.y = y (2D raster)
    # output.grid = grid (2D raster)
    # output.t = time (2D raster)
    
    
    # Get X,Y resolution
    if isinstance(resolution,np.integer):
        xResolution = float(resolution)
        yResolution = float(resolution)
    elif isinstance(resolution,int):
        xResolution = float(resolution)
        yResolution = float(resolution)
    elif isinstance(resolution,float):
        xResolution = resolution
        yResolution = resolution
    elif isinstance(resolution,np.ndarray):
        if len(resolution) == 1:
            xResolution = float(resolution)
            yResolution = float(resolution)
        elif len(resolution) == 2:
            xResolution = float(resolution[0])
            yResolution = float(resolution[1])
        else:
            print("Incorrect resolution input")
    elif isinstance(resolution,list):
        xResolution = float(resolution[0])
        yResolution = float(resolution[1])
    elif isinstance(resolution,str):
        strList = resolution.split(",")
        if len(strList) == 1:
            xResolution = float(resolution)
            yResolution = float(resolution)
        elif len(strList) == 2:
            xResolution = float(strList[0])
            yResolution = float(strList[1])
        else:
            print("Incorrect resolution input")
    else:
        print("Incorrect resolution input")

    # Get grid method
    if(method.lower() == 'min'):
        npOperation = np.nanmin
    elif(method.lower() == 'max'):
        npOperation = np.nanmax
    elif(method.lower() == 'mean'):
        npOperation = np.nanmean
    elif(method.lower() == 'median'):
        npOperation = np.nanmedian
    elif(method.lower() == 'range'):
        npOperation = np.range
    elif(method.lower() == 'std'):
        npOperation = np.nanstd
    elif(method.lower() == 'numel'):
        npOperation = np.size
    else:
        npOperation = np.mean
    # EndIf
            
    # Round all incoming X,Y data
    # if type(origin) == type(None):
    xRnd = (np.round(x/xResolution)*xResolution).astype(int)
    yRnd = (np.round(y/yResolution)*yResolution).astype(int)
    # else:
    #     x_corner, y_corner = origin[0], origin[1]
    #     xRnd = np.floor((x - x_corner) / xResolution).astype(int)
    #     yRnd = np.floor((y - y_corner) / yResolution).astype(int)
    #     # yRnd = np.floor(-(y - y_corner) / res_y).astype(int)
    
    # Get output X,Y grid cells
    if(any(xAllArray) and any(yAllArray)):
        
        xAll = xAllArray
        yAll = np.flipud(yAllArray)
        
    else:
        # Get min,max of rounded X,Y data
        xRndMin = xRnd.min()
        xRndMax = xRnd.max()
        yRndMin = yRnd.min()
        yRndMax = yRnd.max()
        
        # Get all possible grid combinations
        xAll = np.arange(xRndMin, xRndMax + xResolution, xResolution)
        yAll = np.arange(yRndMax, yRndMin - yResolution, -yResolution)

    # endIf
        
    # Get X,Y array of all pts
    xAllArray, yAllArray = np.meshgrid(xAll,yAll)
    
    # Get raster X, Y, Z space
    rasterDataX = xAllArray;
    rasterDataY = yAllArray;
    rasterDataZ = fillValue*np.ones((np.shape(rasterDataX)))
    
    # Get x-rastered, y-rastered, and z data into array
    data = np.column_stack([xRnd, yRnd, z])
    
    # Put array into Pandas dataframe
    df = pd.DataFrame(data, columns=['xRnd', 'yRnd', 'z_hae'])
    
    # Do groupby.agg to get get rastered operation for group
    groupedData = df.groupby(['xRnd', 'yRnd']).agg({'z_hae': [npOperation]})
    groupedData.columns = ['z_agg']
    groupedData = groupedData.reset_index()
    zValsNew = np.array(groupedData['z_agg'])
    
    # Determine new row, column indices to place rastered Z data into
    # if discrete_res:
    df_xRnd_min = np.min(groupedData['xRnd'])
    df_yRnd_min = np.min(groupedData['yRnd'])
    colIndsNew = ((np.array(groupedData['xRnd']) - df_xRnd_min)/xResolution).astype(int)
    rowIndsNew = ((np.array(groupedData['yRnd']) - df_yRnd_min)/yResolution).astype(int)

    # else:
    #     xRndRange = np.arange(min(xRnd), max(xRnd)+1, 1)
    #     yRndRange = np.arange(min(yRnd), max(yRnd)+1, 1)
    #     xRnd0 = xRndRange[0] # min(xRnd)
    #     yRnd0 = yRndRange[0] # min(yRnd)
    #     colIndsNew = np.array(groupedData['xRnd'] - xRnd0).astype(int)
    #     rowIndsNew = np.array(groupedData['yRnd'] - yRnd0).astype(int)

    # Populate rastered Z data into array
    rasterDataZ[rowIndsNew, colIndsNew] = zValsNew
    rasterDataZ = np.flipud(rasterDataZ)
    
    # Grid 'time' array if necessary
    if(any(time)):
        
        # Get x-rastered, y-rastered, and time data into array
        dataTime = np.column_stack([xRnd, yRnd, time])
        
        # Put array into Pandas dataframe
        dfTime = pd.DataFrame(dataTime, columns=['xRnd', 'yRnd', 'time'])
        
        # Do groupby.agg to get get rastered operation for group
        groupedDataTime = dfTime.groupby(['xRnd', 'yRnd']).agg({'time': [npOperation]})
        groupedDataTime.columns = ['time_agg']
        groupedDataTime = groupedDataTime.reset_index()
        tValsNew = np.array(groupedDataTime['time_agg'])
        
        # Populate rastered Z data into array
        rasterDataT = fillValue*np.ones((np.shape(rasterDataX)))
        rasterDataT[rowIndsNew, colIndsNew] = tValsNew
        rasterDataT = np.flipud(rasterDataT)
        
    else:
        
        rasterDataT = []
        
    # EndIf

    # Return output
    return GridStruct(rasterDataX, rasterDataY, rasterDataZ, rasterDataT)

# endDef
    

# Function to write LAZ 1.4 files using pylas (for use with Mega Tiler)
def writeLAZ(df, outfilepath, utm=True):     

    # Create new VLR
    new_vlr = pylas.vlrs.VLR(user_id = 'LASF_Projection',
                             record_id = 2112,
                             description = 'OGC Coordinate System WKT')
    
    # Determine wkt and xx, yy, and zz data
    if(utm):
        zone = df['zone'].iloc[0]
        hemi = df['hemi'].iloc[0]
        if(hemi=='arctic' or hemi=='antarctic'):
            proj = hemi
        else:
            proj = 'utm'
        # endIf
        wkt = selectwkt(proj, zone=zone, hemi=hemi)
        xx = np.array(df['x'])
        yy = np.array(df['y'])   
    else:
        proj = 'wgs84'  
        wkt = selectwkt(proj)
        xx = np.array(df['lon'])
        yy = np.array(df['lat'])
    # endIf
    zz = np.array(df['z_hae'])
    
    # Assign wkt string
    new_vlr.record_data = wkt
    
    # Create new header
    las = pylas.create(point_format_id=6, file_version='1.4')
    
    # Create VLR headers
    las.vlrs.append(new_vlr)
    if(utm):
        
        xmin = np.min(xx)
        ymin = np.min(yy)
        zmin = np.min(zz)
        
        xmax = np.max(xx)
        ymax = np.max(yy)
        zmax = np.max(zz)
    
        # Calculate x, y, and z scale factor
        if(xmax==xmin):
            las.header.x_scale = 1
        else:
            las.header.x_scale = (xmax - xmin) / 100000000;
        # endIf
        
        if(ymax==ymin):
            las.header.y_scale = 1
        else:
            las.header.y_scale = (ymax - ymin) / 100000000;
        # endIf
        
        if(zmax==zmin):
            las.header.z_scale = 1
        else:
            las.header.z_scale = (zmax - zmin) / 100000000;
        # endIf
        
    else:
        
        las.header.global_encoding.wkt = 1
        las.header.x_scale = 0.000001
        las.header.y_scale = 0.000001
        las.header.z_scale = 0.000001
    
    # endIf
        
    # Create extra fields and data types for them
    names=['geoid','beam_number','snow_flag','landcover','yrdoy','ph_h','solar_elevation','signal_conf']
    data_types= ['uint8','uint8','uint8','uint8','uint32','uint8','uint8','uint8']
    descriptions = ['Geoid (m)','Beam Descriptor','Modis Snowcover Flag','Landcover class','Year and Day of Year','Canopy Height (m)','Solar Elevation','Photon Signal Confidence']
    
    # Create descriptions
    for idx, nm in enumerate(names):
        las.add_extra_dim(name=nm, type=data_types[idx], description=descriptions[idx])
    # endFor
    
    # Write dataframe to LAZ object    
    las.x = xx
    las.y = yy
    las.z = zz
    las.intensity = np.array(df['intensity'])
    las.return_number = np.array(df['return_number'])
    las.number_of_returns = np.array(df['number_of_returns'])
    las.synthetic = np.array(df['synthetic'])
    las.key_point = np.array(df['key_point'])
    las.withheld = np.array(df['withheld'])
    las.overlap = np.array(df['overlap'])
    las.scanner_channel = np.array(df['scanner_channel']) 
    las.scan_direction_flag = np.array(df['scan_direction_flag'])
    las.edge_of_flight_line = np.array(df['edge_of_flight_line'])
    las.classification = np.array(df['classification'])
    las.user_data = np.array(df['user_data'])
    las.scan_angle_rank = np.array(df['scan_angle_rank'])
    las.point_source_id = np.array(df['point_source_id'])
    las.gps_time = np.array(df['gps_time'])
    las.geoid = np.array(df['geoid'])
    las.beam_number = np.array(df['beam_number'])
    las.snow_flag = np.array(df['snow_flag'])
    las.landcover = np.array(df['landcover'])
    las.yrdoy = np.array(df['yrdoy'])
    las.ph_h = np.array(df['rel_canopy_height'])
    las.solar_elevation = np.array(df['solar_elev'])
    las.signal_conf = np.array(df['signal_conf'])
    
    # Write LAZ object to LAZ file  
    las.write(outfilepath, do_compress=True)

# endDef
    

# Function to create WKT
def selectwkt(proj,hemi=None,zone=None):
    
    if proj.lower() == 'utm':
        if zone:
            zone = str(zone)
            if hemi.lower() == 'n':
                if zone == '1' or zone == '01':
                    wkt = b'''PROJCS["WGS 84 / UTM zone 1N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-177],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32601"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "2" or zone == "02":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 2N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-171],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32602"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "3" or zone == "03":
                    wkt = b'''ROJCS["WGS 84 / UTM zone 3N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-165],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32603"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "4" or zone == "04":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 4N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-159],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32604"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "5" or zone == "05":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 5N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-153],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32605"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "6" or zone == "06":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 6N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-147],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32606"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "7" or zone == "07":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 7N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-141],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32607"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "8" or zone == "08":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 8N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-135],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32608"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "9" or zone == "09":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 9N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-129],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32609"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "10": 
                    wkt = b'''PROJCS["WGS 84 / UTM zone 10N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-123],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32610"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''                        
                elif zone == "11":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 11N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-117],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32611"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "12":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 12N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-111],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32612"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "13":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 13N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-105],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32613"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "14":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 14N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-99],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32614"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "15":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 15N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-93],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32615"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "16":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 16N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-87],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32616"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "17":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 17N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-81],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32617"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "18":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 18N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-75],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32618"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "19":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 19N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-69],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32619"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "20":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 20N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-63],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32620"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "21":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 21N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-57],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32621"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "22":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 22N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-51],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32622"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "23":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 23N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-45],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32623"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "24":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 24N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-39],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32624"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "25":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 25N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-33],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32625"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "26":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 26N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-27],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32626"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "27":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 27N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-21],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32627"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "28":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 28N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-15],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32628"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "29":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 29N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-9],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32629"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "30":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 30N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-3],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32630"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "31":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 31N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",3],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32631"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "32":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 32N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",9],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32632"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "33":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 33N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",15],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32633"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "34":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 34N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",21],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32634"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "35":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 35N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",27],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32635"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "36":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 36N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",33],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32636"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "37":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 37N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",39],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32637"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "38":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 38N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",45],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32638"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "39":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 39N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",51],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32639"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "40":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 40N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",57],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32640"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "41":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 41N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",63],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32641"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "42":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 42N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",69],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32642"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "43":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 43N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",75],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32643"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "44":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 44N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",81],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32644"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "45":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 45N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",87],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32645"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "46":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 46N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",93],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32646"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "47":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 46N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",93],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32646"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "48":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 48N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",105],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32648"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "49":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 49N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",111],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32649"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "50":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 50N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",117],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32650"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "51":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 51N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",123],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32651"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "52":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 52N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",129],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32652"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "53":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 53N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",135],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32653"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "54":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 54N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",141],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32654"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "55":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 55N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",147],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32655"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "56":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 56N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",153],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32656"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "57":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 57N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",159],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32657"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "58":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 58N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",165],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32658"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "59":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 59N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",171],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32659"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "60":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 60N",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",177],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],AUTHORITY["EPSG","32660"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
            elif hemi.lower() == "s":
                if zone == "1" or zone == "01":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 1S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-177],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32701"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "2" or zone == "02":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 2S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-171],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32702"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "3" or zone == "03":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 3S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-165],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32703"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "4" or zone == "04":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 4S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-159],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32704"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "5" or zone == "05":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 5S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-153],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32705"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "6" or zone == "06":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 6S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-147],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32706"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "7" or zone == "07":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 7S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-141],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32707"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "8" or zone == "08":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 8S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-135],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32708"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "9" or zone == "09":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 9S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-129],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32709"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "10":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 10S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-123],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32710"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "11":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 11S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-117],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32711"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "12":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 12S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-111],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32712"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "13":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 13S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-105],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32713"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "14":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 14S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-99],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32714"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "15":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 15S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-93],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32715"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "16":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 16S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-87],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32716"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "17":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 17S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-81],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32717"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "18":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 18S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-75],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32718"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "19":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 19S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-69],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32719"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "20":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 20S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-63],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32720"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "21":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 21S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-57],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32721"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "22":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 22S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-51],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32722"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "23":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 23S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-45],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32723"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "24":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 24S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-39],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32724"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "25":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 25S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-33],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32725"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "26":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 26S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-27],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32726"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "27":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 27S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-21],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32727"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "28":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 28S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-15],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32728"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "29":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 29S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-9],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32729"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "30":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 30S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-3],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32730"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "31":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 31S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",3],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32731"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "32":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 32S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",9],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32732"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "33":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 33S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",15],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32733"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "34":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 34S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",21],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32734"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "35":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 35S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",27],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32735"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "36":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 36S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",33],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32736"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "37":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 37S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",39],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32737"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "38":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 38S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",45],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32738"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "39":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 39S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",51],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32739"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "40":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 40S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",57],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32740"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "41":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 41S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",63],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32741"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "42":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 42S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",69],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32742"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "43":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 43S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",75],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32743"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "44":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 44S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",81],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32744"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "45":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 45S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",87],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32745"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "46":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 46S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",93],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32746"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "47":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 47S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",99],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32747"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "48":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 48S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",105],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32748"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "49":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 49S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",111],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32749"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "50":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 50S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",117],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32750"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "51":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 51S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",123],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32751"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "52":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 52S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",129],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32752"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "53":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 53S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",135],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32753"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "54":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 54S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",141],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32754"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "55":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 55S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",147],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32755"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "56":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 56S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",153],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32756"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "57":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 57S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",159],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32757"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "58":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 58S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",165],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32758"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "59":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 59S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",171],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32759"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
                elif zone == "60":
                    wkt = b'''PROJCS["WGS 84 / UTM zone 60S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",177],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],AUTHORITY["EPSG","32760"],AXIS["Easting",EAST],AXIS["Northing",NORTH]]'''
    elif proj.lower() == "arctic" or proj.lower() == "north polar":
        wkt = b'''PROJCS["NSIDC Sea Ice Polar Stereographic North",GEOGCS["Unspecified datum based upon the Hughes 1980 ellipsoid",DATUM["Not_specified_based_on_Hughes_1980_ellipsoid",SPHEROID["Hughes 1980",6378273,298.279411123064,AUTHORITY["EPSG","7058"]],AUTHORITY["EPSG","6054"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4054"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Polar_Stereographic"],PARAMETER["latitude_of_origin",70],PARAMETER["central_meridian",-45],PARAMETER["scale_factor",1],PARAMETER["false_easting",0],PARAMETER["false_northing",0],AUTHORITY["EPSG","3411"],AXIS["X",UNKNOWN],AXIS["Y",UNKNOWN]]'''
    elif proj.lower() == "antarctic" or proj.lower() == "south polar": 
        wkt = b'''PROJCS["NSIDC Sea Ice Polar Stereographic South",GEOGCS["Unspecified datum based upon the Hughes 1980 ellipsoid",DATUM["Not_specified_based_on_Hughes_1980_ellipsoid",SPHEROID["Hughes 1980",6378273,298.279411123064,AUTHORITY["EPSG","7058"]],AUTHORITY["EPSG","6054"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4054"]],UNIT["metre",1,AUTHORITY["EPSG","9001"]],PROJECTION["Polar_Stereographic"],PARAMETER["latitude_of_origin",-70],PARAMETER["central_meridian",0],PARAMETER["scale_factor",1],PARAMETER["false_easting",0],PARAMETER["false_northing",0],AUTHORITY["EPSG","3412"],AXIS["X",UNKNOWN],AXIS["Y",UNKNOWN]]'''
    elif proj.lower() == "wgs84":
        wkt = b'''GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]''' 
    else:
        print("No defined Projected Coordinate System Selected")
    return wkt

# endDef
    

# Function to write to KML file
def writeKML(df, kmlName):

    # Suppress warnings that may come from simple kml
    if not sys.warnoptions:
        warnings.simplefilter('ignore')
    # endif
    
    # Open Simple KML
    kml = simplekml.Kml()
    
#    # Open Simple KML style editor
#    style = simplekml.Style()
#    style.iconstyle.icon.href = 'http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png'
#    style.iconstyle.color = 'ff0000ff'  # Red
    
    # Get unique tracks to make kml lines for each ICESat-2 track
    uniqueTracks = list(df.groupby(['yymmdd','track_num','beam_number']).indices.values())
    
    # Loop through each unique file (track)
    for j in range(0,len(uniqueTracks)):
        
        # Get indices for current track
        currTrack = uniqueTracks[j][0]
        currFile = 'ATL03_' + str(df['yymmdd'][currTrack]) + '_' + f"{df['track_num'][currTrack]:04}" + '_beam_' + str(df['beam_number'][currTrack])
        matchingInds = uniqueTracks[j]
    
        # Get all lat/lon coords
        lons_all = df['lon'][matchingInds].to_numpy()
        lats_all = df['lat'][matchingInds].to_numpy()
        
        # Get step value for discrete markers
        latStep = 0.01 # degrees
        latVals = np.arange(np.min(lats_all), np.max(lats_all) + 1, latStep)
    
        # Get closest time values from IceSat MEASURED data
        _, indsToUse = getClosest(lats_all, latVals)
        
        indsToUse = np.unique(indsToUse)
    
        # Reduce lat/lon values to user-specified scale
        lons = lons_all[indsToUse]
        lats = lats_all[indsToUse]
        
        # Plot line for ground track
        latLon = np.column_stack((lons, lats))
        latLonTuple = tuple(latLon)
        ls = kml.newlinestring(name=currFile, coords=latLonTuple)
        ls.extrude = 1
        ls.altitudemode = simplekml.AltitudeMode.clamptoground
        ls.style.linestyle.width = 5
        ls.style.linestyle.color = simplekml.Color.blue
    
#        # Loop through all lat/lon values and make KML markers
#        for i in range(0,len(lons)):
#            
#            # Plot marker points
#            pnt = kml.newpoint(name='', coords=[(lons[i], lats[i])])
#            pnt.style = style
#            
#        # endFor

    # endFor
    
    # Save output KML file
    kml.save(kmlName)
    
# endDef
    
    
# Function to find closest points in an array
def getClosest(inputArray, closestPts):
    
    # Initialize outputs
    minInd = np.zeros(np.shape(closestPts), dtype = int)
    minVal = np.zeros(np.shape(closestPts))
    
    # Loop through closest points array and find closest point to input array
    for i in range(0,len(closestPts)):
        closestPt = closestPts[i]
        arrayDif = np.abs(inputArray - closestPt)
        minInd[i] = np.argmin(arrayDif) 
        minVal[i] = inputArray[minInd[i]]
    # EndFor
        
    # Return outputs
    return minVal, minInd

# endDef
    

# Function to write raster data to .tif file
def writeTIF(df_all, tifName, utm=True, colorByAlt=True):
    
    # Raster UTM data
    rasterRes = 1000
    if(colorByAlt):
        rasterParam = 'z_hae' 
        rasterOp = 'mean'
        rasterData = getRaster(df_all['x'], df_all['y'], df_all[rasterParam], resolution=rasterRes, method=rasterOp, fillValue = np.nan, time = [], xAllArray = [], yAllArray = [], origin=None)
    else:
        rasterParam = 'classification' 
        rasterOp = 'median'
        rasterData = getRaster(df_all['x'], df_all['y'], df_all[rasterParam], resolution=rasterRes, method=rasterOp, fillValue = np.nan, time = [], xAllArray = [], yAllArray = [], origin=None)
        rasterData.grid = np.ceil(rasterData.grid)
    # endIf
    
    # Get zone, hemi
    hemi = df_all['hemi'].iloc[0]
    zone = df_all['zone'].iloc[0]
        
    # Store x,y data (if lat/lon, then covert utm data)
    if(utm):
        xx = rasterData.x
        yy = rasterData.y
        x_pixel = 1
        y_pixel = 1
    else:
        lats, lons = getUTM2LatLon(rasterData.x, rasterData.y, zone, hemi)
        xx = lons
        yy = lats
        x_pixel = 0.000001
        y_pixel = 0.000001
    # endIf
    zz = rasterData.grid

    # Get EPSG code
    if(utm):
        if(len(zone)==1):
            zone = '0' + zone
        # endIf
        if(hemi=='N'):
            epsg = '326' + zone
        elif(hemi=='S'):
            epsg = '327' + zone
        elif(hemi=='arctic'):
            epsg = '3413'
        elif(hemi=='antarctic'):
            epsg = '3976'
        # endif
    else:
        epsg = '4326'
    # endIf
    epsg = int(epsg)
    
    # Get offsets
    x_offset = 0.0
    y_offset = 0.0

    # Set pixels
    x_pixels = xx.shape[1]
    y_pixels = yy.shape[0]
    
    # Set driver
    driver = gdal.GetDriverByName('GTiff')
    
    # Set spatial reference system
    srs = ogr.osr.SpatialReference()
    srs.ImportFromEPSG(epsg)
    
    # Set dataset
    x_min = np.min(xx)
    y_max = np.max(yy)
    dataset = driver.Create(tifName, x_pixels, y_pixels, 1, gdal.GDT_Float64)
    dataset.SetGeoTransform((x_min + x_offset, x_pixel, 0,                      
                             y_max + y_offset, 0, -y_pixel))  
    dataset.SetProjection(srs.ExportToWkt())
    dataset.GetRasterBand(1).WriteArray(zz)
    dataset.FlushCache()

# endDef
    
    
# Function to print stock output to screen after processing/loading
def print_standard_results(df_all):
    
    # Get default stats on dataframe
    df_stats = df_all.describe().transpose().apply(lambda s: s.apply('{0:.5f}'.format))
    df_stats_short = df_stats[['min','max','mean','std']]
    
    # Print results
    print()
    print('Dataframe Quicklook:')
    print('--------------------')
    print(df_all)
    print()
    print('Dataframe Stats:')
    print('----------------')
    print(df_stats_short)
    print()
            
# endDef
    

# Function to filter data
def filterData(df_all, columnName, filterOperation, filterValue):
    
    # Filter data based on user's input
    eval_str = "df_all[df_all['" + columnName + "']" + filterOperation + filterValue + ']'
    df_all = eval(eval_str)
    
    # Print default results
    print_standard_results(df_all)
    
    return df_all
    
# endDef
    

# Function to get LAZ file names into array
def getFileNamesInLatLon(inpPath, latMin, latMax, lonMin, lonMax):
    
    # Get all LAZ files in input directory
    files = np.array(glob.glob(inpPath + '\\*.laz'))
    
    # Get lat/lon lower/upper bounds from file name
    latUpperBoundInt = np.array([int(os.path.basename(x).split('_')[4].split('N')[0]) for x in files])
    latLowerBoundInt = np.array([int(os.path.basename(x).split('_')[5].split('N')[0]) for x in files])
    lonUpperBoundInt = np.array([int(os.path.basename(x).split('_')[4].split('N')[1][0:-1]) for x in files])
    lonLowerBoundInt = np.array([int(os.path.basename(x).split('_')[5].split('N')[1][0:-1]) for x in files])
    
    # Find Lat/Lon in bounds
    latMinTF = np.logical_and(latMin >= latLowerBoundInt, latMin <= latUpperBoundInt)
    latMaxTF = np.logical_and(latMax >= latLowerBoundInt, latMax <= latUpperBoundInt)
    lonMinTF = np.logical_and(lonMin >= lonLowerBoundInt, lonMin <= lonUpperBoundInt)
    lonMaxTF = np.logical_and(lonMax >= lonLowerBoundInt, lonMax <= lonUpperBoundInt)
    
    # Find files TF that are in bounds
    filesTF = np.logical_and.reduce((latMinTF, latMaxTF, lonMinTF, lonMaxTF))
    
    # Get matching files
    filesMatch = files[filesTF].tolist()
    
    return filesMatch
    
# endDef