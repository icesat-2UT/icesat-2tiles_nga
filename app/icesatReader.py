#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script that contains ATL03 and ATL 08 H5 Reader functions for PhoREAL

Copyright 2019 Applied Research Laboratories, University of Texas at Austin

This package is free software; the copyright holder gives unlimited
permission to copy and/or distribute, with or without modification, as
long as this notice is preserved.

Authors:
    Eric Guenther
    Mike Alonzo
    
Date: February 27, 2019
"""
import pickle
pickle.HIGHEST_PROTOCOL = 4
import pandas as pd
import numpy as np
import os
import h5py

import geopandas as gpd


from icesatUtils import getH5Keys
from icesatIO import readAtl03H5
from icesatIO import readAtlH5
from icesatIO import readAtl03DataMapping
from icesatIO import readAtl08DataMapping
from icesatUtils import getAtl08Mapping
from icesatUtils import wgs84_to_utm_find_and_transform
from icesatUtils import wgs84_to_epsg_transform
from icesatUtils import getCoordRotFwd
from icesatUtils import getNameParts
from icesatUtils import get_h5_meta
from icesatIO import readTruthRegionsTxtFile
from icesatUtils import identify_hemi_zone
from icesatIO import writeLas
from icesatIO import readHeaderMatFile                           
from getAtlMeasuredSwath_auto import atl03Struct as Atl03StructLegacy



class AtlRotationStruct:
    
    # Define class with designated fields
    def __init__(self, R_mat, xRotPt, yRotPt, desiredAngle, phi):
        
        self.R_mat = R_mat
        self.xRotPt = xRotPt
        self.yRotPt = yRotPt
        self.desiredAngle = desiredAngle
        self.phi = phi


class AtlStruct:
        
    # Define class with designated fields
    def __init__(self, df, gtNum, epsg, zone, hemi, kmlRegionName, 
                 headerFilePath, truthFilePath, atlFilePath, atlFileName, 
                 trackDirection, atlProduct, alth5Info, dataIsMapped, 
                 rotation_data, ancillary=None, orbit_info=None):
            
        self.df = df
        self.gtNum = gtNum
        self.epsg = epsg
        self.zone = zone
        self.hemi = hemi
        self.kmlRegionName = kmlRegionName
        self.headerFilePath = headerFilePath
        self.truthFilePath = truthFilePath
        self.atlFilePath = atlFilePath
        self.atlFileName = atlFileName
        self.trackDirection = trackDirection
        self.atlProduct = atlProduct
        self.atlVersion = alth5Info.atlVersion
        self.year = alth5Info.year
        self.month = alth5Info.month
        self.day = alth5Info.day
        self.hour = alth5Info.hour
        self.minute = alth5Info.minute
        self.second = alth5Info.second
        self.trackNum = alth5Info.trackNum
        self.unknown = alth5Info.unknown
        self.releaseNum = alth5Info.releaseNum
        self.incrementNum = alth5Info.incrementNum
        self.dataIsMapped = dataIsMapped
        self.rotationData = rotation_data
        self.ancillary = ancillary
        self.orbit_info = orbit_info
                

# Read ATL03 Heights, put in Pandas DF
def read_atl03_heights_data(atl03filepath, gt):
    # Iterate through keys for "Heights"
    keys = getH5Keys(atl03filepath,gt + '/heights')
    keys.append('geoid') #Include the geoid correction for future post-processing
    keys.append('gps_time_offset')
    # Read each key, put it in pandas df
    for idx, key in enumerate(keys):
        if key =='geoid':
            data = readAtl03H5(atl03filepath, '/geophys_corr/' + key, gt)
        # elif key =='gps_time_offset':
        #     data = readAtl03H5(atl03filepath, '/METADATA/' + key,
        else:
            data = readAtl03H5(atl03filepath, '/heights/' + key, gt)
            
        if idx == 0:
            df = pd.DataFrame(data,columns=[key])
        else:
            df = pd.concat([df,pd.DataFrame(data,columns=[key])],axis=1)
    df['beam'] = gt
    return df

def read_atl03_geolocation(atl03filepath, gt):
    # Iterate through keys for "Heights"
    keys = getH5Keys(atl03filepath,gt + '/geolocation')
    key_info = get_H5_keys_info(atl03filepath, gt + '/geolocation')
    
    # Read each key, put it in pandas df
    for idx, key in enumerate(keys):
        data = readAtl03H5(atl03filepath, '/geolocation/' + key, gt)
        if key_info[idx][1] != 'Group':
            if idx == 0:
                df = pd.DataFrame(data,columns=[key.split('/')[-1]])
            else:
                if len(data.shape) == 2:
                    cols = data.shape[1]
                    for idx2 in range(0,cols):
                        df = pd.concat([df,pd.DataFrame(data[:,idx2],columns=\
                                                        [key.split('/')[-1] +\
                                                         '_' + str(idx2)])],
                                       axis=1)                        
                else:
                    df = pd.concat([df,pd.DataFrame(data,columns=\
                                                    [key.split('/')[-1]])],
                                   axis=1)
    return df

    
# Read ATL03 Heights, put in Pandas DF
def read_atl08_land_segments(atl08filepath, gt):
    # Iterate through keys for "Land Segments"
    keys = getH5Keys(atl08filepath,gt + '/land_segments')
    #/signal_photons/ph_h/
    key_info = get_H5_keys_info(atl08filepath,gt + '/land_segments')
    # Read each key, put it in pandas df
    for idx, key in enumerate(keys):
        data = readAtl03H5(atl08filepath, '/land_segments/' + key, gt)
        if key_info[idx][1] != 'Group':
            if idx == 0:
                df = pd.DataFrame(data,columns=[key.split('/')[-1]])
            else:
                if len(data.shape) == 2:
                    cols = data.shape[1]
                    for idx2 in range(0,cols):
                        df = pd.concat([df,pd.DataFrame(data[:,idx2],columns=\
                                                        [key.split('/')[-1] +\
                                                         '_' + str(idx2)])],
                                       axis=1)                        
                else:
                    df = pd.concat([df,pd.DataFrame(data,columns=\
                                                    [key.split('/')[-1]])],
                                   axis=1)
    return df

def read_atl08_canopy_h(atl08filepath, gt):
    # Iterate through keys for "Land Segments"
    keys = getH5Keys(atl08filepath,gt + '/signal_photons')#signal_photons/ph_h/
    #/signal_photons/ph_h/
    key_info = get_H5_keys_info(atl08filepath,gt + '/signal_photons')
    # Read each key, put it in pandas df
    for idx, key in enumerate(keys):
        data = readAtl03H5(atl08filepath, '/signal_photons/' + key, gt)
        if key_info[idx][1] != 'Group':
            if idx == 0:
                df = pd.DataFrame(data,columns=[key.split('/')[-1]])
            else:
                if len(data.shape) == 2:
                    cols = data.shape[1]
                    for idx2 in range(0,cols):
                        df = pd.concat([df,pd.DataFrame(data[:,idx2],columns=\
                                                        [key.split('/')[-1] +\
                                                         '_' + str(idx2)])],
                                       axis=1)                        
                else:
                    df = pd.concat([df,pd.DataFrame(data,columns=\
                                                    [key.split('/')[-1]])],
                                   axis=1)
    return df

# Read ATL03 Heights, put in Pandas DF
def read_atl09_hr_profile(atl09filepath, gt):
    # Iterate through keys for "Land Segments"
    subgroup = 'profile_' + gt[2] + '/high_rate/'
    keys = getH5Keys(atl09filepath, subgroup)
    key_info = get_H5_keys_info(atl09filepath, subgroup)
    # Read each key, put it in pandas df
    for idx, key in enumerate(keys):
        data = readAtlH5(atl09filepath, subgroup + '/' + key + '/')
        if key == 'ds_layers':
            ds_layers = np.array(data)
        elif key == 'ds_va_bin_h':
            ds_va_bin_h = np.array(data)
        elif key == 'cab_prof':
            cab_prof = np.array(data)
        elif key == 'density_pass1':
            density_pass1 = np.array(data)
        elif key == 'density_pass2':
            density_pass2 = np.array(data)
        else:
            
            if key_info[idx][1] != 'Group':
                if idx == 0:
                    df = pd.DataFrame(data,columns=[key.split('/')[-1]])
                else:
                    if len(data.shape) == 2:
                        cols = data.shape[1]
                        for idx2 in range(0,cols):
                            df = pd.concat(
                                [df,pd.DataFrame(data[:,idx2],columns=\
                                                 [key.split('/')[-1] +'_' + 
                                                  str(idx2)])],axis=1)                        
                    else:
                        df = pd.concat(
                            [df,pd.DataFrame(data,columns=\
                                             [key.split('/')[-1]])],axis=1)
        
    return df, ds_layers, ds_va_bin_h, cab_prof, density_pass1, density_pass2

def read_atl09_ancillary_data(atl09filepath):
    # Iterate through keys for ancillary data
    subgroup = '/ancillary_data/'
    keys = getH5Keys(atl09filepath, subgroup)
    key_info = get_H5_keys_info(atl09filepath, subgroup)
    
    byte_encoded = ['control', 'data_start_utc', 'data_end_utc',
                    'granule_start_utc', 'granule_end_utc']
    
    # Read each key, put it in pandas df
    for idx, key in enumerate(keys):
        data = readAtlH5(atl09filepath, subgroup + '/' + key + '/')
        
        if np.isin(key, ['release', 'version']):
            data = int(data)
            
        if np.isin(key, byte_encoded):
            data = data[0]             
            data = data.decode('utf-8')
        
        if key_info[idx][1] != 'Group':
            if len(key.split('/')) == 1:
                if idx == 0:
                    df = pd.Series(data,
                                   index=[key.split('/')[-1]],
                                   dtype=object)
                else:
                    df = pd.concat(
                        [df, pd.Series(data, 
                                       index=[key.split('/')[-1]],
                                       dtype=object)])
            
    return df

def read_atl09_orbit_info(atl09filepath):

    subgroup = '/orbit_info/'
    keys = getH5Keys(atl09filepath, subgroup)
    key_info = get_H5_keys_info(atl09filepath, subgroup)

    # Read each key, put it in pandas df
    for idx, key in enumerate(keys):
        data = readAtlH5(atl09filepath, subgroup + '/' + key + '/')
        if idx == 0:
            df = pd.Series(data, index=[key], dtype=object)
        else:
            df_key = pd.Series(data, index=[key], dtype=object)
            df = pd.concat([df, df_key])
    
    return df

# Map classifications from ATL08, map back to ATL03 Photons
def get_atl03_classification(atl03filepath, atl08filepath, df, gt):
    # Read ATL03 metrics for class mapping
    atl03_ph_index_beg, atl03_segment_id = \
    readAtl03DataMapping(atl03filepath,gt)
    
    # Read ATL08 for class mapping
    atl08_classed_pc_indx, atl08_classed_pc_flag, atl08_segment_id = \
    readAtl08DataMapping(atl08filepath, gt)
    
    # Map ATL08 classifications to ATL03 Photons
    allph_classed = getAtl08Mapping(atl03_ph_index_beg, atl03_segment_id, 
                                    atl08_classed_pc_indx, 
                                    atl08_classed_pc_flag, 
                                    atl08_segment_id)
    
    # Add classifications to ATL03 DF
    df = pd.concat([df,pd.DataFrame(allph_classed,
                                    columns=['classification'])],axis=1)
    
    # Replace nan with -1 (unclassified)
    df.replace({'classification' : np.nan}, -1)
    return df

# Calculate alongtrack time
def get_atl_time(df):
    delta_time = np.array(df['delta_time'])
    min_detla_time = np.min(delta_time[np.nonzero(delta_time)])
    time = delta_time - min_detla_time
    
    df = pd.concat([df,pd.DataFrame(time,
                                    columns=['time'])],axis=1)
    return df
    
# Calcualte Easting/Northing
def get_atl_coords(df, epsg = None):
    columns = list(df.columns)
    
    if 'lon_ph' in columns:
        lon = np.array(df['lon_ph'])
        lat = np.array(df['lat_ph'])
    elif 'longitude' in columns:
        lon = np.array(df['longitude'])
        lat = np.array(df['latitude'])
    elif 'reference_photon_lon' in columns:
        lon = np.array(df['reference_photon_lon'])
        lat = np.array(df['reference_photon_lat'])        
    
    # Specify EPSG Code or automatically find zone
    if epsg:
        xcoord, ycoord = wgs84_to_epsg_transform(epsg, lon, lat)
    else:
        xcoord, ycoord, epsg = wgs84_to_utm_find_and_transform(lon, lat)
        
    if 'easting' not in columns:
        df = pd.concat([df,pd.DataFrame(xcoord,
                                        columns=['easting'])],axis=1)
        df = pd.concat([df,pd.DataFrame(ycoord,
                                        columns=['northing'])],axis=1)
    else:
        print('Warning: Overwritting Existing Coordinates')
        df = df.drop(columns = ['easting'])
        df = df.drop(columns = ['northing'])

        df = pd.concat([df,pd.DataFrame(xcoord,
                                        columns=['easting'])],axis=1)
        df = pd.concat([df,pd.DataFrame(ycoord,
                                        columns=['northing'])],axis=1)        
    
    return df, epsg

def get_atl_alongtrack(df, atl03struct = None):
    easting = np.array(df['easting'])
    northing = np.array(df['northing'])
    
    if atl03struct:
        R_mat = atl03struct.rotationData.R_mat
        xRotPt = atl03struct.rotationData.xRotPt
        yRotPt = atl03struct.rotationData.yRotPt
        desiredAngle = 90
        crossTrack, alongTrack, R_mat, xRotPt, yRotPt, phi = \
        getCoordRotFwd(easting, northing, R_mat, xRotPt, yRotPt, [])    
    else:
        desiredAngle = 90
        crossTrack, alongTrack, R_mat, xRotPt, yRotPt, phi = \
        getCoordRotFwd(easting, northing, [], [], [], desiredAngle)

    if 'crosstrack' not in list(df.columns):
        df = pd.concat([df,pd.DataFrame(crossTrack,
                                    columns=['crosstrack'])],axis=1)
        df = pd.concat([df,pd.DataFrame(alongTrack,
                                        columns=['alongtrack'])],axis=1)
    else:
        print('Warning: Overwritting Existing Alongtrack/Crosstrack')
        df = df.drop(columns = ['crosstrack'])
        df = df.drop(columns = ['alongtrack'])
        df = pd.concat([df,pd.DataFrame(crossTrack,
                                        columns=['crosstrack'])],axis=1)
        df = pd.concat([df,pd.DataFrame(alongTrack,
                                        columns=['alongtrack'])],axis=1)  

    
    rotation_data = AtlRotationStruct(R_mat, xRotPt, yRotPt, desiredAngle, phi)
    
    return df, rotation_data
    
def get_atl03_df(atl03filepath, atl08filepath, gt, epsg = None):
    df = read_atl03_heights_data(atl03filepath, gt)
    df = get_atl03_classification(atl03filepath, atl08filepath, df, gt)
    df = get_atl_time(df)
    df, epsg = get_atl_coords(df, epsg = None)
    df, rotationData = get_atl_alongtrack(df)
    return df

def get_direction(lat):
    if(np.abs(lat[-1]) > np.abs(lat[0])):
        track_direction = 'Ascending'
    else:
        track_direction = 'Descending'
    return track_direction

def get_file_name(filepath):
    filepath = os.path.normpath(os.path.abspath(filepath))
    filename = os.path.splitext(os.path.basename(filepath))[0]
    return filename
  
def get_kml_region(lat,lon, kml_bounds_txt):
    # Determine if ATL03 track goes over a lidar truth region
    kmlBoundsTextFile = kml_bounds_txt
    kmlRegionName = False
    headerFilePath = False
    truthFilePath = False
    
    if kmlBoundsTextFile and (os.path.exists(kmlBoundsTextFile)):
        
        # Message to user
        print('   Finding Truth Region...')
        
        # Read kmlBounds.txt file and get contents
        kmlInfo = readTruthRegionsTxtFile(kmlBoundsTextFile)
        
        # Loop through kmlBounds.txt and find matching TRUTH area
        maxCounter = len(kmlInfo.regionName) - 1
        counter = 0
        while(not kmlRegionName):
            latInFile = (lat >= kmlInfo.latMin[counter]) & \
                (lat <= kmlInfo.latMax[counter])
            lonInFile = (lon >= kmlInfo.lonMin[counter]) & \
                (lon <= kmlInfo.lonMax[counter])
            trackInRegion = any(latInFile & lonInFile)
            if(trackInRegion):
                
                # Get truth region info
                kmlRegionName = kmlInfo.regionName[counter]
                headerFilePath = \
                    os.path.normpath(kmlInfo.headerFilePath[counter])
                truthFilePath = \
                    os.path.normpath(kmlInfo.truthFilePath[counter])
                
                # Print truth region
                print('   Truth File Region: %s' % kmlRegionName)
            
            if(counter >= maxCounter):
                print('   No Truth File Region Found in kmlBounds.txt')
                break
            counter += 1           
    else:
        kmlRegionName = None
        headerFilePath = None
        truthFilePath = None
        
            # Could not read kmlBounds.txt file
        
    return kmlRegionName, headerFilePath, truthFilePath
    
def write_atl03_las(atlstruct, outpath):
    xx = np.array(atlstruct.df.easting)
    yy = np.array(atlstruct.df.northing)
    zz = np.array(atlstruct.df.h_ph)
    cc = np.array(atlstruct.df.classification)
    ii = np.array(atlstruct.df.signal_conf_ph)
    sigconf = np.array(atlstruct.df.signal_conf_ph)
    hemi = atlstruct.hemi
    zone = atlstruct.zone
    
    print('   Writing ATL03 .las file...', end = " ")
    try:
        outname = atlstruct.atl03FileName + '_' + atlstruct.gtNum + '.las'
    except AttributeError:
        # occasionally atl03FileName is not in the atlstruct
        outname = atlstruct.atlFileName + '_' + atlstruct.gtNum + '.las'

    outfile = os.path.normpath(outpath + '/' + outname)
    
    if(not os.path.exists(os.path.normpath(outpath))):
        os.mkdir(os.path.normpath(outpath))
    
    # Get projection
    if(atlstruct.zone=='3413' or atlstruct.zone=='3976'):
        # 3413 = Arctic, 3976 = Antarctic
        lasProjection = atlstruct.hemi
        # Write .las file
        writeLas(xx,yy,zz,lasProjection,outfile,cc,ii,sigconf)

    else:
        # Write .las file for UTM projection case
        writeLas(xx,yy,zz,'utm',outfile,cc,ii,sigconf,hemi,zone)

    print('Complete') 

def get_H5_keys_info(atl08filepath,gt):
    keys = getH5Keys(atl08filepath, gt)
    h = h5py.File(atl08filepath, 'r')
    key_name = []
    key_type = []
    key_len = []
    for key in keys:
        try:
            data = h[gt + '/' + key]
            kname = str(key)
            ktype = str(data.dtype)
            klen = int(len(data))
            key_name.append(kname)
            key_type.append(ktype)
            key_len.append(klen)
        except:
            kname = str(key)
            ktype = 'Group'
            klen = 0
            key_name.append(kname)
            key_type.append(ktype)
            key_len.append(klen)
    key_info = [list(a) for a in zip(key_name, key_type, key_len)]
    return key_info

def match_atl_to_atl03(df, atl03struct):
    # Calculate Time
    delta_time03 = np.array(atl03struct.df.delta_time)
    time = np.array(df.delta_time) -\
        np.min(delta_time03[np.nonzero(delta_time03)])
    df = pd.concat([df,pd.DataFrame(time,
                            columns=['time'])],axis=1)
    # Calculate Projected Coordinates
    df, epsg = get_atl_coords(df, atl03struct.epsg)
    # Calculate Along track
    df, rotation_data = get_atl_alongtrack(df, atl03struct)        
    return df, rotation_data, epsg
        
def get_atl03_struct(atl03filepath, gt, atl08filepath = None, epsg = None, 
                     kml_bounds_txt = None, header_file_path = None):
    df = read_atl03_heights_data(atl03filepath, gt)
    if atl08filepath:
        try:
            df = get_atl03_classification(atl03filepath, atl08filepath, df, gt)
            dataIsMapped = True
        except:
            dataIsMapped = False
    else:
        dataIsMapped = False
    df = get_atl_time(df)
    df, epsg = get_atl_coords(df, epsg)
    df, rotation_data = get_atl_alongtrack(df)
    track_direction = get_direction(np.array(df.lat_ph))
    atl03filename = get_file_name(atl03filepath)
    atl03_info = getNameParts(atl03filename)
    hemi, zone = identify_hemi_zone(epsg)

    kml_region_name, header_file_path_out, truth_file_path = \
        get_kml_region(np.array(df.lat_ph),np.array(df.lon_ph), kml_bounds_txt)

    if header_file_path:
        del header_file_path_out
    else:
        header_file_path = header_file_path_out

    if dataIsMapped:
        # for some reason, need to reset nan to -1,
        # evemn though it's done in get_atl03_classification..
        c = np.array(df.classification)
        nan_index = np.where(np.isnan(c))
        c[nan_index] = -1 # assign nan to -1
        c = c.astype(int)
        df.classification = c
        # df.replace({'classification' : np.nan}, -1)

    # Assign everything to the struct
    atl03Struct = AtlStruct(df, gt, epsg, zone, hemi,kml_region_name, 
                              header_file_path, truth_file_path, 
                              atl03filepath, atl03filename, track_direction,
                              'ATL03', atl03_info, dataIsMapped, rotation_data)
    
    if header_file_path:
        headerData = readHeaderMatFile(header_file_path)
    else:
        headerData = None
    setattr(atl03Struct,'headerData',headerData)

        
    return atl03Struct

def get_atl08_struct(atl08filepath, gt, atl03struct = None, epsg = None, 
                     kml_bounds_txt = None):
    df = read_atl08_land_segments(atl08filepath, gt)
    
    # If ATL03 Struct is available
    if atl03struct:
        # Calculate Time
        df, rotation_data, epsg = match_atl_to_atl03(df, atl03struct) 
        kml_region_name = atl03struct.kmlRegionName
        header_file_path = atl03struct.headerFilePath
        truth_file_path = atl03struct.truthFilePath

    # If ATL03 Struct is not available
    else:
        # Calculate Time
        df = get_atl_time(df)
        # Calculate Projected Coordinates
        df, epsg = get_atl_coords(df, epsg)    
        # Calculate Along Track    
        df, rotation_data = get_atl_alongtrack(df)
        kml_region_name, header_file_path, truth_file_path = \
            get_kml_region(np.array(df.longitude), np.array(df.latitude),
                           kml_bounds_txt)
    track_direction = get_direction(np.array(df.latitude))
    atl08filename = get_file_name(atl08filepath)
    atl08_info = getNameParts(atl08filename)
    hemi, zone = identify_hemi_zone(epsg)
    dataIsMapped = True

    # Assign everything to the struct
    atl08Struct = AtlStruct(df, gt, epsg, zone, hemi,kml_region_name, 
                              header_file_path, truth_file_path, 
                              atl08filepath, atl08filename, track_direction,
                              'ATL08', atl08_info, dataIsMapped, rotation_data)
    return atl08Struct

def get_atl08canopy_struct(atl08filepath, gt, atl03struct = None, epsg = None, 
                     kml_bounds_txt = None):
    df = read_atl08_canopy_h(atl08filepath,gt)
    # If ATL03 Struct is available
    if atl03struct:
        # Calculate Time
        df, rotation_data, epsg = match_atl_to_atl03(df, atl03struct) 
        print('echo')
        kml_region_name = atl03struct.kmlRegionName
        header_file_path = atl03struct.headerFilePath
        truth_file_path = atl03struct.truthFilePath
        
    # If ATL03 Struct is not available
    else:
        # Calculate Time
        df = get_atl_time(df)
        print('echo')
        # Calculate Projected Coordinates
        df, epsg = get_atl_coords(df, epsg)    
        # Calculate Along Track    
        df, rotation_data = get_atl_alongtrack(df)
        kml_region_name, header_file_path, truth_file_path = \
            get_kml_region(np.array(df.longitude), np.array(df.latitude),
                           kml_bounds_txt)
    track_direction = get_direction(np.array(df.latitude))
    atl08filename = get_file_name(atl08filepath)
    atl08_info = getNameParts(atl08filename)
    hemi, zone = identify_hemi_zone(epsg)
    dataIsMapped = True
   

    # Assign everything to the struct
    atl08Struct = AtlStruct(df, gt, epsg, zone, hemi,kml_region_name, 
                              header_file_path, truth_file_path, 
                              atl08filepath, atl08filename, track_direction,
                              'ATL08', atl08_info, dataIsMapped, rotation_data)
    return atl08Struct


def get_atl09_struct(atl09filepath, gt, atl03struct = None, epsg = None, 
                     kml_bounds_txt = None):
    df, ds_layers, ds_va_bin_h, cab_prof, density_pass1, density_pass2 =\
        read_atl09_hr_profile(atl09filepath, gt)
    
    # If ATL03 Struct is available
    if atl03struct:
        df, rotation_data, epsg = match_atl_to_atl03(df, atl03struct) 
        kml_region_name = atl03struct.kmlRegionName
        header_file_path = atl03struct.headerFilePath
        truth_file_path = atl03struct.truthFilePath

    # If ATL03 Struct is not available
    else:
        # Calculate Time
        df = get_atl_time(df)
        # Calculate Projected Coordinates
        df, epsg = get_atl_coords(df, epsg)    
        # Calculate Along Track    
        df, rotation_data = get_atl_alongtrack(df)
        kml_region_name, header_file_path, truth_file_path = \
            get_kml_region(np.array(df.longitude), np.array(df.latitude),
                           kml_bounds_txt)
    track_direction = get_direction(np.array(df.latitude))
    atl09filename = get_file_name(atl09filepath)
    atl09_info = getNameParts(atl09filename)
    hemi, zone = identify_hemi_zone(epsg)
    dataIsMapped = True
    ancillary = read_atl09_ancillary_data(atl09filepath)
    orbit_info =  read_atl09_orbit_info(atl09filepath)
    # Assign everything to the struct
    atl09Struct = AtlStruct(df, gt, epsg, zone, hemi,kml_region_name, 
                              header_file_path, truth_file_path, 
                              atl09filepath, atl09filename, track_direction,
                              'ATL09', atl09_info, dataIsMapped, rotation_data,
                              ancillary, orbit_info)
    setattr(atl09Struct,'ds_layers', ds_layers)
    setattr(atl09Struct,'ds_va_bin_h', ds_va_bin_h)
    setattr(atl09Struct,'cab_prof', cab_prof)
    setattr(atl09Struct,'density_pass1',density_pass1)
    setattr(atl09Struct,'density_pass2',density_pass2)

    return atl09Struct

def get_geolocation_mapping(height_len, ph_index_beg, segment_ph_cnt, target):
    data = np.zeros(height_len)
    for i_id in range(0,len(target)):
        data[ph_index_beg[i_id]-1:\
             ph_index_beg[i_id]-1 + segment_ph_cnt[i_id]] =\
             np.full((segment_ph_cnt[i_id]), target[i_id]) 
    return data

def append_atl03_geolocation(heights, geolocation, fields = ['segment_id']):
    height_len = len(heights)
    ph_index_beg = np.array(geolocation.ph_index_beg)
    segment_ph_cnt = np.array(geolocation.segment_ph_cnt)
    for field in fields:
        target = np.array(geolocation[field])
        data = get_geolocation_mapping(height_len, ph_index_beg, 
                                       segment_ph_cnt, target)
        heights = pd.concat([heights,pd.DataFrame(data,
                                columns=[field])],axis=1)
    
    return heights
    
def convert_atl03_to_legacy(atl03):
    intensity = np.zeros(len(atl03.df))
    atl03h5Info = getNameParts(atl03.atlFileName)
    atl03legacy = Atl03StructLegacy(atl03.df.lat_ph, atl03.df.lon_ph, 
                              atl03.df.easting, atl03.df.northing, 
                              atl03.df.crosstrack, atl03.df.alongtrack, 
                              atl03.df.h_ph, atl03.df.time, 
                              atl03.df.delta_time, atl03.df.signal_conf_ph, 
                              atl03.df.classification, intensity, atl03.gtNum,
                              atl03.zone, atl03.hemi, atl03.kmlRegionName, 
                              atl03.headerFilePath, atl03.truthFilePath, 
                              atl03.atlFileName, atl03.atlFileName, 
                              atl03.trackDirection, atl03h5Info, 
                              atl03.dataIsMapped)
    rotationData = atl03.rotationData
    headerData = atl03.headerData
    return atl03legacy, rotationData, headerData

def write_pickle(data, filename):
    import pickle
    fp = open(filename, 'wb')
    pickle.dump(data, fp)
    fp.close()

def read_pickle(filename):
    import pickle
    fp = open(filename, 'rb')
    data = pickle.load(fp)
    fp.close()
    return data

def convert_df_to_mat(df,outfilename):
    from scipy import io
    comps =  outfilename.split('.')
    if comps[-1] != 'mat':
        outfilename = outfilename + ".mat"
    # scipy.io.savemat(outfilename, {'struct':df.to_dict("list")})
    io.savemat(outfilename, {'struct':df.to_dict("list")})

def lightweight_atl_h5(heights, basepath, atl_light_name):
    save_light_location = basepath + atl_light_name
    heights_light = heights.loc[heights['classification']>0]
    heights_light.to_hdf(save_light_location, 
                                       key = 'heights_light')
    sll = save_light_location[:-3]
    convert_df_to_mat(heights_light, sll)
    return save_light_location, heights_light

def lightweight_atl_csv(heights, basepath, atl_light_name):
    save_light_location = basepath + atl_light_name
    heights_light = heights.loc[heights['classification']>0]
    sll = save_light_location[:-3]

    heights_light.to_csv(sll)
    # sll = save_light_location[:-3]
    #convert_df_to_mat(heights_light, sll)
    return save_light_location, heights_light

def geoid_corrector(heights):
    heights_corr = heights
    heights_corr[heights_corr.geoid>10e10].geoid = np.nan
    return heights_corr

def rolling_variance_filter(dataframe, window_size, variance_threshold_stringency):
    canopy_data = dataframe[dataframe.classification == 3].h_ph
    canopy_mov_std = canopy_data.rolling(window_size, win_type='boxcar', 
                                          center = True).std()
    variance_threshold = variance_threshold_stringency * np.median(canopy_mov_std.dropna())
    cf = dataframe.loc[canopy_mov_std<variance_threshold]
    time_std_index = dataframe['delta_time'].loc[canopy_mov_std<variance_threshold]
    return cf, canopy_mov_std, time_std_index

def ground_neighbor_filter(dataframe, height_diff_threshold):
    segs = np.unique(dataframe.segment_id)#.astype('int')
    segs_bin = pd.DataFrame(index = segs)
    for ii in segs:
        bool_seg = dataframe.segment_id == ii
        segs_bin.loc[ii,'segment_id'] = ii
        segs_bin.loc[ii,'mean_h'] = np.nanmean(dataframe[bool_seg & (dataframe.classification==1)].h_ph)
    
    dataframe_merged = pd.merge(dataframe, segs_bin, how='outer', left_on='segment_id', right_on='segment_id')
    dataframe_merged['ground_diff'] = dataframe_merged.h_ph- dataframe_merged.mean_h
    dataframe_merged.classification.loc[(dataframe_merged.classification == 2) & (dataframe_merged.ground_diff>height_diff_threshold)] = -2   #dummy value to change in future versions
    return dataframe_merged

def append_atl03_modislandcover(heights_i, segment_landcover):
    # segsf= lambda x : np.arange(x.segment_id_beg, x.segment_id_end+1)
    # segsary=segment_landcover.apply(segsf, axis=1)
    # segsary.name='segsary'
    # segment_landcover['segment_id'] = segsary   
    segment_landcover['segment_id'] = segment_landcover.apply(lambda x:  np.arange(x.segment_id_beg, x.segment_id_end+1), axis=1)
    expanded_df = segment_landcover.explode('segment_id')
    heights_expanded = pd.merge(heights_i, expanded_df, how='outer', left_on='segment_id', right_on='segment_id')
    
    return heights_expanded

def legacy_name(df):
    df_legacy_names = df.rename(columns={'Delta Time (sec)':'delta_time', # This is sloppy, but needed
                        'Height (m HAE)':'h_ph', 'Latitude (deg)':'lat_ph',
                        'Longitude (deg)':'lon_ph', 'Classification':'classification',
                        'Time (sec)':'time', 'UTM Easting (m)':'easting',
                          'UTM Northing (m)':'northing','Geoid (m)':'geoid',
                          'refDEM (m)':'dem_h', 'Segment ID':'segment_id'})
    return df_legacy_names

def gui_name(df):
    df_gui_names = df.rename(columns={'delta_time':'Delta Time (sec)',
                        'h_ph':'Height (m HAE)', 'lat_ph':'Latitude (deg)',
                        'lon_ph':'Longitude (deg)', 'classification':'Classification',
                        'time':'Time (sec)', 'easting':'UTM Easting (m)',
                          'northing':'UTM Northing (m)','geoid':'Geoid (m)',
                          'dem_h':'refDEM (m)', 'segment_id':'Segment ID'})
    return df_gui_names

def get_attribute_info(atlfilepath, gt):
    # add year/doy, sc_orient, beam_number/type to 08 dataframe
    year, doy = get_h5_meta(atlfilepath, meta='date', rtn_doy=True)

    with h5py.File(atlfilepath, 'r') as fp:
        try:
            fp_a = fp[gt].attrs
            description = (fp_a['Description']).decode()
            beam_type = (fp_a['atlas_beam_type']).decode()
            atlas_pce = (fp_a['atlas_pce']).decode()
            spot_number = (fp_a['atlas_spot_number']).decode()
            atmosphere_profile = (fp_a['atmosphere_profile']).decode()
            groundtrack_id = (fp_a['groundtrack_id']).decode().lower()
            sc_orient = (fp_a['sc_orientation']).decode().lower()
        except:
            description = ''
            beam_type = ''
            atlas_pce = ''
            spot_number = ''
            atmosphere_profile = ''
            groundtrack_id = ''
            sc_orient = ''
    info_dict = {
        "description" : description,
        "atlas_beam_type" : beam_type,
        "atlas_pce" : atlas_pce,
        "atlas_spot_number" : spot_number,
        'atmosphere_profile' : atmosphere_profile,
        "groundtrack_id" : groundtrack_id,
        "sc_orientation" : sc_orient,
        "year" : year,
        "doy" : doy
        # "gps_time_offset" : gps_time_offset
        }
    
    return info_dict
