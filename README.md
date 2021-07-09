# ICESat-2 Tiles
## A Convenient Processed Format for ICESat-2 Data

This module provides the framework for ingesting ATL03 and ATL08 data into a tiled format in compressed LAZ format. It assumes local access to ATL03 and ATL08 data products from [ICESat-2](https://icesat-2.gsfc.nasa.gov/). The code is built off the sourcecode for [PhoReal](https://github.com/icesat-2UT/PhoREAL).

### Getting Started

A working environment can be found in /envs/requirements.txt, though there are several packages that may need to be installed manually. Even after installing the requirements.txt, we ***strongly recommend*** that you manually install pylas with the command:
> pip install pylas[lazrs]
> 
This ensures the correct pylas version and LAZ-compression backend is built. Most of the compression errors we run into can be attributed to a build issue with pylas.
We also recommend installing python-geohash with care. It can have compatibility issues with Microsoft Visual C++. 
Presently, we recommend running the tile-creating program in a separate python environment than the tile-viewer. Some of the geospatial software dependancies are fragile. 

### Running the code 
To run locally, the user must change the example_local.inp file to reflect their internal directory structure.

The code can be executed by navigating to the app/ directory and running:
> python local_tiling.py
> 
The code will take input variables from the example_local.inp and begin processing ATL03 and ATL08 files into the specified tile size. 

#### Notes on the verbose output:
If the .h5 data is successfully processed into tiles, the script outputs: 
>Success on [.h5 file] beam [beam name] index: [index of file in directory]

The notable function execution times are also printed, as well as the names of the successfully created tiles.

Failures are typically a result of the heavy filtering of photon data before tiling. Empty data will cause the script to fail, particularly in the geohashing functions that rely on geolocated photons.
Another failure mode is if the ALT03 does not have a corresponding ATL08 file. The script will simply skip it, and move on to the next file.

### MegaTile Data Format
The data within the “MegaTiles” are stored in LAZv 1.4, point format 6. This format is self-describing, with formal fields intrinsic to the version. However, one can also append extra fields for each point, making this version far more adaptable. By adopting LAZ 1.4, the bulk data size of the ICESat-2 data is reduced by a factor of 10 from the initial h5 format to LAZ. ICESat-2 collects ~100Tb annually and we expect the tiled global dataset to be approximately an order of magnitude smaller. Individual file size can be changed to user needs, though the maximum recommended tile size is 5 degrees by 5 degree. At that size, the files are large enough that a typical desktop would have difficulty reading into memory. Nominally, a single 5x5 LAZ file is on the order of 150Mb, but can be more than 500Mb if the tile is densely covered. Currently, we are recommending a tile size of 1x1 degree or 2x2 degrees for desktop processing.

NGA-MegaTile LAZ fields are:

"ClassFlags": 	empty

"Classification"	ATL08 classification of individual photon,

"EdgeOfFlightLine": 	empty

"GpsTime"	empty

"Intensity": 	empty

"NumberOfReturns": 	empty

"PointId": 	Point number in each MegaTile

"PointSourceId": 	empty

"ReturnNumber": 	

"ScanAngleRank": 	

"ScanChannel": 	

"ScanDirectionFlag": 	

"UserData": 	

"X": 	WGS84 Longitude.  Longitude can be reprojected into UTM coordinates using PhoREAL MTtools.

"Y": 	WGS84 latitude. Latitude can be reprojected into UTM coordinates using PhoREAL MTtools.

"YYMMDD": 	Year, Month, Day (YYMMDD),

"Z": 	Elevation (m) above the WGS84 ellipsoid,

"beam_number": 	ATLAS beam number.  Beams 1, 3, 5 are strong beams; Beams 2, 4, 6 are the weak beams

"delta_time": 	ATLAS time of each outgoing laser shot. The “delta_time” parameter is consistent across all ICESat-2 data products,

"geoid": 	Geoid height (m) of EGM2008,

"landcover": 	landcover value from MODIS (rel 005 will use Copernicus LC),

"msw_flag": 	Multiple Scatter Warning flag from ATL09. Value of 0 indicate no observed scattering (clear skies). Values greater than 0 indicate observed scattering at different heights in the atmosphere

"ph_h": 	ATL03 Photon height (m) above estimated ground,

"signal_conf":,	ATL03 Signal Confidence

"snow_flag": 	NOAA snow flag,

"solar_elevation": 	Solar Elevation (Negative values indicate sun is below horizon; positive values sun is above horizon)

"track_num": 	ICESat-2 Orbit Number,

"yrdoy": 	Year and Day of Year (YYDDD)

### MegaTile Naming Convention
At present the MegaTile file naming convention is:
ATL_tile_<ICESat-2 RELEASE NUMBER>_<ICESat-2 CYCLE NUMBER>_<NORTH EAST TILE CORNER>_<SOUTH WEST TILE CORNER>_<TILE EDGE LENGTH in Degrees>.laz

For example:
ATL_tile_r004_cycle06_60N30E_55N25E_5.laz
