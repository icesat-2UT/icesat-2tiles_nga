# ICESat-2 Tiles
## A Convenient Processed Format for ICESat-2 Data

This module provides the framework for ingesting ATL03 and ATL08 data into a tiled format in compressed LAZ format. It assumes local access to ATL03 and ATL08 data products from [ICESat-2](https://icesat-2.gsfc.nasa.gov/). The code is built off the sourcecode for [PhoReal](https://github.com/icesat-2UT/PhoREAL).

### Getting Started

A working environment can be found in /envs/requirements.txt, though there are several packages that may need to be installed manually. Even after installing the requirements.txt, we ***strongly recommend*** that you manually install pylas with the command:
> pip install pylas[lazrs]
> 
This ensures the correct pylas version and LAZ-compression backend is built. Most of the compression errors we run into can be attributed to a build issue with pylas.
We also recommend installing python-geohash with care. It can have compatibility issues with Microsoft Visual C++.

### Running the code 
To run locally, the user must change the example_local.inp file to reflect their internal directory structure.

The code can be executed by navigating to the app/ directory and running:
> python local_tiling.py
> 
The code will take global variables from the example_local.inp and begin processing ATL03 and ATL08 files into the specified tile size. 

Notes on the verbose output:
If the .h5 data is successfully processed into tiles, the script outputs: 
>Success on [.h5 file] beam [beam name] index: [index of file in directory]

The notable function execution times are also printed, as well as the names of the successfully created tiles.

Failures are typically a result of the heavy filtering of photon data before tiling. Empty data will cause the script to fail, particularly in the geohashing functions that rely on geolocated photons.
Another failure mode is if the ALT03 does not have a corresponding ATL08 file. The script will simply skip it, and move on to the next file.
