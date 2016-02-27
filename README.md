# nomads

This repository contains python code that runs an agent-based model of mobile pastoralists as well as relevant input files.

### What is this repository for? ###

* Quick summary: This repository was established for developing and disseminating a model for simulating the spatial behaviour of mobile pastoralists who live in arid and semi-arid environments. This model explicitly simulates the complex movement patterns of pastoralists and computes the resultant natural resource access that they can achieve under different circumstances. The repository also contains a series of empirical spatial data (GeoTIFF format) on monthly vegetation distribution in northeast Nigeria, so you can simulate pastoralist behaviour in this particular geographical setting.
* Version: 1.0 (September 24, 2015)

### How do I get set up? ###

* Summary of set up: the model (nomads.py) runs on Python 2.7 extended with other widely-used modules (see below). Before executing it, obtain PyCX simulator developed by [PyCX project](http://pycx.sourceforge.net) and place the implementation file (pycxsimulator.py) in the model's root directory. Then everything gets ready. Move to the model's root directly and enter 'python nomads.py' from the command line. You can control the simulation conditions by manipulating the control panel that appears as well as by altering the parameter values that are specified in the header section of the main source code.
* Dependencies: the model requires SciPy (0.15.1), Matplotlib (1.4.3), and GDAL Python package (2.0.0) for its running along with PyCX simulator motioned above.
* Note on input data: the input data (monthly vegetation distributions in northeast Nigeria) are placed in '/Input/NENGA' directory. These were derived from MODIS satellite imagery. MODIS datasets are publicly available from USGS's [EarthExplorer](http://earthexplorer.usgs.gov). In order to employ these data for simulation, you need to convert the data to GeoTIFF format.

### Known issues ###

### Who do I talk to? ###

* If there is something wrong with this distribution, contact [Takuto Sakamoto](mailto:takutos@mac.com).