# truss-model-updating

## Description

This program was written to show a numerical model updating example on 3D trusses. 

### Features

The program has 2 major features:

#### 1) Truss Solver

General program for solving a truss with the following attributes:

* Geometry (nodes, elements, cross-sections)
* Materials (E)
* Supports
* Forces on the nodes (no moments)

#### 2) Model Updating

Using the `Truss Solver`, considering the given loads, the displacements can be calculated. Of course
if we add measurements, a delta vector can be computed from the calculated and measured displacements.
Using this, by modifying the used numerical model the error vector can be reduced. This part of the code
tries to minimize the error vector by iterating. The resulted model is called updated model.  

## Technical Requirements

* Python 3.6 (or later)
* External packages:
    * numpy
    * pyserial
    * matplotlib
    * mpl_toolkits

## Quick Start

### Normal run
The program can be run the following way:

    python3 truss.py [project_title] [compatibility_mode] [simulation] <input_file.str>
    python3 truss.py -t example -c 1 -s 0 truss.str  

#### Simplest run
*[defaults: t=<input_file's name> c=2 s=0]*

    python3 truss.py bridge
    
The ".str" ending is optional at the arguments.
    
#### Model updating example

* with Arduino serial input:

        python3 truss.py -c 1 -s 0 bridge.str

* with simulation:

        python3 truss.py -c 1 -s 1 bridge.str
        
#### Help

    python3 truss.py -h

## Configuration

Many settings help customization. These settings are collected in predefined profiles as written below.
The following options are available in the configuration process:

* Logging: Prints log data to the console
* Graphics: Creates and saves diagrams/pictures about the structure (uses external libraries)
* NumPy solver: Force uses the Numpy solver instead of the built-in one
* Oslib: OS file action availability
* Updating: Enables model updating
* Arduino: Enables Arduino input stream
* Debug: Speeds up runtime by some tweaks. DO NOT USE 'in production'
* Realistic simulations: Opens the saved input stream and fetches data with realistic delays/timing according to the timestamps. Effective only with simulation ON.

These settings can be found in the truss_framework.py file below the following lines 

    if self.compatibility_mode == 0:
        ### User defined ###
        # Modify as needed #

If needed, please edit the User defined mode only.

### Compatibility Modes

Mode | Mode's name | Logging | Graphics | Numpy solver | OSlib | Updating | Arduino | Debug | Realistic simulation
:-------: | :---------: | :-----: | :------: | :------: |:------:| :------: | :------: | :------: | :------:
**0** | User defined* | ✔ | ✔ | ✔ | ✔ |   |   | ✔ |   
**1** | Informative | ✔ | ✔ | ✔ | ✔ | ✔ | ✔ |   | ✔ 
**2** | Maximum compatibility |  |  |  |  |  |  |  | ✔
**3** | Android mode** | ✔ |  |  | ✔ |  |  |  | ✔  

*User defined settings may vary according to the local configurations.

**Code can be run on Android devices with 3rd party applications, like [QPython3](https://play.google.com/store/apps/details?id=org.qpython.qpy3) or [	
QPy3.6.](https://play.google.com/store/apps/details?id=org.qpython.qpy36)

### Simulations

Generally, model updating is possible in this application by supplying measurement data by an Arduino device.
This approach has a few cons:

* Updating method can be run in real time only
* Arduino connection is required
* Experiment cannot be repeated, since the input stream is unique

As a workaround, the input stream can be recorded and saved. This backup stream is capable to be
reused by the simulation anytime later.

Turning on simulation opens the backup file and loads the measurement stream. This can be done in two different ways:
* fast-forward (default): processing the data as fast as possible
* realistic-mode: each measurement is given to the program in a scheduled way, using the timestamps. 
 