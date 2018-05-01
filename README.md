# truss-model-updating

## Description

This program was written to show a numerical model updating example on 3D trusses.

## Requirements

* Python ^3.6
* External packages:
    * numpy
    * pyserial
    * matplotlib
    * mpl_toolkits

## Quick start

The program can be run the following way:

    python3 truss.py <input_file.str> [project_title] [compatibility_mode] [simulation]
    python3 truss.py truss.str -t example -c 3 -s 0  

Simplest run [defaults: -t structure -c 3 -s 0]:

    python3 truss.py bridge.str
    
For help:

    python 3 truss.py -h

## Configuration

Many settings help customization. These settings are collected in predefined profiles as written below.
The following options are available in the configuration process:

* Logging: Puts log data to the console
* Graphics: Creates and saves diagrams/pictures about the structure (uses external libraries)
* Numpy solver: Force uses the Numpy solver instead of the built-in one
* Oslib: OS file action availability
* Updating: Enables model updating
* Arduino: Enables Arduino input stream
* Debug: Speeds up runtime by some tweaks. DO NOT USE 'in production'
* Realistic simulations: Opens the backed up input stream and fetches data with realistic delays/timing according to the timestamps.

### Compatibility modes

Mode | Mode's name | Logging | Graphics | Numpy solver | OSlib | Updating | Arduino | Debug | Realistic simulation*
:-------: | :---------: | :-----: | :------: | :------: |:------:| :------: | :------: | :------: | :------:
**0** | User defined | ✔ | ✔ | ✔ | ✔ |   |   | ✔ |   
**1** | Informative | ✔ | ✔ | ✔ | ✔ | ✔ | ✔ |   | ✔ 
**2** | Maximum compatibility |  |  |  |  |  |  |  | ✔
**3** | Android mode | ✔ |  |  | ✔ |  |  |  | ✔  


