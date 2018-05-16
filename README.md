# flagUtils

This repo holds all necessary scripts to reduce, calibrate, and perform basic analysis on FLAG spectral line data. This README provides basic info on each of the contained scripts as well as syntactically correct examples. These scripts are specifically set up to be used within the GBO computing environment. 


## Getting Started

Again, many of the constituent scripts are designed to work specifically within the GBO computing environment. Much of the directory paths are hard coded, so one must change these manually. To get set up, clone the repo within your desired directory:
git clone https://github.com/nipingel/flagUtils.git

### Prerequisites

These scripts are written for python, GBTIDL, and bash. 

### Installing

only a simple 'git clone' should be necessary

## Authors

* **Nickolas Pingel** - *Initial work* - [nipingel](https://github.com/nipingel)

## Scripts

### calcSysFlux_Grid.pro

This GBTIDL script averages the OFF power from the six offs associated with a calibration grid observation, while also computing the maximum power in associated ON scans of the grid (Poff and Pon, respectively). The maximum power value is taken to be the ON associated power. The system equivalent flux density is computed by using Equation 1 from Perley and Butler 2017 based on the user provided calSource name and is given in Jy. The equation is: 
Ssys = Ssrc*Poff/(Pon-Poff)
