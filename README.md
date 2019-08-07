#  Project SCALLOPS

## Introduction

SCALLOPS (Simulation of Crystal Amplification of Laser Light and Optical Pulse Shaping) is a software to help simulate the development of laser pulse via amplification in crystal. The project is developed as a part of my work as a LCLS summer intern in 2019. The primary purpose is to use this simulation for long pulse (ns) and short pulse (fs) laser systems at Matter in Extreme Conditions (MEC), one of the experimental hutches at Linac Coherent Light Source (LCLS) at SLAC National Laboratory. The simulation could also be used in future projects, such as Laser CARP and potentially MEC's novel peta watt laser system.

## Installation

The code utilizes the following python packages:
- [Numpy](https://numpy.org/)
- [Scipy](https://www.scipy.org/) 
- [Pandas](https://pandas.pydata.org/) 
- [Pint](https://pint.readthedocs.io/en/0.9/)

I recommend install [Anaconda](https://www.anaconda.com/) distribution which contains everything except Pint.

## Examples

### Laser CARP

Laser CARP (Collecter of At-Risk Photons) is MEC's novel project to further ensure safety of short-pulse laser system. Additional laser "CARP" will be installed to MPA 2 of short-pulse laser to take off leftover energy inside MPA2's Ti:Sapphire crystal after the main laser pulse got amplified. This simulation will help determine the plausibility of CARP.

First, amplify the main laser pulse. Evolution of energy inside pulse and energy in crystal are plotted below. We can see that there are significant energy left inside the crystal after.

![laser amplification in MPA2](https://github.com/patinkaew/Project-SCALLOPS/blob/master/pics/laser_MPA2.png "laser amplification in MPA2")

Suppose that 1% of this amplified laser pulse got reflected back and travel back to MPA2. Due to the setup of MPA2, it will experience 3 passses of amplification. Evolution of energy inside pulse and energy in crystal are plotted below.

![reflected light amplification in MPA2](https://github.com/patinkaew/Project-SCALLOPS/blob/master/pics/reflected_MPA2.png "reflected light amplification in MPA2")

The plot shows that there is still ~ 13 gain inside the crystal. 

Now, instead of leave the energy stored in Ti:Sapphire, we fire another laser "CARP" after main laser pulse exits MPA2. CARP will take leftover energy from Ti:Sapphire crystal.

### Long-pulse laser system

## Acknowledgement
