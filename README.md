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

Now, instead of leave the energy stored in Ti:Sapphire, we fire another laser "CARP" after main laser pulse exits MPA2. CARP will take leftover energy from Ti:Sapphire crystal. CARP is an OPO laser with pulse energy ~ 100 mJ and pulse duration ~ 6 ns. We allow CARP to pass through Ti:Sapphire 5 times. Below is the plot of amplification of CARP inside MPA2.

![CARP amplification in MPA2](https://github.com/patinkaew/Project-SCALLOPS/blob/master/pics/carp_MPA2.png)

Amplified CARP light will be dumped after. Then, the reflected light (1% from main light) arrive back to MPA2. The plot shows, with laser CARP, the gain decrease to < 1.5.

![reflected amplification in MPA2 after CARP](https://github.com/patinkaew/Project-SCALLOPS/blob/master/pics/reflected_carp_MPA2.png)

Since one pass will take about 10 ns for CARP to go through the crystal, for 5 passes of CARP, we know that it would take ~ 50 ns to finish taking off the energy. However, it takes ~ 100 ns for reflected light to get back to MPA2. Thus, this setup is plausible.

### Long-pulse laser system

At MEC, long-pulse laser consists of multiple amplification steps. The front end laser produce initial light and then it gets amplified in Nd:Glass (diameter 25 mm). Then, amplified laser light is splitted into 4 beams. Each beam gets further amplified in Nd:Glass (diameter 50 mm) again. Lastly, amplified four beams recombine into one beam. Code's sequence resembles these processes and the simple flowchart is shown below.

![long pulse flowchart](https://github.com/patinkaew/Project-SCALLOPS/blob/master/pics/long_pulse_flowchart.png width = 20)

First, we will demonstrate 'gain shifting'. Suppose that we seed with gaussian pulse. After amplifications, we could see that the peaks get shifted earlier because the earlier parts of the pulse will take some energy stored in crystal. Earlier parts will see more energy stored in crystal, compared to later parts. Thus, earlier parts will get more amplified.

![gain shifting](https://github.com/patinkaew/Project-SCALLOPS/blob/master/pics/long_gain_shift.png)

For this reason, suppose we want the output to be a square pulse, if we seed with a square pulse, we would not get square pulse back.

![long pulse square input](https://github.com/patinkaew/Project-SCALLOPS/blob/master/pics/long_output1.png)

Now, we could try seed with different pulse shape to make the output look like sqaure pulse.

![long pulse square output](https://github.com/patinkaew/Project-SCALLOPS/blob/master/pics/long_input1.png)

### Short-pulse laser system

![short pulse flowchart](https://github.com/patinkaew/Project-SCALLOPS/blob/master/pics/short_pulse_flowchart.png)

## Acknowledgements
