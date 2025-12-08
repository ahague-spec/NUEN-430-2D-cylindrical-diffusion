#start of code, import some stuff
import numpy as np

#input data Rod and clad dimensions from: https://www.researchgate.net/figure/BWR-fuel-assembly-dimensions-Fensin-2004-Mueller-et-al-2013a_tbl1_323487844
Cylinder_Radius = 5.13e-3  #In meters, done this way to result in lower numbers to avoid potential overflows. This represents fuel + gap + clad radius
Fuel_ Height = 3.81 #In meters
#Everything past the cladding is assumed to be coolant
235_Abs_CS= "PLaceholder" #Absorption cross section of 235 isotope UO2. PLaceholder for now
238_Abs_CS= "Placeholder" #Asked the professor for reference material, he might be able to direct us to something.
#Will be some matrices here for scattering cross section

#Also gonna get data for fission cross sections, 
"""
It might actually be a good idea to create code in a spyder file, then paste it here.
Data from this source: https://www.oecd-nea.org/science/wprs/eg3drtb/NEA-C5G7MOX.PDF
There are 7 total energy groups, with tables of data provided for each of them that I will port over Monday evening
Also, this data homogenizes the fuel and clad
"""







