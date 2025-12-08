#start of code, import some stuff
import numpy as np

#input data Rod and clad dimensions from: https://www.researchgate.net/figure/BWR-fuel-assembly-dimensions-Fensin-2004-Mueller-et-al-2013a_tbl1_323487844
Fuel_Radius = 4.38e-3  #In meters, done this way to result in lower numbers because I want to 
Fuel_ Height = 3.81 #In meters
Gap_Thickness = 9e-5 #in meters, this region we can assume is empty space for purposes of neutron interactions
Clad_Thickness = 6.6e-4 #in meters
#Everything past the cladding is assumed to be coolant
Enrichment = 0.05 #Fraction of UO2 that's U-235. Remaining parts are U-238
235_Abs_CS= "PLaceholder" #Absorption cross section of 235 isotope UO2. PLaceholder for now
238_Abs_CS= "Placeholder" #Asked the professor for reference material, he might be able to direct us to something.


#Also gonna get data for fission cross sections, 
"""
It might actually be a good idea to create code in a spyder file, then paste it here.
Data from this source: https://www.oecd-nea.org/science/wprs/eg3drtb/NEA-C5G7MOX.PDF
There are 7 total energy groups, and we can assume that all neutrons are born fast (in group 1)
Also, this data homogenizes the fuel and clad
"""





