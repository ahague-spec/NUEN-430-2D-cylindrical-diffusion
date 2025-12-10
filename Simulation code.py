#start of code, import some stuff
import numpy as np
import matplotlib.pyplot as plt

#input data Rod and clad dimensions from: https://www.researchgate.net/figure/BWR-fuel-assembly-dimensions-Fensin-2004-Mueller-et-al-2013a_tbl1_323487844
#All units are at cm-level
One_Group_Toggle=True #Toggle whether or not we're using one group or seven, for analytical testing
Mesh_Points=100 #How many mesh points are there going to be?
if One_Group_Toggle==True: #Set number of groups. Code could theoreticall accomodate any number, but we're only using 1 and 7
    Num_Groups= 1
else:
    Num_Groups = 7
    
Cylinder_Radius = 0.513 #Fuel and Clad are homogenized, this is in cm
"""Fuel_Height = 381 #Height of overall rod, doesn't matter due to 1D assumption
#Everything past the cladding is assumed to be coolant, fuel and clad are homogenized
#Properties for fuel/clad and moderator. Acquired from this source: https://www.oecd-nea.org/science/wprs/eg3drtb/NEA-C5G7MOX.PDF"""

""""#Diffusion coefficients are derived from transport cross section (1/(3*(Tr_CS + Abs_CS)))"""
  #Fuel/Clad data
Fuel_Clad_Scat_Matr=[[1.27537E-01, 4.23780E-02, 9.43740E-06, 5.51630E-09, 0.00000E+00, 0.00000E+00, 0.00000E+00],
[0.00000E+00, 3.24456E-01, 1.63140E-03, 3.14270E-09, 0.00000E+00, 0.00000E+00, 0.00000E+00],
[0.00000E+00, 0.00000E+00, 4.50940E-01, 2.67920E-03, 0.00000E+00, 0.00000E+00, 0.00000E+00],
[0.00000E+00, 0.00000E+00, 0.00000E+00, 4.52565E-01, 5.56640E-03, 0.00000E+00, 0.00000E+00],
[0.00000E+00, 0.00000E+00, 0.00000E+00, 1.25250E-04, 2.71401E-01, 1.02550E-02, 1.00210E-08],
[0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.29680E-03, 2.65802E-01, 1.68090E-02],
[0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 8.54580E-03, 2.73080E-01]]

Fuel_Clad_Fission_Spectrum=[0.5879,0.4118,3.391e-4,1.176e-7,0,0,0]
Fuel_Clad_Abs_CS=[8.025e-3,3.717e-3,2.677e-2,9.624e-2,3.002e-2,0.1113,0.2828]
Fuel_Clad_Total_CS=[2.12450E-01,3.55470E-01,4.85540E-01,5.59400E-01,3.18030E-01,4.01460E-01,5.70610E-01 ]
Fuel_Clad_Fission_CS=[7.212e-3,8.193e-4,6.453e-3,1.8565e-2,1.781e-2,8.303e-2,0.216]
Fuel_Clad_Diff_Coeffs=[1.7924,0.9994,0.6573,0.5123,0.9752,0.6582,0.3935]

  #Moderator data
Mod_Scat_Matr = [[4.44777E-02, 1.13400E-01, 7.23470E-04, 3.74990E-06, 5.31840E-08, 0.00000E+00, 0.00000E+00],
[0.00000E+00, 2.82334E-01, 1.29940E-01, 6.23400E-04, 4.80020E-05, 7.44860E-06, 1.04550E-06],
[0.00000E+00, 0.00000E+00, 3.45256E-01, 2.24570E-01, 1.69990E-02, 2.64430E-03, 5.03440E-04],
[0.00000E+00, 0.00000E+00, 0.00000E+00, 9.10284E-02, 4.15510E-01, 6.37320E-02, 1.21390E-02],
[0.00000E+00, 0.00000E+00, 0.00000E+00, 7.14370E-05, 1.39138E-01, 5.11820E-01, 6.12290E-02],
[00.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 1.32440E-01, 2.48070E+00]]
    #No fission in the moderator (duh)
Mod_Diff_Coeffs =[2.0858,0.8071,0.5644,0.5686,0.4606,0.2627,0.1240] #units are cm^-2
Mod_Abs_CS= [6.0105e-4,1.5793e-5,3.37160E-04,1.94060E-03,5.74160E-03,1.50010E-02,3.72390E-02]
Mod_Total_CS= [2.30070E-01,7.76460E-01,1.4842,1.5052,1.5592,2.0254,3.3057]

"""#Scattering is from column group to row group 
#(i.e 0.2823 is  group 2 -> group 2, 0.1324 is 7->6 etc.)
#Coordinates in calls are [rows down][columns right] from top left (Which is 0,0)
This is just reference material for me, cuz I tend to forget python does x,y in the opposite way it usually is"""

Greatest_Extents=np.empty(Num_Groups) #based on zero-edge assumptuion, find the zero-point for each energy group
for i in range(Num_Groups): #Loop populates extrapolation distance array. Contains positions for zero-flux assumption
    Greatest_Extents[i]=Cylinder_Radius+2*(Mod_Diff_Coeffs[i]) 
    
X_Points=np.linspace(0,max(Greatest_Extents),Mesh_Points) #10 points between centerline and greatest extent. Not all groups will reach this far
Flux_Array=np.zeros((Num_Groups,len(X_Points)))#Create empty array to store flux values at each point
#Loop just defines values for graph testing,
for i in range(Num_Groups): #Run this loop
    for j in range(Mesh_Points):
        if X_Points[j]>Greatest_Extents[i]:
            Flux_Array[i][j]=0
        else:
            Flux_Array[i][j]=i**2 + j**2
        
#Print out flux values for each group
plt.plot(X_Points,Flux_Array[0],color="red",label="Group 1 Flux")
if (One_Group_Toggle == False): #Only try printing out other groups if multi-group assumption happens
    plt.plot(X_Points,Flux_Array[1],color="orange",label="Group 2 Flux")
    plt.plot(X_Points,Flux_Array[2],color="yellow",label="Group 3 Flux")
    plt.plot(X_Points,Flux_Array[3],color="green",label="Group 4 Flux")
    plt.plot(X_Points,Flux_Array[4],color="blue",label="Group 5 Flux")
    plt.plot(X_Points,Flux_Array[5],color="indigo",label="Group 6 Flux")
    plt.plot(X_Points,Flux_Array[6],color="violet",label="Group 7 Flux")
plt.axvline(x=Cylinder_Radius, color='black',linestyle='--', label='Edge of Fuel/Clad') #Indicate edge of the cylinder
plt.xlabel("Radial distance from center (cm)")
plt.ylabel("Group Flux (n / (cm\u00b2 * s)")
plt.legend()
plt.show()
#End of code
