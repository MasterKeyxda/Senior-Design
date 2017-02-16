import 127_3Dflow

#WING SETTINGS

#Wing Name (list, first entry is plot title name, second is file name)
wingname = ['Rectangular Test Wing', 'rect']

#Airfoil name 
    #(polar file must be located in 'Data' folder and named 'name.polar.dat)
    #Must be given as list, allowing for blended airfoils
foil = ['naca1412']
#Number of Spanwise Stations
ny = 101
#Wingspan
b = 1.0
#Aspect Ratio
AR = 20.0
#Wing Area
S = b ** 2 / AR
#Rectangular Wing Chord Distribution
c = b/AR * np.ones(ny)
#Untwisted Chord Distribution
t = np.zeros(ny)

#Global Angle of Attack
alpha = 0
#Freestream Velocity
Vinf = 1.0
#Freestream Density
rhoinf = 1.0
# Reynold's Number
Re = 1.0e5

#Get Data From XFOIL
cmd = ["%s" % foil[0],\
       "oper",\
       "VISC %0.5f" % Re,\
       "ITER %i" % 1000,\
       "PACC",\
       "%s.polar.dat" % foil[0],\
       "",\
       "ASEQ %4.2f %4.2f %4.2f" % (-10.0, 30.0, 0.1),\
       "PACC"]

xfoil(cmd)

#Move Data to "Data" Folder
shutil.move(os.getcwd() + "\\%s.polar.dat" % foil[0], os.getcwd() + "\\Data\\%s.polar.dat" % foil[0])

#Iteration Parameters
maxiter, maxres, damping = 300, 1e-4, 0.05   #Anderson's conditions

#CALL MAIN FUNCTION
testwing = LiftingLineMain(wingname, foil, S, b, c, t, alpha, Vinf, rhoinf,
                            maxiter, maxres, damping)