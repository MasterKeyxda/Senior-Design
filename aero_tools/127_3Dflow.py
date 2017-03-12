import numpy as np
import scipy as sci
from scipy import integrate
from scipy.interpolate import InterpolatedUnivariateSpline
import matplotlib.pyplot as plt
import pickle

import subprocess as sp
import os
import shutil
import sys
import string

#DEFAULT FIGURE SIZING
factor = 1 #scaling factor
plt.rcParams['figure.figsize'] = 6*factor, 6*factor #square plot

#HOW TO DISABLE AUTOSCROLL OF RESULTS
from IPython.display import display, Javascript
disable_autoscroll = """
IPython.OutputArea.prototype._should_scroll = function(lines) {
    return false;
}
"""
display(Javascript(disable_autoscroll))

#MY CUSTOM PLOTTING SETTINGS

#dot, start, x, tri-line, plus
smallmarkers = ['.', '*', 'd', '1', '+']
bigmarkers = ['o', 'v', 'd', 's', '*', 'D', 'p', '>', 'H', '8']
scattermarkers = ['o', 'v', 'd', 's', 'p']

#DEFAULT FONT STYLING
Ttl = 32
Lbl = 32
Box = 28
Leg = 28
Tck = 22
params = {
          'axes.labelsize' : Lbl,
          'axes.titlesize' : Ttl,
#           'text.fontsize'  : Box,
          'font.size'      : Box,
          'legend.fontsize': Leg,
          'xtick.labelsize': Tck,
          'ytick.labelsize': Tck,
          'font.family': 'helvetica'
}
import matplotlib
matplotlib.rcParams.update(params)

#DEFAULT PLOT COLORS
xkcdcolors = ["windows blue", "dusty purple", "leaf green", "macaroni and cheese",  "cherry" , 
              "greyish", "charcoal", "salmon pink", "sandstone",      "tangerine",]
xkcdhex =    ['#3778bf',      '#825f87',      '#5ca904',    '#efb435',              '#cf0234', 
              '#a8a495', "#343837" , "fe7b7c"     , "#c9ae74"  ,      "ff9408"   ,]

import seaborn as sns
#No Background fill, legend font scale, frame on legend
sns.set(style='whitegrid', font_scale=1.5, rc={'legend.frameon': True})
#Mark ticks with border on all four sides (overrides 'whitegrid')
sns.set_style('ticks')
#ticks point in
sns.set_style({"xtick.direction": "in","ytick.direction": "in"})
#Nice blue, purple, green
sns.set_palette(sns.xkcd_palette(xkcdcolors))
#FIX INVISIBLE MARKER BUG
sns.set_context(rc={'lines.markeredgewidth': 0.1})

#PLOTTING PARAMETERS
savetype = '.png'
savedir = 'results/'

def MakeOutputDir(savedir):
    """make results output directory if it does not already exist.
    instring --> directory path from script containing folder
    """
    #split individual directories
    splitstring = savedir.split('/')
    prestring = ''
    for string in splitstring:
        prestring += string + '/'
        try:
            os.mkdir(prestring)
        except Exception:
            pass

MakeOutputDir(savedir)

def MakeFileName(S,b,alpha,name,ext='dat'):
    """ Make filename for saving wing objects with pickle. Savename of form:
    'savedir/Sxxxbxxxaxxx_name.filetype' 
    name --> additional save name text 
    ext --> file extension 
    """ 
    return '{}/S{}b{}a{}_{}.{}'.format(savedir, S, b, alpha, name, ext)
	
class Wing:
    """Stores information and solutions for a given wing"""
    def __init__(self, loadfile=''):
        """Initialize wing class.  Load results from file if filename specified
        """
        if loadfile != '':
            self.LoadData(loadfile)

    def init(self, name, Vinf, rho, a, S, b, c, polars, t=0, yfoil=0):
        """Create class from scratch (dont load data)
        name --> wing name list (1st element is title name, 2nd is filename)
        Vinf --> freestream velocity (for non-dimensionalization)
        rho --> freestream density (for non-dimensionalization)
        a --> geometric alpha [deg] of wing
        b --> wingspan
        c --> chord distribution (# points is # spanwise stations, use odd #)
        polars --> polars associated with airfoils in wing.  List of polar
                    classess if blended wing, or single polar class if
                    constant airfoil
        t --> wing geometric twist distribution (default=0)
        yfoil --> list that specifies spanwise location of different airfoils.
                    Each entry indicates the length of the airfoil distribution
                    (i.e. first airfoil: y/b/2=0 --> yfoil[0],
                    2nd foil y/b/2=0.5 --> yfoil[1])
        """
        #GENERAL WING INFO
        self.name = name
        #freestream velocity
        self.Vinf = Vinf
        #freestream density
        self.rho = rho
        #angle of attack
        self.SetAoA(a)

        #WING GEOMETRY
        #wing span
        self.b = b
        #chord distribution
        self.c = c
        #points in spanwise direction (same and number of chord points)
        self.ny = len(self.c)
        if self.ny%2 == 0:
            print('\n\nWARNING: NUMBER OF SPANWISE STATIONS MUST BE ODD\n\n')
        #root chord index (for odd # of points)
        self.iroot = int( (self.ny - 1) / 2 + 1 - 1 )
        #non-dimensional wingspan vector
        self.y = np.linspace(-self.b/2., self.b/2., self.ny)
        self.dy = self.y[1] - self.y[0]
        #wing area
        self.S = S
        #twist distribution (positive rotates wingtip upward)
        if type(t) == int or type(t) == float:
            self.t = np.ones(self.ny) * t
        else:
            self.t = t

        #Airfoil dististribution (blending)
        if yfoil == 0:
            #only one airfoil in wing, make list of zeros (1st index)
            self.yfoil = np.zeros(self.ny)
        else:
            #multiple airfoils
            #list same length as y, contains index of polar of each airfoil used
            #as each y location

            #Start in center of wing, moving right (mirror later)
            self.yfoil = []
            ind = 0 #index of current airfoil polar
            for y in self.y[self.iroot:]:
                #if current locaiton is new airfoil, change foil index
                if y >= yfoil[ind+1]:
                    ind += 1 #toggle airfoil polar index
                self.yfoil.append( ind )
            #Mirror yfoil for left side of wing
            self.yfoil = np.array(list(self.yfoil[::-1][:-1])
                                + list(self.yfoil))

        #AERODYNAMC PROPERTIES
        #polars of each airfoil used (list)
        if type(polars) != list:
            polars = [polars]
        self.polars = polars

        #Initial root circulation guess (elliptic)
        Gam0 = np.pi * self.c[self.iroot] / (1 + np.pi * self.c[self.iroot]/2)
        #Initial circulation distribution (elliptic)
        self.Gam = Gam0 * np.sqrt(1 - (2*self.y / self.b)**2)

        #INITIALIZE ITERATION PARAMETERS
        #old circulation dist in iteration
        self.Gamold = np.zeros(self.ny)
        #vorticity distribution
        self.dGdy = np.zeros(self.ny)
        #number of iterations in solution
        self.iter = 0
        #induced AoA dist
        self.ai = np.zeros(self.ny)
        #effective AoA dist
        self.aeff = np.zeros(self.ny)
        #downwash dist
        self.w = np.zeros(self.ny)
        #section Cl dist
        self.Cl = np.zeros(self.ny)
        #section lift dist
        self.Lpr = np.zeros(self.ny)

        #3D TERMS
        self.CL = 0
        self.L = 0
        self.CDi = 0
        self.Di = 0
        
        
    def FileName(self, ext='.obj'):
        """Make savename
        name --> additional save name text
        ext --> file extension
        """
        return MakeFileName(int(self.S), self.b, self.a_deg, self.name[1], ext)
        # return '{}/S{}b{}a{}_{}{}'.format(savedir, int(self.S), self.b,
        #                                             self.a_deg, name, ext)

    def SaveData(self):
        """Save wing data to file
        """
        savename = MakeFileName(int(self.S), self.b, self.a_deg, self.name[1], '.obj')
        f = open(savename, 'wb')
        pickle.dump(self.__dict__, f, pickle.HIGHEST_PROTOCOL)
        f.close()

    def LoadData(self, savename):
        """load wing data from file
        """
        print('loading {}'.format(savename))
        f = open(savename, 'rb')
        tmpdic = pickle.load(f)
        f.close()
        self.__dict__.update(tmpdic)

    def SetAoA(self, alpha):
        """Change geometric angle of attack"""
        #geometric AoA [deg]
        self.a_deg = alpha
        #geometric AoA [rad]
        self.a = self.a_deg*np.pi/180

    def InducedAoA(self):
        """Integrate circulation along wingspan to calculate induced AoA
        Used induced AoA to calculate effective AoA and downwash.
        Use dummy variable integration to avoid singularities
        """

        def wInt(i, j):
            """The integral in the downwash equation. i is index of point where
            downwash is being calculated, j is index of trailing vortex whose
            effect we are calculating
            """
            return self.dGdy[j] / (self.y[i] - self.y[j])

        
        #DUMMY VARIABLE INTEGRATION

        #dummy variable of integration (halfway between real wing points)
        xi = [y + self.dy/2 for y in self.y]

        #d/dy(Gamma) (forward diff in y is central diff in xi)
        for i in range(0, self.ny-1):
            self.dGdy[i] = (self.Gam[i+1] - self.Gam[i]) / (self.dy)

        #For each spanwise location, sum effects of each horseshoe vortex
        for i in range(0, self.ny):
            #each horseshoe's effect is an integral
            integ = np.zeros(self.ny)
            for j, (y, g) in enumerate(zip(xi, self.dGdy)):
                integ[j] = g / (self.y[i] - y)
            #integrate effect of each horseshoe at i point
                #exclude endpoint, which is beyond right wingtip
            J = np.trapz(integ[:-1], xi[:-1])

            #induced AoA
            self.ai[i] = 1/(4*np.pi*self.Vinf) * J
            #effective AoA
            self.aeff[i] = self.a - self.ai[i] + self.t[i]
            #downwash dist
            self.w[i] = -self.Vinf * self.ai[i]

    def LiftDist(self):
        """Get section lift distribution from airfoil polar and
        effective AoA distribution.  Get sectional lift force from
        Kutta-Joukowski Lift Theorem"""
        #save old solution
        self.Gamold = np.array(self.Gam)
        #for each effective AoA, find Cl from sectional data
        for i, aeff in enumerate(self.aeff):
            #INTERPOLATE SECTION CL
            self.Cl[i] = self.polars[ int(self.yfoil[i]) ].Interp(aeff)

            #ZERO CIRCULAION AT TIPS (DUMMY INTEGRATION REQUIRES THIS)
            self.Cl[0], self.Cl[-1] = 0, 0

            #SECTION LIFT
            self.Lpr[i] = 0.5 * self.rho * self.Vinf**2 * self.c[i] * self.Cl[i]
            #NEW CIRCULATION DISTRIBUTION (KJ THM)
            self.Gam[i] = 0.5 * self.Vinf * self.c[i] * self.Cl[i]

    def IterateWing(self, maxiter, maxres, damp, showfig=False):
        """Iterate wing circulation distribution until the residual is within
        acceptable bounds.  Use Successive Under Relaxation (SUR) iterative
        method.
        maxiter --> maximum number of iterations
        maxres --> maximum size of residual that will be accepted
        damp --> iteration damping factor
        showfig --> plot each iterative result for circulation on same plot
        """

        #RESET CIRCULATION GUESS
        #Initial root circulation guess (elliptic)
        Gam0 = np.pi * self.c[self.iroot] / (1 + np.pi * self.c[self.iroot]/2)
        #Initial circulation distribution (elliptic)
        self.Gam = Gam0 * np.sqrt(1 - (2*self.y / self.b)**2)
        #self.Gam = 1 - np.linspace(0, 1, self.ny) #linear guess


        print('Beginning Iteration...')


        if showfig:
            plt.figure()

        iter = 1
        res = maxres+1
        j  = 0
        while iter < maxiter and j < 5:
            #CALC NEW CIRCULATION DITRIBUTION
            self.InducedAoA()
            self.LiftDist()
            #CALC RESIDUAL (difference between current and prev solution)
            R = self.Gam - self.Gamold

            #mid = 20
            #R = np.array(list(self.Gam[:self.iroot-mid] - self.Gamold[:self.iroot-mid])
            #              + list(self.Gam[self.iroot+mid:] - self.Gamold[self.iroot+mid:]))

            res = max(abs(R))


            #stop iteration when residual has been less than maxres for 5 iters
            if res < maxres:
                j += 1
                #j = 5 #skip the 5 steps after residual

            #SAVE NEW CIRCULATION DISTRIBUTION (use SUR)
            self.Gam = self.Gamold + damp * (self.Gam - self.Gamold)

            #increment iteration count
            iter += 1

            if showfig:
                xx = self.y/(self.b/2)
                yy = self.Gam/max(self.Gam) #normalized
                yy = self.Gam               #non-normalized
                #plt.plot(self.y/(self.b/2), self.Gam/max(self.Gam))#normalized
                if iter == 2:
                    plt.plot(xx, yy, color='black', linewidth=3)
                else:
                    plt.plot(xx, yy)



        self.iter = iter - 1 #save final iteration count
        ires = np.where( R == max(abs(R)) )
        self.maxresloc = self.y[ires]
        print('Iteration Complete! (iters={}, res={}, maxresloc={})'.format(
                    self.iter,res, self.maxresloc/(self.b/2)) )

        if showfig:
            plt.title('Circulation iter=%s'%(self.iter), fontdict=font_ttl)
            plt.xlabel(r'$\frac{y}{b/2}$', fontdict=font_lbl)
            plt.ylabel(r'$\frac{\Gamma}{\Gamma_0}$', fontdict=font_lbl)
            plt.axis([0, 1, 0, 1])
            plt.show()

    def GetLift(self):
        """Calculate 3D lift of wing from ciruclation distribution.
        Only integrate inboard points, exclude wingtips
        """
        self.CL = (2 / (self.Vinf * self.S)
                    * np.trapz(self.Gam, self.y))
        self.L = self.CL * 0.5 * self.rho * self.Vinf ** 2

    def GetDrag(self):
        """Calculate induced drag of wing from ciruclation distribution
        """
        self.CDi = (2 / (self.Vinf * self.S)
                   * np.trapz( (self.Gam * self.ai), self.y))
        self.Di = self.CDi * 0.5 * self.rho * self.Vinf ** 2

    def GetLiftDist(self, maxiter=300, maxres=1e-4, damping=0.05, showfig=False):
        """Calculate lift/circulation distribution of given wing and given
        angle of attack
        """

        #Converge wing for current AoA
        self.IterateWing(maxiter, maxres, damping, showfig=showfig)

        #CALCULATE AERODYANAMIC FORCES
        self.GetLift()
        self.GetDrag()

    def AlphaSweep(self, alphas=np.linspace(0,20,10), overwrite=True,
                    maxiter=300, maxres=1e-4, damping=0.05):
        """Solve flow over a wing at various angles of attack.
        If an angle of attack for this wing has already been solved and saved
        as a pickled object, load this data rather than calculating it.
        self --> wing object
        alphas --> list of geometric angles of attack to run
        overwrite --> overwrite existing wing objects
        """
        #STORE 3D POLAR VARIABLES AS LISTS IN DICTIONARY
        polar = {}
        keys = ['alpha', 'Cl', 'Cdi']
        for key in keys:
            polar[key] = []

        #SOLVE WING FOR EACH ALPHA
        for alpha in alphas:
            #Set current AoA to find coefficients for
            self.SetAoA(alpha)
            #CHECK IF ALREADY SOLVED FOR
            wingname = self.FileName() #current object file name
            if not os.path.isfile(wingname) or overwrite:
                #SOLVE WING FOR GIVEN AOA IF NO SAVED FILE
                #Reiterate lift distribution
                self.GetLiftDist(maxiter=maxiter, maxres=maxres, damping=damping)
            else:
                #LOAD WING DATA FROM FILE
                self.LoadData(wingname)

            #SAVE POLAR DATA
            polar['alpha'].append( alpha    )
            polar['Cl'].append(    self.CL  )
            polar['Cdi'].append(   self.CDi )

        #SAVE DATA
        savename = '{}/{}_polar.dat'.format(savedir, self.name[1])
        np.savetxt(savename, np.c_[polar['alpha'],
                                   polar['Cl'],
                                   polar['Cdi']] )

        return polar

    def PlotWing(self, color='black'):
        """Plot wing planform geometry.
        Assumes chord distribution is symmetric about x=0.
        color --> line color
        """
        #plot leading edge
        handle, = plt.plot(self.y, self.c/2.0, color=color)
        #plot trailing edge
        plt.plot(self.y, -self.c/2.0, color=color)
        #plot left wing end
        plt.plot([self.y[0], self.y[0]], [-self.c[0]/2, self.c[0]/2],
                                                            color=color)
        #plt right wing end
        plt.plot([self.y[-1], self.y[-1]], [-self.c[-1]/2, self.c[-1]/2],
                                                            color=color)
        return handle

    def PlotDistProps(self, showplot=False):
        """Plot wingspan distribution properties (induced alpha, circulation,
            downwash, etc.
        """
        #Plot Wing Dist Properties
        f, axarr = plt.subplots(5, sharex=True, figsize=[5, 10])
        #PLOT CIRCULATION DISTRIBUTION
        #axarr[0].set_title(r'$\Gamma$')
        axarr[0].plot(self.y,self.Gam, label=r'$\Gamma$', marker='.',)
        axarr[0].set_ylabel(r'$\Gamma$')
        #PLOT VORTICITY DIST
        #axarr[1].set_title(r'$d\Gamma/dy$', fontdict=font_ttl)
        axarr[1].plot(self.y,self.dGdy, label=r'$\frac{d\Gamma}{dy}$'
                            , marker='.',)
        axarr[1].set_ylabel(r'$d\Gamma/dy$')
        #PLOT DOWNWASH
        #axarr[2].set_title('Downwash', fontdict=font_ttl)
        axarr[2].plot(self.y,self.w, label = 'w', marker='.',)
        axarr[2].set_ylabel('Downwash')
        #PLOT INDUCED ANGLE OF ATTACK
        #axarr[3].set_title(r'$\alpha_i$', fontdict=font_ttl)
        axarr[3].plot(self.y,self.ai*180/np.pi, label=r'$\alpha_i$',
                            marker='.',)
        axarr[3].set_ylabel(r'$\alpha_i$')
        #PLOT EFFECTIVE ANGLE OF ATTACK
        #axarr[4].set_title(r'$\alpha_{eff}$', fontdict=font_ttl)
        axarr[4].plot(self.y,self.aeff*180/np.pi, label=r'$\alpha_{eff}$',
                            marker='.',)
        axarr[4].set_ylabel(r'$\alpha_{eff}$')
        #axarr[4].legend(bbox_to_anchor=(1.01, 0.5), loc='center left',
                    #fontsize='large',
                    ##fancybox=True, frameon=True,
                    #framealpha=0.75,
                    #numpoints=1, scatterpoints=1,
                    #borderpad=0.1, borderaxespad=0.1, handletextpad=0.2,
                    #handlelength=1.0, labelspacing=0,
                    #)
        axarr[4].set_xlabel(r'$\frac{y}{b/2}$')
		
class Polar:
    """Stores polar data for airfoils"""

    def __init__(self, name, alpha=0, Cl=0):
        """Can manually create polar data by providing alpha and Cl arrays.
        Otherwise, enter default and load XFOIL data with later function
        name --> polar name
        alpha --> list or single of AoA in degrees (gets converted to radians)
        Cl    --> list or single of lift coefficient
        """
        self.name = name
        self.alpha_deg = alpha
        self.alpha = alpha*np.pi/180
        self.Cl = Cl

    def ReadXfoilPolar(self, ifile, Cd=0):
        """Reads XFOIL polar output. Alpha and Cl ready automatically.
        Other info is optional and read if specified in function input
        ifile --> path of input file (string)
        Cd    --> read drag coefficient if for Cd=1
        """
        data = np.loadtxt(ifile, skiprows=12, unpack=False)
        self.alpha_deg = data[:,0]
        self.alpha = self.alpha_deg * np.pi / 180
        self.Cl = data[:,1]
        if Cd != 0:
            self.Cd = data[:,2]

    def Interp(self, alpha):
        """Interpolate Cl for given alpha
        """
        return np.interp(alpha, self.alpha, self.Cl)

    def LiftCurve(self):
        """Plot lift curve of polar"""
        plt.figure()
        plt.title(self.name[0] + ' Lift Curve', fontdict=font_ttl)
        plt.xlabel(r'$\alpha [deg]$', fontdict=font_lbl)
        plt.ylabel(r'$C_l$', fontdict=font_lbl)
        plt.plot(self.alpha*180/np.pi, self.Cl)
        plt.show()
		
def LiftingLineMain(name, airfoils, S, b, c, t=0, alpha=0, Vinf=1, rhoinf=1,
                        maxiter=300, maxres=1e-4, damping=0.05):
    """Simulate flow over a finite wing of given parameters.
    name --> name of wing (used for file name and plotting label purposes)
    airfoils --> list of airfoil section file names to use
    S --> wing area
    b --> wing span
    c --> spanwise chord distribution, give number of spanwise stations
    t --> twist distribution
    alpha --> angle of attack to compute lifting line calcs at
    Vinf --> used in non-dimensionalizing coefficients
    rhoinf --> used in non-dimensionalizing coefficients
    maxiter --> iteration limit for lifting line convergence
    maxres --> maximum allowable residual for convergence
    damping -- SUR damping ratio
    """

    #READ POLAR DATA FOR WING
    polars = []
    for i, foil in enumerate(airfoils):
        #initialize each polar with its name
        polars.append( Polar(foil) )
        #read xfoil data
        polars[-1].ReadXfoilPolar( 'Data/{}.polar.dat'.format(foil), Cd=1)

    #INITIALIZE WING OBJECT
    wing = Wing()
    wing.init(name, Vinf, rhoinf, alpha, S, b, c, polars, t=t, yfoil=0)

    #SOLVE LIFT DISTRIBUTION OVER WING
    wing.GetLiftDist(maxiter, maxres, damping, showfig=False)

    return wing
    
def xfoil(cmd):
    # Open XFOIL on command line
    ps = sp.Popen(['xfoil.exe'],
                  stdin=sp.PIPE,
                  stdout=None,
                  stderr=None)

    # Pass commands into XFOIL command line
    res = ps.communicate( bytes('\n'.join(cmd), 'utf=8' ))
    print (res)