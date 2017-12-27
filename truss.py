# -*- coding: utf-8 -*-
"""
Created on Sat Dec 30 18:49:40 2018

3D truss model updater program created by Máté Szedlák (2016-2018).
Copyright MIT, Máté Szedlák 2016-2018.
"""
### Code for "Processing" only
#size(640, 360)
#background(126)

import math
import itertools
from copy import deepcopy
from extra_math import mat_vec_mult as mat_vec_mult
from extra_math import invert as invert
#from extra_math import check_for_all_zeros as check_for_all_zeros
#from extra_math import swap_row as swap_row
from extra_math import swap_col as swap_col
#from extra_math import make_identity as make_identity
from config import Configuration



# COMPATIBILITY MODES:
    # 0: User defined
    # 1: Processing3
    # 2: Android
    # 3: Most information (with numpy)
    # 4: Maximum compatibility

_COMPATIBLE_MODE = 0
_SIMULATION = 1                     # Simulating measurements based on input file

Conf = Configuration(_COMPATIBLE_MODE, _SIMULATION)


PORTS = ['COM1', 'COM2', 'COM3']      # List of possible communication ports
PORTNUMBER = 0                      # Applied communication port

if _COMPATIBLE_MODE == 0*0:
    ### User defined ###
    # Modify as needed #
    Conf.mode_name = "User defined"
    Conf.log = 1                 # Logging time
    Conf.graphics = 1            # Graphical features
    Conf.solver = 1              # 0: Basic solver, 1: NumPy solver
    Conf.OSlib = 1  # Basic OS file features (e.g. file size)
    Conf.updating = 0            # Model Updating: On/ Off
    Conf.arduino = 0             # Arduino input: On/Off
    Conf.debug = 0               # Debugging mode
    Conf.realistic_simulation = 0 # Wait as long as it was originally. Only valid with _SIMULATION = 1

if Conf.OSlib:
    import os
    try:
        if not os.path.exists("./Results/"):
            os.makedirs("./Results")
    except FileNotFoundError:
        print("Error: Please manually create the 'Results' folder" 
              "in the root of the project")
        

if Conf.simulation or not Conf.updating:
    _ARDUINO = 0

if _COMPATIBLE_MODE == 2:
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

if Conf.log:
    import time
    import datetime
    TIC = time.time()
    print('------------------------------------')
    print('Truss calculational program')
    print('Created by Máté Szedlák (23/11/2016)')
    print('Compatibility mode: ' + Conf.mode_name)
    if Conf.solver == 0:
        print('- Solver is set to default')
    elif Conf.solver == 1:
        print('- Solver is set to NumPy')
    else:
        raise Exception("Solver settings are invalid!")
    if Conf.updating:
        print('+ Model updating is turned ON')
        if Conf.simulation:
            print('Input data is SIMULATED!')
    else:
        print('- Model updating is turned OFF')
    print('------------------------------------')


if Conf.arduino:
    SER = 0
    try:
        import serial
    except ImportError:
        print("You tried to import \'serial\' in Windows mode without installing \'pySerial\'.")
        print("Please first install pySerial: http://playground.arduino.cc/Interfacing/Python")
        raise Exception('Android mode denied: pyserial not found')

    PORTNUMBER -= 1
    while SER == 0:
        PORTNUMBER += 1
        if PORTNUMBER >= len(PORTS):
            PORTNUMBER = 0
        time.sleep(0.6)
        print('Opening serial at port ' + str(PORTS[PORTNUMBER]))
        try:
            SER.close()
        except Exception:
            pass
        try:
            SER = serial.Serial(PORTS[PORTNUMBER], 9600, timeout=0)
        except serial.SerialException:
            Exception(PORTS[PORTNUMBER] + ' port is busy. It might be occupied by this program or another one :/ Be careful or try resetting this program')
            SER = 0
            try:
                SER.close()
            except Exception:
                pass
        except Exception:
            SER = 0
            try:
                SER.close()
            except Exception:
                pass


#if _ARDUINO or _SIMULATION:
#    try:
#        mappingfile = 'arduino_mapping.txt'
#        with open(mappingfile, "r") as textfile:
#            line = textfile.readline().strip()
#            arduino_mapping = line.upper().split(',')
#
#    except IOError:
#        raise Exception('File not found: ' + mappingfile)


if Conf.solver:
    # NumPy library for solving linear equations in another way
    import numpy as np

if Conf.graphics:
    # library for drawing
    from graphics import Arrow3D


# XXX DEPRECATED !
def endoffile(givenfile, line):
    """
    Check if end of file is reached. Implemented due to compatibility reasons.
    """
    if Conf.OSlib:
        return givenfile.tell() < os.fstat(givenfile.fileno()).st_size
    else:
        return not line == "EOF"

def logtime(prev_time, title):
    """
    Calculating and printing the time consumption of tasks

    Should be called with the previously saved part-time and the name of the actual task
    At the first call, should be called with TIC value. The input argument
    should be overwritten by this funvtion's return value.
    """
    if Conf.log:
        new_time = time.time()
        print(title)
        print('Time: ' + str("{:10.3f}".format(new_time - prev_time)))
        print('------------------------------------')
        return new_time
    else:
        return 0

def error(delta):
    """
    Error function using least-square method
    """
    sumerr = 0
    for deltaelement in delta:
        sumerr += deltaelement**2
    sumerr = math.sqrt(sumerr)

    return sumerr

class Truss(object):
    """
    General structure class
    """
    def __init__(self, name):
        self.name = name              # Name of structure
        self.known_f_a = []           # Nodes without supports
        self.known_f_notzero = []     # Nodes with loads
        self.dof = 3                  # Truss's degree of freedom
        self.node = []                # Element's end nodes
        self.constraint = []          # Supports
        self.force = []               # Force
        self.nodalcoord = []          # Coordinates of nodes
        self.nodalcoord_def = []      # Coordinates after deformations
        self.area = []                # Cross-sectional areas
        self.el_mod = []                   # Material data
        self.nodenum = 0              # Number of nodes
        self.elenum = 0               # Number of elements
        self.eledof = []              # Mapping between DOF and node
        self.stiffness = []           # Global stiffness matrix
        self.mod_stiffnesses = []     # Modified stiffnesses in a hyper-matrix
        self.ele_length = []                   # Length of the elements
        self._norm_stiff = []                  # E/L
        self._cx = []
        self._cy = []
        self._cz = []
        self._s_loc = []
        self._loc_stiff = []          # Local stiffnes matrix
        self.dis_new = []
        self.force_new = []
        self.stiff_new = []
        self.displacement = []        # Relative displacements
        self._stiffisfresh = 0
        self._postprocessed = 0
        self.init_disp = []
        self.stresscolor = []         # Color mapping for stresses
        self.known_dis_a = []
        self.stress = []              # Element's stresses
        self._io_origin = 0           # Array's first element number during IO. Default is 0.
        self.analysis = {}
        self._mod_stiffisfresh = 0
        self.mod_displacements = []
        self.keypoint = []
        self.keypnum = 0
        self.effect = []
        self.toteffect = []
        self.sortedeff = []
        self.specdof_inputstring = ''
        self.tresshold = 0.1
        self.effectratio = []
        self.processeddata = []          # To store last input line
        self.modifications = []          # Storing modifications for model updating
        self.readelements = [0]*9
        self.arduino_mapping = []
        self.errorlimit = 0.5
        self.modificationlimit = 0.6
        self.unitmodification = 0.05
        self.measurement = [0.]
        self.numofupdates = [0, 0, 0]        # [#Successfully updated model, #Updates with overflow exit, #Updates where there were no more modification option]
        self.iterationlimit = 20

    def read(self, filename):
        """
        Input file for TRUSS.py program
        All commands must be written with uppercase characters
        *** The values MUST be written in the exact following line of a command
        Only lines with the command and nothing more counts.
        Everything else will be neglected. Even hastags are useless :)
        The order of the commands are indifferent.

        Commands and their format (example):
            DOF - Degree of freedom: 3
            ELEMENTS - Elements given by end-nodes: 0, 1 | 0, 2 ...
            COORDINATES - Nodal coordinates: 0., 0., 0., | 0., 3., 0. ...
            CROSS-SECTIONS - This data will be evaluated in Python: 3.0*(10**(-4)), 5.0*(10**(-4)) ...
            MATERIALS - This data will be evaluated in Python: 70.0*(10**9), 100.0*(10**9) ...
            FORCES - Selected DOF + Force: 11, +1000000.0 | 12, +1000000.0 ...
            SUPPORTS - Selected DOF + Prescribed displacement: 0, 0.0 | 1, 0.0 ...
            SPECDOF - Selected node's DOF will be analysed during Model Updating: 1, xyz | 3 y | 10 xz ...

            EOF - For compatibility reasons EOF should be placed after the commands
        """
        self._io_origin = 0
        readelementnames = ["Origin", "DOF", "Elements", "Coordinates", "Cross-sections", "Materials", "Forces", "Supports", "Measured DOFs"]

        with open(filename, "r") as sourcefile:
            sourceline = ""
            while sourceline != "EOF":
                sourceline = sourcefile.readline().strip()

                if sourceline.upper() == "_ORIGIN":
                    sourceline = sourcefile.readline().strip()
                    self._io_origin = int(sourceline)
                    self.readelements[0] = 1

                if sourceline.upper() == "DOF":
                    sourceline = sourcefile.readline().strip()
                    self.setdof(int(sourceline))
                    self.readelements[1] = 1

                if sourceline.upper() == "ELEMENTS":
                    sourceline = sourcefile.readline().strip()
                    inpstr = []
                    inpnum = []
                    inpstr = [x.split(',') for x in sourceline.split('|')]
                    if len(inpstr[0]) == 1:
                        inpstr = [x.split(';') for x in sourceline.split('|')]
                    if [''] in inpstr:
                        inpstr.remove([''])
                    inpnum = [[int(x[0]) - self._io_origin, int(x[1]) - self._io_origin] for x in inpstr]
                    self.setelements(inpnum)
                    self.readelements[2] = 1

                if sourceline.upper() == "COORDINATES":
                    sourceline = sourcefile.readline().strip()
                    inpstr = []
                    inpnum = []
                    inpstr = [x.split(',') for x in sourceline.split('|')]
                    if len(inpstr[0]) == 1:
                        inpstr = [x.split(';') for x in sourceline.split('|')]
                    if [''] in inpstr:
                        inpstr.remove([''])
                    if self.dof == 3:
                        inpnum = [[float(x[0]), float(x[1]), float(x[2])] for x in inpstr]
                    elif self.dof == 2:
                        inpnum = [[float(x[0]), float(x[1]), 0.] for x in inpstr]
                    self.setcoordinates(inpnum)
                    self.readelements[3] = 1

                if sourceline.upper() == "CROSS-SECTIONS":
                    sourceline = sourcefile.readline().strip()
                    inpstr = []
                    inpnum = []
                    inpstr = sourceline.split(',')
                    if len(inpstr) == 1:
                        inpstr = sourceline.split(';')
                    if '' in inpstr:
                        inpstr.remove('')
                    inpnum = [float(eval(x)) for x in inpstr]
                    self.setcrosssections(inpnum)
                    self.readelements[4] = 1

                if sourceline.upper() == "MATERIALS":
                    sourceline = sourcefile.readline().strip()
                    inpstr = []
                    inpnum = []
                    inpstr = sourceline.split(',')
                    if len(inpstr) == 1:
                        inpstr = sourceline.split(';')
                    if '' in inpstr:
                        inpstr.remove('')
                    inpnum = [float(eval(x)) for x in inpstr]
                    self.setmaterials(inpnum)
                    self.readelements[5] = 1

                if sourceline.upper() == "FORCES":
                    sourceline = sourcefile.readline().strip()
                    inpstr = []
                    inpnum = []
                    inpstr = [x.split(',') for x in sourceline.split('|')]
                    if len(inpstr[0]) == 1:
                        inpstr = [x.split(';') for x in sourceline.split('|')]
                    if [''] in inpstr:
                        inpstr.remove([''])
                    inpnum = [[int(x[0]) - self._io_origin, float(x[1])] for x in inpstr]
                    self.setforces(sorted(inpnum))
                    self.readelements[6] = 1

                if sourceline.upper() == "SUPPORTS":
                    sourceline = sourcefile.readline().strip()
                    inpstr = []
                    inpnum = []
                    inpstr = [x.split(',') for x in sourceline.split('|')]
                    if len(inpstr[0]) == 1:
                        inpstr = [x.split(';') for x in sourceline.split('|')]
                    if [''] in inpstr:
                        inpstr.remove([''])
                    inpnum = [[int(x[0]) - self._io_origin, float(x[1])] for x in inpstr]
                    self.setsupports(sorted(inpnum))
                    self.readelements[7] = 1

                if sourceline.upper() == "MEASUREMENTS":
                    sourceline = sourcefile.readline().strip()
                    self.specdof_inputstring = sourceline
                    inpstr = []
                    self.arduino_mapping = sourceline.split(',')
                    self.setspecdofs(self.arduino_mapping)
                    self.readelements[8] = 1

        terminate = False
        for i, value in enumerate(self.readelements):
            if i > 0 and (i < 8 or Conf.updating):
            #if i > 0:
                if value == 0:
                    print("The following was not found: " + readelementnames[i])
                    terminate = True
        if terminate:
            raise Exception

    def plot(self, showorig, showresult, showsupports, showforces, \
             showreactions, scaledisplacement, scaleforce, scalez, saveplot):
        """
        Plot function of the Truss class
        This method calls the more general plotstructure() method.

            Plot settings:
                O: Original D: Deformed S: Supports F: Forces R: Reactions
                ScD: Scale displacments (Z-axis) (def:1.0) ScF: Scale forces (def:1.0)
                ScS: Scale Support signs (Z-axis) (def:1.0)
                Save: Save plot to file

                plot(O, D, S, F, R, ScD, ScF, ScS, Save)
        """
        _showvalues = 1     # Show values of forces

        if self._postprocessed == 0:
            print('Postprocess is needed before plotting structure!')
        else:
            if scaledisplacement == 0:
                scaledisplacement = 1.0           # Scale drwaing of displacements
            if scaleforce == 0:
                scaleforce = 1.0                  # Scale force sign
            if scalez == 0:
                scalez = 0.3                      # Scale z-axis

            Arrow3D.plotstructure(self, showorig, showresult, showsupports, showforces, showreactions, \
                 scaledisplacement, scaleforce, scalez, _showvalues, saveplot, Conf.log)

    def __checkcoordinates(self, ignorable):
        """
        Checking coordinates for repeating elements.

            ignorable: [True | False] If the warning is ignorable, only message apperas and the input becomes neglected.
                If the warning is not ignorable, then exceptions will be raised.
            return: [0 | 1] 1 f no error found, otherwise 0.
        """
        if len(self.nodalcoord) != len(list(k for k, _ in itertools.groupby(sorted(self.nodalcoord)))):
            if ignorable == 0:
                raise Exception('Coordinate list has repeating items. Calculation is terminated')
            else:
                print("This node already exists. Input is ignored.")
            return 0
        else:
            return 1

    def setdof(self, dof):
        """
        Setting problem's degree of freedom

            dof: [2 | 3] Model's Degree Of Freedom.
        """
        self.dof = dof
        if self.dof != 2 and self.dof != 3:
            raise Exception('DOF must be 2 or 3.')
        self._stiffisfresh = 0
        self._mod_stiffisfresh = 0
        self._postprocessed = 0

    def setelements(self, node):
        """
        Setting elements (nodal connections) in bulk mode
        """
        self.node = node
        self.nodenum = len(set(list(itertools.chain.from_iterable(sorted(self.node)))))
        self.elenum = len(self.node)
        self._stiffisfresh = 0
        self._mod_stiffisfresh = 0
        self._postprocessed = 0

        # Creating mapping tool for elements
        for node in self.node:
            self.eledof.append([node[0]*3, node[0]*3+1, node[0]*3+2, \
                                node[1]*3, node[1]*3+1, node[1]*3+2])

        # Initialazing matrix for all matrices
        self.init_disp = [0.]*(3*self.nodenum)
        self.force = [0.]*(3*self.nodenum)
        self.stiffness = [0.]*(3*self.nodenum)
        self.known_f_a = []
        self.known_f_notzero = []

    def setcoordinates(self, coordinates):
        """
        Setting coordinates in bulk mode
        """
        self.nodalcoord = coordinates
        self._stiffisfresh = 0
        self._mod_stiffisfresh = 0
        if self.nodenum > len(self.nodalcoord):
            raise Exception('More coordinates are needed')
        elif self.node == []:
            raise Exception('Nodes must be set before defining elements')
        self.__checkcoordinates(False)

    def modcoordinate(self, node, coordinate):
        """
        Modify coordinate
        """
        if self.__checkcoordinates(True):
            self.nodalcoord[node] = coordinate
        self._stiffisfresh = 0
        self._mod_stiffisfresh = 0

    def setcrosssections(self, area):
        """
        Setting cross-sections in bulk mode
        """
        self.area = area
        self._stiffisfresh = 0
        self._mod_stiffisfresh = 0
        self._postprocessed = 0

    def modcrosssection(self, element, area):
        """
        Modifying cross-sections by elements
        """
        self.area[element] = area
        self._stiffisfresh = 0
        self._mod_stiffisfresh = 0
        self._postprocessed = 0

    def setmaterials(self, el_mod):
        """
        Setting material data in bulk mode
        """
        self.el_mod = el_mod
        self._stiffisfresh = 0
        self._mod_stiffisfresh = 0
        self._postprocessed = 0

    def modmaterial(self, element, el_mod):
        """
        Modifying material data by elements
        """
        self.el_mod[element] = el_mod
        self._stiffisfresh = 0
        self._mod_stiffisfresh = 0
        self._postprocessed = 0

    def setforces(self, forces):
        """
        Set forces
        """
        for fdof, force in forces:
            if self.dof == 3:
                self.force[fdof] = force
            elif self.dof == 2:
                self.force[fdof + (fdof//2)] = force
        self._postprocessed = 0

    def modforce(self, element, force):
        """
        Modifying forces by each
        """
        self.force[element] = force
        self._postprocessed = 0

    def setsupports(self, constraints):
        """
        Set supports
        """
        for cdof, constraint in constraints:
            if self.dof == 3:
                self.constraint.append([cdof, constraint])
            elif self.dof == 2:
                self.constraint.append([cdof + (cdof // 2), constraint])
        self._stiffisfresh = 0
        self._mod_stiffisfresh = 0
        self._postprocessed = 0

    def setspecdofs(self, specdofs):
        """
        Set special nodal DOFs
        """
        self.analysis = {}

        for dofname in specdofs:
            node = int(dofname[:len(dofname)-1])-self._io_origin
            if 'X' in dofname:
                self.analysis[dofname] = node*3+0
                self.keypoint.append(node*3+0)
            if 'Y' in dofname:
                self.analysis[dofname] = node*3+1
                self.keypoint.append(node*3+1)

            if 'Z' in dofname:
                if self.dof == 3:
                    self.analysis[dofname] = node*3+2
                    self.keypoint.append(node*3+2)
                else:
                    print("Z-direction is not allowed in 2D structures. Please check the \'MEASUREMENTS\' section in the input file.")
                    raise Exception

        self.keypnum = len(self.analysis)
        if self.keypnum == 0 and Conf.updating:
            print("There is no valid measured DOF. Please check the \'MEASUREMENTS\' section in the input file.")
            raise Exception


    def calcstiffness(self):
        """
        Stiffness matrix compilation
        """
        self._postprocessed = 0

        if self.dof == 2:
            for zdof in range(self.nodenum):
                self.constraint.append([int(zdof*3+2), 0.])
        self.constraint = list(k for k, _ in itertools.groupby(sorted(self.constraint)))

        #Setting known forces
        for dofloc in range(3*self.nodenum):
            self.known_f_a.append(dofloc)
            if self.force[dofloc] != 0:
                self.known_f_notzero.append(dofloc)

        self.known_dis_a = []
        for constr in self.constraint:
            self.init_disp[constr[0]] = constr[1]
            self.known_dis_a.append(constr[0])
            try:
                self.known_f_a.remove(constr[0])
                self.known_f_notzero.remove(constr[0])
            except ValueError:
                pass

        ele_length = [0.]*self.elenum
        self._norm_stiff = [0.]*self.elenum
        self._cx = [0.]*self.elenum
        self._cy = [0.]*self.elenum
        self._cz = [0.]*self.elenum
        self._s_loc = [0.]*self.elenum
        self._loc_stiff = [0.]*self.elenum
        self.stress = [0.]*self.elenum

        self.stiffness = [[0.]*(len(self.nodalcoord)*3)]*(len(self.nodalcoord)*3)

        for i in range(self.elenum):
            ele_length[i] = math.sqrt((self.nodalcoord[self.node[i][1]][0]-self.nodalcoord[self.node[i][0]][0])**2+ \
                (self.nodalcoord[self.node[i][1]][1]-self.nodalcoord[self.node[i][0]][1])**2 + \
                (self.nodalcoord[self.node[i][1]][2]-self.nodalcoord[self.node[i][0]][2])**2)

            self._cx[i] = (self.nodalcoord[self.node[i][1]][0]-self.nodalcoord[self.node[i][0]][0])/ele_length[i]
            self._cy[i] = (self.nodalcoord[self.node[i][1]][1]-self.nodalcoord[self.node[i][0]][1])/ele_length[i]
            self._cz[i] = (self.nodalcoord[self.node[i][1]][2]-self.nodalcoord[self.node[i][0]][2])/ele_length[i]
            self._norm_stiff[i] = self.el_mod[i]/ele_length[i]
            self._s_loc[i] = [[self._cx[i]**2, self._cx[i]*self._cy[i], self._cx[i]*self._cz[i], -self._cx[i]**2, -self._cx[i]*self._cy[i], -self._cx[i]*self._cz[i]], \
                [self._cx[i]*self._cy[i], self._cy[i]**2, self._cy[i]*self._cz[i], -self._cx[i]*self._cy[i], -self._cy[i]**2, -self._cy[i]*self._cz[i]], \
                [self._cx[i]*self._cz[i], self._cy[i]*self._cz[i], self._cz[i]**2, -self._cx[i]*self._cz[i], -self._cy[i]*self._cz[i], -self._cz[i]**2], \
                [-self._cx[i]**2, -self._cx[i]*self._cy[i], -self._cx[i]*self._cz[i], self._cx[i]**2, self._cx[i]*self._cy[i], self._cx[i]*self._cz[i]], \
                [-self._cx[i]*self._cy[i], -self._cy[i]**2, -self._cy[i]*self._cz[i], self._cx[i]*self._cy[i], self._cy[i]**2, self._cy[i]*self._cz[i]], \
                [-self._cx[i]*self._cz[i], -self._cy[i]*self._cz[i], -self._cz[i]**2, self._cx[i]*self._cz[i], self._cy[i]*self._cz[i], self._cz[i]**2]]
            self._loc_stiff[i] = [[y* self.area[i]* self._norm_stiff[i] for y in x] for x in self._s_loc[i]]
            ele_dof_vec = self.eledof[i]

            stiffincrement = [0.]*(len(self.nodalcoord)*3)

            for j in range(3*2):
                for k in range(3*2):
                    stiffincrement[ele_dof_vec[k]] = self._loc_stiff[i][j][k]
                self.stiffness[ele_dof_vec[j]] = [x + y for x, y in zip(self.stiffness[ele_dof_vec[j]], stiffincrement)]
        self._stiffisfresh = 1

    def calcmodstiffness(self, index, magnitude):
        """
        Convergency step in stiffness matrix modification
        """
        if self.mod_stiffnesses == []:
            self.mod_stiffnesses = [0.]*(self.elenum+1)

        #for loopindex in range(self.elenum):
        _mod_stiffnesses_temp = [[0.]*(len(self.nodalcoord)*3)]*(len(self.nodalcoord)*3)

        for i in range(self.elenum):
            if i == index:
                _mod_norm_stiff = self._norm_stiff[i] * (1.0 + self.modifications[i] + magnitude) #self.el_mod[i]/ele_length[i]
            else:
                _mod_norm_stiff = self._norm_stiff[i] * (1.0 + self.modifications[i]) #self.el_mod[i]/ele_length[i]

            _mod_loc_stiff = [[y*self.area[i]*_mod_norm_stiff for y in x] for x in self._s_loc[i]]

            ele_dof_vec = self.eledof[i]

            stiffincrement = [0.]*(len(self.nodalcoord)*3)

            for j in range(3*2):
                for k in range(3*2):
                    stiffincrement[ele_dof_vec[k]] = _mod_loc_stiff[j][k]
                _mod_stiffnesses_temp[ele_dof_vec[j]] = [x + y for x, y in zip(_mod_stiffnesses_temp[ele_dof_vec[j]], stiffincrement)]

        self.mod_stiffnesses[index] = _mod_stiffnesses_temp


    def solve(self):
        """
        Main solver of the code
        """
        if self._stiffisfresh == 0:
            if Conf.log:
                print('Stiffness matrix is recalculated')
            self.calcstiffness()

        self.dis_new = [0.]*(self.nodenum*3-len(self.constraint))
        self.force_new = [0.]*(self.nodenum*3-len(self.constraint))
        self.stiff_new = [[0.]*(self.nodenum*3-len(self.constraint))]*(self.nodenum*3-len(self.constraint))

        # known force array
        for i, known_f_a in enumerate(self.known_f_a):
            self.force_new[i] = self.force[known_f_a]

        stiffincrement = [0.]*(self.nodenum*3-len(self.constraint))
        for i, kfai in enumerate(self.known_f_a):
            for j, kfaj in enumerate(self.known_f_a):
                stiffincrement[j] = self.stiffness[kfai][kfaj]
            self.stiff_new[i] = [x + y for x, y in zip(self.stiff_new[i], stiffincrement)]

        # SOLVING THE STRUCTURE
        if Conf.solver == 0:
            if Conf.log:
                print('Built-in solver')
            self.dis_new = mat_vec_mult(invert(self.stiff_new), self.force_new)
        else:
            if Conf.log:
                print('NumPy solver')
            self.dis_new = np.linalg.solve(np.array(self.stiff_new), np.array(self.force_new))

        self.displacement = deepcopy(self.init_disp)

        for i, known_f_a in enumerate(self.known_f_a):
            self.displacement[known_f_a] = self.dis_new[i]

        # Deformed shape
        self.nodalcoord_def = []
        for i in range(self.nodenum):
            self.nodalcoord_def.append([self.nodalcoord[i][0]+ self.displacement[i*3+0], \
                self.nodalcoord[i][1]+ self.displacement[i*3+1], self.nodalcoord[i][2]+ self.displacement[i*3+2]])

        # Postrpocesses
        self.postprocess()

        self.mod_displacements = [0.]*(self.elenum+1)


    def solvemodstruct(self, index):
        """
        Solver for the modified structures. 'Index' shows the actual modification number.
        """

        self.mod_displacements[index] = [0.]*(self.nodenum*3)

        dis_new = [0.]*(self.nodenum*3-len(self.constraint))
        stiff_new = [[0.]*(self.nodenum*3-len(self.constraint))]*(self.nodenum*3-len(self.constraint))

        stiffincrement = [0.]*(self.nodenum*3-len(self.constraint))
        for i, kfai in enumerate(self.known_f_a):
            for j, kfaj in enumerate(self.known_f_a):
                stiffincrement[j] = self.mod_stiffnesses[index][kfai][kfaj]
            stiff_new[i] = [x + y for x, y in zip(stiff_new[i], stiffincrement)]

        # SOLVING THE MODIFIED STRUCTURE
        if Conf.solver == 0:
            dis_new = mat_vec_mult(invert(stiff_new), self.force_new)
        else:
            dis_new = np.linalg.solve(np.array(stiff_new), np.array(self.force_new))

        mod_displacement_temp = deepcopy(self.init_disp)

        for i, kfa in enumerate(self.known_f_a):
            mod_displacement_temp[kfa] = dis_new[i] - self.dis_new[i]

        self.mod_displacements[index] = [x + y for x, y in zip(self.mod_displacements[index], mod_displacement_temp)]



    def evaluate(self):
        """
        Calculates the relative displacement of each individual available unit-modification
        compared to the measured differencies (delta).

        delta: [DOF number, difference]

        return effect: [efffect on 1. point, effect on 2. point, ..., modification number]
                       where each line number shows the corresponding modification number
        """
        self.effect = [[0.]*(self.keypnum + 2)]*self.elenum
        self.toteffect = [0.]*self.keypnum
        self.sortedeff = [[[0.]*(self.keypnum + 2)]*self.elenum]*self.keypnum

        effect_temp = [0.]*(self.keypnum + 2)

        for modnum in range(self.elenum):
            effect_temp[self.keypnum] = int(modnum)
            for j, dofnum in enumerate(self.keypoint):
                try:
                    effect_temp[j] = self.mod_displacements[modnum][dofnum]
                    self.effect[modnum] = [x for x in effect_temp]
                    self.toteffect[j] += abs(self.effect[modnum][j])
                except IndexError:
                    print("Maybe the mapping data is invalid.")
                    print("Please check the \'arduino_mapping.txt\' input whether the given DOFs are correct or not.")
                    raise IndexError

        self.effectratio = deepcopy(self.effect)
        for i in range(self.elenum):
            for j in range(self.keypnum):
                if self.toteffect[j] > 0:
                    self.effectratio[i][j] = abs(self.effectratio[i][j]/self.toteffect[j])
                else:
                    self.effectratio[i][j] = 0

        #print("   \'effectratio\' is not used yet")

        # Sort by effectiveness
        for i in range(self.keypnum):
            self.sortedeff[i] = deepcopy(self.effect)

            # Check sign of the effect
            for ktemp in range(self.elenum):
                if self.sortedeff[i][ktemp][i] < 0:
                    for jtemp in range(self.keypnum):
                        self.sortedeff[i][ktemp][jtemp] = abs(self.sortedeff[i][ktemp][jtemp])
                        self.sortedeff[i][ktemp][self.keypnum +1] = -1
                else:
                    self.sortedeff[i][ktemp][self.keypnum +1] = +1

            for j in range(self.keypnum):
                if i != j and j != 0:
                    self.sortedeff[i] = swap_col(sorted(swap_col(self.sortedeff[i], 0, j), reverse=True), 0, j)
            if i != 0:
                self.sortedeff[i] = swap_col(sorted(swap_col(self.sortedeff[i], 0, i), reverse=True), 0, i)
            else:
                self.sortedeff[i] = sorted(self.sortedeff[i], reverse=True)

    def difference(self, num_displ, measurement):
        """
        Calculate the difference between the Numerical solution and Real-life measurement.
        The Real-life measurement should be given the following format:

            MEASUREMENT: [[13X, -2.154], [16Y, 5.256], ...]
        """
        #print(nodenumber option should be added! <XXX>)

        delta = []

        for loc, measured in measurement:
            try:
                dof = self.analysis[loc.upper()]
            except KeyError:
                print('The given measurement location cannot be matched with the input data.')
                print('The available nodes are: {\'NAMES\': mapping addresses}')
                print(self.analysis)
                SER.close()
                raise Exception('Watchpoint name error')

            delta.append(measured - num_displ[dof])

        return delta

    def optimize(self, delta):
        """
        Modell updating - core function
        """
        #modnum = min(10, self.elenum)
        modnum = self.elenum
        self.modifications = [0.0]*self.elenum



        if not _SIMULATION:
            appendix = ""
        else:
            appendix = " - SIMULATED"

        newdelta = delta
        j = 0

        print("-----")
        print("Step: 0/"+ str(self.iterationlimit))

        while (error(newdelta) > self.errorlimit and j <= self.iterationlimit and (self.capable() or j <= 1)):     # Optimization loop
            j += 1

            print("Error: " + str(error(newdelta)))
            print("-----")
            print("Step: " + str(j) + "/"+ str(self.iterationlimit))

            ratio = [0.]*modnum
            unit = 0

            prevmodifications = self.modifications

            for index in range(self.elenum):
                self.modifications[index] = min(abs(self.modifications[index] - self.unitmodification), self.modificationlimit) *math.copysign(1, self.modifications[index]- self.unitmodification)
                self.calcmodstiffness(index, self.modifications[index])
                self.solvemodstruct(index)
            self.evaluate()

            self.calcmodstiffness(self.elenum, 0)
            self.solvemodstruct(self.elenum)

            newdelta = self.difference(self.mod_displacements[self.elenum], self.measurement)

            for i, effect in enumerate(self.toteffect):
                if effect == 0.0:
                    print("None of the variables has effect on " + str(self.arduino_mapping[i]))
                    print("Model updating has no solution.")
                    raise Exception

            for i in range(self.elenum):
                modificationnumber = self.sortedeff[0][i][1]
                ratio[modificationnumber] = abs(self.sortedeff[0][i][0] / self.toteffect[0])*math.copysign(1, self.sortedeff[0][i][2])
                unit += abs(ratio[modificationnumber]*self.sortedeff[0][i][0])

            scale = newdelta[0]/unit
            for i in range(self.elenum):
                modificationnumber = self.sortedeff[0][i][1]
                self.modifications[modificationnumber] = min(abs(prevmodifications[modificationnumber] - self.unitmodification*ratio[modificationnumber]), self.modificationlimit)\
                    *math.copysign(1, prevmodifications[modificationnumber] - self.unitmodification*ratio[modificationnumber])
                    # the last part is already the sign itself without the sign function

            print("Ratio: " + str(scale))

        print("Final error: " + str(error(newdelta)))

        if not self.capable() and j > 1:
            print("Optimization could not be finished successfully.")
            print("The remaining error is: " + str(error(newdelta)))


        with open(self.name + ' - UpdateResults'+ appendix +'.txt', 'a') as outfile:
            if j > 1:
                if j <= self.iterationlimit and self.capable():
                    self.numofupdates[0] += 1
                    outfile.write("Update state: SUCCESSFUL\n")
                if not j <= self.iterationlimit:
                    self.numofupdates[1] += 1
                    outfile.write("Update state: Run out of iteration limit\n")
                if not self.capable() and j > 1:
                    self.numofupdates[2] += 1
                    outfile.write("Update state: No more possible modification\n")
            else:
                outfile.write("Update state: Optimization was skipped\n")
            outfile.write("Requiered iterations: " + str(j) + "\n")
            outfile.write("Measurement: " + str(self.measurement) + "\n")
            outfile.write("Original delta: " + str(delta) + "\n")
            outfile.write("New delta: " + str(newdelta) + " (limit: " + str(self.errorlimit) +")\n")
            outfile.write("Final error: " + str(error(newdelta)) + "\n")
            outfile.write("Modifications [%]: \n")
            outfile.write(str(self.modifications) + "\n")
            outfile.write("Original displacements: \n")
            outfile.write(str(self.displacement) + "\n")
            if j > 1:
                outfile.write("New displacements: \n")
                outfile.write(str(self.mod_displacements[self.elenum]) + "\n")
            outfile.write("----------------------\n")


    def capable(self):
        """
        Function telling whether there are more options to modify
        """
        capable = False
        for variable in self.modifications:
            if abs(variable) <= 0.95*self.modificationlimit and abs(variable) > 0.01:
                capable = True
        return capable

    def seterrorlimit(self, errorlimit):
        """
        Setting general stop parameter for model updating
        """
        if errorlimit > 0.0:
            self.errorlimit = errorlimit
        else:
            print("The error limit must be a positive number")
            raise Exception

    def setmodificationlimit(self, modificationlimit):
        """
        Setting modification limit for members (model updating)
        """
        if modificationlimit > 0.0 and modificationlimit < 1.0:
            self.modificationlimit = modificationlimit
        else:
            print("The modification limit must be higher than 0.0 and lower than 1.0")
            raise Exception

    def setunitmodification(self, unitmodification):
        """
        Setting modification step (model updating)
        """
        if abs(unitmodification) >= 0.01 and abs(unitmodification) < 0.5:
            self.unitmodification = unitmodification
        else:
            print("The absolut value of the unit modification must be minimum 0.01 and maximum 0.5")
            raise Exception

    def setiterationlimit(self, iterationlimit):
        """
        Setting maximum number of iterations (model updating)
        """
        if int(iterationlimit) > 1 and int(iterationlimit) <= math.pow(10, 4):
            self.iterationlimit = int(iterationlimit)
        else:
            print("The iterationlimit must be between 2 and 10.000")
            raise Exception

    def readarduino(self, base, saveinput):
        """
        Read data from Arduino
        """
        # Read data from Arduino

        maxdifference = 0.8          # Maximum input difference treshold in mm
        arduinovalues = []
        data = [0.]*len(self.arduino_mapping)
        newdata = False
        bigdifference = False
        readerror = False

        try:
            arduinoline = SER.readline()
            if len(arduinoline) > 0:
                arduinovalues = arduinoline.split(',')
                try:
                    if arduinovalues[0][len(arduinovalues)-1] == '.':
                        arduinovalues[0] = arduinovalues[0][:len(arduinovalues[0])-2]
                    else:
                        del arduinovalues[len(arduinovalues)-1]
                except IndexError:
                    print("Index Error... continuing")
                if len(arduinovalues) == len(self.arduino_mapping):
                    try:
                        for i in range(len(self.arduino_mapping)):
                            data[i] = float(arduinovalues[i]) - float(base[i][1])

                            if abs(data[i] - self.processeddata[i]) > maxdifference:
                                bigdifference = True
                            if abs(float(arduinovalues[i])) < 2.0:
                                readerror = True

                        self.processeddata = data
                        newdata = True

                    except ValueError:
                        print("Value error... continuing")
                        SER.flushInput()
                        time.sleep(0.5)
                    except Exception:
                        print("Type error: " + str(arduinovalues) + "... continuing")
                        SER.flushInput()
                        time.sleep(0.5)

            SER.flushInput()

        except serial.SerialTimeoutException:
            print("Data could not be read... continuing")
            SER.flushInput()
            time.sleep(0.5)


        if newdata and not bigdifference and not readerror:
            self.measurement = zip(self.arduino_mapping, data)

            saveinput.write(str(data) +', '+ str(time.time()) + "\n")
            # Calculate differences
            delta = self.difference(self.displacement, self.measurement)

            print("Delta: " + str(delta))

            newdata = False
            bigdifference = False
            readerror = False
            return delta

        newdata = False
        bigdifference = False
        readerror = False

    def simulatearduino(self, arduinoline, prevline):
        """
        Simulate data, based on previous measurement
        """
        arduinovalues = []
        data = [0.]*len(self.arduino_mapping)

        skip = 0
        sleeptime = 0.
        try:
            try:
                prevreadtime = float(str(prevline.split(']')[1]).split(',')[1])
                nowreadtime = float(str(arduinoline.split(']')[1]).split(',')[1])
                try:
                    if Conf.realistic_simulation:
                        sleeptime = nowreadtime - prevreadtime
                except Exception:
                    pass
            except Exception:
                skip = 1
                sleeptime = 0.

            if not skip:
                if not sleeptime > 0:
                    sleeptime = 0.
                arduinoline = str(arduinoline.split(']')[0])+"]"
                arduinovalues = eval(arduinoline)
                try:
                    for i in range(len(self.arduino_mapping)):
                        data[i] = float(arduinovalues[i])
                    self.processeddata = data
                except Exception:
                    print("Type error: " + str(arduinovalues) + "... continuing")

                self.measurement = zip(self.arduino_mapping, data)

                # Calculate differences
                delta = self.difference(self.displacement, self.measurement)
                time.sleep(sleeptime)
                print(delta)
                return delta

        except IndexError:
            print("IndexError")
            #pass
        except Exception:
            print("Exception in simulation data")
            #pass

    def updatemodel(self):
        """
        General function to manage model updatin procedure.
        """
        self.processeddata = [0.]*len(self.arduino_mapping)

        if not _SIMULATION:
            base = self.calibrate()
            filemode = 'a'
        else:
            base = ['SIMULATION']
            try:
                os.remove(self.name + ' - UpdateResults - SIMULATED.txt')
            except Exception:
                pass
            filemode = 'r'

        with open(self.name + ' - Input Data.txt', filemode) as inputfile:
            # Saving input data
            if not _SIMULATION:
                inputfile.write('Input data of \'' + self.name + '\':\n\n')
                inputfile.write('Start Time: ' + str(datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')) + "\n")
                inputfile.write('Base: ' + str(base) + '\n')

            inputline = "[],0.0"
            for i in range(1000):
                if not _SIMULATION:
                    delta = self.readarduino(base, inputfile)
                    self.optimize(delta)
                else:
                    try:
                        delta = None
                        previnputline = inputline
                        inputline = inputfile.readline()
                        if not inputline == '':
                            delta = self.simulatearduino(inputline, previnputline)
                    except Exception:
                        pass
                    if not delta is None:
                        self.optimize(delta)
        print("Update statistics:")
        print("Totally updated models: " + str(TRUSS.numofupdates[0] + TRUSS.numofupdates[1]++ TRUSS.numofupdates[2]))
        print("  Successfully updated models: " + str(TRUSS.numofupdates[0]))
        print("  Updates with running out of possibilities: " + str(TRUSS.numofupdates[2]))
        print("  Updates did not finshed: " + str(TRUSS.numofupdates[1]))

    def calibrate(self):
        """
        Calibration for Arduino measurement. All the measurements will describe the dispalecements from the claibration-state.
        """
        answer_1 = '0'
        restart = '0'
        accept = '0'
        arduinovalues = []
        print("Before starting the model updating, the measuring tools must be calibrated.")
        print("The calibration should be done in load-free state.")
        while (answer_1 not in ['Y', 'N']):
            answer_1 = input('Can we start the calibration? (y/n) ').upper()
        if answer_1 == 'N':
            SER.close()
            raise Exception('Calibration is terminated')
        else:
            try:
                SER.flushInput()
                #time.sleep(0.2)
                arduinoline = '' #SER.readline()
                while len(arduinoline) == 0:
                    time.sleep(0.2)
                    arduinoline = SER.readline()

                if len(arduinoline) > 0:
                    arduinovalues = arduinoline.split(',')
                    del arduinovalues[len(arduinovalues)-1]               # if needed!!!
                    if len(arduinovalues) == len(self.arduino_mapping):
                        measurement = zip(self.arduino_mapping, arduinovalues)
                        print("Calibration result:")
                        print(measurement)
                        while (accept not in ['Y', 'N']):
                            accept = input('Ok? Can we start the main part? Put on the loads! (y/n) ').upper()
                            if accept == 'N':
                                restart = 'Y'
                    else:
                        print("Data error. Calibartion is restarting.")
                        print("Arduino values:" + str(arduinovalues))
                        restart = 'Y'
                else:
                    print('The calibration cannot be done: no data')
                    while (restart not in ['Y', 'N']):
                        restart = input('Do you want to restart calibration? (y/n) ').upper()
            except Exception:
                print('The calibration cannot be done: exception was raised')
                while (restart not in ['Y', 'N']):
                    restart = input('Do you want to restart calibration? (y/n) ').upper()

        if restart == 'Y':
            print("Restarting calibration")
            SER.flushInput()
            return self.calibrate()
        elif restart == 'N':
            SER.close()
            raise Exception('Calibration is terminated')
        if accept == 'Y':
            return measurement


    def postprocess(self):
        """
        Calculates reaction forces and stresses
        """
        self._reactions()
        self._stresses()
        self._postprocessed = 1

    def _reactions(self):
        """
        Calculates reaction forces
        """
        for i in self.known_dis_a:
            self.force[i] = 0
            for j, displ in enumerate(self.displacement):
                self.force[i] += self.stiffness[i][j]*displ

    def _stresses(self):
        """
        Calculates stress in elements

        Last part: Coloring elements for graphical output
        """
        self.stress = [0.]*self.elenum
        for element in range(self.elenum):
            locstiff = [-self._cx[element], -self._cy[element], -self._cz[element], \
                         self._cx[element], self._cy[element], self._cz[element]]
            for i in range(3*2):
                self.stress[element] += locstiff[i]*self.displacement[self.eledof[element][i]]
            self.stress[element] = self.stress[element]*self._norm_stiff[element]

        smax = max([abs(min(self.stress)), max(self.stress), 0.000000001])
        self.stresscolor = [float(x)/float(smax) for x in self.stress]

    def postprocessed(self):
        """
        Tells if the structure's postprocess part is already calcuated
        """
        return self._postprocessed

    def writeresults(self, fname):
        """
        Writing results to file.
        """
        out_element = ''
        for i in self.node:
            out_element += str(i[0] + self._io_origin) + ', ' + str(i[1] + self._io_origin) + ' | '
        out_coords = ''
        for i in self.nodalcoord:
            out_coords += str(i[0]) + ', ' + str(i[1]) + ', ' + str(i[2]) + ' | '
        out_crsect = ''
        for i in self.area:
            out_crsect += str(i) + ', '
        out_materials = ''
        for i in self.el_mod:
            out_materials += str(i) + ', '
        out_forces = ''
        for forcedof in self.known_f_notzero:
            if self.dof == 3:
                out_forces += str(forcedof + self._io_origin) + ', ' + str(self.force[forcedof]) + ' | '
            elif self.dof == 2 and i % 3 != 2:
                out_forces += str(forcedof - forcedof//3 + self._io_origin) + ', ' + str(self.force[forcedof]) + ' | '
        out_supports = ''
        for i in self.constraint:
            if self.dof == 3:
                out_supports += str(i[0] + self._io_origin) + ', ' + str(i[1]) + ' | '
            elif i[0] % 3 != 2:
                out_supports += str(i[0] - i[0]//3 + self._io_origin) + ', ' + str(i[1]) + ' | '
        # Not elegant solution
        out_specdofs = self.specdof_inputstring
        try:
            with open(fname, 'w') as outfile:
                # Writing data
                outfile.write('Calculation of \'' + self.name + '\':\n\n')
    
                outfile.write('Reactions\n')
                #for i in range(len(self.force)//3):
                prev = -1
                for i in self.known_dis_a:
                    if self.dof == 3 or i%3 != 2:
                        if i//3 != prev:
                            if i < 100:
                                outfile.write(' ')
                                if i < 9:
                                    outfile.write(' ')
                            nodalforce = ''
                            if (i//3)*3+0 in self.known_dis_a:
                                nodalforce += "{:10.2f}".format(self.force[(i//3)*3+0]) + ', '
                            else:
                                nodalforce += '            '
                            if (i//3)*3+1 in self.known_dis_a:
                                nodalforce += "{:10.2f}".format(self.force[(i//3)*3+1]) + ', '
                            else:
                                nodalforce += '            '
                            if self.dof != 2 and (i//3)*3+2 in self.known_dis_a:
                                nodalforce += "{:10.2f}".format(self.force[(i//3)*3+2]) + '\n'
                            else:
                                nodalforce += '          \n'
                            if nodalforce != '                                  \n':
                                outfile.write(str(i//3 + self._io_origin) + ', ' +  nodalforce)
                        prev = i//3
                outfile.write('\n')
    
                outfile.write('Displacements\n')
                for i in range(len(self.displacement)//3):
                    if i < 100:
                        outfile.write(' ')
                        if i < 9:
                            outfile.write(' ')
                    outfile.write(str(i + self._io_origin) + ', ' +  "{:10.3f}".format(self.displacement[i*3 +0]) + ', ' \
                          + "{:10.3f}".format(self.displacement[i*3 +1]) + ', ' + "{:10.3f}".format(self.displacement[i*3 +2]) + ', ' + '\n')
                outfile.write('\n')
    
                outfile.write('Stresses\n')
                for i, stress in enumerate(self.stress):
                    if i < 100:
                        outfile.write(' ')
                        if i < 9:
                            outfile.write(' ')
                    outfile.write(str(i + self._io_origin) + ', ' +  "{:10.3f}".format(stress) + '\n')
                outfile.write('\n')
    
                # Saving original input
                outfile.write('----- Original input: -----\n\n')
                outfile.write('_ORIGIN\n')
                outfile.write(str(self._io_origin) + '\n\n')
                outfile.write('DOF\n')
                outfile.write(str(self.dof) + '\n\n')
                outfile.write('ELEMENTS\n')
                outfile.write(out_element + '\n\n')
                outfile.write('COORDINATES\n')
                outfile.write(out_coords + '\n\n')
                outfile.write('CROSS-SECTIONS\n')
                outfile.write(out_crsect + '\n\n')
                outfile.write('MATERIALS\n')
                outfile.write(out_materials + '\n\n')
                outfile.write('FORCES\n')
                outfile.write(out_forces + '\n\n')
                outfile.write('SUPPORTS\n')
                outfile.write(out_supports + '\n\n')
                outfile.write('SPECDOF\n')
                outfile.write(out_specdofs + '\n\n')
                outfile.write('EOF\n')
        except FileNotFoundError:
            print("Error: Please manually create the 'Results' folder" 
              "in the root of the project")
            

##################################
#   BEGINNING OF THE MAIN PART   #
##################################

PARTTIME = logtime(TIC, "Initialization")

#  Define new truss
TRUSS = Truss('bridge')
if not Conf.debug:
    TRUSS.name = input('Test name: ')
else:
    print("*** Debug mode ***")
    print("*** The following file will be opened: " + TRUSS.name + ".str")

# Read input file
#TRUSS.read('lab_01.txt')

try:
    TRUSS.read(TRUSS.name + ".str")
except IOError:
    print("The following file could not be opened: " + TRUSS.name + ".str")
    print("Please make sure that the structural data is available for the program in the running directory.")
    raise IOError

#if _ARDUINO or _SIMULATION:                # deprecated
#    TRUSS.setspecdofs(arduino_mapping)

PARTTIME = logtime(PARTTIME, "Setting up structure")

# Calculate stiffness-matrix
TRUSS.calcstiffness()
#TRUSS.calcstiffness_plate()

PARTTIME = logtime(PARTTIME, "Calculating Stiffness Matrix")

#Solve structure
TRUSS.solve()
#TRUSS.solve_plate()

PARTTIME = logtime(PARTTIME, "Solving")

if Conf.updating:
    TRUSS.setunitmodification(0.05)
    TRUSS.seterrorlimit(1.2)
    TRUSS.setmodificationlimit(0.7)
    TRUSS.setiterationlimit(100)
    TRUSS.updatemodel()

PARTTIME = logtime(PARTTIME, "Updating numerical model")

if Conf.graphics:
    # Plot settings:
    # O: Original D: Deformed S: Supports F: Forces R: Reactions
    # ScD: Scale displacments (Z-axis) (def:1.0) ScF: Scale forces (def:1.0)
    # ScS: Scale Support signs (Z-axis) (def:1.0)
    # Save: Save plot to file

    #     plot(O, D, S, F, R, ScD, ScF, ScS, Save)
    TRUSS.plot(1, 0, 1, 1, 0, 1.0, 0.0, 0.0, True)
    TRUSS.plot(1, 1, 1, 0, 0, 1.0, 0.0, 0.0, True)
    TRUSS.plot(0, 1, 1, 1, 1, 2.0, 0.0, 0.0, True)
    #pass

PARTTIME = logtime(PARTTIME, "Plotting")

# Write results to file
TRUSS.writeresults("./Results/" + TRUSS.name + ' - Results.txt')

PARTTIME = logtime(PARTTIME, "Writing results to the output file")

if _ARDUINO:
    # Closing Arduino port
    SER.close()

if Conf.log:
    TAC = time.time()
    TOTALTIME = TAC-TIC
    if Conf.updating:
        print("Update statistics:")
        print("Totally updated models: " + str(TRUSS.numofupdates[0] + TRUSS.numofupdates[1]++ TRUSS.numofupdates[2]))
        print("  Successfully updated models: " + str(TRUSS.numofupdates[0]))
        print("  Updates with running out of possibilities: " + str(TRUSS.numofupdates[2]))
        print("  Updates did not finshed: " + str(TRUSS.numofupdates[1]))
    print('Total time: ' + str("{:10.3f}".format(TOTALTIME)))
