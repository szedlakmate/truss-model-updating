# -*- coding: utf-8 -*-
"""
Nightly version of truss.py code created by Máté Szedlák (2016). Copyright and other notes are the same as in the main version.

This code should be used for simulations. Please remember, that the simulation should be configured at:
"def calcstiffness_plate(self):"
"""
#Code for Processing only
#size(640, 360)
#background(126)

import math
import itertools
from copy import deepcopy


# COMPATIBILITY MODES:
    # 0: User defined
    # 1: Processing3
    # 2: Android
    # 3: Most information (with numpy)
    # 4: Maximum compatibility

_COMPATIBLE_MODE = 0

_SIMULATION = 0                     # Simulating measurements

PORT = 'COM1'

fakedata = [['13X', -2.6], ['13Y', +1.0]]

if _COMPATIBLE_MODE == 00:
    ### User defined ###
    # Modify as needed #
    _MODE_NAME = "User defined"
    _LOG = 1                # Logging time
    _GRAPHICS = 0           # Graphical features
    _SOLVER = 1             # 0: Basic solver, 1: NumPy solver
    _OSLIBRARYAVAILABLE = 1 # Basic OS file features (e.g. file size)
    _UPDATING = 1           # Model Updating: On/ Off
    _ARDUINO = 1            # Arduino input: On/Off

##############################################
           ###  DO NOT MODIFY ###            #
                                             #
elif _COMPATIBLE_MODE == 01:                 #
    ### "Processing 3" mode ###              #
    _MODE_NAME = "Processing 3"              #
    _LOG = 1*1                               #
    _GRAPHICS = 0*0                          #
    _SOLVER = 0*0                            #
    _OSLIBRARYAVAILABLE = 0*0                #
    _UPDATING = 1*1                          #
    _ARDUINO = 1*1                           #
                                             #
elif _COMPATIBLE_MODE == 02:                 #
    ### Android mode ###                     #
    # DO NOT MODIFY                          #
    _MODE_NAME = "Android"                   #
    _LOG = 1*1                               #
    _GRAPHICS = 0*0                          #
    _SOLVER = 0*0                            #
    _OSLIBRARYAVAILABLE = 1*1                #
    _UPDATING = 0*0                          #
    _ARDUINO = 0*0                           #
                                             #
elif _COMPATIBLE_MODE == 03:                 #
    ### Most information ###                 #
    # DO NOT MODIFY                          #
    _MODE_NAME = "Most information"          #
    _LOG = 1*1                               #
    _GRAPHICS = 1*1                          #
    _SOLVER = 1*1                            #
    _OSLIBRARYAVAILABLE = 1*1                #
    _UPDATING = 1*1                          #
    _ARDUINO = 1*1                           #
                                             #
else:                                        #
    ### Maximum compatibility ###            #
    # DO NOT MODIFY                          #
    _MODE_NAME = "Maximum compatibility"     #
    _LOG = 0*0                               #
    _GRAPHICS = 0*0                          #
    _SOLVER = 0*0                            #
    _OSLIBRARYAVAILABLE = 0*0                #
    _UPDATING = 0*0                          #
    _ARDUINO = 0*0                           #
                                             #
           ###  DO NOT MODIFY ###            #
##############################################


if _OSLIBRARYAVAILABLE:
    import os

if _SIMULATION:
    _ARDUINO = 0

if _COMPATIBLE_MODE == 2:
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

if _LOG:
    import time
    import datetime
    TIC = time.time()
    print '------------------------------------'
    print 'Truss calculational program'
    print 'Created by Máté Szedlák (23/11/2016)'
    print 'Compatibility mode: ' + _MODE_NAME
    if _SOLVER == 0:
        print '- Solver is set to default'
    elif _SOLVER == 1:
        print '- Solver is set to NumPy'
    else:
        raise Exception("Solver settings are invalid!")
    if _UPDATING:
        print '+ Model updating is turned ON'
        if _SIMULATION:
            print 'Input data is SIMULATED!'
    else:
        print '- Model updating is turned OFF'
    print '------------------------------------'
    
    
if _ARDUINO:
    SER = 0
    try:
        import serial
    except ImportError:
        print "You tried to import \'serial\' in Windows mode without installing \'pySerial\'."
        print "Please first install pySerial. Try this: http://playground.arduino.cc/Interfacing/Python"
        raise Exception('Android mode denied: pyserial not found')

    while (SER == 0):
        time.sleep(0.6)
        print 'Opening serial at port ' + str(PORT)
        try:
            SER.close()
        except Exception:
            pass
        try:
            SER = serial.Serial(PORT, 9600, timeout=0)
        except serial.SerialException:
            Exception(PORT + ' port is busy. It might be occupied by this program or another one :/ Be careful or try resetting this program')
            SER = 0
        except Exception:
            SER = 0

            
if _ARDUINO or _SIMULATION:            
    try:
        mappingfile = 'arduino_mapping.txt'
        with open(mappingfile, "r") as textfile:
            line = textfile.readline().strip()
            arduino_mapping = line.upper().split(',')
    except IOError:
        raise Exception('File not found: ' + mappingfile)


if _SOLVER:
    # NumPy library for solving linear equations in another way
    import numpy as np

if _GRAPHICS:
    # libraries for drawing
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch
    from mpl_toolkits.mplot3d import proj3d

# From here:
def mat_vec_mult(mat_a, vec_b):
    """
    Multiplying matrix with a vector, giving the result as a vector

    Source:
    https://stackoverflow.com/questions/10508021/matrix-multiplication-in-python
    """
    vec_c = [0.]*len(mat_a)
    for i, row in enumerate(mat_a):
        for j, elem in enumerate(vec_b):
            vec_c[i] += row[j]*elem
    return vec_c

def invert(mat_x):
    """
    Invert a matrix X according to gauss-jordan elimination
    In gauss-jordan elimination, we perform basic row operations to turn a matrix into
    row-echelon form.  If we concatenate an identity matrix to our input
    matrix during this process, we will turn the identity matrix into our inverse.
    X - input list of lists where each list is a matrix row
    output - inverse of X

    Source:
    http://www.vikparuchuri.com/blog/inverting-your-very-own-matrix/
    """
    #copy X to avoid altering input
    mat_x = deepcopy(mat_x)

    #Get dimensions of X
    rows = len(mat_x)
    cols = len(mat_x[0])

    #Get the identity matrix and append it to the right of mat_x
    #This is done because our row operations will make the identity into the inverse
    identity = make_identity(rows, cols)
    for i in xrange(0, rows):
        mat_x[i] += identity[i]

    i = 0
    for j in xrange(0, cols):
        #print("On col {0} and row {1}".format(j, i))
        #Check to see if there are any nonzero values below the current row in the current column
        zero_sum, first_non_zero = check_for_all_zeros(mat_x, i, j)
        #If everything is zero, increment the columns
        if zero_sum == 0:
            if j == cols:
                return mat_x
            raise Exception("Matrix is singular")
        #If mat_x[i][j] is 0, and there is a nonzero value below it, swap the two rows
        if first_non_zero != i:
            mat_x = swap_row(mat_x, i, first_non_zero)
        #Divide mat_x[i] by mat_x[i][j] to make mat_x[i][j] equal 1
        mat_x[i] = [m/mat_x[i][j] for m in mat_x[i]]

        #Rescale all other rows to make their values 0 below mat_x[i][j]
        for k in xrange(0, rows):
            if k != i:
                scaled_row = [mat_x[k][j] * m for m in mat_x[i]]
                mat_x[k] = [mat_x[k][m] - scaled_row[m] for m in xrange(0, len(scaled_row))]
        #If either of these is true, we have iterated through the matrix, and are done
        if i == rows or j == cols:
            break
        i += 1

    #Get just the right hand matrix, which is now our inverse
    for i in xrange(0, rows):
        mat_x[i] = mat_x[i][cols:len(mat_x[i])]
    return mat_x

def check_for_all_zeros(mat_x, i, j):
    """
    Check matrix mat_x to see if only zeros exist at or below row i in column j
    mat_x - a list of lists
    i - row index
    j - column index
    returns -
        zero_sum - the count of non zero entries
        first_non_zero - index of the first non value
    """
    non_zeros = []
    first_non_zero = -1
    for k in xrange(i, len(mat_x)):
        non_zero = mat_x[k][j] != 0
        non_zeros.append(non_zero)
        if first_non_zero == -1 and non_zero:
            first_non_zero = k
    zero_sum = sum(non_zeros)
    return zero_sum, first_non_zero

def swap_row(mat_x, i, j):
    """
    Swap row i and row j in a list of lists
    mat_x - list of lists
    i - row index
    j - row index
    returns- modified matrix
    """
    mat_x[j], mat_x[i] = mat_x[i], mat_x[j]
    return mat_x

def swap_col(mat_x, i, j):
    """
    Swap colum i and column j in a list of lists
    mat_x - list of lists
    i - column index
    j - column index
    returns- modified matrix
    """
    for item in mat_x:
        item[i], item[j] = item[j], item[i]
    return mat_x

def make_identity(row_num, col_num):
    """
    Make an identity matrix with dimensions rxc
    row_num - number of rows
    col_num - number of columns
    returns - list of lists corresponding to  the identity matrix
    """
    identity = []
    for i in xrange(0, row_num):
        row = []
        for j in xrange(0, col_num):
            elem = 0
            if i == j:
                elem = 1
            row.append(elem)
        identity.append(row)
    return identity

if _GRAPHICS:
    class Arrow3D(FancyArrowPatch):
        """
        Vector drawer module from the internet
        """
        def __init__(self, xs, ys, zs, *args, **kwargs):
            FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
            self._verts3d = xs, ys, zs

        def draw(self, renderer):
            _xs3d, _ys3d, _zs3d = self._verts3d
            _xs, _ys, _zs = proj3d.proj_transform(_xs3d, _ys3d, _zs3d, renderer.M)
            self.set_positions((_xs[0], _ys[0]), (_xs[1], _ys[1]))
            FancyArrowPatch.draw(self, renderer)

    def plotstructure(struct, showorig, showresult, showsupports, \
            showforces, showreactions, scaledisp, scale_f, z_corr, showvalues, saveplot):
        """
        General plotting method for structures

            scaledisp: Scale drwaing of displacements
            scale_f:   Scale force sign
            z_corr:    Scale z-axis
        """
        plotname = struct.name
        plot_width = 18.0              # Plot width in inches
        xframe = 0                     # Frame width at X direction
        yframe = 0                     # Frame width at Y direction
        zframe = 0                     # Frame width at Z direction
        scale_sup = scale_f*0.3        # Scale support sign  # All the others are input parameters

        # Stress coloring settings [R G B] - Examples considering pressure:
        # 0: [1, 0, 0]              Plain red
        # 1: [x, 0, 0]              Red to Black
        # 2: [1, 1-x, 0-x]          Red to White
        # 3: [1, (1-x)/2, (1-x)/2]  Red to MildRed - Distincts pressure and tension
        # 4: [x, 1-x, 0]            Red to Green
        _coloring = 3 # € [0, 1, 2, 3, 4]

        fig = plt.figure()
        _ax = fig.add_subplot(111, projection='3d')

        if struct.dof == 2:
            _ax.view_init(elev=90., azim=-90.)
            _ax.w_zaxis.line.set_lw(0.)
            _ax.set_zticklabels([])

        xmin = min(list(struct.nodalcoord[x][0] for x in range(struct.nodenum)))
        xmax = max(list(struct.nodalcoord[x][0] for x in range(struct.nodenum)))
        ymin = min(list(struct.nodalcoord[x][1] for x in range(struct.nodenum)))
        ymax = max(list(struct.nodalcoord[x][1] for x in range(struct.nodenum)))
        zmin = min(list(struct.nodalcoord[x][2] for x in range(struct.nodenum)))
        zmax = max(list(struct.nodalcoord[x][2] for x in range(struct.nodenum)))

        deltax = xmax - xmin
        deltay = ymax - ymin
        xframe = max(deltax * 0.05, 2)
        yframe = max(deltay * 1.5, 2)
        plot_height = plot_width * ((deltay + yframe*2)/max(300, (deltax + xframe*2)))
        fig.set_size_inches(plot_width, plot_height)

        _ax.set_xlim3d(xmin - xframe, xmax + xframe)
        _ax.set_ylim3d(ymin - yframe, ymax + yframe)
        _ax.set_zlim3d(zmin - zframe, zmax + zframe)

        if showorig == showresult:
            _coloring = 0

        # Giving plot names
        if showorig == 1 and showresult == 0 and showsupports == 1 and showreactions == 0:
            plotname += ' - Initial structure'
            if showforces:
                plotname += ' with forces'
        elif showorig == 1 and showresult == 1:
            plotname += ' - Deformation'
            if showreactions == 0:
                plotname += ' with reactions'
        elif showorig == 0 and showresult == 1:
            plotname += ' - Stresses'
            if showreactions == 0:
                plotname += ' with reactions'
        else:
            plotname += ' - Unnamed'

        print plotname + ": "
        if showresult:
            dipslaydisplacement = deepcopy(struct.nodalcoord_def)
            if scaledisp != 1.0:
                if _LOG:
                    print('Displacements are scaled with factor: ') + str(scaledisp)
                for i in range(struct.nodenum):
                    for j in range(3):
                        dipslaydisplacement[i][j] = (struct.nodalcoord_def[i][j] -\
                        struct.nodalcoord[i][j]) * scaledisp + struct.nodalcoord[i][j]

        for i in range(struct.elenum):
            # Plot undeformed structure
            if showorig:
                _ax.plot([struct.nodalcoord[struct.node[i][1]][0], struct.nodalcoord[struct.node[i][0]][0]], \
                    [struct.nodalcoord[struct.node[i][1]][1], struct.nodalcoord[struct.node[i][0]][1]], \
                    zs=[struct.nodalcoord[struct.node[i][1]][2], struct.nodalcoord[struct.node[i][0]][2]], color='b')
            # Plot deformed structure
            if showresult:
                if struct.postprocessed():
                    if struct.stresscolor[i] > 0:
                        if _coloring == 1:
                            rgb_col = [0, 0, abs(struct.stresscolor[i])]
                        elif _coloring == 2:
                            rgb_col = [1-abs(struct.stresscolor[i]), \
                                       1-abs(struct.stresscolor[i]), 1]
                        elif _coloring == 3:
                            rgb_col = [(1-abs(struct.stresscolor[i]))/2, \
                                       (1-abs(struct.stresscolor[i]))/2, 1]
                        elif _coloring == 4:
                            rgb_col = [0, 1-abs(struct.stresscolor[i]), \
                                       abs(struct.stresscolor[i])]
                        else:
                            rgb_col = [1, 0, 0]
                    else:
                        if _coloring == 1:
                            rgb_col = [abs(struct.stresscolor[i]), 0, 0]
                        elif _coloring == 2:
                            rgb_col = [1, 1-abs(struct.stresscolor[i]), \
                                       1-abs(struct.stresscolor[i])]
                        elif _coloring == 3:
                            rgb_col = [1, (1-abs(struct.stresscolor[i]))/2, \
                                       (1-abs(struct.stresscolor[i]))/2]
                        elif _coloring == 4:
                            rgb_col = [abs(struct.stresscolor[i]), \
                                       1-abs(struct.stresscolor[i]), 0]
                        else:
                            rgb_col = [1, 0, 0]
                else:
                    print 'Stresses are not calculated'
                    rgb_col = [1, 0, 0]
                _ax.plot([dipslaydisplacement[struct.node[i][1]][0], dipslaydisplacement[struct.node[i][0]][0]], \
                        [dipslaydisplacement[struct.node[i][1]][1], dipslaydisplacement[struct.node[i][0]][1]], \
                        zs=[dipslaydisplacement[struct.node[i][1]][2], dipslaydisplacement[struct.node[i][0]][2]], color=rgb_col)

        if showforces:
            for i in struct.known_f_notzero:
                if struct.force[i] < 0:
                    value = -1.0
                else:
                    value = 1.0
                if i % 3 == 0:
                    f_dir = [value*scale_f, 0., 0.]
                elif i % 3 == 1:
                    f_dir = [0., value*scale_f, 0.]
                else:
                    f_dir = [0., 0., value*scale_f*z_corr]
                f_arrow = Arrow3D([struct.nodalcoord[i//3][0], struct.nodalcoord[i//3][0] + f_dir[0]], \
                                  [struct.nodalcoord[i//3][1], struct.nodalcoord[i//3][1] + f_dir[1]], \
                                  [struct.nodalcoord[i//3][2], struct.nodalcoord[i//3][2] + f_dir[2]], \
                                  mutation_scale=20, lw=1, arrowstyle="-|>", color="k")
                _ax.add_artist(f_arrow)

        if showreactions:
            e_previous = -100
            for i in struct.known_dis_a:
                value = 0.0             # BIZTOSAN KELL IDE? <XXX>
                if struct.force[i] < 0:
                    value = -1.0
                elif struct.force[i] > 0:
                    value = 1.0
                if i % 3 == 0:
                    f_dir = [value*scale_f, 0., 0.]
                elif i % 3 == 1:
                    f_dir = [0., value*scale_f, 0.]
                else:
                    f_dir = [0., 0., value*scale_f*z_corr]
                if abs(struct.force[i]) > 0:
                    f_arrow = Arrow3D([struct.nodalcoord[i//3][0], struct.nodalcoord[i//3][0] + f_dir[0]], \
                                      [struct.nodalcoord[i//3][1], struct.nodalcoord[i//3][1] + f_dir[1]], \
                                      [struct.nodalcoord[i//3][2], struct.nodalcoord[i//3][2] + f_dir[2]], \
                                      mutation_scale=20, lw=1, arrowstyle="-|>", color="darkolivegreen")
                    _ax.add_artist(f_arrow)
                    if showvalues:
                        _ax.set_xticklabels([])
                        _ax.set_yticklabels([])
                        _ax.set_zticklabels([])
                        if not i//3 == e_previous//3:
                            if struct.dof == 3:
                                _ax.text(struct.nodalcoord[i//3][0], \
                                    struct.nodalcoord[i//3][1], \
                                    struct.nodalcoord[i//3][2], \
                                    "{:10.2f}".format(struct.force[(i//3)*3+0])+'\n'+\
                                    "{:10.2f}".format(struct.force[(i//3)*3+1])+'\n'+\
                                    "{:10.2f}".format(struct.force[(i//3)*3+2]),\
                                    fontsize=12, horizontalalignment='right')
                            elif struct.dof == 2:
                                _ax.text(struct.nodalcoord[i//3][0], \
                                    struct.nodalcoord[i//3][1], \
                                    struct.nodalcoord[i//3][2], \
                                    "{:10.2f}".format(struct.force[(i//3)*3+0])+'\n'+\
                                    "{:10.2f}".format(struct.force[(i//3)*3+1]),\
                                    fontsize=12, horizontalalignment='right')
                e_previous = i

        if showsupports:
            for i in struct.known_dis_a:
                if i % 3 == 0:
                    f_dir = [-1.0 * scale_sup, 0., 0.]
                    col = 'g'
                elif i % 3 == 1:
                    f_dir = [0., -1.0 * scale_sup, 0.]
                    col = 'y'
                else:
                    f_dir = [0., 0., -1.0 * scale_sup * z_corr]
                    col = 'brown'
                if i % 3 != 2 or struct.dof == 3:
                    _ax.plot([struct.nodalcoord[i//3][0], struct.nodalcoord[i//3][0]+f_dir[0]], \
                        [struct.nodalcoord[i//3][1], struct.nodalcoord[i//3][1]+f_dir[1]], \
                        zs=[struct.nodalcoord[i//3][2], struct.nodalcoord[i//3][2]+f_dir[2]], \
                        color=col, linewidth=4.0)
        plt.show()
        if saveplot:
            fig.savefig(plotname + '.png')
            print '\'' + plotname +'.png\' is saved.'
            print '------------------------------------'
        return

def endoffile(textfile, line):
    """
    Check if end of file is reached. Implemented due to compatibility reasons.
    """
    if _OSLIBRARYAVAILABLE:
        return textfile.tell() < os.fstat(textfile.fileno()).st_size
    else:
        return not line == "EOF"

def logtime(prev_time, title):
    """
    Calculating and printing the time consumption of tasks

    Should be called with the previously saved part-time and the name of the actual task
    At the first call, should be called with TIC value. The input argument
    should be overwritten by this funvtion's return value.
    """
    if _LOG:
        new_time = time.time()
        print title
        print 'Time: ' + str("{:10.3f}".format(new_time - prev_time))
        print '------------------------------------'
        return new_time
    else:
        return 0

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
            CROSS-SECTIONS - This data will be evaulated in Python: 3.0*(10**(-4)), 5.0*(10**(-4)) ...
            MATERIALS - This data will be evaulated in Python: 70.0*(10**9), 100.0*(10**9) ...
            FORCES - Selected DOF + Force: 11, +1000000.0 | 12, +1000000.0 ...
            SUPPORTS - Selected DOF + Prescribed displacement: 0, 0.0 | 1, 0.0 ...
            SPECDOF - Selected node's DOF will be analysed during Model Updating: 1, xyz | 3 y | 10 xz ...

            EOF - For compatibility reasons EOF should be placed after the commands
        """
        self._io_origin = 0

        with open(filename, "r") as textfile:
            line = ""
            while endoffile(textfile, line):
                line = textfile.readline().strip()

                if line.upper() == "_ORIGIN":
                    line = textfile.readline().strip()
                    self._io_origin = int(line)

                if line.upper() == "DOF":
                    line = textfile.readline().strip()
                    self.setdof(int(line))

                if line.upper() == "ELEMENTS":
                    line = textfile.readline().strip()
                    inpstr = []
                    inpnum = []
                    inpstr = [x.split(',') for x in line.split('|')]
                    if len(inpstr[0]) == 1:
                        inpstr = [x.split(';') for x in line.split('|')]
                    if [''] in inpstr:
                        inpstr.remove([''])
                    inpnum = [[int(x[0]) - self._io_origin, int(x[1]) - self._io_origin] for x in inpstr]
                    self.setelements(inpnum)

                if line.upper() == "COORDINATES":
                    line = textfile.readline().strip()
                    inpstr = []
                    inpnum = []
                    inpstr = [x.split(',') for x in line.split('|')]
                    if len(inpstr[0]) == 1:
                        inpstr = [x.split(';') for x in line.split('|')]
                    if [''] in inpstr:
                        inpstr.remove([''])
                    if self.dof == 3:
                        inpnum = [[float(x[0]), float(x[1]), float(x[2])] for x in inpstr]
                    elif self.dof == 2:
                        inpnum = [[float(x[0]), float(x[1]), 0.] for x in inpstr]
                    self.setcoordinates(inpnum)

                if line.upper() == "CROSS-SECTIONS":
                    line = textfile.readline().strip()
                    inpstr = []
                    inpnum = []
                    inpstr = line.split(',')
                    if len(inpstr) == 1:
                        inpstr = line.split(';')
                    if '' in inpstr:
                        inpstr.remove('')
                    inpnum = [float(eval(x)) for x in inpstr]
                    self.setcrosssections(inpnum)

                if line.upper() == "MATERIALS":
                    line = textfile.readline().strip()
                    inpstr = []
                    inpnum = []
                    inpstr = line.split(',')
                    if len(inpstr) == 1:
                        inpstr = line.split(';')
                    if '' in inpstr:
                        inpstr.remove('')
                    inpnum = [float(eval(x)) for x in inpstr]
                    self.setmaterials(inpnum)

                if line.upper() == "FORCES":
                    line = textfile.readline().strip()
                    inpstr = []
                    inpnum = []
                    inpstr = [x.split(',') for x in line.split('|')]
                    if len(inpstr[0]) == 1:
                        inpstr = [x.split(';') for x in line.split('|')]
                    if [''] in inpstr:
                        inpstr.remove([''])
                    inpnum = [[int(x[0]) - self._io_origin, float(x[1])] for x in inpstr]
                    self.setforces(sorted(inpnum))

                if line.upper() == "SUPPORTS":
                    line = textfile.readline().strip()
                    inpstr = []
                    inpnum = []
                    inpstr = [x.split(',') for x in line.split('|')]
                    if len(inpstr[0]) == 1:
                        inpstr = [x.split(';') for x in line.split('|')]
                    if [''] in inpstr:
                        inpstr.remove([''])
                    inpnum = [[int(x[0]) - self._io_origin, float(x[1])] for x in inpstr]
                    self.setsupports(sorted(inpnum))

                if line.upper() == "SPECDOF":
                    line = textfile.readline().strip()
                    self.specdof_inputstring = line
                    inpstr = []
                    inpnum = []
                    inpstr = [x.split(',') for x in line.split('|')]
                    if len(inpstr[0]) == 1:
                        inpstr = [x.split(';') for x in line.split('|')]
                    if [''] in inpstr:
                        inpstr.remove([''])
                    inpnum = [[int(x[0]) - self._io_origin, str(x[1]).upper()] for x in inpstr]
                    self.setspecdofs(sorted(inpnum))

    def plot(self, showorig, showresult, showsupports, showforces, \
             showreactions, scaledisplacement, scaleforce, scalez, saveplot):
        """
        Plot function of the Truss class
        This method calls the more general plotstructure() method.
        """
        _showvalues = 1     # Show values of forces

        if self._postprocessed == 0:
            print 'Postprocess is needed before plotting structure!'
        else:
            if scaledisplacement == 0:
                scaledisplacement = 1.0           # Scale drwaing of displacements
            if scaleforce == 0:
                scaleforce = 1.0                  # Scale force sign
            if scalez == 0:
                scalez = 0.3                      # Scale z-axis

            plotstructure(self, showorig, showresult, showsupports, showforces, showreactions, \
                 scaledisplacement, scaleforce, scalez, _showvalues, saveplot)

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
                print "This node already exists. Input is ignored."
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
        self.known_f_a = []                   # <XXX> These won't be good here during updating the structure
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

        # <XXX> MODSUPPORT?????????????????????????????????????????????????????????????????????

    def setspecdofs(self, specdofs):
        """
        Set special nodal DOFs
        """
        self.analysis = {}

        for node, dofdir in specdofs:
            if 'X' in dofdir:
                self.analysis[str(node + self._io_origin)+'X'] = node*3+0
                self.keypoint.append(node*3+0)
            if 'Y' in dofdir:
                self.analysis[str(node + self._io_origin)+'Y'] = node*3+1
                self.keypoint.append(node*3+1)
            if self.dof == 3:
                if 'Z' in dofdir:
                    self.analysis[str(node + self._io_origin)+'Z'] = node*3+2
                    self.keypoint.append(node*3+2)

        self.keypnum = len(self.analysis)


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
        #if _UPDATING == 1:
        #    self.calcmodstiffnesses()


    def calcmodstiffnesses___(self):         #<XXX> Delete function
        """
        Unity-modificated stiffness matrixes
        """
        self.mod_stiffnesses = []

        for loopindex in range(self.elenum):
            _mod_stiffnesses_temp = [[0.]*(len(self.nodalcoord)*3)]*(len(self.nodalcoord)*3)

            for i in range(self.elenum):
                if i == loopindex:
                    _mod_norm_stiff = self._norm_stiff[i] *0.99 #self.el_mod[i]/ele_length[i]
                else:
                    _mod_norm_stiff = self._norm_stiff[i] #self.el_mod[i]/ele_length[i]

                _mod_loc_stiff = [[y*self.area[i]*_mod_norm_stiff for y in x] for x in self._s_loc[i]]

                ele_dof_vec = self.eledof[i]

                stiffincrement = [0.]*(len(self.nodalcoord)*3)

                for j in range(3*2):
                    for k in range(3*2):
                        stiffincrement[ele_dof_vec[k]] = _mod_loc_stiff[j][k]
                    _mod_stiffnesses_temp[ele_dof_vec[j]] = [x + y for x, y in zip(_mod_stiffnesses_temp[ele_dof_vec[j]], stiffincrement)]

            self.mod_stiffnesses.append(_mod_stiffnesses_temp)
        self._mod_stiffisfresh = 1

    def calcmodstiffness(self, index, change):
        """
        Convergency step in stiffness matrix modification
        """
        #self.mod_stiffnesses = []

        #for loopindex in range(self.elenum):
        _mod_stiffnesses_temp = [[0.]*(len(self.nodalcoord)*3)]*(len(self.nodalcoord)*3)

        for i in range(self.elenum):
            if i == index:
                _mod_norm_stiff = self._norm_stiff[i] * (1 + change) #self.el_mod[i]/ele_length[i]
            else:
                _mod_norm_stiff = self._norm_stiff[i] #self.el_mod[i]/ele_length[i]

            _mod_loc_stiff = [[y*self.area[i]*_mod_norm_stiff for y in x] for x in self._s_loc[i]]

            ele_dof_vec = self.eledof[i]

            stiffincrement = [0.]*(len(self.nodalcoord)*3)

            for j in range(3*2):
                for k in range(3*2):
                    stiffincrement[ele_dof_vec[k]] = _mod_loc_stiff[j][k]
                _mod_stiffnesses_temp[ele_dof_vec[j]] = [x + y for x, y in zip(_mod_stiffnesses_temp[ele_dof_vec[j]], stiffincrement)]

        self.mod_stiffnesses[index] = _mod_stiffnesses_temp
        #self._mod_stiffisfresh = 1

    def solve(self):
        """
        Main solver of the code
        """
        if self._stiffisfresh == 0:
            if _LOG:
                print 'Stiffness matrix is recalculated'
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
        if _SOLVER == 0:
            if _LOG:
                print 'Built-in solver'
            self.dis_new = mat_vec_mult(invert(self.stiff_new), self.force_new)
        else:
            if _LOG:
                print 'NumPy solver'
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
        if _UPDATING:
            for index in range(self.elenum):
                self.solvemodstruct(index, -0.05)
                self.evaluate(index)

    def solvemodstruct(self, index, change):
        """
        <XXX> Main solver of the code
        """
        if self._mod_stiffisfresh == 0:
            if _LOG:
                print 'Modified stiffness matrices are recalculated'

            self.mod_stiffnesses = [0.] * range(self.elenum)
            #for index in range(self.elenum):
            self.calcmodstiffness(index, change)

        self.mod_displacements = [[0.]*(self.nodenum*3)]*self.elenum

        #for loop in range(self.elenum):

        dis_new = [0.]*(self.nodenum*3-len(self.constraint))
        stiff_new = [[0.]*(self.nodenum*3-len(self.constraint))]*(self.nodenum*3-len(self.constraint))

        stiffincrement = [0.]*(self.nodenum*3-len(self.constraint))
        for i, kfai in enumerate(self.known_f_a):
            for j, kfaj in enumerate(self.known_f_a):
                stiffincrement[j] = self.mod_stiffnesses[index][kfai][kfaj]
            stiff_new[i] = [x + y for x, y in zip(stiff_new[i], stiffincrement)]

        # SOLVING THE MODIFIED STRUCTURE
        if _SOLVER == 0:
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
                effect_temp[j] = self.mod_displacements[modnum][dofnum]
                self.effect[modnum] = [x for x in effect_temp]
                self.toteffect[j] += abs(self.effect[modnum][j])

        self.effectratio = deepcopy(self.effect)
        for i in range(self.elenum):
            for j in range(self.keypnum):
                if self.toteffect[j] > 0:
                    self.effectratio[i][j] = abs(self.effectratio[i][j]/self.toteffect[j])
                else:
                    self.effectratio[i][j] = 0

        print "   \'effectratio\' is not used yet"

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
        #print "Print nodenumber option must be added!" #<XXX>

        delta = []

        for loc, measured in measurement:
            try:
                dof = self.analysis[loc.upper()]
            except KeyError:
                print 'The given measurement location cannot be matched with the input data.'
                print 'The available nodes are: {\'NAMES\': mapping addresses}'
                print self.analysis
                SER.close()
                raise Exception('Watchpoint name error')

            delta.append(measured - num_displ[dof])

        return delta

    def calcstiffness_plate(self):
	"""
	self.L, self.F, and self.clamping defines the configuration of the simulation.
	Do not forget to set them right before running the code.
	"""
        self.E = 210000.            # N/mm^2
        self.I = 110.634
		#self.L = 665.              # mm
        self.L = 893.               # mm
        self.x = self.L
        self.analysis = {}
        self.analysis['0Z'] = 0
        #self.F = 4.905               # N, 500 g
        #self.F = 5.886               # N, 600 g
        self.F = 11.772               # N, 1200 g
        self.clamping =  115.      # Nm/rad
        #self.clamping =  5000000000000000000000000000000000      # Nm/rad

    def solve_plate(self):
        self.displacement = [0.]
        self.displacement[0] += -(self.F*self.L*self.x*self.x/(2*self.E*self.I) - \
                self.F*self.x*self.x*self.x/(6*self.E*self.I))
        self.displacement[0] += -(self.F*self.L/self.clamping)
        print "Deflection: " + str(self.displacement)

        self.evaulate_plate()

    def evaulate_plate(self):
        self.effect = [0.]*2

        self.effect[0] = -(self.F*self.L*self.x*self.x/(2*self.E*self.I*0.95) - \
                self.F*self.x*self.x*self.x/(6*self.E*self.I*0.95) + self.F*self.L/self.clamping)
        self.effect[1] = -(self.F*self.L*self.x*self.x/(2*self.E*self.I) - \
                self.F*self.x*self.x*self.x/(6*self.E*self.I) + self.F*self.L/(self.clamping*0.9))

        # self.effect describes that x% of weakening causes self.effect big change in the displacements.
        self.effect[0] = self.effect[0] - self.displacement[0]
        self.effect[1] = self.effect[1] - self.displacement[0]

        self.toteffect = [abs(self.effect[0]) + abs(self.effect[1])]

        print self.effect

    def optimize(self, delta):
        """
        Modell updating - core function
        """
        #print "Hejhó :) <XXX>"
        # print "optimize_plate is running!" <XXX>
        self.update = [0., 0.]
        ratio = [0., 0.]
        
        if not _SIMULATION:
            appendix = ''
        else:
            appendix = ' - SIMULATED'
        
        
        
        ratio[0] = abs(self.effect[0] / self.toteffect[0])
        ratio[1] = abs(self.effect[1] / self.toteffect[0])

        unit = abs(ratio[0]*self.effect[0]) + abs(ratio[1]*self.effect[1])

        scale = delta[0]/unit
        print "Delta: " + str(delta[0]) + " Unit: " + str(unit)

        self.update[0] = (1 - (math.copysign(1, self.effect[0]) *0.05 * ratio[0]) *scale) *self.I
        self.update[1] = (1 - (math.copysign(1, self.effect[1]) *0.10 * ratio[1]) *scale) *self.clamping

        print "Updated I: " + str(self.update[0]) + " Clamp.: " + str(self.update[1])
        with open(self.name + ' - UpdateResults'+ appendix +'.txt', 'a') as outfile:
            outfile.write(str([scale, self.update[0], self.update[1]]) + "\n")  # <XXX> MUST BE UPDATED HERE!!!!

    def readarduino(self, base, saveinput):
        """
        Read data from Arduino
        """
        # Read data from Arduino

        maxdifference = 0.8          # Maximum input difference treshold in mm
        arduinovalues = []
        data = [0.]*len(arduino_mapping)
        newdata = False
        bigdifference = False
        readerror = False

        try:
            arduinoline = SER.readline()
            if len(arduinoline) > 0:
                arduinovalues = arduinoline.split(',')
                try:
                    if arduinovalues[0][len(arduinovalues)-1]=='.':
                        arduinovalues[0] = arduinovalues[0][:len(arduinovalues[0])-2]
                    else:
                        del arduinovalues[len(arduinovalues)-1]
                except IndexError:
                    print "Index Error... continuing"
                if len(arduinovalues) == len(arduino_mapping):
                    try:
                        for i in range(len(arduino_mapping)):
                            data[i] = float(arduinovalues[i]) - float(base[i][1])

                            if (abs(data[i] - self.processeddata[i]) > maxdifference):
                                bigdifference = True
                            if abs(float(arduinovalues[i])) < 2.0:
                                readerror = True

                        self.processeddata = data
                        newdata = True

                    except ValueError:
                        print('Value error... continuing')
                        SER.flushInput()
                        time.sleep(0.5)
                    except Exception:
                        print('Type error: ' + str(arduinovalues) + "... continuing")
                        SER.flushInput()
                        time.sleep(0.5)

            SER.flushInput()

        except serial.SerialTimeoutException:
            print('Data could not be read... continuing')
            SER.flushInput()
            time.sleep(0.5)

        #if bigdifference:
        #    print "Big difference: "
        #    print "Data :" + str(data)
        #if readerror:
        #    print "Reading error: "
        #    print data

        if newdata and not bigdifference and not readerror:
            measurement =  zip(arduino_mapping, data)

            saveinput.write(str(data) +', '+ str(time.time()) + "\n")
            # Calculate differences
            delta = self.difference(self.displacement, measurement)

            print delta

            # Optimize structure
            #print '*********** HJACKED *************'          # Debugging tool
            #self.displacement[0] = -1.
            #self.optimize([-1])
            self.optimize(delta)

        newdata = False
        bigdifference = False
        readerror = False
        time.sleep(0.5)             #<XXX> To be deleted!!!
        
        

    def simulatearduino(self, simulateinput):
        """
        Simulate data, based on previous measurement
        """
        arduinovalues = []
        arduinoline = ' '
        data = [0.]*len(arduino_mapping)
        
        try:
            arduinoline = simulateinput.readline()
            
            if arduinoline[0] == '[':
                arduinoline = arduinoline[1:len(arduinoline)-2]
                arduinovalues = arduinoline.split(',')
    
                try:
                    for i in range(len(arduino_mapping)):
                        data[i] = float(arduinovalues[i])
                    self.processeddata = data
                except Exception:
                    print('Type error: ' + str(arduinovalues) + "... continuing")
    
                measurement =  zip(arduino_mapping, data)
                
                # Calculate differences
                delta = self.difference(self.displacement, measurement)
    
                print delta
    
                # Optimize structure
                self.optimize(delta) 
                
        except IndexError:
            pass
        except Exception:
            pass

    def updatemodel(self):        
        self.processeddata = [0.]*len(arduino_mapping)
        
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

            for i in range(1000):
                if not _SIMULATION:
                    self.readarduino(base, inputfile)
                else:
                    self.simulatearduino(inputfile)

    def calibrate(self):

        c = '0'
        r = '0'
        accept = '0'
        arduinovalues = []
        print "Before starting the model updating, the measuring tools must be calibrated."
        print "The calibration should be done in load-free state."
        while (c not in ['Y','N']):
            c = raw_input('Can we start the calibration? (y/n) ').upper()
        if c == 'N':
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
                    if len(arduinovalues) == len(arduino_mapping):
                        measurement =  zip(arduino_mapping, arduinovalues)
                        print "Calibration result:"
                        print measurement
                        while (accept not in ['Y','N']):
                            accept = raw_input('Ok? Can we start the main part? Put on the loads! (y/n) ').upper()
                            if accept == 'N':
                                r = 'Y'
                    else:
                        print "Data error. Calibartion is restarting."
                        print "Arduino values:" + str(arduinovalues)
                        r = 'Y'
                else:
                    print 'The calibration cannot be done: no data'
                    while (r not in ['Y','N']):
                        r = raw_input('Do you want to restart calibration? (y/n) ').upper()
            except Exception:
                print 'The calibration cannot be done: exception was raised'
                while (r not in ['Y','N']):
                    r = raw_input('Do you want to restart calibration? (y/n) ').upper()

        if r == 'Y':
            print "Restarting calibration"
            SER.flushInput()
            self.calibrate()
        elif r == 'N':
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

##################################
#   BEGINNING OF THE MAIN PART   #
##################################

PARTTIME = logtime(TIC, "Initialization")

#  Define new truss
TRUSS = Truss('Debug')
TRUSS.name = raw_input('Test name: ')

# Read input file
#TRUSS.read('lab_01.txt')

PARTTIME = logtime(PARTTIME, "Setting up structure")

# Calculate stiffness-matrix
#TRUSS.calcstiffness()
TRUSS.calcstiffness_plate()

PARTTIME = logtime(PARTTIME, "Calculating Stiffness Matrix")

#Solve structure
#TRUSS.solve()

TRUSS.solve_plate()

PARTTIME = logtime(PARTTIME, "Solving")

if _UPDATING:
    TRUSS.updatemodel()

PARTTIME = logtime(PARTTIME, "Updating numerical model")

if _GRAPHICS:
    # Plot settings:
    # O: Original D: Deformed S: Supports F: Forces R: Reactions
    # ScD: Scale displacments (Z-axis) (def:1.0) ScF: Scale forces (def:1.0)
    # ScS: Scale Support signs (Z-axis) (def:1.0)
    # Save: Save plot to file

    #     plot(O, D, S, F, R, ScD, ScF, ScS, Save)
    #TRUSS.plot(1, 0, 1, 1, 0, 1.0, 0.0, 0.0, True)
    #TRUSS.plot(1, 1, 1, 0, 0, 1.0, 0.0, 0.0, True)
    #TRUSS.plot(0, 1, 1, 1, 1, 2.0, 0.0, 0.0, True)
    pass

PARTTIME = logtime(PARTTIME, "Plotting")

# Write results to file
#TRUSS.writeresults(TRUSS.name + ' - Results.txt')

PARTTIME = logtime(PARTTIME, "Writing results to the output file")

if _ARDUINO:
    # Closing Arduino port
    SER.close()

if _LOG:
    TAC = time.time()
    TOTALTIME = TAC-TIC
    print 'Total time: ' + str("{:10.3f}".format(TOTALTIME))
