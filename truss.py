# -*- coding: utf-8 -*-
"""
Created on April 30 18:49:40 2018

3D truss model updater program created by Máté Szedlák.
Copyright MIT, Máté Szedlák 2016-2018.
"""
import os
import itertools
import math
import time
import datetime
import numpy as np
try:
    import serial
except ImportError:
    print("You tried to import \'serial\' without installing \'pySerial\'.")
    print("Please first install pySerial: http://playground.arduino.cc/Interfacing/Python")
    raise Exception('pyserial package not found')
from copy import deepcopy
try:
    from truss_extras import TrussFramework
    from graphics import Arrow3D
    from config import Configuration
    from extra_math import invert as invert
    from extra_math import mat_vec_mult as mat_vec_mult
    from extra_math import swap_col as swap_col
except ImportError:
    print("Input data is missing")
    print("Please check the truss_extras.py, config.py, graphics.py and extra_math.py files.")
    raise Exception('Requirement is missing')


"""
 COMPATIBILITY MODES:
     0: User defined
     1: DEPRECATED
     2: Android
     3: Most information (with numpy)
     4: Maximum compatibility
"""

_COMPATIBLE_MODE = 0
_SIMULATION = 1                     # Simulating measurements based on input file

Conf = Configuration(_COMPATIBLE_MODE, _SIMULATION)


if _COMPATIBLE_MODE == 0*0:
    #** User defined ***
    # Modify as needed #
    Conf.mode_name = "User defined"
    Conf.log = 1                 # Logging time
    Conf.graphics = 1            # Graphical features
    Conf.solver = 1              # 0: Basic solver, 1: NumPy solver
    Conf.OSlib = 1  # Basic OS file features (e.g. file size)
    Conf.updating = 0            # Model Updating: On/ Off
    Conf.arduino = 0             # Arduino input: On/Off
    Conf.debug = 0               # Debugging mode
    Conf.realistic_simulation = 0  # Wait as long as it was originally. Only valid with _SIMULATION = 1

if Conf.OSlib:
    try:
        if not os.path.exists("./Structures/"):
            os.makedirs("./Structures")
    except FileNotFoundError:
        print("Error: Please manually create the 'Structures' folder"
              "in the root of the project")


if Conf.simulation or not Conf.updating:
    _ARDUINO = 0

if _COMPATIBLE_MODE == 2:
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

if Conf.log:
    TIC = time.time()
    print('------------------------------------')
    print('Truss calculation program')
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


"""
if _ARDUINO or _SIMULATION:
    try:
        mapping_file = 'arduino_mapping.txt'
        with open(mapping_file, "r") as textfile:
            line = textfile.readline().strip()
            arduino_mapping = line.upper().split(',')

    except IOError:
        raise Exception('File not found: ' + mapping_file)
"""


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
    sum_of_errors = 0
    for delta_element in delta:
        sum_of_errors += delta_element**2
        sum_of_errors = math.sqrt(sum_of_errors)

    return sum_of_errors


class Truss(TrussFramework):
    """
    General structure class
    """

    def __checkcoordinates(self, ignorable):
        """
        Checking coordinates for repeating elements.

            ignorable: [True | False] If the warning is ignorable, only message apperas and the input becomes neglected.
                If the warning is not ignorable, then exceptions will be raised.
            return: [0 | 1] 1 f no error found, otherwise 0.
        """
        if len(self.nodal_coord) != len(list(k for k, _ in itertools.groupby(sorted(self.nodal_coord)))):
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
        self._stiff_is_fresh = 0
        self._mod_stiff_is_fresh = 0
        self._post_processed = 0

    def setelements(self, node):
        """
        Setting elements (nodal connections) in bulk mode
        """
        self.node = node
        self.node_num = len(set(list(itertools.chain.from_iterable(sorted(self.node)))))
        self.element_num = len(self.node)
        self._stiff_is_fresh = 0
        self._mod_stiff_is_fresh = 0
        self._post_processed = 0

        # Creating mapping tool for elements
        for node in self.node:
            self.element_DOF.append([node[0]*3, node[0]*3+1, node[0]*3+2,
                                node[1]*3, node[1]*3+1, node[1]*3+2])

        # Initialazing matrix for all matrices
        self._init_displacement = [0.]*(3*self.node_num)
        self.force = [0.]*(3*self.node_num)
        self.stiffness = [0.]*(3*self.node_num)
        self.known_f_a = []
        self.known_f_not_zero = []

    def setcoordinates(self, coordinates):
        """
        Setting coordinates in bulk mode
        """
        self.nodal_coord = coordinates
        self._stiff_is_fresh = 0
        self._mod_stiff_is_fresh = 0
        if self.node_num > len(self.nodal_coord):
            raise Exception('More coordinates are needed')
        elif not self.node:
            raise Exception('Nodes must be set before defining elements')
        self.__checkcoordinates(False)

    def modcoordinate(self, node, coordinate):
        """
        Modify coordinate
        """
        if self.__checkcoordinates(True):
            self.nodal_coord[node] = coordinate
        self._stiff_is_fresh = 0
        self._mod_stiff_is_fresh = 0

    def setcrosssections(self, area):
        """
        Setting cross-sections in bulk mode
        """
        self.area = area
        self._stiff_is_fresh = 0
        self._mod_stiff_is_fresh = 0
        self._post_processed = 0

    def modcrosssection(self, element, area):
        """
        Modifying cross-sections by elements
        """
        self.area[element] = area
        self._stiff_is_fresh = 0
        self._mod_stiff_is_fresh = 0
        self._post_processed = 0

    def setmaterials(self, el_mod):
        """
        Setting material data in bulk mode
        """
        self.elastic_modulo = el_mod
        self._stiff_is_fresh = 0
        self._mod_stiff_is_fresh = 0
        self._post_processed = 0

    def modmaterial(self, element, el_mod):
        """
        Modifying material data by elements
        """
        self.elastic_modulo[element] = el_mod
        self._stiff_is_fresh = 0
        self._mod_stiff_is_fresh = 0
        self._post_processed = 0

    def setforces(self, forces):
        """
        Set forces
        """
        for fdof, force in forces:
            if self.dof == 3:
                self.force[fdof] = force
            elif self.dof == 2:
                self.force[fdof + (fdof//2)] = force
        self._post_processed = 0

    def modforce(self, element, force):
        """
        Modifying forces by each
        """
        self.force[element] = force
        self._post_processed = 0

    def setsupports(self, constraints):
        """
        Set supports
        """
        for cdof, constraint in constraints:
            if self.dof == 3:
                self.constraint.append([cdof, constraint])
            elif self.dof == 2:
                self.constraint.append([cdof + (cdof // 2), constraint])
        self._stiff_is_fresh = 0
        self._mod_stiff_is_fresh = 0
        self._post_processed = 0

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
                    print("Z-direction is not allowed in 2D structures. "
                          "Please check the 'MEASUREMENTS' section in the input file.")
                    raise Exception

        self.keypoint_num = len(self.analysis)
        if self.keypoint_num == 0 and Conf.updating:
            print("There is no valid measured DOF. Please check the \'MEASUREMENTS\' section in the input file.")
            raise Exception

    def calcstiffness(self):
        """
        Stiffness matrix compilation
        """
        self._post_processed = 0

        if self.dof == 2:
            for zdof in range(self.node_num):
                self.constraint.append([int(zdof*3+2), 0.])
        self.constraint = list(k for k, _ in itertools.groupby(sorted(self.constraint)))

        # Setting known forces
        for dofloc in range(3*self.node_num):
            self.known_f_a.append(dofloc)
            if self.force[dofloc] != 0:
                self.known_f_not_zero.append(dofloc)

        self.known_dis_a = []
        for constr in self.constraint:
            self._init_displacement[constr[0]] = constr[1]
            self.known_dis_a.append(constr[0])
            try:
                self.known_f_a.remove(constr[0])
                self.known_f_not_zero.remove(constr[0])
            except ValueError:
                pass

        ele_length = [0.]*self.element_num
        self._norm_stiff = [0.]*self.element_num
        self._cx = [0.]*self.element_num
        self._cy = [0.]*self.element_num
        self._cz = [0.]*self.element_num
        self._s_loc = [0.]*self.element_num
        self._loc_stiff = [0.]*self.element_num
        self.stress = [0.]*self.element_num

        self.stiffness = [[0.]*(len(self.nodal_coord)*3)]*(len(self.nodal_coord)*3)

        for i in range(self.element_num):
            ele_length[i] = math.sqrt((self.nodal_coord[self.node[i][1]][0]-self.nodal_coord[self.node[i][0]][0])**2 +
                                      (self.nodal_coord[self.node[i][1]][1]-self.nodal_coord[self.node[i][0]][1])**2 +
                                      (self.nodal_coord[self.node[i][1]][2]-self.nodal_coord[self.node[i][0]][2])**2)

            self._cx[i] = (self.nodal_coord[self.node[i][1]][0]-self.nodal_coord[self.node[i][0]][0])/ele_length[i]
            self._cy[i] = (self.nodal_coord[self.node[i][1]][1]-self.nodal_coord[self.node[i][0]][1])/ele_length[i]
            self._cz[i] = (self.nodal_coord[self.node[i][1]][2]-self.nodal_coord[self.node[i][0]][2])/ele_length[i]
            self._norm_stiff[i] = self.elastic_modulo[i]/ele_length[i]
            self._s_loc[i] = [[self._cx[i]**2, self._cx[i]*self._cy[i], self._cx[i]*self._cz[i], -self._cx[i]**2, -self._cx[i]*self._cy[i], -self._cx[i]*self._cz[i]],
                              [self._cx[i]*self._cy[i], self._cy[i]**2, self._cy[i]*self._cz[i], -self._cx[i]*self._cy[i], -self._cy[i]**2, -self._cy[i]*self._cz[i]],
                              [self._cx[i]*self._cz[i], self._cy[i]*self._cz[i], self._cz[i]**2, -self._cx[i]*self._cz[i], -self._cy[i]*self._cz[i], -self._cz[i]**2],
                              [-self._cx[i]**2, -self._cx[i]*self._cy[i], -self._cx[i]*self._cz[i], self._cx[i]**2, self._cx[i]*self._cy[i], self._cx[i]*self._cz[i]],
                              [-self._cx[i]*self._cy[i], -self._cy[i]**2, -self._cy[i]*self._cz[i], self._cx[i]*self._cy[i], self._cy[i]**2, self._cy[i]*self._cz[i]],
                              [-self._cx[i]*self._cz[i], -self._cy[i]*self._cz[i], -self._cz[i]**2, self._cx[i]*self._cz[i], self._cy[i]*self._cz[i], self._cz[i]**2]]
            self._loc_stiff[i] = [[y * self.area[i] * self._norm_stiff[i] for y in x] for x in self._s_loc[i]]
            ele_dof_vec = self.element_DOF[i]

            stiffincrement = [0.]*(len(self.nodal_coord)*3)

            for j in range(3*2):
                for k in range(3*2):
                    stiffincrement[ele_dof_vec[k]] = self._loc_stiff[i][j][k]
                self.stiffness[ele_dof_vec[j]] = [x + y for x, y in zip(self.stiffness[ele_dof_vec[j]], stiffincrement)]
        self._stiff_is_fresh = 1

    def calcmodstiffness(self, index, magnitude):
        """
        Convergency step in stiffness matrix modification
        """
        if not self.mod_stiffnesses:
            self.mod_stiffnesses = [0.]*(self.element_num+1)

        # for loopindex in range(self.element_num):
        _mod_stiffnesses_temp = [[0.]*(len(self.nodal_coord)*3)]*(len(self.nodal_coord)*3)

        for i in range(self.element_num):
            if i == index:
                _mod_norm_stiff = self._norm_stiff[i] * (1.0 + self.modifications[i] + magnitude)  # self.elastic_modulo[i]/ele_length[i]
            else:
                _mod_norm_stiff = self._norm_stiff[i] * (1.0 + self.modifications[i])  # self.elastic_modulo[i]/ele_length[i]

            _mod_loc_stiff = [[y*self.area[i]*_mod_norm_stiff for y in x] for x in self._s_loc[i]]

            ele_dof_vec = self.element_DOF[i]

            stiffincrement = [0.]*(len(self.nodal_coord)*3)

            for j in range(3*2):
                for k in range(3*2):
                    stiffincrement[ele_dof_vec[k]] = _mod_loc_stiff[j][k]
                _mod_stiffnesses_temp[ele_dof_vec[j]] = [x + y for x, y in
                                                         zip(_mod_stiffnesses_temp[ele_dof_vec[j]], stiffincrement)]

        self.mod_stiffnesses[index] = _mod_stiffnesses_temp

    def solve(self):
        """
        Main solver of the code
        """
        if self._stiff_is_fresh == 0:
            if Conf.log:
                print('Stiffness matrix is recalculated')
            self.calcstiffness()

        self.dis_new = [0.]*(self.node_num*3-len(self.constraint))
        self.force_new = [0.]*(self.node_num*3-len(self.constraint))
        self.stiff_new = [[0.]*(self.node_num*3-len(self.constraint))]*(self.node_num*3-len(self.constraint))

        # known force array
        for i, known_f_a in enumerate(self.known_f_a):
            self.force_new[i] = self.force[known_f_a]

        stiffincrement = [0.]*(self.node_num*3-len(self.constraint))
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

        self.displacement = deepcopy(self._init_displacement)

        for i, known_f_a in enumerate(self.known_f_a):
            self.displacement[known_f_a] = self.dis_new[i]

        # Deformed shape
        self.nodal_coord_def = []
        for i in range(self.node_num):
            self.nodal_coord_def.append([self.nodal_coord[i][0] + self.displacement[i*3+0],
                                        self.nodal_coord[i][1] + self.displacement[i*3+1], self.nodal_coord[i][2] + self.displacement[i*3+2]])

        # Postrpocesses
        self.postprocess()

        self.mod_displacements = [0.]*(self.element_num+1)

    def solvemodstruct(self, index):
        """
        Solver for the modified structures. 'Index' shows the actual modification number.
        """

        self.mod_displacements[index] = [0.]*(self.node_num*3)

        dis_new = [0.]*(self.node_num*3-len(self.constraint))
        stiff_new = [[0.]*(self.node_num*3-len(self.constraint))]*(self.node_num*3-len(self.constraint))

        stiffincrement = [0.]*(self.node_num*3-len(self.constraint))
        for i, kfai in enumerate(self.known_f_a):
            for j, kfaj in enumerate(self.known_f_a):
                stiffincrement[j] = self.mod_stiffnesses[index][kfai][kfaj]
            stiff_new[i] = [x + y for x, y in zip(stiff_new[i], stiffincrement)]

        # SOLVING THE MODIFIED STRUCTURE
        if Conf.solver == 0:
            dis_new = mat_vec_mult(invert(stiff_new), self.force_new)
        else:
            dis_new = np.linalg.solve(np.array(stiff_new), np.array(self.force_new))

        mod_displacement_temp = deepcopy(self._init_displacement)

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
        self.effect = [[0.]*(self.keypoint_num + 2)]*self.element_num
        self.total_effect = [0.]*self.keypoint_num
        self.sorted_effect = [[[0.]*(self.keypoint_num + 2)]*self.element_num]*self.keypoint_num

        effect_temp = [0.]*(self.keypoint_num + 2)

        for modnum in range(self.element_num):
            effect_temp[self.keypoint_num] = int(modnum)
            for j, dofnum in enumerate(self.keypoint):
                try:
                    effect_temp[j] = self.mod_displacements[modnum][dofnum]
                    self.effect[modnum] = [x for x in effect_temp]
                    self.total_effect[j] += abs(self.effect[modnum][j])
                except IndexError:
                    print("Maybe the mapping data is invalid.")
                    print("Please check the \'arduino_mapping.txt\' input whether the given DOFs are correct or not.")
                    raise IndexError

        self.effect_ratio = deepcopy(self.effect)
        for i in range(self.element_num):
            for j in range(self.keypoint_num):
                if self.total_effect[j] > 0:
                    self.effect_ratio[i][j] = abs(self.effect_ratio[i][j]/self.total_effect[j])
                else:
                    self.effect_ratio[i][j] = 0

        # print("   \'effectratio\' is not used yet")

        # Sort by effectiveness
        for i in range(self.keypoint_num):
            self.sorted_effect[i] = deepcopy(self.effect)

            # Check sign of the effect
            for ktemp in range(self.element_num):
                if self.sorted_effect[i][ktemp][i] < 0:
                    for jtemp in range(self.keypoint_num):
                        self.sorted_effect[i][ktemp][jtemp] = abs(self.sorted_effect[i][ktemp][jtemp])
                        self.sorted_effect[i][ktemp][self.keypoint_num + 1] = -1
                else:
                    self.sorted_effect[i][ktemp][self.keypoint_num + 1] = +1

            for j in range(self.keypoint_num):
                if i != j and j != 0:
                    self.sorted_effect[i] = swap_col(sorted(swap_col(self.sorted_effect[i], 0, j), reverse=True), 0, j)
            if i != 0:
                self.sorted_effect[i] = swap_col(sorted(swap_col(self.sorted_effect[i], 0, i), reverse=True), 0, i)
            else:
                self.sorted_effect[i] = sorted(self.sorted_effect[i], reverse=True)

    def difference(self, num_displ, measurement):
        """
        Calculate the difference between the Numerical solution and Real-life measurement.
        The Real-life measurement should be given the following format:

            MEASUREMENT: [[13X, -2.154], [16Y, 5.256], ...]
        """
        # print(nodenumber option should be added! <XXX>)

        delta = []

        for loc, measured in measurement:
            try:
                dof = self.analysis[loc.upper()]
            except KeyError:
                print('The given measurement location cannot be matched with the input data.')
                print('The available nodes are: {\'NAMES\': mapping addresses}')
                print(self.analysis)
                self.serial.close()
                raise Exception('Watchpoint name error')

            delta.append(measured - num_displ[dof])

        return delta

    def optimize(self, delta):
        """
        Modell updating - core function
        """
        # modnum = min(10, self.element_num)
        modnum = self.element_num
        self.modifications = [0.0]*self.element_num

        if not _SIMULATION:
            appendix = ""
        else:
            appendix = " - SIMULATED"

        newdelta = delta
        j = 0

        print("-----")
        print("Step: 0/" + str(self.iteration_limit))

        # Optimization loop
        while error(newdelta) > self.error_limit and j <= self.iteration_limit and (self.capable() or j <= 1):
            j += 1

            print("Error: " + str(error(newdelta)))
            print("-----")
            print("Step: " + str(j) + "/" + str(self.iteration_limit))

            ratio = [0.]*modnum
            unit = 0

            prevmodifications = self.modifications

            for index in range(self.element_num):
                self.modifications[index] = min(abs(self.modifications[index] - self.unit_modification),
                                                self.modification_limit) * \
                                            math.copysign(1, self.modifications[index] - self.unit_modification)
                self.calcmodstiffness(index, self.modifications[index])
                self.solvemodstruct(index)
            self.evaluate()

            self.calcmodstiffness(self.element_num, 0)
            self.solvemodstruct(self.element_num)

            newdelta = self.difference(self.mod_displacements[self.element_num], self.measurement)

            for i, effect in enumerate(self.total_effect):
                if effect == 0.0:
                    print("None of the variables has effect on " + str(self.arduino_mapping[i]))
                    print("Model updating has no solution.")
                    raise Exception

            for i in range(self.element_num):
                modificationnumber = self.sorted_effect[0][i][1]
                ratio[modificationnumber] = abs(self.sorted_effect[0][i][0] / self.total_effect[0]) * \
                                            math.copysign(1, self.sorted_effect[0][i][2])
                unit += abs(ratio[modificationnumber]*self.sorted_effect[0][i][0])

            scale = newdelta[0]/unit
            for i in range(self.element_num):
                modificationnumber = self.sorted_effect[0][i][1]
                self.modifications[modificationnumber] = min(abs(prevmodifications[modificationnumber] -
                                                                 self.unit_modification*ratio[modificationnumber]),
                                                             self.modification_limit) \
                                                         * math.copysign(1, prevmodifications[modificationnumber] -
                                                                         self.unit_modification*ratio[modificationnumber])
                # the last part is already the sign itself without the sign function

            print("Ratio: " + str(scale))

        print("Final error: " + str(error(newdelta)))

        if not self.capable() and j > 1:
            print("Optimization could not be finished successfully.")
            print("The remaining error is: " + str(error(newdelta)))

        with open("./Structures/" + self.name + ' - UpdateResults' + appendix + '.txt', 'a') as outfile:
            if j > 1:
                if j <= self.iteration_limit and self.capable():
                    self.number_of_updates[0] += 1
                    outfile.write("Update state: SUCCESSFUL\n")
                if not j <= self.iteration_limit:
                    self.number_of_updates[1] += 1
                    outfile.write("Update state: Run out of iteration limit\n")
                if not self.capable() and j > 1:
                    self.number_of_updates[2] += 1
                    outfile.write("Update state: No more possible modification\n")
            else:
                outfile.write("Update state: Optimization was skipped\n")
            outfile.write("Requiered iterations: " + str(j) + "\n")
            outfile.write("Measurement: " + str(self.measurement) + "\n")
            outfile.write("Original delta: " + str(delta) + "\n")
            outfile.write("New delta: " + str(newdelta) + " (limit: " + str(self.error_limit) + ")\n")
            outfile.write("Final error: " + str(error(newdelta)) + "\n")
            outfile.write("Modifications [%]: \n")
            outfile.write(str(self.modifications) + "\n")
            outfile.write("Original displacements: \n")
            outfile.write(str(self.displacement) + "\n")
            if j > 1:
                outfile.write("New displacements: \n")
                outfile.write(str(self.mod_displacements[self.element_num]) + "\n")
            outfile.write("----------------------\n")

    def capable(self):
        """
        Function telling whether there are more options to modify
        """
        for variable in self.modifications:
            if 0.01 < abs(variable) <= 0.95*self.modification_limit:
                capable = True
        return capable

    def seterrorlimit(self, errorlimit):
        """
        Setting general stop parameter for model updating
        """
        if errorlimit > 0.0:
            self.error_limit = errorlimit
        else:
            print("The error limit must be a positive number")
            raise Exception

    def setmodificationlimit(self, modificationlimit):
        """
        Setting modification limit for members (model updating)
        """
        if 0.0 < modificationlimit < 1.0:
            self.modification_limit = modificationlimit
        else:
            print("The modification limit must be higher than 0.0 and lower than 1.0")
            raise Exception

    def setunitmodification(self, unitmodification):
        """
        Setting modification step (model updating)
        """
        if 0.01 <= abs(unitmodification) < 0.5:
            self.unit_modification = unitmodification
        else:
            print("The absolut value of the unit modification must be minimum 0.01 and maximum 0.5")
            raise Exception

    def setiterationlimit(self, iterationlimit):
        """
        Setting maximum number of iterations (model updating)
        """
        if 1 < int(iterationlimit) <= math.pow(10, 4):
            self.iteration_limit = int(iterationlimit)
        else:
            print("The iterationlimit must be between 2 and 10.000")
            raise Exception

    def updatemodel(self):
        """
        General function to manage model updatin procedure.
        """
        self.processed_data = [0.]*len(self.arduino_mapping)

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

        with open("./Structures/" + self.name + ' - Input Data.txt', filemode) as inputfile:
            # Saving input data
            if not _SIMULATION:
                inputfile.write('Input data of \'' + self.name + '\':\n\n')
                inputfile.write('Start Time: ' +
                                str(datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')) + "\n")
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
                    if delta:
                        self.optimize(delta)
        print("Update statistics:")
        print("Totally updated models: " + str(TRUSS.numofupdates[0] + TRUSS.numofupdates[1]++ TRUSS.numofupdates[2]))
        print("  Successfully updated models: " + str(TRUSS.numofupdates[0]))
        print("  Updates with running out of possibilities: " + str(TRUSS.numofupdates[2]))
        print("  Updates did not finshed: " + str(TRUSS.numofupdates[1]))

    def postprocess(self):
        """
        Calculates reaction forces and stresses
        """
        self._reactions()
        self._stresses()
        self._post_processed = 1

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
        self.stress = [0.]*self.element_num
        for element in range(self.element_num):
            locstiff = [-self._cx[element], -self._cy[element], -self._cz[element],
                        self._cx[element], self._cy[element], self._cz[element]]
            for i in range(3*2):
                self.stress[element] += locstiff[i]*self.displacement[self.element_DOF[element][i]]
            self.stress[element] = self.stress[element]*self._norm_stiff[element]

        smax = max([abs(min(self.stress)), max(self.stress), 0.000000001])
        self.stress_color = [float(x)/float(smax) for x in self.stress]

    def postprocessed(self):
        """
        Tells if the structure's postprocess part is already calcuated
        """
        return self._post_processed

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
# TRUSS.read('lab_01.txt')

try:
    TRUSS.read(TRUSS.name + ".str")
except IOError:
    print("The following file could not be opened: " + "./Structures/" + TRUSS.name + ".str")
    print("Please make sure that the structural data is available for the program in the running directory.")
    raise IOError

# if _ARDUINO or _SIMULATION:                # deprecated
#    TRUSS.setspecdofs(arduino_mapping)

PARTTIME = logtime(PARTTIME, "Setting up structure")

# Calculate stiffness-matrix
TRUSS.calcstiffness()
# TRUSS.calcstiffness_plate()

PARTTIME = logtime(PARTTIME, "Calculating Stiffness Matrix")

# Solve structure
TRUSS.solve()
# TRUSS.solve_plate()

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
    # pass

PARTTIME = logtime(PARTTIME, "Plotting")

# Write results to file
TRUSS.writeresults("./Structures/" + TRUSS.name + ' - Results.txt')

PARTTIME = logtime(PARTTIME, "Writing results to the output file")

if Conf.arduino:
    # Closing Arduino port
    TRUSS.close_serial()

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
