# -*- coding: utf-8 -*-
"""
Created on April 30 18:49:40 2018

3D truss model updater program created by Máté Szedlák.
For more information please check the README
Copyright MIT, Máté Szedlák 2016-2018.
"""
import os
import argparse
import itertools
import math
import time
import datetime
import numpy as np
from copy import deepcopy
try:
    from truss_framework import TrussFramework
    from truss_graphics import Arrow3D
    from extra_math import invert as invert
    from extra_math import mat_vec_mult as multiply_matrix_vector
    from extra_math import swap_col as swap_columns
except ImportError:
    print("Input data is missing")
    print("Please check the truss_framework.py, truss_graphics.py and extra_math.py files.")
    raise Exception('External file is missing')


def error(delta):
    """
    Error function using least-square method
    :param delta: error vector
    :return: sum of errors using least-square method
    """
    sum_of_errors = 0
    for delta_element in delta:
        sum_of_errors += delta_element**2
        sum_of_errors = math.sqrt(sum_of_errors)

    return sum_of_errors


class Truss(TrussFramework):
    """
    Object for solving a truss
    """

    def _check_coordinates(self, ignorable):
        """
        Checking coordinates for repeated elements during input file reading.

        :param ignorable: [True|False] If the warning is ignorable, only message appears and the input is neglected.
                If the warning is not ignorable, exceptions will be raised.
        :return: [True|False] True if no error found, otherwise False.
        """
        if len(self.nodal_coord) != len(list(k for k, _ in itertools.groupby(sorted(self.nodal_coord)))):
            if ignorable == 0:
                raise Exception('Coordinate list has repeating items. Calculation is terminated')
            else:
                print("This node already exists. Input is ignored.")
            return False
        else:
            return True

    def set_model_DOF(self, DOF):
        """
        Setting problem's degree of freedom
        :param DOF: [2 | 3] Model's Degree Of Freedom
        :return: None
        """
        self.DOF = DOF
        if self.DOF not in [2, 3]:
            raise Exception('DOF must be 2 or 3.')
        # Set freshness flags after geometry modification
        self._stiff_is_fresh = 0
        self._mod_stiff_is_fresh = 0
        self._post_processed = 0

    def bulk_set_elements(self, nodal_connection_list):
        """
        Setting elements (nodal connections) in bulk mode

        :param nodal_connection_list: an array of arrays, where each sub-array is a pair of integers, namely [i, j].
        i, j are the ID's of the nodes and an element is i -> j.
        :return: None
        """
        # Set attribute
        self.nodal_connections = nodal_connection_list
        # Recalculate pending variables
        self.number_of_nodes = len(set(list(itertools.chain.from_iterable(sorted(self.nodal_connections)))))
        self.element_num = len(self.nodal_connections)
        # Set freshness flags after geometry modification
        self._stiff_is_fresh = 0
        self._mod_stiff_is_fresh = 0
        self._post_processed = 0

        # Creating mapping tool for elements
        for node in self.nodal_connections:
            self.element_DOF.append([node[0]*3, node[0]*3+1, node[0]*3+2, node[1]*3, node[1]*3+1, node[1]*3+2])

        # Initializing defaults for all matrices
        self._init_displacement = [0.]*(3*self.number_of_nodes)
        self.force = [0.]*(3*self.number_of_nodes)
        self.stiffness = [0.]*(3*self.number_of_nodes)
        self.known_f_a = []
        self.known_f_not_zero = []

    def bulk_set_coordinates(self, coordinate_list):
        """
        Setting coordinates in bulk mode

        :param coordinate_list: An array of coordinate arrays enlisting ALL the nodal coordinates
        :return: None
        """
        # Set attribute
        self.nodal_coord = coordinate_list
        # Set freshness flags after geometry modification
        self._stiff_is_fresh = 0
        self._mod_stiff_is_fresh = 0
        # Validity check
        if self.number_of_nodes > len(self.nodal_coord):
            raise Exception('More coordinates are needed')
        elif not self.nodal_connections:
            raise Exception('Nodes must be set before defining elements')
        self._check_coordinates(False)

    def modify_coordinate(self, node_ID, new_coordinate):
        """
        Modify one node's coordinates

        :param node_ID: the ID of the changeable node
        :param new_coordinate: new coordinates wrapped by an array
        :return: Bool - True uf the modification was applied, False otherwise.

        @example: self.modify_coordinate(5, [1.2, 3.6])  # 2D
        """
        if self._check_coordinates(True):
            self.nodal_coord[node_ID] = new_coordinate
            self._stiff_is_fresh = 0
            self._mod_stiff_is_fresh = 0
            return True
        else:
            return False

    def bulk_set_cross_sections(self, area_list):
        """
        Setting cross-sections in bulk mode

        :param area_list: an array enlisting the cross-sectional data
        :return: None
        """
        self.cross_sectional_area_list = area_list
        self._stiff_is_fresh = 0
        self._mod_stiff_is_fresh = 0
        self._post_processed = 0

    def modify_cross_section(self, element_ID, new_area):
        """
        Modifying cross-sections by elements

        :param element_ID: the ID of the changeable element
        :param new_area: the new cross-sectional area
        :return: None
        """
        self.cross_sectional_area_list[element_ID] = new_area
        self._stiff_is_fresh = 0
        self._mod_stiff_is_fresh = 0
        self._post_processed = 0

    def bulk_set_materials(self, el_mod):
        """
        Setting material data in bulk mode
        """
        self.elastic_modulo = el_mod
        self._stiff_is_fresh = 0
        self._mod_stiff_is_fresh = 0
        self._post_processed = 0

    def modify_material(self, element, el_mod):
        """
        Modifying material data by elements
        """
        self.elastic_modulo[element] = el_mod
        self._stiff_is_fresh = 0
        self._mod_stiff_is_fresh = 0
        self._post_processed = 0

    def bulk_set_forces(self, forces):
        """
        Set forces
        """
        for location, force in forces:
            if self.DOF == 3:
                self.force[location] = force
            elif self.DOF == 2:
                self.force[location + (location//2)] = force
        self._post_processed = 0

    def modify_force(self, element, force):
        """
        Modifying forces by each
        """
        self.force[element] = force
        self._post_processed = 0

    def bulk_set_supports(self, constraints):
        """
        Set supports
        """
        for location, constraint in constraints:
            if self.DOF == 3:
                self.constraint.append([location, constraint])
            elif self.DOF == 2:
                self.constraint.append([location + (location // 2), constraint])
        self._stiff_is_fresh = 0
        self._mod_stiff_is_fresh = 0
        self._post_processed = 0

    def bulk_set_measurement_points(self, measurement_points):
        """
        Set special nodal DOFs
        """
        self.analysis = {}
        print(measurement_points)

        for location in measurement_points:
            node = int(location[:len(location)-1])-self._io_origin
            if 'X' in location:
                self.analysis[location] = node*3+0
                self.keypoint.append(node*3+0)
            if 'Y' in location:
                self.analysis[location] = node*3+1
                self.keypoint.append(node*3+1)

            if 'Z' in location:
                if self.DOF == 3:
                    self.analysis[location] = node*3+2
                    self.keypoint.append(node*3+2)
                else:
                    print("Z-direction is not allowed in 2D structures. "
                          "Please check the 'MEASUREMENTS' section in the input file.")
                    raise Exception

        self.number_of_keypoints = len(self.analysis)
        if self.number_of_keypoints == 0 and self.configuration.updating:
            print("There is no valid measured DOF. Please check the \'MEASUREMENTS\' section in the input file.")
            raise Exception

    def calculate_stiffness_matrix(self):
        """
        Stiffness matrix compilation
        """
        self._post_processed = 0

        if self.DOF == 2:
            for zdof in range(self.number_of_nodes):
                self.constraint.append([int(zdof*3+2), 0.])
        self.constraint = list(k for k, _ in itertools.groupby(sorted(self.constraint)))

        # Setting known forces
        for dofloc in range(3*self.number_of_nodes):
            self.known_f_a.append(dofloc)
            if self.force[dofloc] != 0:
                self.known_f_not_zero.append(dofloc)

        self.known_displacement_a = []
        for constr in self.constraint:
            self._init_displacement[constr[0]] = constr[1]
            self.known_displacement_a.append(constr[0])
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
            ele_length[i] = math.sqrt((self.nodal_coord[self.nodal_connections[i][1]][0]-self.nodal_coord[self.nodal_connections[i][0]][0])**2 +
                                      (self.nodal_coord[self.nodal_connections[i][1]][1]-self.nodal_coord[self.nodal_connections[i][0]][1])**2 +
                                      (self.nodal_coord[self.nodal_connections[i][1]][2]-self.nodal_coord[self.nodal_connections[i][0]][2])**2)

            self._cx[i] = (self.nodal_coord[self.nodal_connections[i][1]][0]-self.nodal_coord[self.nodal_connections[i][0]][0])/ele_length[i]
            self._cy[i] = (self.nodal_coord[self.nodal_connections[i][1]][1]-self.nodal_coord[self.nodal_connections[i][0]][1])/ele_length[i]
            self._cz[i] = (self.nodal_coord[self.nodal_connections[i][1]][2]-self.nodal_coord[self.nodal_connections[i][0]][2])/ele_length[i]
            self._norm_stiff[i] = self.elastic_modulo[i]/ele_length[i]
            self._s_loc[i] = [[self._cx[i]**2, self._cx[i]*self._cy[i], self._cx[i]*self._cz[i], -self._cx[i]**2, -self._cx[i]*self._cy[i], -self._cx[i]*self._cz[i]],
                              [self._cx[i]*self._cy[i], self._cy[i]**2, self._cy[i]*self._cz[i], -self._cx[i]*self._cy[i], -self._cy[i]**2, -self._cy[i]*self._cz[i]],
                              [self._cx[i]*self._cz[i], self._cy[i]*self._cz[i], self._cz[i]**2, -self._cx[i]*self._cz[i], -self._cy[i]*self._cz[i], -self._cz[i]**2],
                              [-self._cx[i]**2, -self._cx[i]*self._cy[i], -self._cx[i]*self._cz[i], self._cx[i]**2, self._cx[i]*self._cy[i], self._cx[i]*self._cz[i]],
                              [-self._cx[i]*self._cy[i], -self._cy[i]**2, -self._cy[i]*self._cz[i], self._cx[i]*self._cy[i], self._cy[i]**2, self._cy[i]*self._cz[i]],
                              [-self._cx[i]*self._cz[i], -self._cy[i]*self._cz[i], -self._cz[i]**2, self._cx[i]*self._cz[i], self._cy[i]*self._cz[i], self._cz[i]**2]]
            self._loc_stiff[i] = [[y * self.cross_sectional_area_list[i] * self._norm_stiff[i] for y in x] for x in self._s_loc[i]]
            ele_dof_vec = self.element_DOF[i]

            stiffness_increment = [0.]*(len(self.nodal_coord)*3)

            for j in range(3*2):
                for k in range(3*2):
                    stiffness_increment[ele_dof_vec[k]] = self._loc_stiff[i][j][k]
                self.stiffness[ele_dof_vec[j]] = [x + y for x, y in zip(self.stiffness[ele_dof_vec[j]], stiffness_increment)]
        self._stiff_is_fresh = 1

    def calculate_modified_stiffness_matrix(self, index, magnitude):
        """
        Convergence step in stiffness matrix modification
        """
        if not self.mod_stiffnesses:
            self.mod_stiffnesses = [0.]*(self.element_num+1)

        # for loopindex in range(self.element_num):
        _mod_stiffnesses_temp = [[0.]*(len(self.nodal_coord)*3)]*(len(self.nodal_coord)*3)

        for i in range(self.element_num):
            if i == index:
                _mod_norm_stiff = self._norm_stiff[i] * (1.0 + self.modifications[i] + magnitude)  # E[i]/L[i]
            else:
                _mod_norm_stiff = self._norm_stiff[i] * (1.0 + self.modifications[i])  # E[i]/L[i]

            _mod_loc_stiff = [[y*self.cross_sectional_area_list[i]*_mod_norm_stiff for y in x] for x in self._s_loc[i]]

            ele_dof_vec = self.element_DOF[i]

            stiffness_increment = [0.]*(len(self.nodal_coord)*3)

            for j in range(3*2):
                for k in range(3*2):
                    stiffness_increment[ele_dof_vec[k]] = _mod_loc_stiff[j][k]
                _mod_stiffnesses_temp[ele_dof_vec[j]] = [x + y for x, y in
                                                         zip(_mod_stiffnesses_temp[ele_dof_vec[j]], stiffness_increment)]

        self.mod_stiffnesses[index] = _mod_stiffnesses_temp

    def solve(self):
        """
        Main solver of the code
        """
        if self._stiff_is_fresh == 0:
            if self.configuration.log:
                print('Stiffness matrix is recalculated')
            self.calculate_stiffness_matrix()

        self.dis_new = [0.]*(self.number_of_nodes*3-len(self.constraint))
        self.force_new = [0.]*(self.number_of_nodes*3-len(self.constraint))
        self.stiff_new = [[0.]*(self.number_of_nodes*3-len(self.constraint))]*(self.number_of_nodes*3-len(self.constraint))

        # known force array
        for i, known_f_a in enumerate(self.known_f_a):
            self.force_new[i] = self.force[known_f_a]

        stiffness_increment = [0.]*(self.number_of_nodes*3-len(self.constraint))
        for i, kfai in enumerate(self.known_f_a):
            for j, kfaj in enumerate(self.known_f_a):
                stiffness_increment[j] = self.stiffness[kfai][kfaj]
            self.stiff_new[i] = [x + y for x, y in zip(self.stiff_new[i], stiffness_increment)]

        # SOLVING THE STRUCTURE
        if self.configuration.solver == 0:
            if self.configuration.log:
                print('Built-in solver')
            self.dis_new = multiply_matrix_vector(invert(self.stiff_new), self.force_new)
        else:
            if self.configuration.log:
                print('NumPy solver')
            self.dis_new = np.linalg.solve(np.array(self.stiff_new), np.array(self.force_new))

        self.displacement = deepcopy(self._init_displacement)

        for i, known_f_a in enumerate(self.known_f_a):
            self.displacement[known_f_a] = self.dis_new[i]

        # Deformed shape
        self.nodal_coord_def = []
        for i in range(self.number_of_nodes):
            self.nodal_coord_def.append([self.nodal_coord[i][0] + self.displacement[i*3+0],
                                        self.nodal_coord[i][1] + self.displacement[i*3+1], self.nodal_coord[i][2] + self.displacement[i*3+2]])

        # Postrpocesses
        self.post_process()

        self.modified_displacements = [0.]*(self.element_num+1)

    def solve_modified_structure(self, index):
        """
        Solver for the modified structures. 'Index' shows the actual modification number.
        """

        self.modified_displacements[index] = [0.]*(self.number_of_nodes*3)

        dis_new = [0.]*(self.number_of_nodes*3-len(self.constraint))
        stiff_new = [[0.]*(self.number_of_nodes*3-len(self.constraint))]*(self.number_of_nodes*3-len(self.constraint))

        stiffness_increment = [0.]*(self.number_of_nodes*3-len(self.constraint))
        for i, kfai in enumerate(self.known_f_a):
            for j, kfaj in enumerate(self.known_f_a):
                stiffness_increment[j] = self.mod_stiffnesses[index][kfai][kfaj]
            stiff_new[i] = [x + y for x, y in zip(stiff_new[i], stiffness_increment)]

        # SOLVING THE MODIFIED STRUCTURE
        if self.configuration.solver == 0:
            dis_new = multiply_matrix_vector(invert(stiff_new), self.force_new)
        else:
            dis_new = np.linalg.solve(np.array(stiff_new), np.array(self.force_new))

        mod_displacement_temp = deepcopy(self._init_displacement)

        for i, kfa in enumerate(self.known_f_a):
            mod_displacement_temp[kfa] = dis_new[i] - self.dis_new[i]

        self.modified_displacements[index] = [x + y for x, y in zip(self.modified_displacements[index], mod_displacement_temp)]

    def evaluate(self):
        """
        Calculates the relative displacement of each individual available unit-modification
        compared to the measured differences (delta).

        delta: [DOF number, difference]

        return effect: [effect on 1. point, effect on 2. point, ..., modification number]
                       where each line number shows the corresponding modification number
        """
        self.effect = [[0.]*(self.number_of_keypoints + 2)]*self.element_num
        self.total_effect = [0.]*self.number_of_keypoints
        self.sorted_effect = [[[0.]*(self.number_of_keypoints + 2)]*self.element_num]*self.number_of_keypoints

        effect_temp = [0.]*(self.number_of_keypoints + 2)

        for modnum in range(self.element_num):
            effect_temp[self.number_of_keypoints] = int(modnum)
            for j, dofnum in enumerate(self.keypoint):
                try:
                    effect_temp[j] = self.modified_displacements[modnum][dofnum]
                    self.effect[modnum] = [x for x in effect_temp]
                    self.total_effect[j] += abs(self.effect[modnum][j])
                except IndexError:
                    print("Maybe the mapping data is invalid.")
                    print("Please check the \'arduino_mapping.txt\' input whether the given DOFs are correct or not.")
                    raise IndexError

        self.effect_ratio = deepcopy(self.effect)
        for i in range(self.element_num):
            for j in range(self.number_of_keypoints):
                if self.total_effect[j] > 0:
                    self.effect_ratio[i][j] = abs(self.effect_ratio[i][j]/self.total_effect[j])
                else:
                    self.effect_ratio[i][j] = 0

        # print("   \'effectratio\' is not used yet")

        # Sort by effectiveness
        for i in range(self.number_of_keypoints):
            self.sorted_effect[i] = deepcopy(self.effect)

            # Check sign of the effect
            for ktemp in range(self.element_num):
                if self.sorted_effect[i][ktemp][i] < 0:
                    for jtemp in range(self.number_of_keypoints):
                        self.sorted_effect[i][ktemp][jtemp] = abs(self.sorted_effect[i][ktemp][jtemp])
                        self.sorted_effect[i][ktemp][self.number_of_keypoints + 1] = -1
                else:
                    self.sorted_effect[i][ktemp][self.number_of_keypoints + 1] = +1

            for j in range(self.number_of_keypoints):
                if i != j and j != 0:
                    self.sorted_effect[i] = swap_columns(sorted(swap_columns(self.sorted_effect[i], 0, j), reverse=True), 0, j)
            if i != 0:
                self.sorted_effect[i] = swap_columns(sorted(swap_columns(self.sorted_effect[i], 0, i), reverse=True), 0, i)
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
                self.configuration.disconnect()
                raise Exception('Watchpoint name error')

            delta.append(measured - num_displ[dof])

        return delta

    def optimize(self, delta):
        """
        Model updating - core function
        """
        # modnum = min(10, self.element_num)
        modnum = self.element_num
        self.modifications = [0.0]*self.element_num

        if not self.configuration.simulation:
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
                self.calculate_modified_stiffness_matrix(index, self.modifications[index])
                self.solve_modified_structure(index)
            self.evaluate()

            self.calculate_modified_stiffness_matrix(self.element_num, 0)
            self.solve_modified_structure(self.element_num)

            newdelta = self.difference(self.modified_displacements[self.element_num], self.measurement)

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

        with open("./Structures/" + self.title + ' - UpdateResults' + appendix + '.txt', 'a') as outfile:
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
            outfile.write("Required iterations: " + str(j) + "\n")
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
                outfile.write(str(self.modified_displacements[self.element_num]) + "\n")
            outfile.write("----------------------\n")

    def capable(self):
        """
        Function telling whether there are more options to modify
        """
        capable = False
        for variable in self.modifications:
            if 0.01 < abs(variable) <= 0.95*self.modification_limit:
                capable = True
        return capable

    def set_error_limit(self, errorlimit):
        """
        Setting general stop parameter for model updating
        """
        if errorlimit > 0.0:
            self.error_limit = errorlimit
        else:
            print("The error limit must be a positive number")
            raise Exception

    def set_modification_limit(self, modification_limit):
        """
        Setting modification limit for members (model updating)
        """
        if 0.0 < modification_limit < 1.0:
            self.modification_limit = modification_limit
        else:
            print("The modification limit must be higher than 0.0 and lower than 1.0")
            raise Exception

    def set_unit_modification(self, unit_modification):
        """
        Setting modification step (model updating)
        """
        if 0.01 <= abs(unit_modification) < 0.5:
            self.unit_modification = unit_modification
        else:
            print("The absolute value of the unit modification must be minimum 0.01 and maximum 0.5")
            raise Exception

    def set_iteration_limit(self, iteration_limit):
        """
        Setting maximum number of iterations (model updating)
        """
        if 1 < int(iteration_limit) <= math.pow(10, 4):
            self.iteration_limit = int(iteration_limit)
        else:
            print("The iteration limit must be between 2 and 10.000")
            raise Exception

    def update_model(self):
        """
        General function to manage model updating procedure.
        """
        self.processed_data = [0.]*len(self.arduino_mapping)

        if not self.configuration.simulation:
            base = self.calibrate()
            filemode = 'a'
        else:
            base = ['SIMULATION']
            try:
                os.remove(str(self.title) + ' - UpdateResults - SIMULATED.txt')
            except Exception:
                pass
            filemode = 'r'

        with open("./Structures/" + str(self.title) + ' - Input Data.txt', filemode) as input_file:
            # Saving input data
            if not self.configuration.simulation:
                input_file.write('Input data of \'' + self.title + '\':\n\n')
                input_file.write('Start Time: ' +
                                str(datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')) + "\n")
                input_file.write('Base: ' + str(base) + '\n')

            new_line = "[],0.0"
            for i in range(1000):
                if not self.configuration.simulation:
                    delta = self.readarduino(base, input_file)
                    self.optimize(delta)
                else:
                    delta = None
                    try:
                        previous_line = new_line
                        new_line = input_file.readline()
                        if not new_line == '':
                            delta = self.simulate_arduino(new_line, previous_line)
                    except Exception:
                        pass
                    if delta:
                        self.optimize(delta)
        print("Update statistics:")
        print("Totally updated models: " + str(TRUSS.num_of_updates[0] + TRUSS.num_of_updates[1] + TRUSS.num_of_updates[2]))
        print("  Successfully updated models: " + str(TRUSS.num_of_updates[0]))
        print("  Updates with running out of possibilities: " + str(TRUSS.num_of_updates[2]))
        print("  Updates did not finished: " + str(TRUSS.num_of_updates[1]))

    def post_process(self):
        """
        Calculates reaction forces and stresses

        :return None
        """
        self._reactions()
        self._stresses()
        self._post_processed = 1

    def _reactions(self):
        """
        Calculates reaction forces
        """
        for i in self.known_displacement_a:
            self.force[i] = 0
            for j, displacement in enumerate(self.displacement):
                self.force[i] += self.stiffness[i][j]*displacement

    def _stresses(self):
        """
        Calculates stress in elements

        Last part: Coloring elements for graphical output
        """
        self.stress = [0.]*self.element_num
        for element in range(self.element_num):
            local_stiffness = [-self._cx[element], -self._cy[element], -self._cz[element],
                        self._cx[element], self._cy[element], self._cz[element]]
            for i in range(3*2):
                self.stress[element] += local_stiffness[i]*self.displacement[self.element_DOF[element][i]]
            self.stress[element] = self.stress[element]*self._norm_stiff[element]

        s_max = max([abs(min(self.stress)), max(self.stress), 0.000000001])
        self.stress_color = [float(x)/float(s_max) for x in self.stress]


##################################
#   BEGINNING OF THE MAIN PART   #
##################################


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Input file, stored in the ./Structure folder [*.str]", )
    parser.add_argument("-t", "--title", help="Title of the project", default="structure")
    parser.add_argument("-c", "--compatibility", metavar='int', type=int, choices=range(4), default=2,
                        help="0: User defined, 1: Most information (with numpy), "
                             "2: Maximum compatibility mode, 3: Android")
    parser.add_argument("-s", "--simulation", metavar='int', type=int, choices=range(2), default=0,
                        help="0: No / 1: Yes")
    args = parser.parse_args()

    # Define new structure
    TRUSS = Truss(str(args.input), str(args.title), args.compatibility, args.simulation)

    if TRUSS.configuration.log:
        TRUSS.configuration.start_logging()

    TRUSS.read()
    TRUSS.configuration.part_time("Setting up structure")

    # Calculate stiffness-matrix
    TRUSS.calculate_stiffness_matrix()
    TRUSS.configuration.part_time("Calculating Stiffness Matrix")

    # Solve structure
    TRUSS.solve()
    TRUSS.configuration.part_time("Solving")

    # Update iteration
    if TRUSS.configuration.updating:
        TRUSS.set_unit_modification(0.05)
        TRUSS.set_error_limit(1.2)
        TRUSS.set_modification_limit(0.7)
        TRUSS.set_iteration_limit(100)
        TRUSS.update_model()
        TRUSS.configuration.part_time("Updating numerical model")

    # Plotting
    if TRUSS.configuration.graphics:
        # Plot settings:
        # O: Original D: Deformed S: Supports F: Forces R: Reactions
        # ScD: Scale displacements (Z-axis) (def:1.0) ScF: Scale forces (def:1.0)
        # ScS: Scale Support signs (Z-axis) (def:1.0)
        # Save: Save plot to file

        #     plot(O, D, S, F, R, ScD, ScF, ScS, Save)
        TRUSS.plot(1, 0, 1, 1, 0, 1.0, 0.0, 0.0, True)
        TRUSS.plot(1, 1, 1, 0, 0, 1.0, 0.0, 0.0, True)
        TRUSS.plot(0, 1, 1, 1, 1, 2.0, 0.0, 0.0, True)
        TRUSS.configuration.part_time("Plotting")

    # Write results to file
    TRUSS.write_results("Structures/" + TRUSS.title + ' - Results.txt')
    TRUSS.configuration.part_time("Wrote results to the output file")

    # Closing Arduino port
    if TRUSS.configuration.arduino:
        TRUSS.close_serial()

    # End logging
    if TRUSS.configuration.log:
        TRUSS.end_logging()
