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
try:
    import numpy
except ImportError:
    print("NumPy could not be loaded")
from copy import deepcopy
try:
    from truss_framework import TrussFramework
    from extra_math import invert
    from extra_math import mat_vec_mult as multiply_matrix_vector
    from extra_math import swap_col as swap_columns
except ImportError:
    print("Input is missing")
    print("Please check the truss_framework.py, truss_graphics.py and extra_math.py files.")
    print("Graphical libraries could not be loaded. GUI can not be used.")


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
        self.number_of_elements = len(self.nodal_connections)
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
        self.stiffness_matrix = [0.]*(3*self.number_of_nodes)
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

    def bulk_set_measurement_points(self, measurement_points=['13Y']):
        """
        Set special nodal DOFs
        """
        self.analysis = {}

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
            for dof_z in range(self.number_of_nodes):
                self.constraint.append([int(dof_z*3+2), 0.])
        self.constraint = list(k for k, _ in itertools.groupby(sorted(self.constraint)))

        # Setting known forces
        for location in range(3*self.number_of_nodes):
            self.known_f_a.append(location)
            if self.force[location] != 0:
                self.known_f_not_zero.append(location)

        self.known_displacement_a = []
        for constraint in self.constraint:
            self._init_displacement[constraint[0]] = constraint[1]
            self.known_displacement_a.append(constraint[0])
            try:
                self.known_f_a.remove(constraint[0])
                self.known_f_not_zero.remove(constraint[0])
            except ValueError:
                pass

        elements_lengths = [0.]*self.number_of_elements
        self._norm_stiff = [0.]*self.number_of_elements
        self._cx = [0.]*self.number_of_elements
        self._cy = [0.]*self.number_of_elements
        self._cz = [0.]*self.number_of_elements
        self._s_loc = [0.]*self.number_of_elements
        local_stiffness_matrix = [0.]*self.number_of_elements
        self.stress = [0.]*self.number_of_elements

        self.stiffness_matrix = [[0.]*(self.number_of_nodes*3)]*(self.number_of_nodes*3)

        for i in range(self.number_of_elements):
            elements_lengths[i] =\
                math.sqrt(sum([(j-i)**2 for j, i
                               in zip(self.nodal_coord[self.nodal_connections[i][1]],
                                      self.nodal_coord[self.nodal_connections[i][0]])]))

            self._cx[i] = (self.nodal_coord[self.nodal_connections[i][1]][0]-self.nodal_coord[self.nodal_connections[i][0]][0])/elements_lengths[i]
            self._cy[i] = (self.nodal_coord[self.nodal_connections[i][1]][1]-self.nodal_coord[self.nodal_connections[i][0]][1])/elements_lengths[i]
            self._cz[i] = (self.nodal_coord[self.nodal_connections[i][1]][2]-self.nodal_coord[self.nodal_connections[i][0]][2])/elements_lengths[i]
            self._norm_stiff[i] = self.elastic_modulo[i]/elements_lengths[i]

            # local stiffness matrix calculation
            self._s_loc[i] = [[self._cx[i]**2, self._cx[i]*self._cy[i], self._cx[i]*self._cz[i], -self._cx[i]**2, -self._cx[i]*self._cy[i], -self._cx[i]*self._cz[i]],
                              [self._cx[i]*self._cy[i], self._cy[i]**2, self._cy[i]*self._cz[i], -self._cx[i]*self._cy[i], -self._cy[i]**2, -self._cy[i]*self._cz[i]],
                              [self._cx[i]*self._cz[i], self._cy[i]*self._cz[i], self._cz[i]**2, -self._cx[i]*self._cz[i], -self._cy[i]*self._cz[i], -self._cz[i]**2],
                              [-self._cx[i]**2, -self._cx[i]*self._cy[i], -self._cx[i]*self._cz[i], self._cx[i]**2, self._cx[i]*self._cy[i], self._cx[i]*self._cz[i]],
                              [-self._cx[i]*self._cy[i], -self._cy[i]**2, -self._cy[i]*self._cz[i], self._cx[i]*self._cy[i], self._cy[i]**2, self._cy[i]*self._cz[i]],
                              [-self._cx[i]*self._cz[i], -self._cy[i]*self._cz[i], -self._cz[i]**2, self._cx[i]*self._cz[i], self._cy[i]*self._cz[i], self._cz[i]**2]]
            local_stiffness_matrix[i] = [[y * self.cross_sectional_area_list[i] * self._norm_stiff[i] for y in x] for x in self._s_loc[i]]
            ele_dof_vec = self.element_DOF[i]

            stiffness_increment = [0.]*(self.number_of_nodes*3)

            for j in range(3*2):
                for k in range(3*2):
                    stiffness_increment[ele_dof_vec[k]] = local_stiffness_matrix[i][j][k]
                self.stiffness_matrix[ele_dof_vec[j]] = [x + y for x, y in zip(self.stiffness_matrix[ele_dof_vec[j]], stiffness_increment)]
        self._stiff_is_fresh = 1

    def calculate_modified_stiffness_matrix(self, index, magnitude):
        """
        Convergence step in stiffness matrix modification
        """
        if not self.modified_stiffness_matrices:
            self.modified_stiffness_matrices = [0.]*(self.number_of_elements+1)

        # for loopindex in range(self.number_of_elements):
        modified_stiffness_matrices = [[0.]*(self.number_of_nodes*3)]*(self.number_of_nodes*3)

        for i in range(self.number_of_elements):
            if i == index:
                _mod_norm_stiff = self._norm_stiff[i] * (1.0 + self.modifications[i] + magnitude)  # E[i]/L[i]
            else:
                _mod_norm_stiff = self._norm_stiff[i] * (1.0 + self.modifications[i])  # E[i]/L[i]

            _mod_loc_stiff = [[y*self.cross_sectional_area_list[i]*_mod_norm_stiff for y in x] for x in self._s_loc[i]]

            ele_dof_vec = self.element_DOF[i]

            stiffness_increment = [0.]*(self.number_of_nodes*3)

            for j in range(3*2):
                for k in range(3*2):
                    stiffness_increment[ele_dof_vec[k]] = _mod_loc_stiff[j][k]
                modified_stiffness_matrices[ele_dof_vec[j]] = [x + y for x, y in
                                                         zip(modified_stiffness_matrices[ele_dof_vec[j]], stiffness_increment)]

        self.modified_stiffness_matrices[index] = modified_stiffness_matrices

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
                stiffness_increment[j] = self.stiffness_matrix[kfai][kfaj]
            self.stiff_new[i] = [x + y for x, y in zip(self.stiff_new[i], stiffness_increment)]

        # SOLVING THE STRUCTURE
        if self.configuration.solver == 0:
            if self.configuration.log:
                print('Built-in solver')
            self.dis_new = multiply_matrix_vector(invert(self.stiff_new), self.force_new)
        else:
            if self.configuration.log:
                print('NumPy solver')
            self.dis_new = numpy.linalg.solve(numpy.array(self.stiff_new), numpy.array(self.force_new))

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

        self.updating_container.modified_displacements = [0.]*(self.number_of_elements+1)

    def solve_modified_structure(self, index):
        """
        Solver for the modified structures. 'Index' shows the actual modification number.
        """

        self.updating_container.modified_displacements[index] = [0.]*(self.number_of_nodes*3)

        dis_new = [0.]*(self.number_of_nodes*3-len(self.constraint))
        stiff_new = [[0.]*(self.number_of_nodes*3-len(self.constraint))]*(self.number_of_nodes*3-len(self.constraint))

        stiffness_increment = [0.]*(self.number_of_nodes*3-len(self.constraint))
        for i, kfai in enumerate(self.known_f_a):
            for j, kfaj in enumerate(self.known_f_a):
                stiffness_increment[j] = self.modified_stiffness_matrices[index][kfai][kfaj]
            stiff_new[i] = [x + y for x, y in zip(stiff_new[i], stiffness_increment)]

        # SOLVING THE MODIFIED STRUCTURE
        if self.configuration.solver == 0:
            dis_new = multiply_matrix_vector(invert(stiff_new), self.force_new)
        else:
            dis_new = numpy.linalg.solve(numpy.array(stiff_new), numpy.array(self.force_new))

        mod_displacement_temp = deepcopy(self._init_displacement)

        for i, kfa in enumerate(self.known_f_a):
            mod_displacement_temp[kfa] = dis_new[i] - self.dis_new[i]

        self.updating_container.modified_displacements[index] = [x + y for x, y in zip(self.updating_container.modified_displacements[index], mod_displacement_temp)]

    def evaluate(self):
        """
        Calculates the relative displacement of each individual available unit-modification
        compared to the measured differences (delta).

        delta: [DOF number, difference]

        return effect: [effect on 1. point, effect on 2. point, ..., modification number]
                       where each line number shows the corresponding modification number
        """
        self.updating_container.effect = [[0.]*(self.number_of_keypoints + 2)]*self.number_of_elements
        self.updating_container.total_effect = [0.]*self.number_of_keypoints
        self.updating_container.sorted_effect = [[[0.]*(self.number_of_keypoints + 2)]*self.number_of_elements]*self.number_of_keypoints

        effect_temp = [0.]*(self.number_of_keypoints + 2)

        for modnum in range(self.number_of_elements):
            effect_temp[self.number_of_keypoints] = int(modnum)
            for j, dofnum in enumerate(self.keypoint):
                try:
                    effect_temp[j] = self.updating_container.modified_displacements[modnum][dofnum]
                    self.updating_container.effect[modnum] = [x for x in effect_temp]
                    self.updating_container.total_effect[j] += abs(self.updating_container.effect[modnum][j])
                except IndexError:
                    print("Maybe the mapping data is invalid.")
                    print("Please check the \'arduino_mapping.txt\' input whether the given DOFs are correct or not.")
                    raise IndexError

        self.effect_ratio = deepcopy(self.updating_container.effect)
        for i in range(self.number_of_elements):
            for j in range(self.number_of_keypoints):
                if self.updating_container.total_effect[j] > 0:
                    self.effect_ratio[i][j] = abs(self.effect_ratio[i][j]/self.updating_container.total_effect[j])
                else:
                    self.effect_ratio[i][j] = 0

        # print("   \'effectratio\' is not used yet")

        # Sort by effectiveness
        for i in range(self.number_of_keypoints):
            self.updating_container.sorted_effect[i] = deepcopy(self.updating_container.effect)

            # Check sign of the effect
            for ktemp in range(self.number_of_elements):
                if self.updating_container.sorted_effect[i][ktemp][i] < 0:
                    for jtemp in range(self.number_of_keypoints):
                        self.updating_container.sorted_effect[i][ktemp][jtemp] = abs(self.updating_container.sorted_effect[i][ktemp][jtemp])
                        self.updating_container.sorted_effect[i][ktemp][self.number_of_keypoints + 1] = -1
                else:
                    self.updating_container.sorted_effect[i][ktemp][self.number_of_keypoints + 1] = +1

            for j in range(self.number_of_keypoints):
                if i != j and j != 0:
                    self.updating_container.sorted_effect[i] = swap_columns(sorted(swap_columns(self.updating_container.sorted_effect[i], 0, j), reverse=True), 0, j)
            if i != 0:
                self.updating_container.sorted_effect[i] = swap_columns(sorted(swap_columns(self.updating_container.sorted_effect[i], 0, i), reverse=True), 0, i)
            else:
                self.updating_container.sorted_effect[i] = sorted(self.updating_container.sorted_effect[i], reverse=True)

    def difference(self, num_displ, measurement):
        """
        Calculate the difference between the Numerical solution and Real-life measurement.
        The Real-life measurement should be given the following format:

            MEASUREMENT: [[13X, -2.154], [16Y, 5.256], ...]
        """
        # print(nodenumber option should be added! <XXX>)

        delta = []

        for measured in measurement:
            #loc, measured = package[0], package[1]
            loc = '13Y'
            try:
                dof = self.analysis[loc.upper()]
            except KeyError:
                print('The given measurement location cannot be matched with the input data.')
                print('The available nodes are: {\'NAMES\': mapping addresses}')
                print(self.analysis)
                self.configuration.disconnect()
                raise Exception('Watchpoint name error')

            delta.append(measured - num_displ[dof])

        print('delta: ' + str(delta))
        return delta

    def optimize(self, delta):
        """
        Model updating - core function
        """

        modnum = self.number_of_elements
        self.modifications = [0.0]*self.number_of_elements

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

            for index in range(self.number_of_elements):
                self.modifications[index] = min(abs(self.modifications[index] - self.unit_modification),
                                                self.modification_limit) * \
                                            math.copysign(1, self.modifications[index] - self.unit_modification)
                self.calculate_modified_stiffness_matrix(index, self.modifications[index])
                self.solve_modified_structure(index)
            self.evaluate()

            self.calculate_modified_stiffness_matrix(self.number_of_elements, 0)
            self.solve_modified_structure(self.number_of_elements)

            newdelta = self.difference(self.updating_container.modified_displacements[self.number_of_elements],
                                       self.updating_container.measurement)

            for i, effect in enumerate(self.updating_container.total_effect):
                if effect == 0.0:
                    print("None of the variables has effect on " + str(self.arduino_mapping[i]))
                    print("Model updating has no solution.")
                    raise Exception

            for i in range(self.number_of_elements):
                modificationnumber = self.updating_container.sorted_effect[0][i][1]
                ratio[modificationnumber] = abs(self.updating_container.sorted_effect[0][i][0] / self.updating_container.total_effect[0]) * \
                                            math.copysign(1, self.updating_container.sorted_effect[0][i][2])
                unit += abs(ratio[modificationnumber]*self.updating_container.sorted_effect[0][i][0])
            print(newdelta)
            scale = newdelta[0]/unit
            for i in range(self.number_of_elements):
                modificationnumber = self.updating_container.sorted_effect[0][i][1]
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

        self.write_output_stream(j, appendix)

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
        self.processed_data = [0.]*1 #*len(self.arduino_mapping)

        if not self.configuration.simulation:
            base = self.set_base()
            filemode = 'a'
        else:
            #self.bulk_set_measurement_points()
            base = ['SIMULATION']
            filemode = 'r'
            with open("./Simulation/" + str(self.title) + ' - Input Data.txt', filemode) as input_file:
                new_line = "[],0.0"
                for i in range(1000):
                        delta = None
                        previous_line = new_line
                        new_line = input_file.readline()
                        if not new_line == '':
                            delta = self.mock_delta(new_line, previous_line)
                        if delta:
                            self.optimize(delta)

        print("Update statistics:")
        print("Totally updated models: " + str(self.updating_container.number_of_updates[0] + self.updating_container.number_of_updates[1] + self.updating_container.number_of_updates[2]))
        print("  Successfully updated models: " + str(self.updating_container.number_of_updates[0]))
        print("  Updates with running out of possibilities: " + str(self.updating_container.number_of_updates[2]))
        print("  Updates did not finished: " + str(self.updating_container.number_of_updates[1]))

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
                self.force[i] += self.stiffness_matrix[i][j]*displacement

    def _stresses(self):
        """
        Calculates stress in elements

        Last part: Coloring elements for graphical output
        """
        self.stress = [0.]*self.number_of_elements
        for element in range(self.number_of_elements):
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

# Setup console run
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--title", metavar='str', type=str,
                        help="Manually label project. By default it comes from the input file's name.", default='')
    parser.add_argument("-c", "--compatibility", metavar='int', type=int, choices=range(4), default=1,
                        help="0: User defined, 1: Most information (with numpy), "
                             "2: Maximum compatibility mode, 3: Android")
    parser.add_argument("-s", "--simulation", metavar='int', type=int, choices=range(2), default=1, help="0: No|1: Yes")
    parser.add_argument("-i", "--input", metavar='str', type=str, default="", help="Input file, stored in the ./Structure folder [*.str]")
    args = parser.parse_args()

    if args.input != "":
        input = args.input
    else:
        input = "bridge"
        print("DEMO MODE")
        print("You started the program in demo mode.\nThe following example is based on the bridge.str file which "
              "is located in the Structures folder.")


    # Define new structure
    TRUSS = Truss(input_file=input, title=args.title, compatibility_mode=args.compatibility,
                  simulation=args.simulation)

    if TRUSS.configuration.log:
        TRUSS.configuration.start_logging()

    TRUSS.read()
    TRUSS.configuration.part_time("Setting up structure")

    # Calculate stiffness-matrix
    TRUSS.calculate_stiffness_matrix()
    TRUSS.configuration.part_time("Calculating stiffness matrix")

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
    TRUSS.configuration.part_time("Wrote results to the output file : Structures/{} - Results.txt".format(TRUSS.title))

    # Closing Arduino port
    if TRUSS.configuration.arduino:
        TRUSS.close_serial()

    # End logging
    if TRUSS.configuration.log:
        TRUSS.end_logging()

