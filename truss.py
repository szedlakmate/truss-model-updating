# -*- coding: utf-8 -*-
"""
Created on April 30 18:49:40 2018

3D truss model updater program created by Máté Szedlák.
For more information please check the README
Copyright MIT, Máté Szedlák 2016-2018.
"""
import argparse
import itertools
import math
try:
    import numpy
except ImportError:
    print("NumPy could not be loaded")
from copy import deepcopy
try:
    from truss_framework import TrussFramework
    from truss_framework import error
    from truss_framework import plot
    from truss_framework import animate
    from extra_math import invert
    from extra_math import mat_vec_mult as multiply_matrix_vector
    from extra_math import swap_col as swap_columns
except ImportError:
    print("Input is missing")
    print("Please check the truss_framework.py and extra_math.py files.")
    print("Graphical libraries could not be loaded. GUI can not be used.")


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
        if len(self.truss.nodal_coord) != len(list(k for k, _ in itertools.groupby(sorted(self.truss.nodal_coord)))):
            if ignorable == 0:
                raise Exception('Coordinate list has repeating items. Calculation is terminated')
            else:
                print("This node already exists. Input is ignored.")
            return False
        else:
            return True





    def modify_coordinate(self, node_ID, new_coordinate):
        """
        Modify one node's coordinates

        :param node_ID: the ID of the changeable node
        :param new_coordinate: new coordinates wrapped by an array
        :return: Bool - True uf the modification was applied, False otherwise.

        @example: self.modify_coordinate(5, [1.2, 3.6])  # 2D
        """
        if self._check_coordinates(True):
            self.truss.nodal_coord[node_ID] = new_coordinate
            self.invalidate_stiffness_matrices('all')
            return True
        else:
            return False


    def modify_cross_section(self, element_ID, new_area):
        """
        Modifying cross-sections by elements

        :param element_ID: the ID of the changeable element
        :param new_area: the new cross-sectional area
        :return: None
        """
        self.truss.cross_sectional_area_list[element_ID] = new_area
        self.invalidate_stiffness_matrices('all')


    def modify_material(self, element_ID, E):
        """
        Modifying material data by elements

        :param element_ID: ID of the element ehich should be modified
        :param E: {number} Elastic modulo
        :return: None
        """
        self.truss.elastic_modulo[element_ID] = E
        self.invalidate_stiffness_matrices('all')


    def modify_force(self, element_ID, force):
        """
        Modifying specific force

        :param element_ID: the ID of the element which should be mofidifed
        :param force: new force value
        :return: None
        """
        self.truss.force[element_ID] = force
        self.set_postprocess_needed_flag()


    def calculate_stiffness_matrix(self):
        """
        Stiffness matrix compilation
        """
        self.set_postprocess_needed_flag()

        if self.DOF == 2:
            for dof_z in range(self.truss.number_of_nodes()):
                self.truss.constraint.append([int(dof_z*3+2), 0.])
        self.truss.constraint = list(k for k, _ in itertools.groupby(sorted(self.truss.constraint)))

        # Setting known forces
        for location in range(3*self.truss.number_of_nodes()):
            self.truss.known_f_a.append(location)
            if self.truss.force[location] != 0:
                self.truss.known_f_not_zero.append(location)

        self.truss.known_displacement_a = []
        for constraint in self.truss.constraint:
            self._init_displacement[constraint[0]] = constraint[1]
            self.truss.known_displacement_a.append(constraint[0])
            try:
                self.truss.known_f_a.remove(constraint[0])
                self.truss.known_f_not_zero.remove(constraint[0])
            except ValueError:
                pass

        elements_lengths = [0.]*self.truss.number_of_elements()
        self._norm_stiff = [0.]*self.truss.number_of_elements()
        self._cx = [0.]*self.truss.number_of_elements()
        self._cy = [0.]*self.truss.number_of_elements()
        self._cz = [0.]*self.truss.number_of_elements()
        self._s_loc = [0.]*self.truss.number_of_elements()
        local_stiffness_matrix = [0.]*self.truss.number_of_elements()
        self.truss.stress = [0.]*self.truss.number_of_elements()

        self.stiffness_matrix = [[0.]*(self.truss.number_of_nodes()*3)]*(self.truss.number_of_nodes()*3)

        for i in range(self.truss.number_of_elements()):
            elements_lengths[i] =\
                math.sqrt(sum([(j-i)**2 for j, i
                               in zip(self.truss.nodal_coord[self.truss.nodal_connections[i][1]],
                                      self.truss.nodal_coord[self.truss.nodal_connections[i][0]])]))

            self._cx[i] = (self.truss.nodal_coord[self.truss.nodal_connections[i][1]][0]-self.truss.nodal_coord[self.truss.nodal_connections[i][0]][0])/elements_lengths[i]
            self._cy[i] = (self.truss.nodal_coord[self.truss.nodal_connections[i][1]][1]-self.truss.nodal_coord[self.truss.nodal_connections[i][0]][1])/elements_lengths[i]
            self._cz[i] = (self.truss.nodal_coord[self.truss.nodal_connections[i][1]][2]-self.truss.nodal_coord[self.truss.nodal_connections[i][0]][2])/elements_lengths[i]
            self._norm_stiff[i] = self.truss.elastic_modulo[i]/elements_lengths[i]

            # local stiffness matrix calculation
            self._s_loc[i] = [[self._cx[i]**2, self._cx[i]*self._cy[i], self._cx[i]*self._cz[i], -self._cx[i]**2, -self._cx[i]*self._cy[i], -self._cx[i]*self._cz[i]],
                              [self._cx[i]*self._cy[i], self._cy[i]**2, self._cy[i]*self._cz[i], -self._cx[i]*self._cy[i], -self._cy[i]**2, -self._cy[i]*self._cz[i]],
                              [self._cx[i]*self._cz[i], self._cy[i]*self._cz[i], self._cz[i]**2, -self._cx[i]*self._cz[i], -self._cy[i]*self._cz[i], -self._cz[i]**2],
                              [-self._cx[i]**2, -self._cx[i]*self._cy[i], -self._cx[i]*self._cz[i], self._cx[i]**2, self._cx[i]*self._cy[i], self._cx[i]*self._cz[i]],
                              [-self._cx[i]*self._cy[i], -self._cy[i]**2, -self._cy[i]*self._cz[i], self._cx[i]*self._cy[i], self._cy[i]**2, self._cy[i]*self._cz[i]],
                              [-self._cx[i]*self._cz[i], -self._cy[i]*self._cz[i], -self._cz[i]**2, self._cx[i]*self._cz[i], self._cy[i]*self._cz[i], self._cz[i]**2]]
            local_stiffness_matrix[i] = [[y * self.cross_sectional_area_list[i] * self._norm_stiff[i] for y in x] for x in self._s_loc[i]]
            ele_dof_vec = self.element_DOF[i]

            stiffness_increment = [0.]*(self.truss.number_of_nodes()*3)

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
            self.modified_stiffness_matrices = [0.]*(self.truss.number_of_elements()+1)

        # for loopindex in range(self.truss.number_of_elements()):
        modified_stiffness_matrices = [[0.]*(self.truss.number_of_nodes()*3)]*(self.truss.number_of_nodes()*3)

        for i in range(self.truss.number_of_elements()):
            if i == index:
                _mod_norm_stiff = self._norm_stiff[i] * (1.0 + self.updating_container.modifications[i] + magnitude)  # E[i]/L[i]
            else:
                _mod_norm_stiff = self._norm_stiff[i] * (1.0 + self.updating_container.modifications[i])  # E[i]/L[i]

            _mod_loc_stiff = [[y*self.cross_sectional_area_list[i]*_mod_norm_stiff for y in x] for x in self._s_loc[i]]

            ele_dof_vec = self.element_DOF[i]

            stiffness_increment = [0.]*(self.truss.number_of_nodes()*3)

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

        self.dis_new = [0.]*(self.truss.number_of_nodes()*3-len(self.truss.constraint))
        self.truss.force_new = [0.]*(self.truss.number_of_nodes()*3-len(self.truss.constraint))
        self.stiff_new = [[0.]*(self.truss.number_of_nodes()*3-len(self.truss.constraint))]*(self.truss.number_of_nodes()*3-len(self.truss.constraint))

        # known force array
        for i, known_f_a in enumerate(self.truss.known_f_a):
            self.truss.force_new[i] = self.truss.force[known_f_a]

        stiffness_increment = [0.]*(self.truss.number_of_nodes()*3-len(self.truss.constraint))
        for i, kfai in enumerate(self.truss.known_f_a):
            for j, kfaj in enumerate(self.truss.known_f_a):
                stiffness_increment[j] = self.stiffness_matrix[kfai][kfaj]
            self.stiff_new[i] = [x + y for x, y in zip(self.stiff_new[i], stiffness_increment)]

        # SOLVING THE STRUCTURE
        if self.configuration.solver == 0:
            if self.configuration.log:
                print('Built-in solver')
            self.dis_new = multiply_matrix_vector(invert(self.stiff_new), self.truss.force_new)
        else:
            if self.configuration.log:
                print('NumPy solver')
            self.dis_new = numpy.linalg.solve(numpy.array(self.stiff_new), numpy.array(self.truss.force_new))

        self.displacement = deepcopy(self._init_displacement)

        for i, known_f_a in enumerate(self.truss.known_f_a):
            self.displacement[known_f_a] = self.dis_new[i]

        # Deformed shape
        self.truss.nodal_coord_def = []
        for i in range(self.truss.number_of_nodes()):
            self.truss.nodal_coord_def.append([self.truss.nodal_coord[i][0] + self.displacement[i*3+0],
                                        self.truss.nodal_coord[i][1] + self.displacement[i*3+1], self.truss.nodal_coord[i][2] + self.displacement[i*3+2]])

        # Postrpocesses
        self.post_process()

        self.updating_container.modified_displacements = [0.]*(self.truss.number_of_elements()+1)

    def solve_modified_structure(self, index):
        """
        Solver for the modified structures. 'Index' shows the actual modification number.
        """

        self.updating_container.modified_displacements[index] = [0.]*(self.truss.number_of_nodes()*3)

        dis_new = [0.]*(self.truss.number_of_nodes()*3-len(self.truss.constraint))
        stiff_new = [[0.]*(self.truss.number_of_nodes()*3-len(self.truss.constraint))]*(self.truss.number_of_nodes()*3-len(self.truss.constraint))

        stiffness_increment = [0.]*(self.truss.number_of_nodes()*3-len(self.truss.constraint))
        for i, kfai in enumerate(self.truss.known_f_a):
            for j, kfaj in enumerate(self.truss.known_f_a):
                stiffness_increment[j] = self.modified_stiffness_matrices[index][kfai][kfaj]
            stiff_new[i] = [x + y for x, y in zip(stiff_new[i], stiffness_increment)]

        # SOLVING THE MODIFIED STRUCTURE
        if self.configuration.solver == 0:
            dis_new = multiply_matrix_vector(invert(stiff_new), self.truss.force_new)
        else:
            dis_new = numpy.linalg.solve(numpy.array(stiff_new), numpy.array(self.truss.force_new))

        mod_displacement_temp = deepcopy(self._init_displacement)

        for i, kfa in enumerate(self.truss.known_f_a):
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
        self.updating_container.effect = [[0.]*(self.number_of_keypoints + 2)]*self.truss.number_of_elements()
        self.updating_container.total_effect = [0.]*self.number_of_keypoints
        self.updating_container.sorted_effect = [[[0.]*(self.number_of_keypoints + 2)]*self.truss.number_of_elements()]*self.number_of_keypoints

        effect_temp = [0.]*(self.number_of_keypoints + 2)

        for modnum in range(self.truss.number_of_elements()):
            effect_temp[self.number_of_keypoints] = int(modnum)
            for j, dofnum in enumerate(self.truss.keypoints):
                try:
                    effect_temp[j] = self.updating_container.modified_displacements[modnum][dofnum]
                    self.updating_container.effect[modnum] = [x for x in effect_temp]
                    self.updating_container.total_effect[j] += abs(self.updating_container.effect[modnum][j])
                except IndexError:
                    print("Maybe the mapping data is invalid.")
                    print("Please check the \'arduino_mapping.txt\' input whether the given DOFs are correct or not.")
                    raise IndexError

        self.effect_ratio = deepcopy(self.updating_container.effect)
        for i in range(self.truss.number_of_elements()):
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
            for ktemp in range(self.truss.number_of_elements()):
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


    def optimize(self, delta):
        """
        Model updating - core function

        :param delta: difference between the latest measurement and the calculated model 
        :return: True - the optimization step was successful | False - Optimization failed
        """

        modnum = self.truss.number_of_elements()
        self.updating_container.modifications = [0.0]*self.truss.number_of_elements()

        if not self.configuration.simulation:
            appendix = ""
        else:
            appendix = " - SIMULATED"

        new_delta = delta
        j = 0

        print("-----")
        print("Step: 0/" + str(self.updating_container.iteration_limit))

        # Optimization loop
        while error(new_delta) > self.updating_container.error_limit and (self.capable() or j <= 1):
            j += 1

            print("Error: " + str(error(new_delta)))
            print("-----")
            print("Step: " + str(j) + "/" + str(self.updating_container.iteration_limit))

            ratio = [0.]*modnum
            unit = 0

            previous_modifications = self.updating_container.modifications
            """
            for index in range(self.truss.number_of_elements()):
                modification = self.updating_container.modifications[index] - self.updating_container.unit_modification
                self.updating_container.modifications[index] =\
                    min(abs(modification), self.updating_container.modification_limit) * math.copysign(1, modification)
                self.calculate_modified_stiffness_matrix(index, self.updating_container.modifications[index])
                self.solve_modified_structure(index)
            """
# TODO: ERROR: self.solve_modified_structure doubled!
            self.evaluate()

            self.calculate_modified_stiffness_matrix(self.truss.number_of_elements(), 0)
            self.solve_modified_structure(self.truss.number_of_elements())

            new_delta = self.difference(self.updating_container.modified_displacements[self.truss.number_of_elements()],
                                       self.updating_container.measurement)

            for i, effect in enumerate(self.updating_container.total_effect):
                if effect == 0.0:
                    print("None of the variables has effect on " + str(self.updating_container.arduino_mapping[i]))
                    print("Model updating has no solution.")
                    self.updating_container.number_of_updates[2] += 1
                    raise Exception

            #for i in range(self.truss.number_of_elements()):
            modified_node_ID = self.updating_container.sorted_effect[0][i][1]

            ratio[modified_node_ID] = min(abs(self.updating_container.total_effect[0]/self.updating_container.sorted_effect[0][i][0]),self.updating_container.modification_limit) * math.copysign(1, self.updating_container.sorted_effect[0][i][2])

            unit += abs(ratio[modified_node_ID]*self.updating_container.sorted_effect[0][i][0])
            print(new_delta)
            scale = new_delta[0]/unit

            #for i in range(self.truss.number_of_elements()):
            #modified_node_ID = self.updating_container.sorted_effect[0][i][1]
            self.updating_container.modifications[modified_node_ID] =\
                min(abs(previous_modifications[modified_node_ID] - self.updating_container.unit_modification*ratio[modified_node_ID]),
                self.updating_container.modification_limit) *\
                math.copysign(1, previous_modifications[modified_node_ID] - self.updating_container.unit_modification*ratio[modified_node_ID])
                # the last part is already the sign itself without the sign function

            print("Ratio: " + str(scale))

            if j >= self.updating_container.iteration_limit:
                print("Iteration limit is reached.")
                print("The remaining error is: " + str(error(new_delta)))
                self.updating_container.number_of_updates[1] += 1
                return False

        print("Final error: " + str(error(new_delta)))

        if not self.capable() and j > 1:
            print("Optimization could not be finished successfully.")
            print("The remaining error is: " + str(error(new_delta)))
            self.updating_container.number_of_updates[2] += 1
            return False

        self.write_output_stream(j, appendix)
        self.updating_container.number_of_updates[0] += 1
        return True

    def capable(self):
        """
        Function telling whether there are more options to modify
        """
        capable = False
        for variable in self.updating_container.modifications:
            if abs(variable) <= 0.95 * self.updating_container.modification_limit:
                capable = True
                break
            else:
                pass
                # TODO: not to reach this part of the code
        return capable
        
    def start_model_updating(self, unit_modification=0.05, error_limit=1.2, modification_limit=0.7, iteration_limit=100):
        """
        Configuring the solver. If the iterations do not converge, settings shall be tuned here

        :param unit_modification: Setting modification step (model updating)
        :param error_limit: Setting general stop parameter for model updating
        :param modification_limit: Setting modification limit for members (model updating)
        :param iteration_limit: Setting maximum number of iterations (model updating)
        :return: None
        """
        if 0.01 <= abs(unit_modification) < 0.5:
            self.updating_container.unit_modification = unit_modification
        else:
            print("The absolute value of the unit modification must be minimum 0.01 and maximum 0.5")
            raise Exception
        if error_limit > 0.0:
            self.updating_container.error_limit = error_limit
        else:
            print("The error limit must be a positive number")
            raise Exception

        if 0.0 < modification_limit < 1.0:
            self.updating_container.modification_limit = modification_limit
        else:
            print("The modification limit must be higher than 0.0 and lower than 1.0")
            raise Exception

        if 1 < int(iteration_limit) <= math.pow(10, 4):
            self.updating_container.iteration_limit = int(iteration_limit)
        else:
            print("The iteration limit must be between 2 and 10.000")
            raise Exception

        # TODO: This call starts the optimization loop in a very unlucky way. This should be refactored!
        self.update_model()

    def update_model(self):
        """
        General function to manage model updating procedure.
        """
        self.processed_data = [0.]*1 #*len(self.arduino_mapping)

        if not self.configuration.simulation:
            base = self.set_base()
            filemode = 'a'
            # TODO: complete this process
        else:
            #self.bulk_set_measurement_points()
            base = ['SIMULATION']
            filemode = 'r'
            self.check_folder('Simulation')
            with open("./Simulation/" + str(self.title) + ' - Input Data.txt', filemode) as input_file:
                new_line = "[],0.0"
                for i in range(1000):
                        delta = None
                        previous_line = new_line
                        new_line = input_file.readline()
                        if not new_line == '':
                            delta = self.mock_delta(new_line, previous_line)
                        if delta:
                            if self.updating_container.original_delta == []:
                                self.updating_container.original_delta = delta
                                self.updating_container.latest_delta = delta
                            if not self.optimize(delta):
                                break


    def post_process(self):
        """
        Calculates reaction forces and stresses

        :return None
        """
        self._reactions()
        self._stresses(self.truss)
        self._post_processed = 1

    def _reactions(self):
        """
        Calculates reaction forces
        """
        for i in self.truss.known_displacement_a:
            self.truss.force[i] = 0
            for j, displacement in enumerate(self.displacement):
                self.truss.force[i] += self.stiffness_matrix[i][j]*displacement


    # TODO: should get a TrussModelData object
    def _stresses(self, truss):
        """
        Calculates stress in elements

        Last part: Coloring elements for graphical output
        """
        truss.stress = [0.]*truss.number_of_elements()
        # TODO: Should be calculated based on inner data (displacements, E, A)

        for element in range(truss.number_of_elements()):

            """
            local_stiffness = [-truss._cx[element], -truss._cy[element], -truss._cz[element],
                        truss._cx[element], truss._cy[element], truss._cz[element]]
            for i in range(3*2):
                truss.stress[element] += local_stiffness[i]*truss.displacement[truss.element_DOF[element][i]]
            truss.stress[element] = truss.stress[element]*truss._norm_stiff[element]
        """
            truss.stress[element] = \
                (truss.element_length(element, 'deformed') - truss.element_length(element, 'original')) / truss.element_length(element, 'original') * truss.elastic_modulo[element]

        s_max = max([abs(min(truss.stress)), max(truss.stress), 0.000000001])


        # TODO: Should be written out using a function
        truss.stress_color = [float(x)/float(s_max) for x in truss.stress]


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

    # Plotting
    if TRUSS.configuration.graphics:
        plot(TRUSS.truss, original=True, result=False, supports=True, forces=True, reactions=False, label_classes='load')
        plot(TRUSS.truss, original=True, result=True, supports=True, forces=False, reactions=False)
        plot(TRUSS.truss, original=False, result=True, supports=True, forces=True, reactions=True, label_classes='load')
        animate(TRUSS.truss, 'load')
        TRUSS.configuration.part_time("Plotting")

    # Write results to file
    TRUSS.write_results()
    TRUSS.configuration.part_time("Wrote results to the output file : Structures/<title> - Results.txt".format(TRUSS.title))

    # Update iteration
    if TRUSS.configuration.updating:
        TRUSS.start_model_updating(unit_modification=0.05, error_limit=1.2, modification_limit=0.7, iteration_limit=100)
        TRUSS.configuration.part_time("Updating numerical model")

    # TODO: plot the updated model

    # Closing Arduino port
    if TRUSS.configuration.arduino:
        TRUSS.configuration.disconnect()

    # End logging
    if TRUSS.configuration.log:
        TRUSS.end_logging()
