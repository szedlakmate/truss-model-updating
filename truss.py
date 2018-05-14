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
    from truss_framework import ModelUpdatingContainer
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

    """
   
    def calculate_modified_stiffness_matrix(self, index, magnitude):
      
        if not self.updating_container.trusses:
            self.updating_container.trusses = [0.]*(self.truss.number_of_elements()+1)

        # for loopindex in range(self.truss.number_of_elements()):
        modified_stiffness_matrices = [[0.]*(self.truss.number_of_nodes()*3)]*(self.truss.number_of_nodes()*3)

        for i in range(self.truss.number_of_elements()):
            if i == index:
                _mod_norm_stiff = self.truss._norm_stiff[i] * (1.0 + self.updating_container.modifications[i] + magnitude)  # E[i]/L[i]
            else: 
                _mod_norm_stiff = self.truss._norm_stiff[i] * (1.0 + self.updating_container.modifications[i])  # E[i]/L[i]

            _mod_loc_stiff = [[y*self.truss.cross_sectional_area_list[i]*_mod_norm_stiff for y in x] for x in self._s_loc[i]]

            ele_dof_vec = self.truss.element_DOF[i]

            stiffness_increment = [0.]*(self.truss.number_of_nodes()*3)

            for j in range(3*2):
                for k in range(3*2):
                    stiffness_increment[ele_dof_vec[k]] = _mod_loc_stiff[j][k]
                modified_stiffness_matrices[ele_dof_vec[j]] = [x + y for x, y in
                                                               zip(modified_stiffness_matrices[ele_dof_vec[j]], stiffness_increment)]

        self.updating_container.trusses[index] = modified_stiffness_matrices
       Convergence step in stiffness matrix modification
          """
    """
    def solve_modified_structure(self, index):
        
        Solver for the modified structures. 'Index' shows the actual modification number.
        

        self.updating_container.modified_displacements[index] = [0.]*(self.truss.number_of_nodes()*3)

        dis_new = [0.]*(self.truss.number_of_nodes()*3-len(self.truss.constraint))
        stiff_new = [[0.]*(self.truss.number_of_nodes()*3-len(self.truss.constraint))]*(self.truss.number_of_nodes()*3-len(self.truss.constraint))

        stiffness_increment = [0.]*(self.truss.number_of_nodes()*3-len(self.truss.constraint))
        for i, kfai in enumerate(self.truss.known_f_a):
            for j, kfaj in enumerate(self.truss.known_f_a):
                stiffness_increment[j] = self.updating_container.trusses[index].stiffness_matrix[kfai][kfaj]
            stiff_new[i] = [x + y for x, y in zip(stiff_new[i], stiffness_increment)]

        # SOLVING THE MODIFIED STRUCTURE
        if self.configuration.solver == 0:
            dis_new = multiply_matrix_vector(invert(stiff_new), self.truss.force_new)
        else:
            dis_new = numpy.linalg.solve(numpy.array(stiff_new), numpy.array(self.truss.force_new))

        mod_displacement_temp = deepcopy(self.truss._init_displacement)

        for i, kfa in enumerate(self.truss.known_f_a):
            mod_displacement_temp[kfa] = dis_new[i] - self.dis_new[i]

        self.updating_container.modified_displacements[index] = [x + y for x, y in zip(self.updating_container.modified_displacements[index], mod_displacement_temp)]
"""
    def evaluate_updates(self):
        """
        Calculates the relative displacement of each individual available unit-modification
        compared to the measured differences (delta).

        delta: [DOF number, difference]

        return effect: [effect on 1. point, effect on 2. point, ..., modification number]
                       where each line number shows the corresponding modification number
        """
        number_of_keypoints = self.truss.number_of_keypoints()
        number_of_elements = self.truss.number_of_elements()
        
        self.updating_container.effect = [[0.]*(number_of_keypoints + 2)]*number_of_elements
        self.updating_container.total_effect = [0.]*number_of_keypoints
        self.updating_container.sorted_effect = [[[0.]*(number_of_keypoints + 2)]*number_of_elements]*number_of_keypoints

        effect_temp = [0.]*(number_of_keypoints + 2)

        for element_ID in range(number_of_elements):
            effect_temp[number_of_keypoints] = int(element_ID)
            for j, dofnum in enumerate(self.truss.keypoints):
                try:
                    effect_temp[j] = self.updating_container[element_ID].displacement[dofnum]
                    self.updating_container.effect[element_ID] = [x for x in effect_temp]
                    self.updating_container.total_effect[j] += abs(self.updating_container.effect[element_ID][j])
                except IndexError:
                    print("Maybe the mapping data is invalid.")
                    print("Please check the \'arduino_mapping.txt\' input whether the given DOFs are correct or not.")
                    raise IndexError

        self.updating_container.effect_ratio = deepcopy(self.updating_container.effect)
        for i in range(number_of_elements):
            for j in range(number_of_keypoints):
                if self.updating_container.total_effect[j] > 0:
                    self.updating_container.effect_ratio[i][j] = abs(self.updating_container.effect_ratio[i][j]/self.updating_container.total_effect[j])
                else:
                    self.updating_container.effect_ratio[i][j] = 0

        # print("   \'effectratio\' is not used yet")

        # Sort by effectiveness
        for i in range(number_of_keypoints):
            self.updating_container.sorted_effect[i] = deepcopy(self.updating_container.effect)

            # Check sign of the effect
            for ktemp in range(number_of_elements):
                if self.updating_container.sorted_effect[i][ktemp][i] < 0:
                    for jtemp in range(number_of_keypoints):
                        self.updating_container.sorted_effect[i][ktemp][jtemp] = abs(self.updating_container.sorted_effect[i][ktemp][jtemp])
                        self.updating_container.sorted_effect[i][ktemp][number_of_keypoints + 1] = -1
                else:
                    self.updating_container.sorted_effect[i][ktemp][number_of_keypoints + 1] = +1

            for j in range(number_of_keypoints):
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

        element_ID = self.truss.number_of_elements()
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

            ratio = [0.]*element_ID
            unit = 0

            previous_modifications = self.updating_container.modifications
            # TODO: apply Truss classes!!!!!!!
            for index in range(self.truss.number_of_elements()):
                modification = self.updating_container.modifications[index] - self.updating_container.unit_modification

                self.updating_container.modifications[index] = \
                    min(abs(modification), self.updating_container.modification_limit) * math.copysign(1, modification)

                #self.calculate_modified_stiffness_matrix(index, self.updating_container.modifications[index])
                self.updating_container.trusses[index].calculate_stiffness_matrix()
                self.updating_container.trusses[index].solve(self.configuration)

            # TODO: ERROR: self.solve_modified_structure doubled!
            self.evaluate_updates()

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
            self.updating_container.modifications[modified_node_ID] = \
                min(abs(previous_modifications[modified_node_ID] - self.updating_container.unit_modification*ratio[modified_node_ID]),
                    self.updating_container.modification_limit) * \
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
        self.updating_container.trusses = [self.truss] * self.truss.number_of_elements()

        if 0.01 <= abs(unit_modification) < 0.5:
            self.updating_container.unit_modification = unit_modification
        else:
            print("The absolute value of the unit modification must be minimum 0.01 and maximum 0.5\nGiven: %s" % unit_modification)
            raise Exception
        if error_limit > 0.0:
            self.updating_container.error_limit = error_limit
        else:
            print("The error limit must be a positive number\nGiven: %s" % error_limit)
            raise Exception

        if 0.0 < modification_limit < 1.0:
            self.updating_container.modification_limit = modification_limit
        else:
            print("The modification limit must be higher than 0.0 and lower than 1.0\nGiven: %s" % modification_limit)
            raise Exception

        if 1 < int(iteration_limit) <= math.pow(10, 4):
            self.updating_container.iteration_limit = int(iteration_limit)
        else:
            print("The iteration limit must be between 2 and 10.000\nGiven: %s" % iteration_limit)
            raise Exception

            # TODO: This call starts the optimization loop in a very unlucky way. This should be refactored!
        self.update_model()

    def update_model(self):
        """
        General function to manage model updating procedure.
        """
        self.processed_data = [0.]*1 #*len(self.arduino_mapping)
        #self.updating_container.trusses = [T] * self.truss.number_of_elements()

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
    TRUSS.truss.calculate_stiffness_matrix()
    TRUSS.configuration.part_time("Calculating stiffness matrix")

    # Solve structure
    TRUSS.truss.solve(TRUSS.configuration)
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
