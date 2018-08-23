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
from truss_framework import TrussFramework
from truss_framework import error
from truss_framework import plot
from truss_framework import animate
from extra_math import swap_col as swap_columns


class Truss(TrussFramework):
    """
    Object for solving a truss
    """

    def _check_coordinates(self, ignorable):
        """
        Checking coordinates for repeated elements during input file reading.

        :param ignorable: {Boolean} If the warning is ignorable, only message appears and the input is neglected.
                If the warning is not ignorable, exceptions will be raised.
        :return: {Boolean} True if no error found, otherwise False.
        """
        if len(self.truss.nodal_coord) != len(list(k for k, _ in itertools.groupby(sorted(self.truss.nodal_coord)))):
            if ignorable is False:
                raise Exception('Coordinate list has recurring items. Calculation is terminated')
            else:
                print("This node already exists. Input is ignored.")
            return False
        else:
            return True

    def modify_coordinate(self, node_id, new_coordinate):
        """
        Modify one node's coordinates

        :param node_id: {number} the Id of the changeable node
        :param new_coordinate: {Array<number>} new coordinates wrapped by an array
        :return: {Boolean} True uf the modification was applied, False otherwise.

        @example: self.modify_coordinate(5, [1.2, 3.6])  # 2D
        """
        if self._check_coordinates(True):
            self.truss.nodal_coord[node_id] = new_coordinate
            self.truss.invalidate_stiffness_matrices()
            return True
        else:
            return False

    def modify_cross_section(self, element_id, new_area):
        """
        Modifying cross-sections by elements

        :param element_id: {number} the Id of the changeable element
        :param new_area: {number} the new cross-sectional area
        :return: None
        """
        self.truss.cross_sectional_area_list[element_id] = new_area
        self.truss.invalidate_stiffness_matrices()

    def modify_material(self, element_id, elastic_modulo):
        """
        Modifying material data by elements

        :param element_id: {number} Id of the element which should be modified
        :param elastic_modulo: {number} Elastic modulo
        :return: None
        """
        self.truss.elastic_modulo[element_id] = elastic_modulo
        self.truss.invalidate_stiffness_matrices()

    def modify_force(self, element_id, force):
        """
        Modifying specific force

        :param element_id: {number} the Id of the element which should be modified
        :param force: {number} new force value
        :return: None
        """
        self.truss.force[element_id] = force
        self.truss.set_postprocess_needed_flag()

    def evaluate_updates(self):
        """
        Calculates the relative displacements of each individual available unit-modification
        compared to the measured differences (delta).

        delta: [dof number, difference]

        return effect: [effect on 1. point, effect on 2. point, ..., modification number]
                       where each line number shows the corresponding modification number

        :return: None
        """
        number_of_keypoints = self.truss.number_of_keypoints()
        number_of_elements = self.truss.number_of_elements()

        # effect = [[statistic at index], [statistic at index], ...]
        self.updating_container.effect = [[0] * number_of_keypoints] * number_of_elements

        # sum of effects
        self.updating_container.total_effect = [0] * number_of_keypoints

        self.updating_container.sorted_effect = \
            [[[0] * number_of_keypoints] * number_of_elements]*number_of_keypoints

        effect_temp = [0] * (number_of_keypoints + 1)

        # Determine the effects of each modification
        for element_id in range(number_of_elements):

            # Modified element's Id
            effect_temp[number_of_keypoints] = element_id

            # Determine displacement vector corresponding to a modification
            for keypoint_count, dof_Id in enumerate(self.truss.keypoints):
                try:
                    # Pick displacement vector at the measurement point
                    effect_temp[keypoint_count] = self.updating_container.trusses[element_id].displacements[dof_Id]

                    # Copy value
                    self.updating_container.effect[element_id] = [x for x in effect_temp]

                    # Add effect of index-th modification
                    # TODO: this comment might be false
                    self.updating_container.total_effect[keypoint_count] += \
                        self.updating_container.effect[element_id][keypoint_count]

                except IndexError:
                    print("The mapping data is probably invalid.")
                    print("Please check the \'arduino_mapping.txt\' input whether the given dofs are correct or not.")
                    raise IndexError

        self.updating_container.effect_ratio = deepcopy(self.updating_container.effect)
        for i in range(number_of_elements):
            for j in range(number_of_keypoints):
                if self.updating_container.total_effect[j] > 0:
                    self.updating_container.effect_ratio[i][j] = \
                        abs(self.updating_container.effect_ratio[i][j] / self.updating_container.total_effect[j])
                else:
                    self.updating_container.effect_ratio[i][j] = 0

        self.updating_container.sorted_effect_sign = [0]*number_of_keypoints

        # Sort by effectiveness
        for i in range(number_of_keypoints):
            self.updating_container.sorted_effect[i] = deepcopy(self.updating_container.effect)
            self.updating_container.sorted_effect_sign[i] = [0]*number_of_elements

            # Check sign of the effect
            for k_temp in range(number_of_elements):
                if self.updating_container.sorted_effect[i][k_temp][i] < 0:
                    for j_temp in range(number_of_keypoints):
                        self.updating_container.sorted_effect[i][k_temp][j_temp] = \
                            abs(self.updating_container.sorted_effect[i][k_temp][j_temp])

                        self.updating_container.sorted_effect_sign[i][k_temp] = -1
                else:
                    self.updating_container.sorted_effect_sign[i][k_temp] = +1

            for j in range(number_of_keypoints):
                if i != j and j != 0:
                    self.updating_container.sorted_effect[i] = \
                        swap_columns(sorted(swap_columns(
                            self.updating_container.sorted_effect[i], 0, j), reverse=True), 0, j)
            if i != 0:
                self.updating_container.sorted_effect[i] = \
                    swap_columns(sorted(swap_columns(
                        self.updating_container.sorted_effect[i], 0, i), reverse=True), 0, i)
            else:
                self.updating_container.sorted_effect[i] = \
                    sorted(self.updating_container.sorted_effect[i], reverse=True)

    # TODO: negative numbers kill the optimization
    def optimize(self, delta):
        """
        Model updating - core function

        :param delta: difference between the latest measurement and the calculated model 
        :return: {Boolean} Shows whether the optimization step successful
        """

        element_id = self.truss.number_of_elements()
        self.updating_container.modifications = [0.0] * self.truss.number_of_elements()

        new_delta = delta

        j = 0

        print("------------------------------------")
        print("Threshold: %.4f" % self.updating_container.error_limit)
        print("Maximal number of steps: %.0f" % self.updating_container.iteration_limit)
        print("Step: 0/%.0f" % self.updating_container.iteration_limit)

        # Optimization loop
        while error(new_delta) > self.updating_container.error_limit and (self.capable() or j <= 1):
            j += 1

            print("Error: %.2f" % error(new_delta))

            print("-----")
            print("Step: %.0f" % j)

            ratio = [0] * element_id
            unit = 0

            previous_modifications = self.updating_container.modifications

            # Loop though the possible modifications
            for index in range(self.truss.number_of_elements()):
                # Calculate the magnitude of a further unit modification on an existing one
                modification = self.updating_container.modifications[index] - self.updating_container.unit_modification

                # Apply the modification considering the limits
                self.updating_container.modifications[index] = \
                    min(abs(modification), self.updating_container.modification_limit) * math.copysign(1, modification)

                # Apply modification on material
                """
                Apply a modification on the index-th structure's index-th parameter
                considering one modifiable parameter on each element
                """
                self.apply_update(self.updating_container.modifications[index], index)

                # Calculate the modified stiffness matrix
                self.updating_container.trusses[index].calculate_stiffness_matrix()

                # Solve the modified structure
                self.updating_container.trusses[index].solve(self.configuration)

            # Evaluate all the modifications
            self.evaluate_updates()

            print('trusses[0] is static!!')
            new_delta = self.difference(self.updating_container.trusses[0].displacements,
                                        self.updating_container.measurement)

            for i, effect in enumerate(self.updating_container.total_effect):
                if effect == 0.0:
                    print("None of the variables has effect on %.0f" % self.updating_container.arduino_mapping[i])
                    print("Model updating has no solution.")
                    self.updating_container.number_of_updates[2] += 1
                    raise Exception

            for i in range(self.truss.number_of_keypoints()):
                modified_node_id = self.updating_container.sorted_effect[0][i][1]

                ratio[modified_node_id] = \
                    abs(self.updating_container.total_effect[0] / self.updating_container.sorted_effect[0][i][0]) * \
                    math.copysign(1, self.updating_container.sorted_effect_sign[i][0])

                unit += abs(ratio[modified_node_id]*self.updating_container.sorted_effect[0][i][0])

            new_delta_string = ""
            for factor in new_delta:
                new_delta_string += " %.2f" % factor

            print("Delta: [%s ]" % new_delta_string)
            scale = new_delta[0]/unit

            for i in range(self.truss.number_of_elements()):
                modified_node_id = self.updating_container.sorted_effect[0][i][1]

                self.updating_container.modifications[modified_node_id] = \
                    min(abs(previous_modifications[modified_node_id] - self.updating_container.unit_modification *
                            ratio[modified_node_id]),
                        self.updating_container.modification_limit) *\
                    math.copysign(1,
                                  previous_modifications[modified_node_id] - self.updating_container.unit_modification *
                                  ratio[modified_node_id])  # the last part is already the sign itself

            print("Tangent: %.2f" % scale)

            if j >= self.updating_container.iteration_limit:
                print("Iteration limit is reached.")
                print("The remaining error is: %.2f" % error(new_delta))
                self.updating_container.number_of_updates[1] += 1
                return False

        print("Final error: %.2f" % error(new_delta))

        if not self.capable() and j > 1:
            print("Optimization could not be finished successfully.")
            print("The remaining error is: %.2f" % error(new_delta))
            self.updating_container.number_of_updates[2] += 1
            return False

        self.write_output_stream(j)
        self.updating_container.number_of_updates[0] += 1
        return True

    def start_model_updating(self, unit_modification=0.05, error_limit=0.9,
                             modification_limit=0.7, iteration_limit=100):
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
            print("The absolute value of the unit modification must be minimum 0.01 and maximum 0.5\nGiven: %s"
                  % unit_modification)
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

        self.update_model()

    def update_model(self):
        """
        General function to manage model updating procedure.
        """
        self.processed_data = [0]*1  # *len(self.arduino_mapping)


        if not self.configuration.simulation:
            base = self.set_base()
            filemode = 'a'
            # TODO: complete this process
        else:
            # self.bulk_set_measurement_points()
            base = ['SIMULATION']
            filemode = 'r'
            self.check_folder('Simulation')

        self.initiate_output_stream()

        with open("./Simulation/" + str(self.title) + ' - Input Data.txt', filemode) as input:
            for i in range(1000):
                delta = []
                new_line = input.readline()

                if not new_line == '':
                    delta = self.mock_delta(new_line)

                if delta:
                    if self.updating_container.original_delta == []:
                        self.updating_container.original_delta = delta

                    self.updating_container.latest_delta = delta
                    self.optimize(delta)
                    self.updating_container.reset(self.truss)


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
    parser.add_argument("-i", "--input", metavar='str', type=str, default="",
                        help="Input file, stored in the ./Structure folder [*.str]")
    args = parser.parse_args()

    if args.input != "":
        input_file = args.input
    else:
        input_file = "bridge"
        print("DEMO MODE")
        print("You started the program in demo mode.\nThe following example is based on the bridge.str file which "
              "is located in the Structures folder.")

    # Define new structure
    TRUSS = Truss(input_file=input_file, title=args.title, compatibility_mode=args.compatibility,
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
        plot(TRUSS.truss, original=True, result=False, supports=True, forces=True, reactions=False,
             label_classes='load')
        plot(TRUSS.truss, original=True, result=True, supports=True, forces=False, reactions=False)
        plot(TRUSS.truss, original=False, result=True, supports=True, forces=True, reactions=True, label_classes='load')
        animate(TRUSS.truss, 'load')
        TRUSS.configuration.part_time("Plotting")

    # Write results to file
    TRUSS.write_results()
    TRUSS.configuration.part_time("Wrote results to the output file: Structures/%s - Results.txt" % TRUSS.title)

    # Update iteration
    if TRUSS.configuration.updating:
        TRUSS.start_model_updating(unit_modification=0.05, error_limit=0.9, modification_limit=0.7, iteration_limit=100)
        print("------------------------------------")
        TRUSS.configuration.part_time("Updating numerical model")

    # Closing Arduino port
    if TRUSS.configuration.arduino:
        TRUSS.configuration.disconnect()

    # End logging
    if TRUSS.configuration.log:
        TRUSS.end_logging()
