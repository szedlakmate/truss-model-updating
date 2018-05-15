# -*- coding: utf-8 -*-
"""
Created on April 30 18:49:40 2018

TrussFramework created by Máté Szedlák.
Copyright MIT, Máté Szedlák 2016-2018.
"""
import os
import time
import math
import itertools
from copy import deepcopy
try:
    import numpy
except ImportError:
    print("NumPy could not be loaded")
try:
    from extra_math import invert
    from extra_math import mat_vec_mult as multiply_matrix_vector
    from extra_math import swap_col as swap_columns
except ImportError:
    print("Input is missing")
    print("Please check the extra_math.py file.")
try:
    import serial
except ImportError:
    print("You tried to import \'serial\' without installing \'pySerial\'.")
    print("Please first install pySerial: http://playground.arduino.cc/Interfacing/Python")
    print("Arduino input can not be used")


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


class PlotConfiguration(object):
    def __init__(self):
        self.labels = []
        self.paths = []


class TrussConfiguration(object):
    """
    Truss configuration file
    """

    def __init__(self, input_file, compatibility_mode=2, simulation=0):
        """
        :param input_file: Input file, stored in the ./Structure folder [*.str]
        :param compatibility_mode: 0: User defined, 1: Most information (with numpy),
        2: Maximum compatibility mode, 3: Android
        :param simulation: 0: False, 1: True
        """
        static_ON = 1
        static_OFF = 0

        self.input_file = input_file
        self.TIC = time.time()
        self.TAC = 0
        self.TOTAL_TIME = 0
        self.previous_time = self.TIC
        self.compatibility_mode = compatibility_mode
        self.simulation = simulation
        self.updating = {'error_limit': 0.5, 'modification_limit': 0.6,
                         'unit_modification': 0.05, 'iteration_limit': 20}

        if self.compatibility_mode == 0:
            ### User defined ###
            # Modify as needed #
            self.mode_name = "User defined"
            self.log = 1  # Logging time
            self.graphics = 1  # Graphical features
            self.solver = 1  # 0: Basic solver, 1: NumPy solver
            self.OSlib = 1  # Basic OS file features (e.g. file size)
            self.updating = 1  # Model Updating: On/ Off
            self.arduino = 0  # Arduino input: On/Off
            self.debug = 1  # Debugging mode
            self.realistic_simulation = 0  # Wait as long as it was originally. Only valid with _SIMULATION = 1

        elif self.compatibility_mode == 1:
            ### Informative ###
            # DO NOT MODIFY
            self.mode_name = "Informative mode"
            self.log = static_ON
            self.graphics = static_ON
            self.solver = static_ON
            self.OSlib = static_ON
            self.updating = static_ON
            self.arduino = static_ON
            self.debug = static_OFF
            self.realistic_simulation = static_ON

        elif self.compatibility_mode == 2:
            ### Maximum compatibility mode ###
            # DO NOT MODIFY
            self.mode_name = "Maximum compatibility mode"
            self.log = static_OFF
            self.graphics = static_OFF
            self.solver = static_OFF
            self.OSlib = static_OFF
            self.updating = static_OFF
            self.arduino = static_OFF
            self.debug = static_OFF
            self.realistic_simulation = static_ON

        elif self.compatibility_mode == 3:
            ### Android mode ###
            # DO NOT MODIFY
            self.mode_name = "Android"
            self.log = static_ON
            self.graphics = static_OFF
            self.solver = static_OFF
            self.OSlib = static_ON
            self.updating = static_OFF
            self.arduino = static_OFF
            self.debug = static_OFF
            self.realistic_simulation = static_ON

        """
        if self.OSlib:
            try:
                if not os.path.exists("./Structures/"):
                    os.makedirs("./Structures")
            except FileNotFoundError:
                print("Error: Please manually create the 'Structures' folder in the root of the project")
        """

        if self.simulation or not self.updating:
            self.arduino = 0

        if self.compatibility_mode == 2:
            os.chdir(os.path.dirname(os.path.abspath(__file__)))

    def start_logging(self):
        self.part_time("Initialization")

        print('------------------------------------')
        print('Truss calculation program')
        print('Created by Máté Szedlák (https://github.com/szedlakmate/truss-model-updating)')
        print('Compatibility mode: ' + self.mode_name)
        if self.solver == 0:
            print('  Solver is set to default')
        elif self.solver == 1:
            print('  Solver is set to NumPy')
        else:
            raise Exception("Solver settings are invalid!")
        if self.updating:
            print('+ Model updating is turned ON')
            if self.simulation:
                print('Input data is SIMULATED!')
        else:
            print('- Model updating is turned OFF')
        print('------------------------------------')

    def end_logging(self, number_of_updates):
        """
        Solver log + Total time

        :param number_of_updates: Update statistic array
        :return: None
        """

        self.TAC = time.time()
        self.TOTAL_TIME = self.TAC - self.TIC
        if self.updating:
            print("Update statistics:")
            print("Totally updated models: " + str(
                number_of_updates[0] + number_of_updates[1] + number_of_updates[2]))
            print("  Successfully updated models: " + str(number_of_updates[0]))
            print("  Updates with running out of possibilities: " + str(number_of_updates[2]))
            print("  Updates did not finished: " + str(number_of_updates[1]))
            print("\nTotal time: {:10.3f}".format(self.TOTAL_TIME))

        else:
            print("\nTotal time: {:10.3f}".format(self.TOTAL_TIME))

    def part_time(self, message):
        """
        Calculating and printing the time consumption of tasks

        Should be called with the previously saved part-time and the name of the actual task
        At the first call, should be called with TIC value. The input argument
        should be overwritten by this function's return value.

        :param message: message will be printed
        :return: returns actual time to start new lap (sub-timer)
        """
        new_time = time.time()
        print(message)
        print('Time: ' + str("{:10.3f}".format(new_time - self.previous_time)))
        print('------------------------------------')
        self.previous_time = new_time

class TrussModelData(object):
    """
    Data model for structures
    """
    def __init__(self, title):
        # Labels
        self.analysis = {}
        self.title = title
        # Truss data
        self.DOF = 3
        self.keypoints = []
        self.nodal_connections = []  # Element's end nodes
        self.constraint = []  # Supports
        self.force = []  # Force
        self.nodal_coord = []  # Coordinates of nodes
        self.cross_sectional_area_list = []  # Cross-sectional areas
        self.elastic_modulo = []  # Material data
        #Additional data
        self.known_f_a = []  # Nodes without supports
        self.known_f_not_zero = []  # Nodes with loads
        self.known_displacement_a = []
        self.known_f_not_zero = []
        self.stiffness_matrix = []  # Global stiffness matrix
        #self.element_length = []  # Length of the elements
        self.element_DOF = []
        self._norm_stiff = []                  # E/L
        self._cx = []
        self._cy = []
        self._cz = []
        self._s_loc = []
        self.dis_new = []
        self.force_new = []
        self.stiff_new = []
        self._stiff_is_fresh = 0
        self._post_processed = 0
        self._init_displacement = []
        # Results
        self.nodal_coord_def = []  # Coordinates after deformations
        # Secondary variables
        self.displacements = []  # Relative displacements
        self.stress = []  # Element's stresses
        self.stress_color = []  # Color mapping for stresses
        # Plot data
        self.plot_data = PlotConfiguration()

    def number_of_nodes(self):
        # TODO: refactor
        return len(set(list(itertools.chain.from_iterable(sorted(self.nodal_connections)))))

    def number_of_elements(self):
        return len(self.nodal_connections)

    def number_of_keypoints(self):
        return len(self.keypoints)

    def element_length(self, ID, target='original'):
        if str(target.lower()) == 'original':
            return math.sqrt(sum([(j-i)**2 for j, i in zip(self.nodal_coord[self.nodal_connections[ID][1]],
                                                           self.nodal_coord[self.nodal_connections[ID][0]])]))
        else:
            return math.sqrt(sum([(j - i) ** 2 for j, i in zip(self.nodal_coord_def[self.nodal_connections[ID][1]],
                                                               self.nodal_coord_def[self.nodal_connections[ID][0]])]))

    def set_postprocess_needed_flag(self):
        """
        Set portprocess required flag

        :return: None
        """
        self._post_processed = 0

    def invalidate_stiffness_matrices(self):
        """
        Set flags for stiffness matrix recalculation

        :return: None
        """
        self._stiff_is_fresh = 0
        self.set_postprocess_needed_flag()

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
        self.invalidate_stiffness_matrices()


    def bulk_set_measurement_points(self, measurement_points, io_origin):
        """
        Set special nodal DOFs
        """
        for location in measurement_points:
            node = int(location[0:len(location)-1])-io_origin
            if 'X' in location:
                self.analysis[location] = node*3+0
                self.keypoints.append(node*3+0)
            if 'Y' in location:
                self.analysis[location] = node*3+1
                self.keypoints.append(node*3+1)

            if 'Z' in location:
                if self.DOF == 3:
                    self.analysis[location] = node*3+2
                    self.keypoints.append(node*3+2)
                else:
                    print("Z-direction is not allowed in 2D structures. "
                          "Please check the 'MEASUREMENTS' section in the input file.")
                    raise Exception

        if self.number_of_keypoints == 0:
            print("There is no valid measured DOF. Please check the \'MEASUREMENTS\' section in the input file.")

    # TODO: This functions is probably not called although it should be part of the set_base() process
    def bulk_set_elements(self, nodal_connection_list):
        """
        Setting elements (nodal connections) in bulk mode

        :param nodal_connection_list: an array of arrays, where each sub-array is a pair of integers, namely [i, j].
        i, j are the ID's of the nodes and an element is i -> j.
        :return: None
        """
        # Set attribute
        self.nodal_connections = nodal_connection_list

        # Set freshness flags after geometry modification
        self.invalidate_stiffness_matrices()

        # Creating mapping tool for elements
        for node in self.nodal_connections:
            self.element_DOF.append(
                [node[0] * 3, node[0] * 3 + 1, node[0] * 3 + 2, node[1] * 3, node[1] * 3 + 1, node[1] * 3 + 2])

        # Initializing defaults for all matrices
        self._init_displacement = [0] * (3 * self.number_of_nodes())
        self.force = [0.] * (3 * self.number_of_nodes())
        self.stiffness_matrix = [0.] * (3 * self.number_of_nodes())
        self.known_f_a = []
        self.known_f_not_zero = []

    def bulk_set_coordinates(self, coordinate_list):
        """
        Setting coordinates in bulk mode

        :param coordinate_list: An array of coordinate arrays (2/3 elements) enlisting ALL the nodal coordinates
        :return: None
        """
        # Set attribute
        self.nodal_coord = coordinate_list

        # Set freshness flags after geometry modification
        self.invalidate_stiffness_matrices()

        # Validity check
        # TODO: the two side comes fromthe same value. The check should reflect to toher variables
        # if self.number_of_nodes() > len(self.nodal_coord):
        #    raise Exception('More coordinates are needed')
        # elif not self.nodal_connections:
        #    raise Exception('Nodes must be set before defining elements')
        # self._check_coordinates(False)

    def bulk_set_cross_sections(self, area_list):
        """
        Setting cross-sections in bulk mode

        :param area_list: cross-sectional data array according to the elements
        :return: None
        """
        self.cross_sectional_area_list = area_list
        self.invalidate_stiffness_matrices()

    def bulk_set_materials(self, E):
        """
        Setting material data in bulk mode

        :param E: array of elastic modulos according to the elements
        :return: None
        """
        self.elastic_modulo = E
        self.invalidate_stiffness_matrices()

    def bulk_set_forces(self, forces):
        """
        Set forces

        :param forces: matrix of forces in the following pattern: [[location, force], ...]
        :return: None
        """
        # DOF dependent mapping
        for location, force in forces:
            if self.DOF == 3:
                self.force[location] = force
            elif self.DOF == 2:
                self.force[location + (location // 2)] = force
        self.set_postprocess_needed_flag()

    def bulk_set_supports(self, constraints):
        """
        Set supports

        :param constraints: matrix of constraints in the following pattern: [[location, constraint], ...]
        :return:
        """
        for location, constraint in constraints:
            if self.DOF == 3:
                self.constraint.append([location, constraint])
            elif self.DOF == 2:
                self.constraint.append([location + (location // 2), constraint])
        self.invalidate_stiffness_matrices()

        self.set_postprocess_needed_flag()

    def calculate_stiffness_matrix(self):
        """
        Stiffness matrix compilation
        """
        self.set_postprocess_needed_flag()

        if self.DOF == 2:
            for dof_z in range(self.number_of_nodes()):
                self.constraint.append([int(dof_z*3+2), 0.])
        self.constraint = list(k for k, _ in itertools.groupby(sorted(self.constraint)))

        # Setting known forces
        known_f_a = []
        known_f_not_zero = []
        known_displacement_a = []
        for location in range(3*self.number_of_nodes()):
            known_f_a.append(location)
            if self.force[location] != 0:
                known_f_not_zero.append(location)


        for constraint in self.constraint:
            self._init_displacement[constraint[0]] = constraint[1]
            known_displacement_a.append(constraint[0])
            try:
                known_f_a.remove(constraint[0])
                known_f_not_zero.remove(constraint[0])
            except ValueError:
                pass

        self.known_f_a = known_f_a
        self.known_displacement_a = known_displacement_a
        self.known_f_not_zero = known_f_not_zero

        elements_lengths = [0.]*self.number_of_elements()
        self._norm_stiff = [0.]*self.number_of_elements()
        self._cx = [0.]*self.number_of_elements()
        self._cy = [0.]*self.number_of_elements()
        self._cz = [0.]*self.number_of_elements()
        self._s_loc = [0.]*self.number_of_elements()
        local_stiffness_matrix = [0.]*self.number_of_elements()
        self.stress = [0.]*self.number_of_elements()

        self.stiffness_matrix = [[0.]*(self.number_of_nodes()*3)]*(self.number_of_nodes()*3)

        for i in range(self.number_of_elements()):
            elements_lengths[i] = \
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

            stiffness_increment = [0.]*(self.number_of_nodes()*3)

            for j in range(3*2):
                for k in range(3*2):
                    stiffness_increment[ele_dof_vec[k]] = local_stiffness_matrix[i][j][k]
                self.stiffness_matrix[ele_dof_vec[j]] = [x + y for x, y in zip(self.stiffness_matrix[ele_dof_vec[j]], stiffness_increment)]
        self._stiff_is_fresh = 1

    def solve(self, configuration):
        """
        Main solver of the code
        """
        if self._stiff_is_fresh == 0:
            if configuration.log:
                print('Stiffness matrix is recalculated')
            self.calculate_stiffness_matrix()

        self.dis_new = [0.]*(self.number_of_nodes()*3-len(self.constraint))
        self.force_new = [0.]*(self.number_of_nodes()*3-len(self.constraint))
        self.stiff_new = [[0.]*(self.number_of_nodes()*3-len(self.constraint))]*(self.number_of_nodes()*3-len(self.constraint))

        # known force array
        for i, DOF in enumerate(self.known_f_a):
            try:
                self.force_new[i] = self.force[DOF]
            except Exception:
                pass

        stiffness_increment = [0.]*(self.number_of_nodes()*3-len(self.constraint))
        for i, kfai in enumerate(self.known_f_a):
            for j, kfaj in enumerate(self.known_f_a):
                stiffness_increment[j] = self.stiffness_matrix[kfai][kfaj]
            self.stiff_new[i] = [x + y for x, y in zip(self.stiff_new[i], stiffness_increment)]

        # SOLVING THE STRUCTURE
        if configuration.solver == 0:
            if configuration.log:
                print('Built-in solver')
            self.dis_new = multiply_matrix_vector(invert(self.stiff_new), self.force_new)
        else:
            if configuration.log:
                print('NumPy solver')
            self.dis_new = numpy.linalg.solve(numpy.array(self.stiff_new), numpy.array(self.force_new))

        self.displacements = deepcopy(self._init_displacement)

        for i, known_f_a in enumerate(self.known_f_a):
            self.displacements[known_f_a] = self.dis_new[i]

        # Deformed shape
        self.nodal_coord_def = []
        for i in range(self.number_of_nodes()):
            self.nodal_coord_def.append([self.nodal_coord[i][0] + self.displacements[i * 3 + 0],
                                         self.nodal_coord[i][1] + self.displacements[i * 3 + 1], self.nodal_coord[i][2] + self.displacements[i * 3 + 2]])

        # Postrpocesses
        self.post_process()

    def post_process(self):
        """
        Calculates reaction forces and stresses

        :return None
        """
        self.calculate_reactions()
        self.calculate_stresses()
        self._post_processed = 1

    def calculate_reactions(self):
        """
        Calculates reaction forces
        """
        for i in self.known_displacement_a:
            self.force[i] = 0
            for j, displacement in enumerate(self.displacements):
                self.force[i] += self.stiffness_matrix[i][j]*displacement


    def calculate_stresses(self):
        """
        Calculates stress in elements

        Last part: Coloring elements for graphical output
        """
        self.stress = [0.]*self.number_of_elements()

        for element in range(self.number_of_elements()):
            self.stress[element] = \
                (self.element_length(element, 'deformed') - self.element_length(element, 'original')) / self.element_length(element, 'original') * self.elastic_modulo[element]

        s_max = max([abs(min(self.stress)), max(self.stress), 0.000000001])

        # TODO: Should be written out using a function
        self.stress_color = [float(x)/float(s_max) for x in self.stress]


class ModelUpdatingContainer(object):
    """
    Model updating container for isolating the updating parts
    """
    def __init__(self):
        self.modifications = []          # Storing modifications for model updating
        self.read_elements = [0]*9
        self.arduino_mapping = []
        self.measurement = [0.]
        self.number_of_updates = [0, 0, 0]        # [#Successfully updated model, #Updates with overflow exit,
        #self.modified_displacements = []
        self.effect = []
        self.effect_ratio = []
        self.total_effect = []
        self.sorted_effect = []
        self.original_delta = []
        self.latest_delta = []
        # Solver configuration
        self.unit_modification = 0
        self.error_limit = 0
        self.modification_limit = 0
        self.iteration_limit = 0
        #self._mod_stiff_is_fresh = 0
        self.trusses = []


class TrussFramework(object):
    """
    Base class for truss computing, collecting several functions:
    - File operations (input/output)
    - Serial communication
    - Arduino management
    - Plot functions
    """

    def __init__(self, input_file, title, compatibility_mode, simulation):
        """
        Truss Updater framework object.
        :param title: Title of the structure
        """
        # Serial connection
        self.serial_connection = False  # Holds serial connection

        # Project data
        if title == '':
            title = input_file
        self.title = title.replace('.str', '')  # Name of structure
        self.configuration = TrussConfiguration(input_file.replace('.str', '') + '.str', compatibility_mode, simulation)

        # Structure
        self.truss = TrussModelData(self.title)

        # Additional truss data
        self.read_elements = [0.] * 9
        # Model Updating related variables (?)
        self.updating_container = ModelUpdatingContainer()

    def _open_serial(port, baudrate):
        try:
            connection = serial.Serial()
            connection.baudrate = baudrate
            connection.port = port

            if connection.is_open:
                connection.close()

            connection.open()
            return connection
        except serial.SerialException:
            print(str(port) + ' port is busy. It might be occupied by this program or another one :'
                              '/ Be careful or try resetting this program')
            return None

    def connect(self, port='', baudrate=9600):
        self.serial_connection = False
        PORTS = ['COM1', 'COM2', 'COM3']  # List of possible communication ports
        port_number = 0

        if port != '':
            print('FORCED: Opening serial at port ' + str(port))
            self.serial_connection = self._open_serial(port, baudrate)
        else:
            while not self.serial_connection:
                port_number += 1
                if port_number > len(PORTS):
                    port_number = 1
                time.sleep(0.6)
                print('Opening serial at port ' + str(PORTS[port_number]))
                self.serial_connection = self._open_serial(port, baudrate)

    def disconnect(self):
        try:
            self.serial_connection.close()
            return True
        except serial.SerialException:
            print('Serial connection cannot be closed')
            return False
        except AttributeError:
            print('Serial connection cannot be closed because it is not managed by this thread')
            return False

    # TODO: This functions is probably not called although it should be part of the set_base() process

    def read(self):
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
            SUPPORTS - Selected DOF + Prescribed displacements: 0, 0.0 | 1, 0.0 ...
            SPECDOF - Selected node's DOF will be analysed during Model Updating: 1, xyz | 3 y | 10 xz ...

            EOF - For compatibility reasons EOF should be placed after the commands
        """
        self._io_origin = 0
        _read_element_names = ["Origin", "DOF", "Elements", "Coordinates",
                               "Cross-sections", "Materials", "Forces", "Supports", "Measured DOFs"]


        try:
            self.check_folder('Structures')
            with open("./Structures/" + self.configuration.input_file, "r") as sourcefile:
                source_line = ""
                while source_line != "EOF":
                    source_line = sourcefile.readline().strip()

                    if source_line.upper() == "_ORIGIN":
                        source_line = sourcefile.readline().strip()
                        self._io_origin = int(source_line)
                        self.read_elements[0] = 1

                    if source_line.upper() == "DOF":
                        source_line = sourcefile.readline().strip()
                        self.truss.set_model_DOF(int(source_line))
                        self.read_elements[1] = 1

                    if source_line.upper() == "ELEMENTS":
                        source_line = sourcefile.readline().strip()
                        input_string = []
                        input_number = []
                        input_string = [x.split(',') for x in source_line.split('|')]
                        if len(input_string[0]) == 1:
                            input_string = [x.split(';') for x in source_line.split('|')]
                        if [''] in input_string:
                            input_string.remove([''])
                        input_number = [[int(x[0]) - self._io_origin, int(x[1]) - self._io_origin] for x in input_string]
                        self.truss.bulk_set_elements(input_number)
                        self.read_elements[2] = 1

                    if source_line.upper() == "COORDINATES":
                        source_line = sourcefile.readline().strip()
                        input_string = []
                        input_number = []
                        input_string = [x.split(',') for x in source_line.split('|')]
                        if len(input_string[0]) == 1:
                            input_string = [x.split(';') for x in source_line.split('|')]
                        if [''] in input_string:
                            input_string.remove([''])
                        if len(input_string[0]) == 3:
                            input_number = [[float(x[0]), float(x[1]), float(x[2])] for x in input_string]
                        elif len(input_string[0]) == 2:
                            input_number = [[float(x[0]), float(x[1]), 0.] for x in input_string]
                        self.truss.bulk_set_coordinates(input_number)
                        self.read_elements[3] = 1

                    if source_line.upper() == "CROSS-SECTIONS":
                        source_line = sourcefile.readline().strip()
                        input_string = []
                        input_number = []
                        input_string = source_line.split(',')
                        if len(input_string) == 1:
                            input_string = source_line.split(';')
                        if '' in input_string:
                            input_string.remove('')
                        input_number = [float(eval(x)) for x in input_string]
                        self.truss.bulk_set_cross_sections(input_number)
                        self.read_elements[4] = 1

                    if source_line.upper() == "MATERIALS":
                        source_line = sourcefile.readline().strip()
                        input_string = []
                        input_number = []
                        input_string = source_line.split(',')
                        if len(input_string) == 1:
                            input_string = source_line.split(';')
                        if '' in input_string:
                            input_string.remove('')
                        input_number = [float(eval(x)) for x in input_string]
                        self.truss.bulk_set_materials(input_number)
                        self.read_elements[5] = 1

                    if source_line.upper() == "FORCES":
                        source_line = sourcefile.readline().strip()
                        input_string = []
                        input_number = []
                        input_string = [x.split(',') for x in source_line.split('|')]
                        if len(input_string[0]) == 1:
                            input_string = [x.split(';') for x in source_line.split('|')]
                        if [''] in input_string:
                            input_string.remove([''])
                        input_number = [[int(x[0]) - self._io_origin, float(x[1])] for x in input_string]
                        self.truss.bulk_set_forces(sorted(input_number))
                        self.read_elements[6] = 1

                    if source_line.upper() == "SUPPORTS":
                        source_line = sourcefile.readline().strip()
                        input_string = []
                        input_number = []
                        input_string = [x.split(',') for x in source_line.split('|')]
                        if len(input_string[0]) == 1:
                            input_string = [x.split(';') for x in source_line.split('|')]
                        if [''] in input_string:
                            input_string.remove([''])
                        input_number = [[int(x[0]) - self._io_origin, float(x[1])] for x in input_string]
                        self.truss.bulk_set_supports(sorted(input_number))
                        self.read_elements[7] = 1

                    if source_line.upper() == "MEASUREMENTS":
                        source_line = sourcefile.readline().strip()
                        self.special_DOF_input_string = source_line
                        input_string = []
                        self.updating_container.arduino_mapping = source_line.split(',')
                        self.truss.bulk_set_measurement_points(self.updating_container.arduino_mapping, self._io_origin)
                        self.read_elements[8] = 1
        except IOError:
            print("The following file could not be opened: " + "./Structures/" + self.configuration.input_file)
            print("Please make sure that the structural data is available for the program in the run directory.")
            raise IOError

        terminate = False
        for i, value in enumerate(self.read_elements):
            if i > 0 and (i < 8 or self.configuration.updating):  # if i > 0:
                if value == 0:
                    print("The following was not found: " + _read_element_names[i])
                    terminate = True
        if terminate:
            raise Exception

    def readarduino(self, base, save_input):
        """
        Read data from Arduino
        """
        # Read data from Arduino

        maxdifference = 0.8          # Maximum input difference threshold in mm
        arduino_values = []
        data = [0.]*len(self.updating_container.arduino_mapping)
        newdata = False
        bigdifference = False
        readerror = False

        try:
            arduino_line = self.serial_connection.readline()
            if len(arduino_line) > 0:
                arduino_values = arduino_line.split(',')
                try:
                    if arduino_values[0][len(arduino_values)-1] == '.':
                        arduino_values[0] = arduino_values[0][:len(arduino_values[0])-2]
                    else:
                        del arduino_values[len(arduino_values)-1]
                except IndexError:
                    print("Index Error... continuing")
                if len(arduino_values) == len(self.updating_container.arduino_mapping):
                    try:
                        for i in range(len(self.updating_container.arduino_mapping)):
                            data[i] = float(arduino_values[i]) - float(base[i][1])

                            if abs(data[i] - self.processed_data[i]) > maxdifference:
                                bigdifference = True
                            if abs(float(arduino_values[i])) < 2.0:
                                readerror = True

                        self.processed_data = data
                        newdata = True

                    except ValueError:
                        print("Value error... continuing")
                        self.serial_connection.flushInput()
                        time.sleep(0.5)
                    except Exception:
                        print("Type error: " + str(arduino_values) + "... continuing")
                        self.serial_connection.flushInput()
                        time.sleep(0.5)

            self.serial_connection.flushInput()

        except serial.SerialTimeoutException:
            print("Data could not be read... continuing")
            self.serial_connection.flushInput()
            time.sleep(0.5)

        if newdata and not bigdifference and not readerror:
            # self.updating_container.measurement = zip(self.updating_container.arduino_mapping, data)
            self.updating_container.measurement = data

            save_input.write(str(data) + ', ' + str(time.time()) + "\n")
            # Calculate differences
            delta = self.difference(self.truss.displacements, self.updating_container.measurement)

            print("Delta: " + str(delta))

            newdata = False
            bigdifference = False
            readerror = False
            return delta

        newdata = False
        bigdifference = False
        readerror = False


    def difference(self, num_displ, measurement):
        """
        Calculate the difference between the Real-life measurement and the Numerical solution.

        :param num_displ:
        :param measurement: [[13X, -2.154], [16Y, 5.256], ...]
        :return: [measurement - calculated] Negative means higher deflections are measured
        """
        # print(nodenumber option should be added! <XXX>)
        delta = []

        for measured_position in measurement:
            #loc, measured = package[0], package[1]
            loc = '1Y'
            try:
                dof = self.truss.analysis[loc.upper()]
            except KeyError:
                print('The given measurement location cannot be matched with the input data.')
                print('The available nodes are: {\'NAMES\': mapping addresses}')
                print(self.truss.analysis)
                #self.configuration.disconnect()
                raise Exception('Watchpoint name error')
            delta_appendix = float(measured_position) - float(num_displ[dof])
            try:
                delta.append(delta_appendix)
            except Exception:
                pass

        return delta


    def mock_delta(self, arduino_line, previous_line):
        """
        Simulate data, based on previous measurement
        """
        arduino_values = []
        data = [0.]*1 #*len(self.updating_container.arduino_mapping)
        try:
            arduino_line = str(arduino_line.split(']')[0])+"]"
            try:
                arduino_values = eval(arduino_line)
            except SyntaxError:
                print('Error: "' + str(arduino_line) + '"')
            try:
                for i in range(len(self.updating_container.arduino_mapping)):
                    data[i] = float(arduino_values[i])
                self.processed_data = data
            except Exception:
                print("Type error: " + str(arduino_values) + "... continuing")

            self.updating_container.measurement = data

            # Calculate differences
            delta = self.difference(self.truss.displacements, self.updating_container.measurement)
            return delta

        except IndexError:
            print("IndexError")
            # pass
        except IndexError:
            print("Exception in simulation data")
            # pass

    def calibrate(self):
        # TODO: Seriously needs refactoring
        """
        Calibration for Arduino measurement.
        All the measurements will describe the dispalecements from the claibration-state.
        """
        answer_1 = '0'
        restart = '0'
        accept = '0'
        arduino_values = []
        print("Before starting the model updating, the measuring tools must be calibrated.")
        print("The calibration should be done in load-free state.")
        while answer_1 not in ['Y', 'N']:
            answer_1 = input('Can we start the calibration? (y/n) ').upper()
        if answer_1 == 'N':
            try:
                self.serial_connection.close()
            except AttributeError:
                print("Connection is not captured by this thread")
            raise Exception('Calibration is terminated')
        else:
            try:
                self.serial_connection.flushInput()
                # time.sleep(0.2)
                arduino_line = ''  # self.serial_connection.readline()
                while len(arduino_line) == 0:
                    time.sleep(0.2)
                    arduino_line = self.serial_connection.readline()

                if len(arduino_line) > 0:
                    arduino_values = arduino_line.split(',')
                    del arduino_values[len(arduino_values)-1]               # if needed!!!
                    if len(arduino_values) == len(self.updating_container.arduino_mapping):
                        measurement = zip(self.updating_container.arduino_mapping, arduino_values)
                        print("Calibration result:")
                        print(measurement)
                        while accept not in ['Y', 'N']:
                            accept = input('Ok? Can we start the main part? Put on the loads! (y/n) ').upper()
                            if accept == 'N':
                                restart = 'Y'
                    else:
                        print("Data error. Calibartion is restarting.")
                        print("Arduino values:" + str(arduino_values))
                        restart = 'Y'
                else:
                    print('The calibration cannot be done: no data')
                    while restart not in ['Y', 'N']:
                        restart = input('Do you want to restart calibration? (y/n) ').upper()
            except Exception:
                print('The calibration cannot be done: exception was raised')
                while restart not in ['Y', 'N']:
                    restart = input('Do you want to restart calibration? (y/n) ').upper()

        if restart == 'Y':
            print("Restarting calibration")
            self.serial_connection.flushInput()
            return self.calibrate()
        elif restart == 'N':
            self.serial_connection.close()
            raise Exception('Calibration is terminated')
        if accept == 'Y':
            return measurement

    def open_simulation_thread(self):
        self.check_folder('Simulation')
        with open("./Simulation/" + str(self.title) + ' - Input Data.txt', 'r') as simulation_file:
            source_line = ""
            while source_line != "EOF":
                source_line = simulation_file.readline().strip()
                print(source_line)


                if source_line.upper() == "_ORIGIN":
                    source_line = simulation_file.readline().strip()
                    self._io_origin = int(source_line)
                    self.read_elements[0] = 1


        return simulation_file

    def set_base(self):

        if self.configuration.simulation:
            #self.arduino_simulation_thread = self.open_simulation_thread()
            #self.configuration.previous_line

            print("HARD CODED PART 3")
            self.configuration.updating.base = 27

            return 0 #self.mock_delta()
        else:
            return self.calibrate()


    def check_folder(self, directory):
        path = str(os.path.dirname('.')) + '/' + directory.replace('/','').replace('.','') + '/'

        if not os.path.exists(path):
            os.makedirs(path)


    def write_results(self):
        """
        Writing results to file.
        """
        self.check_folder('Results')
        path ="Results/" + self.title + ' - Results.txt'


        out_element = ''
        for i in self.truss.nodal_connections:
            out_element += str(i[0] + self._io_origin) + ', ' + str(i[1] + self._io_origin) + ' | '
        out_coords = ''
        for i in self.truss.nodal_coord:
            out_coords += str(i[0]) + ', ' + str(i[1]) + ', ' + str(i[2]) + ' | '
        out_crsect = ''
        for i in self.truss.cross_sectional_area_list:
            out_crsect += str(i) + ', '
        out_materials = ''
        for i in self.truss.elastic_modulo:
            out_materials += str(i) + ', '
        out_forces = ''
        for forcedof, i in enumerate(self.truss.known_f_not_zero):
            if self.truss.DOF == 3:
                out_forces += str(forcedof + self._io_origin) + ', ' + str(self.truss.force[forcedof]) + ' | '
            elif self.truss.DOF == 2 and i % 3 != 2:
                out_forces += str(forcedof - forcedof//3 + self._io_origin) + ', ' + str(self.truss.force[forcedof]) + ' | '
        out_supports = ''
        for i in self.truss.constraint:
            if self.truss.DOF == 3:
                out_supports += str(i[0] + self._io_origin) + ', ' + str(i[1]) + ' | '
            elif i[0] % 3 != 2:
                out_supports += str(i[0] - i[0]//3 + self._io_origin) + ', ' + str(i[1]) + ' | '
        # Not elegant solution
        out_specdofs = self.special_DOF_input_string
        try:
            with open(path, 'w') as outfile:
                # Writing data
                outfile.write('Calculation of \'' + self.title + '\':\n\n')

                outfile.write('Reactions\n')
                # for i in range(len(self.truss.force)//3):
                prev = -1
                for i in self.truss.known_displacement_a:
                    if self.truss.DOF == 3 or i % 3 != 2:
                        if i//3 != prev:
                            if i < 100:
                                outfile.write(' ')
                                if i < 9:
                                    outfile.write(' ')
                            nodalforce = ''
                            if (i//3)*3+0 in self.truss.known_displacement_a:
                                nodalforce += "{:10.2f}".format(self.truss.force[(i//3)*3+0]) + ', '
                            else:
                                nodalforce += '            '
                            if (i//3)*3+1 in self.truss.known_displacement_a:
                                nodalforce += "{:10.2f}".format(self.truss.force[(i//3)*3+1]) + ', '
                            else:
                                nodalforce += '            '
                            if self.truss.DOF != 2 and (i//3)*3+2 in self.truss.known_displacement_a:
                                nodalforce += "{:10.2f}".format(self.truss.force[(i//3)*3+2]) + '\n'
                            else:
                                nodalforce += '          \n'
                            if nodalforce != '                                  \n':
                                outfile.write(str(i//3 + self._io_origin) + ', ' + nodalforce)
                        prev = i//3
                outfile.write('\n')

                outfile.write('Displacements\n')
                for i in range(len(self.truss.displacements) // 3):
                    if i < 100:
                        outfile.write(' ')
                        if i < 9:
                            outfile.write(' ')
                    outfile.write(str(i + self._io_origin) + ', ' + "{:10.3f}".format(self.truss.displacements[i * 3 + 0]) +
                                  ', ' + "{:10.3f}".format(self.truss.displacements[i * 3 + 1]) +
                                  ', ' + "{:10.3f}".format(self.truss.displacements[i * 3 + 2]) + ', ' + '\n')
                outfile.write('\n')

                outfile.write('Stresses\n')
                for i, stress in enumerate(self.truss.stress):
                    if i < 100:
                        outfile.write(' ')
                        if i < 9:
                            outfile.write(' ')
                    outfile.write(str(i + self._io_origin) + ', ' + "{:10.3f}".format(stress) + '\n')
                outfile.write('\n')

                # Saving original input
                outfile.write('----- Original input: -----\n\n')
                outfile.write('_ORIGIN\n')
                outfile.write(str(self._io_origin) + '\n\n')
                outfile.write('DOF\n')
                outfile.write(str(self.truss.DOF) + '\n\n')
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
            print("Error: Please manually create the 'Structures' folder"
                  "in the root of the project")
            print('Missing:\n'+ str(path))
            raise FileNotFoundError

    def write_output_stream(self, j, appendix):
        self.check_folder('Results')
        with open("./Results/" + self.title + ' - UpdateResults' + appendix + '.txt', 'a') as outfile:
            if j > 1:
                if j <= self.updating_container.iteration_limit and self.capable():
                    self.updating_container.number_of_updates[0] += 1
                    outfile.write("Update state: SUCCESSFUL\n")
                if not j <= self.updating_container.iteration_limit:
                    self.updating_container.number_of_updates[1] += 1
                    outfile.write("Update state: Run out of iteration limit\n")
                if not self.capable() and j > 1:
                    self.updating_container.number_of_updates[2] += 1
                    outfile.write("Update state: No more possible modification\n")
            else:
                outfile.write("Update state: Optimization was skipped\n")
            outfile.write("Required iterations: " + str(j) + "\n")
            outfile.write("Measurement: " + str(self.updating_container.measurement) + "\n")
            outfile.write("Original delta: " + str(self.updating_container.original_delta) + "\n")
            outfile.write("New delta: " + str(self.updating_container.latest_delta) + " (limit: " + str(self.updating_container.error_limit) + ")\n")
            outfile.write("Final error: " + str(error(self.updating_container.latest_delta)) + "\n")
            outfile.write("Modifications [%]: \n")
            outfile.write(str(self.updating_container.modifications) + "\n")
            outfile.write("Original displacements: \n")
            outfile.write(str(self.truss.displacements) + "\n")
            if j > 1:
                outfile.write("New displacements: \n")
                outfile.write(str(self.updating_container.modified_displacements[self.number_of_elements()]) + "\n")
            outfile.write("----------------------\n")

    def end_logging(self):
        self.configuration.end_logging(self.updating_container.number_of_updates)


def plot(truss, original=False, result=False, supports=True, forces=True,
         reactions=True, scale_displacement=1.0, scale_forces=1.0, z_correction=0.3, save=True, label_classes=''):
    """
    Plot function of the Truss class
    This method calls the more general plotstructure() method.

    :param truss: Truss object to be plotted
    :param original: Print original structure?
    :param result: Print calculated structure?
    :param supports: Print supports?
    :param forces: Print forces?
    :param reactions: Print reactions?
    :param scale_displacement: Sign graphical scalign
    :param scale_forces: Sign graphical scalign
    :param z_correction: Scale everything along the Z-axis
    :param save: Save the plot into a file?
    :param label_classes: HTML styled class definition
    :return: None
    """
    # Import is located here due to Android compatibility issues
    from truss_graphics import Arrow3D

    force_values = True     # Show values of forces

    file_path = Arrow3D.plotstructure(truss, original, result, supports, forces, reactions,
                                      scale_displacement, scale_forces, z_correction, force_values, save)

    # Add labels of plots for optional GIF creating
    path_array = label_classes.split(' ')
    for label in path_array:
        truss.plot_data.paths += [file_path]
        truss.plot_data.labels += [label]

def animate(truss, label_selector):
    """
    GIF creator

    :param file_paths: Array, pointing to the source files
    :param gif_name: name of the resulted gif. It will be created in the Results folder.
    :return: None
    """
    import imageio

    file_paths = []

    for (path, applied_label) in zip(truss.plot_data.paths, truss.plot_data.labels):
        if applied_label == label_selector:
            file_paths += [path]

    images = []
    for filename in file_paths:
        images.append(imageio.imread(filename))
    imageio.mimsave('./Results/' + truss.title + ' - ' + label_selector + ' animation.gif', images, duration=0.33*len(file_paths))
