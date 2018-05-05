# -*- coding: utf-8 -*-
"""
Created on April 30 18:49:40 2018

TrussFramework created by Máté Szedlák.
Copyright MIT, Máté Szedlák 2016-2018.
"""
import os
import time
try:
    import serial
except ImportError:
    print("You tried to import \'serial\' without installing \'pySerial\'.")
    print("Please first install pySerial: http://playground.arduino.cc/Interfacing/Python")
    raise Exception('pyserial package not found')
try:
    from truss_graphics import Arrow3D
    from extra_math import invert as invert
    from extra_math import mat_vec_mult as mat_vec_mult
    from extra_math import swap_col as swap_col
except ImportError:
    print("Input data is missing")
    print("Please check the config.py, graphics.py and extra_math.py files.")
    raise Exception('Requirement is missing')


class TrussConfiguration(object):
    """
    Truss configuration file
    """
    def __init__(self, input_file, compatibility_mode=2, simulation=0):
        """
        :param input_file: Input file, stored in the ./Structure folder [*.str]
        :param title: Title of the project
        :param compatibility_mode: 0: User defined, 1: Most information (with numpy),
        2: Maximum compatibility mode, 3: Android
        :param simulation: 0: False, 1: True
        """
        DO_NOT_MODIFY_ON = 1
        DO_NOT_MODIFY_OFF = 0

        self.input_file = input_file
        self.TIC = time.time()
        self.TAC = 0
        self.TOTAL_TIME = 0
        self.previous_time = self.TIC
        self.compatibility_mode = compatibility_mode
        self.simulation = simulation

        if self.compatibility_mode == 0:
            ### User defined ###
            # Modify as needed #
            self.mode_name = "User defined"
            self.log = 1  # Logging time
            self.graphics = 1  # Graphical features
            self.solver = 1  # 0: Basic solver, 1: NumPy solver
            self.OSlib = 1  # Basic OS file features (e.g. file size)
            self.updating = 0  # Model Updating: On/ Off
            self.arduino = 0  # Arduino input: On/Off
            self.debug = 1  # Debugging mode
            self.realistic_simulation = 0  # Wait as long as it was originally. Only valid with _SIMULATION = 1

        elif self.compatibility_mode == 1:
            ### Informative ###
            # DO NOT MODIFY
            self.mode_name = "Informative mode"
            self.log = DO_NOT_MODIFY_ON
            self.graphics = DO_NOT_MODIFY_ON
            self.solver = DO_NOT_MODIFY_ON
            self.OSlib = DO_NOT_MODIFY_ON
            self.updating = DO_NOT_MODIFY_ON
            self.arduino = DO_NOT_MODIFY_ON
            self.debug = DO_NOT_MODIFY_OFF
            self.realistic_simulation = DO_NOT_MODIFY_ON

        elif self.compatibility_mode == 2:
            ### Maximum compatibility mode ###
            # DO NOT MODIFY
            self.mode_name = "Maximum compatibility mode"
            self.log = DO_NOT_MODIFY_OFF
            self.graphics = DO_NOT_MODIFY_OFF
            self.solver = DO_NOT_MODIFY_OFF
            self.OSlib = DO_NOT_MODIFY_OFF
            self.updating = DO_NOT_MODIFY_OFF
            self.arduino = DO_NOT_MODIFY_OFF
            self.debug = DO_NOT_MODIFY_OFF
            self.realistic_simulation = DO_NOT_MODIFY_ON

        elif self.compatibility_mode == 3:
            ### Android mode ###
            # DO NOT MODIFY
            self.mode_name = "Android"
            self.log = DO_NOT_MODIFY_ON
            self.graphics = DO_NOT_MODIFY_OFF
            self.solver = DO_NOT_MODIFY_OFF
            self.OSlib = DO_NOT_MODIFY_ON
            self.updating = DO_NOT_MODIFY_OFF
            self.arduino = DO_NOT_MODIFY_OFF
            self.debug = DO_NOT_MODIFY_OFF
            self.realistic_simulation = DO_NOT_MODIFY_ON

        if self.OSlib:
            try:
                if not os.path.exists("./Structures/"):
                    os.makedirs("./Structures")
            except FileNotFoundError:
                print("Error: Please manually create the 'Structures' folder in the root of the project")

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
            print('- Solver is set to default')
        elif self.solver == 1:
            print('- Solver is set to NumPy')
        else:
            raise Exception("Solver settings are invalid!")
        if self.updating:
            print('+ Model updating is turned ON')
            if self.simulation:
                print('Input data is SIMULATED!')
        else:
            print('- Model updating is turned OFF')
        print('------------------------------------')

    def end_logging(self, truss):
        self.TAC = time.time()
        self.TOTAL_TIME = self.TAC - self.TIC
        if self.updating:
            print("Update statistics:")
            print("Totally updated models: " + str(
                truss.num_of_updates[0] + truss.num_of_updates[1] + truss.num_of_updates[2]))
            print("  Successfully updated models: " + str(truss.num_of_updates[0]))
            print("  Updates with running out of possibilities: " + str(truss.num_of_updates[2]))
            print("  Updates did not finished: " + str(truss.num_of_updates[1]))
        # print('Total time: ' + str("{:10.3f}".format(self.TOTAL_TIME)))
        print("Total time: {:10.3f}".format(self.TOTAL_TIME))

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
        self.title = title              # Name of structure
        self.configuration = TrussConfiguration(input_file, compatibility_mode, simulation)
        # Truss data
        self.known_f_a = []           # Nodes without supports
        self.known_f_not_zero = []     # Nodes with loads
        self.DOF = 3                  # Truss's degree of freedom
        self.nodal_connections = []                # Element's end nodes
        self.constraint = []          # Supports
        self.force = []               # Force
        self.nodal_coord = []          # Coordinates of nodes
        self.nodal_coord_def = []      # Coordinates after deformations
        self.cross_sectional_area_list = []                # Cross-sectional areas
        self.elastic_modulo = []                   # Material data
        self.number_of_nodes = 0              # Number of nodes
        self.number_of_elements = 0               # Number of elements
        self.element_DOF = []              # Mapping between DOF and node
        self.stiffness = []           # Global stiffness matrix
        self.mod_stiffnesses = []     # Modified stiffnesses in a hyper-matrix
        self.element_length = []                   # Length of the elements
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
        self.known_displacement_a = []
        self.modified_displacements = []
        self._stiff_is_fresh = 0
        self._post_processed = 0
        self._init_displacement = []
        self.stress_color = []         # Color mapping for stresses
        self.known_dis_a = []
        self.stress = []              # Element's stresses
        self._io_origin = 0           # Array's first element number during IO. Default is 0.
        self.analysis = {}
        self._mod_stiff_is_fresh = 0
        self.mod_displacements = []
        self.keypoint = []
        self.number_of_keypoints = 0
        self.effect = []
        self.total_effect = []
        self.sorted_effect = []
        self.special_DOF_input_string = ''
        self.threshold = 0.1
        self.effect_ratio = []
        self.processed_data = []          # To store last input line
        self.modifications = []          # Storing modifications for model updating
        self.read_elements = [0]*9
        self.arduino_mapping = []
        self.error_limit = 0.5
        self.modification_limit = 0.6
        self.unit_modification = 0.05
        self.measurement = [0.]
        self.number_of_updates = [0, 0, 0]        # [#Successfully updated model, #Updates with overflow exit,
        # #Updates where there were no more modification option]
        self.iteration_limit = 20

    def _open_serial(port, baudrate):
        try:
            connection = serial.Serial()
            connection.baudrate = baudrate
            connection.port = port

            if connection.connection.is_open:
                connection.connection.close()

            connection.connection.open()
            return connection
        except serial.SerialException:
            print(str(port) + ' port is busy. It might be occupied by this program or another one :'
                              '/ Be careful or try resetting this program')
            return False

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
            SUPPORTS - Selected DOF + Prescribed displacement: 0, 0.0 | 1, 0.0 ...
            SPECDOF - Selected node's DOF will be analysed during Model Updating: 1, xyz | 3 y | 10 xz ...

            EOF - For compatibility reasons EOF should be placed after the commands
        """
        self._io_origin = 0
        _read_element_names = ["Origin", "DOF", "Elements", "Coordinates",
                               "Cross-sections", "Materials", "Forces", "Supports", "Measured DOFs"]

        try:
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
                        self.set_model_DOF(int(source_line))
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
                        self.bulk_set_elements(input_number)
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
                        if self.DOF == 3:
                            input_number = [[float(x[0]), float(x[1]), float(x[2])] for x in input_string]
                        elif self.DOF == 2:
                            input_number = [[float(x[0]), float(x[1]), 0.] for x in input_string]
                        self.bulk_set_coordinates(input_number)
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
                        self.bulk_set_cross_sections(input_number)
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
                        self.bulk_set_materials(input_number)
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
                        self.bulk_set_forces(sorted(input_number))
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
                        self.bulk_set_supports(sorted(input_number))
                        self.read_elements[7] = 1

                    if source_line.upper() == "MEASUREMENTS":
                        source_line = sourcefile.readline().strip()
                        self.special_DOF_input_string = source_line
                        input_string = []
                        self.arduino_mapping = source_line.split(',')
                        self.bulk_set_measurement_points(self.arduino_mapping)
                        self.read_elements[8] = 1
        except IOError:
            print("The following file could not be opened: " + "./Structures/" + self.title + ".str")
            print("Please make sure that the structural data is available for the program in the running directory.")
            raise IOError

        terminate = False
        for i, value in enumerate(self.read_elements):
            if i > 0 and (i < 8 or self.configuration.updating):  # if i > 0:
                if value == 0:
                    print("The following was not found: " + _read_element_names[i])
                    terminate = True
        if terminate:
            raise Exception

    def plot(self, show_original, show_result, show_supports, show_forces,
             show_reactions, scale_displacement, scale_force, scale_Z, save_plot):
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

        if self._post_processed == 0:
            print('Postprocess is needed before plotting structure!')
        else:
            if scale_displacement == 0:
                scale_displacement = 1.0           # Scale drwaing of displacements
            if scale_force == 0:
                scale_force = 1.0                  # Scale force sign
            if scale_Z == 0:
                scale_Z = 0.3                      # Scale z-axis

            Arrow3D.plotstructure(self, show_original, show_result, show_supports, show_forces, show_reactions,
                                  scale_displacement, scale_force, scale_Z, _showvalues, save_plot, self.configuration.log)

    def readarduino(self, base, save_input):
        """
        Read data from Arduino
        """
        # Read data from Arduino

        maxdifference = 0.8          # Maximum input difference threshold in mm
        arduinovalues = []
        data = [0.]*len(self.arduino_mapping)
        newdata = False
        bigdifference = False
        readerror = False

        try:
            arduinoline = self.serial_connection.readline()
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

                            if abs(data[i] - self.processed_data[i]) > maxdifference:
                                bigdifference = True
                            if abs(float(arduinovalues[i])) < 2.0:
                                readerror = True

                        self.processed_data = data
                        newdata = True

                    except ValueError:
                        print("Value error... continuing")
                        self.serial_connection.flushInput()
                        time.sleep(0.5)
                    except Exception:
                        print("Type error: " + str(arduinovalues) + "... continuing")
                        self.serial_connection.flushInput()
                        time.sleep(0.5)

            self.serial_connection.flushInput()

        except serial.SerialTimeoutException:
            print("Data could not be read... continuing")
            self.serial_connection.flushInput()
            time.sleep(0.5)

        if newdata and not bigdifference and not readerror:
            self.measurement = zip(self.arduino_mapping, data)

            save_input.write(str(data) + ', ' + str(time.time()) + "\n")
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

    def simulate_arduino(self, arduinoline, prevline):
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
                    if self.configuration.realistic_simulation:
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
                    self.processed_data = data
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
            # pass
        except Exception:
            print("Exception in simulation data")
            # pass

    def calibrate(self):
        """
        Calibration for Arduino measurement.
        All the measurements will describe the dispalecements from the claibration-state.
        """
        answer_1 = '0'
        restart = '0'
        accept = '0'
        arduinovalues = []
        print("Before starting the model updating, the measuring tools must be calibrated.")
        print("The calibration should be done in load-free state.")
        while answer_1 not in ['Y', 'N']:
            answer_1 = input('Can we start the calibration? (y/n) ').upper()
        if answer_1 == 'N':
            self.serial_connection.close()
            raise Exception('Calibration is terminated')
        else:
            try:
                self.serial_connection.flushInput()
                # time.sleep(0.2)
                arduinoline = ''  # self.serial_connection.readline()
                while len(arduinoline) == 0:
                    time.sleep(0.2)
                    arduinoline = self.serial_connection.readline()

                if len(arduinoline) > 0:
                    arduinovalues = arduinoline.split(',')
                    del arduinovalues[len(arduinovalues)-1]               # if needed!!!
                    if len(arduinovalues) == len(self.arduino_mapping):
                        measurement = zip(self.arduino_mapping, arduinovalues)
                        print("Calibration result:")
                        print(measurement)
                        while accept not in ['Y', 'N']:
                            accept = input('Ok? Can we start the main part? Put on the loads! (y/n) ').upper()
                            if accept == 'N':
                                restart = 'Y'
                    else:
                        print("Data error. Calibartion is restarting.")
                        print("Arduino values:" + str(arduinovalues))
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

    def write_results(self, file_name):
        """
        Writing results to file.
        """
        out_element = ''
        for i in self.nodal_connections:
            out_element += str(i[0] + self._io_origin) + ', ' + str(i[1] + self._io_origin) + ' | '
        out_coords = ''
        for i in self.nodal_coord:
            out_coords += str(i[0]) + ', ' + str(i[1]) + ', ' + str(i[2]) + ' | '
        out_crsect = ''
        for i in self.cross_sectional_area_list:
            out_crsect += str(i) + ', '
        out_materials = ''
        for i in self.elastic_modulo:
            out_materials += str(i) + ', '
        out_forces = ''
        for forcedof in self.known_f_not_zero:
            if self.DOF == 3:
                out_forces += str(forcedof + self._io_origin) + ', ' + str(self.force[forcedof]) + ' | '
            elif self.DOF == 2 and i % 3 != 2:
                out_forces += str(forcedof - forcedof//3 + self._io_origin) + ', ' + str(self.force[forcedof]) + ' | '
        out_supports = ''
        for i in self.constraint:
            if self.DOF == 3:
                out_supports += str(i[0] + self._io_origin) + ', ' + str(i[1]) + ' | '
            elif i[0] % 3 != 2:
                out_supports += str(i[0] - i[0]//3 + self._io_origin) + ', ' + str(i[1]) + ' | '
        # Not elegant solution
        out_specdofs = self.special_DOF_input_string
        try:
            with open("./Structures/" + file_name, 'w') as outfile:
                # Writing data
                outfile.write('Calculation of \'' + self.title + '\':\n\n')

                outfile.write('Reactions\n')
                # for i in range(len(self.force)//3):
                prev = -1
                for i in self.known_dis_a:
                    if self.DOF == 3 or i % 3 != 2:
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
                            if self.DOF != 2 and (i//3)*3+2 in self.known_dis_a:
                                nodalforce += "{:10.2f}".format(self.force[(i//3)*3+2]) + '\n'
                            else:
                                nodalforce += '          \n'
                            if nodalforce != '                                  \n':
                                outfile.write(str(i//3 + self._io_origin) + ', ' + nodalforce)
                        prev = i//3
                outfile.write('\n')

                outfile.write('Displacements\n')
                for i in range(len(self.displacement)//3):
                    if i < 100:
                        outfile.write(' ')
                        if i < 9:
                            outfile.write(' ')
                    outfile.write(str(i + self._io_origin) + ', ' + "{:10.3f}".format(self.displacement[i*3 + 0]) +
                                  ', ' + "{:10.3f}".format(self.displacement[i*3 + 1]) +
                                  ', ' + "{:10.3f}".format(self.displacement[i*3 + 2]) + ', ' + '\n')
                outfile.write('\n')

                outfile.write('Stresses\n')
                for i, stress in enumerate(self.stress):
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
                outfile.write(str(self.DOF) + '\n\n')
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

    def end_logging(self):
        self.configuration.end_logging(self)
