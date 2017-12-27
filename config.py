# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 10:06:25 2017

@author: Pragm
"""
######################################################
                   ###  DO NOT MODIFY ###            #
                                                     #
class Configuration():                               #
    def __init__(self, compatibility, simulation = 0):#
        self.simulation = simulation                 #
        if compatibility == 1:                       #
            ### "Processing 3" mode ###              #
            self.mode_name = "Processing 3"          #
            self.log = 1*1                           #
            self.graphics = 0*0                      #
            self.solver = 0*0                        #
            self.OSlib = 0*0                         #
            self.updating = 1*1                      #
            self.arduino = 1*1                       #
            self.debug = 0*0                         #
            self.realistic_simulation = 1*1          #
                                                     #
        elif compatibility == 2:                     #
            ### Android mode ###                     #
            # DO NOT MODIFY                          #
            self.mode_name = "Android"               #
            self.log = 1*1                           #
            self.graphics = 0*0                      #
            self.solver = 0*0                        #
            self.OSlib = 1*1                         #
            self.updating = 0*0                      #
            self.arduino = 0*0                       #
            self.debug = 0*0                         #
            self.realistic_simulation = 1*1          #
                                                     #
        elif compatibility == 3:                     #
            ### Informative ###                      #
            # DO NOT MODIFY                          #
            self.mode_name = "Informative mode"      #
            self.log = 1*1                           #
            self.graphics = 1*1                      #
            self.solver = 1*1                        #
            self.OSlib = 1*1                         #
            self.updating = 1*1                      #
            self.arduino = 1*1                       #
            self.debug = 0*0                         #
            self.realistic_simulation = 1*1          #
                                                     #
        else:                                        #
            ### Maximum compatibility ###            #
            # DO NOT MODIFY                          #
            self.mode_name = "Maximum compatibility" #
            self.log = 0*0                           #
            self.graphics = 0*0                      #
            self.solver = 0*0                        #
            self.OSlib = 0*0                         #
            self.updating = 0*0                      #
            self.arduino = 0*0                       #
            self.debug = 0*0                         #
            self.realistic_simulation = 1*1          #
                                                     #
                   ###  DO NOT MODIFY ###            #
######################################################