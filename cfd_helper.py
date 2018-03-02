# CFD Helper
"""
Written by Matt Blomquist
Last Update: 2017-12-08 (YYYY-MM-DD)

This python module includes helper functions for the input_parameters.py script
that automatically generates input files for the CFD program.

"""

import numpy as np


class Boundary:
      """A class structure for boundaries used for computational fluid dynamics
      simulations.
      
      Definition of terms:
        type:
          inlet - An inlet type represents a boundary where fluid enters the 
                  problem geometry.
          outlet - An outlet type represents a boundary where fluid exits the 
                  problem geometry.
          wall - A wall type represents a boundary where fluid cannot enter or 
                  exit the boundary and a no-slip condition is in place.
          pressure - A pressure type represents a boundary where a constant, 
                  perscribed pressure condition is applied. 
          symmetry - A symmetric boundary represents a condition where flow 
                  does not cross the boundary and there is no flux.
          
          """

      def __init__(self):

            # Set default type as wall
            self.type = 'wall'

            # Set default perscribed conditions
            self.pressure = None
            self.velocity = np.array([0,0,0]) # Velocity x, y, z = 0 for wall conditions.
            self.temperature = None

            return None

      def set_inlet(self, Re, rho, mu, L, flow_direction, temperature=None):
            """
            Input parameters for setting inlet boundary:
                  Re - Reynolds Number
                  rho - fluid density
                  mu - fluid viscosity
                  L - characteristic length
                  flow_direction - direction normal to boundary
            """

            # Set type to inlet
            self.type = "inlet"

            # Set prescribed pressure to None
            self.pressure = None

            # Set velocity from Reynolds Number for appropriate direction
            if flow_direction == 'x':
                  self.velocity = np.array([Re*mu/rho/L, 0, 0])
            elif flow_direction == 'y':
                  self.velocity = np.array([0, Re*mu/rho/L, 0])
            elif flow_direction == 'z':
                  self.velocity = np.array([0,0, Re*mu/rho/L])

            # Set Temperature if defined
            if temperature != None:
                  self.temperature = temperature

            return None

      def set_outlet(self, Re, rho, mu, L, flow_direction, temperature=None):
            """
            Input parameters for setting outlet boundary:
                  Re - Reynolds Number
                  rho - fluid density
                  mu - fluid viscosity
                  L - characteristic length
                  flow_direction - direction normal to boundary (x, y, z)
            """

            # Set type to outlet
            self.type = "outlet"

            # Set prescribed pressure to None
            self.pressure = None

            # Set velocity from Reynolds Number for appropriate direction
            if flow_direction == 'x':
                  self.velocity = np.array([Re*mu/rho/L, 0, 0])
            elif flow_direction == 'y':
                  self.velocity = np.array([0, Re*mu/rho/L, 0])
            elif flow_direction == 'z':
                  self.velocity = np.array([0,0, Re*mu/rho/L])

            # Set Temperature if defined
            if temperature != None:
                  self.temperature = temperature

            return None

      def set_wall(self, temperature=None):

            self.__init__()

            # Set Temperature if defined
            if temperature != None:
                  self.temperature = temperature

            return None

      def set_pressure(self,pressure, temperature=None):
            """
            Input parameters for setting pressure boundary:
                  Re - Reynolds Number
                  rho - fluid density
                  mu - fluid viscosity
                  L - characteristic length
                  flow_direction - direction normal to boundary (x, y, z)
            """
            # Set type as pressure
            self.type = 'pressure'

            # Set pressure
            self.pressure = pressure

            # Set velocity to None
            self.velocity = None

            # Set Temperature if defined
            if temperature != None:
                  self.temperature = temperature

            return None

      def set_symmmetry(self):
            
            # Set type as symmetry
            self.type = 'symmmetry'

            # Set flux conditions to None
            self.pressure = None
            self.velocity = None
            self.temperature = None

            return None
