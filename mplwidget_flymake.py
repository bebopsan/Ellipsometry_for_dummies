#!/usr/bin/env python3

# Programa simple para aprender a usar qt


import sys
     

import matplotlib
matplotlib.use('Qt4Agg')

from PyQt4 import QtGui, QtCore




from polarization_routines import plot_ellipse
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg \
     import NavigationToolbar2QTAgg as NavigationToolbar

from numpy import pi
import numpy as np

__version__ = '1.0'


def Rotate( element, theta):     
    """Applies a rotation operation on any polarizing element given 
    its Jones matrix representation and a angle of rotation.
    The angle is to be given in radians.  
    
    """
    #Pending to import only the most elevant numpy functions for improving speed
    import numpy as np
    
    assert isinstance(element, np.matrix)
    assert element.shape == (2,2)
    
    rotator = np.matrix([[ np.cos(theta), np.sin(theta)],\
                         [-np.sin(theta), np.cos(theta)]])
    return rotator.transpose()*element*rotator

class MplCanvas(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""
    def __init__(self, figure = 1):

        self.QWP = np.matrix([[ np.exp(1.j*pi/4) , 0],\
                                     [ 0 , np.exp(-1.j*pi/4)]])

        fig_1, fig_2 = plot_ellipse( Rotate(self.QWP,-pi/4)*\
                                                   np.matrix([[np.cos(-3*pi/8)],\
                                                                     [np.sin(-3*pi/8)]],\
                                                                      dtype = 'complex'), \
                                                  show = False, retrieve = True)

        if figure ==1:
            self.fig = fig_1
        elif figure == 2:
            self.fig = fig_2
        else:
            raise ValueError("Oops!  That was no valid number.  Try again...")

        self.fig.hold(False)
    
        FigureCanvas.__init__(self, self.fig)
        FigureCanvas.setSizePolicy(self, 
                QtGui.QSizePolicy.Expanding, 
                QtGui.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def PlotEllipses(self, polAngle, qwpAngle, figure = 1):
        """ Draws the polarization ellipse that results from propagating
             a linear polarized field alligned with angle polAngle trough
             a qwp aligned along qwpAngle"""

            fig_1, fig_2 = plot_ellipse( Rotate(self.QWP, qwpAngle)*\
                                                   np.matrix([[np.cos(polAngle )],\
                                                                     [np.sin(polAngle )]],\
                                                                      dtype = 'complex'), \
                                                  show = False, retrieve = True)

        if figure ==1:
            self.fig = fig_1
        elif figure == 2:
            self.fig = fig_2
        else:
            raise ValueError("Oops!  That was no valid number.  Try again...")
        

class MplWidget(QtGui.QWidget):
        """Widget defined in Qt Designer"""
        def __init__(self, parent = None):
                # initialization of Qt MainWindow widget
                QtGui.QWidget.__init__(self, parent)
                # set the canvas to the Matplotlib widget
                self.canvas_1 = MplCanvas(1)
                self.canvas_2 = MplCanvas(2)
                # create a vertical box layout
                #QSplitter(Qt.Horizontal, self)
                self.hbl = QtGui.QHBoxLayout()
                # add mpl widget to vertical box
                self.hbl.addWidget(self.canvas_2,stretch = 2)
                self.hbl.addWidget(self.canvas_1,stretch = 1)
                # set the layout to th vertical box
                self.setLayout(self.hbl)

        def PlotEllipses(self, polAngle, qwpAngle)

      

