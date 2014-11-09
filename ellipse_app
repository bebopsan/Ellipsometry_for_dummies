#!/usr/bin/env python3

# Programa simple para aprender a usar qt


import sys
     
import PySide
     
from PySide.QtCore import QRect
from PySide.QtGui import (QApplication, QMainWindow, QMessageBox,
                          QIcon, QAction, QWidget, QGridLayout,
                          QTextEdit, QMenuBar, QMenu, QStatusBar)

from ellipse_plot import Ui_Dialog
from polarization_routines import plot_ellipse
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from numpy import pi
import numpy as np

__version__ = '1.0'


def Rotate( element, theta):     
    """Applies a rotation operation on any polarizing element given 
    its Jones matrix representation and a angle of rotation.
    The angle is to be given in radians.  
    
    """
    import numpy as np
    
    assert isinstance(element, np.matrix)
    assert element.shape == (2,2)
    
    rotator = np.matrix([[ np.cos(theta), np.sin(theta)],\
                         [-np.sin(theta), np.cos(theta)]])
    return rotator.transpose()*element*rotator


class MainWindow(QMainWindow, UI_Dialog):
        def __init__(self,parent = None):
                super(MainWindow, self).__init__(parent)
                self.setupUi(self)
                self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
                self.setWindowTitle("application main window")
                self.file_menu = QtGui.QMenu('&File', self)
                self.file_menu.addAction('&Quit', self.fileQuit,
                                         QtCore.Qt.CTRL + QtCore.Qt.Key_Q)
                self.menuBar().addMenu(self.file_menu)

                self.help_menu = QtGui.QMenu('&Help', self)
                self.menuBar().addSeparator()
                self.menuBar().addMenu(self.help_menu)

                self.help_menu.addAction('&About', self.about)

                self.main_widget = QtGui.QWidget(self)

                 H = np.matrix([[ 1 , 0],\
                                 [ 0 , 0]])

                QWP = np.matrix([[ np.exp(1.j*pi/4) , 0],\
                                     [ 0 , np.exp(-1.j*pi/4)]])

                fig, fig2 = plot_ellipse( Rotate(QWP,-pi/4)*np.matrix([[np.cos(-3*pi/8)],[np.sin(-3*pi/8)]],dtype = 'complex'))
        
                canvas = FigureCanvas.__init__(self, fig)

                l = QtGui.QVBoxLayout(self.main_widget)

                l.addWidget(canvas)
                
                self.main_widget.setFocus()
                self.setCentralWidget(self.main_widget)
                self.statusBar().showMessage("All hail matplotlib!", 2000)

                
        def fileQuit(self):
                self.close()
        def closeEvent(self, ce):
                self.fileQuit()
        def about(self):
                QtGui.QMessageBox.about(self, "About",
"""embedding_in_qt4.py example
Copyright 2005 Florent Rougon, 2006 Darren Dale

This program is a simple example of a Qt4 application embedding matplotlib
canvases.

It may be used and modified with no restriction; raw copies as well as
modified versions may be distributed without limitation."""
)
                
if __name__ == '__main__':
    app = QApplication(sys.argv)
    frame = MainWindow()
    frame.show()
    sys.exit(app.exec_())


