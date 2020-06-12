#!/usr/bin/env python3

# Programa simple para aprender a usar Qt


from __future__ import with_statement

import sys

import matplotlib

matplotlib.use('Qt5Agg')
from PyQt5 import QtCore, QtGui, QtWidgets
from ellipse_plot import Ui_MplMainWindow

from polarization_routines import plot_ellipse, getAnglesFromEllipse, getAnglesFromJones

from matplotlib.figure import Figure

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

from numpy import pi

import numpy as np


__version__ = '1.0'

try:
    _encoding = QtWidgets.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtCore.QCoreApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtCore.QCoreApplication.translate(context, text, disambig)

class MainWindow(QtWidgets.QMainWindow, Ui_MplMainWindow):
    def __init__(self, parent = None ):
        super(MainWindow, self).__init__(parent)
        self.setupUi(self)
        self.setWindowTitle("Assistant for polarization transformations")
        self.action_About.triggered.connect(self.about)
        self.action_Close.triggered.connect(self.fileQuit)
        self.action_Documentation.triggered.connect(self.get_to_doc)
         # Generate plots in case the button is clicked
        self.generateButton.clicked.connect(self.GeneratePlots)
        self.jonesGroupBox.clicked.connect(self.jonesChecked)
        self.ellipticityGroupBox.clicked.connect(self.ellipticityChecked)
        self.jonesTranslatePushButton.clicked.connect(self.TranslateJones)
        self.ellipticityTranslatePushButton.clicked.connect(self.TranslateEllipse)
        
    def TranslateJones(self):
        """ This function takes the values of jones vector representation,
         translates them into the ellipse and lab representations and updates those fields.  """

        psi =  eval(self.PsiLineEdit.text())
        delta = eval(self.deltaLineEdit.text())
        angles = getAnglesFromJones(psi, delta )
        self.qwpLineEdit.setText(_translate("MplMainWindow", str(angles['alpha']), None))
        self.polLineEdit.setText(_translate("MplMainWindow", str(angles['psi']), None))
        self.directionLineEdit.setText(_translate("MplMainWindow", str(angles['phi']), None))
        self.eangleLineEdit.setText(_translate("MplMainWindow", str(angles['theta']), None))
        
    def TranslateEllipse(self):
        """ This function takes the values of ellipse parameters representation,
         translates them into the jones vector and lab representations and updates those fields.  """
        
        phi =  eval(self.directionLineEdit.text())
        theta = eval(self.eangleLineEdit.text())
        angles =  getAnglesFromEllipse(phi, theta)
        self.qwpLineEdit.setText(_translate("MplMainWindow", str(angles['alpha']), None))
        self.polLineEdit.setText(_translate("MplMainWindow", str(angles['psi']), None))
        self.PsiLineEdit.setText(_translate("MplMainWindow", str(angles['psi']), None))
        self.deltaLineEdit.setText(_translate("MplMainWindow", str(angles['delta']), None))
        
             
    def jonesChecked(self):
        self.ellipticityGroupBox.setChecked(False)
    def ellipticityChecked(self):
        self.jonesGroupBox.setChecked(False)
        
        
    def GeneratePlots(self):
        """ Reads values of angles and sends them to the mplWidget """
        polAngle = eval(self.polLineEdit.text())
        qwpAngle = eval(self.qwpLineEdit.text())
        self.mpl_1.plot_ellipses(polAngle, qwpAngle)

    def fileQuit(self):
        self.close()

    def keyPressEvent(self, event):
        key = event.key()
        if key == QtCore.Qt.Key_Return:
            self.GeneratePlots()
            
    def about(self):
        QtWidgets.QMessageBox.about(self, "About",
     """
    <b>PolStaPlot 
    <br>Polarization State Plotter
    <br>Platform for plotting and transforming polarization ellipses
    <p>Copyright &copy; 2013-2015 Universidad EAFIT
    <br>Licensed under the terms of the MIT License
    <p>Created by Santiago Echeverri
    <br>Developed and maintained by the 
    <a href="http://www.eafit.edu.co/investigacion/grupos/escuela-ciencias-humanidades/optica-aplicada/Paginas/optica-aplicada.aspx#.VMgNycZuZhQ">Applied Optics Group at EAFIT University</a>
    <p>For further info contact Santiago Echeverri at:
    <ul><li> sechev14@eafit.edu.co 
    </li></ul>
    <p>This program is intended for plotting the polarization ellipse of
    a state of polarization given by an arbitrary Jones vector. It is also
    a tool for easily identifying  the correct translation between different
    representations of polarization states.
    <br>
    <p>For details about use and installation please visit the 
    <a href="https://github.com/bebopsan/Ellipsometry_for_dummies.git">Ellipsometry for Dummies repository</a>.
    """ )
    def get_to_doc(self):
        from PyQt5.QtGui import QDesktopServices
        from PyQt5.QtCore import QUrl

        url = QUrl()
        url.setUrl('https://github.com/bebopsan/Ellipsometry_for_dummies.git')
        QDesktopServices.openUrl(url)
                
if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    frame = MainWindow()
    frame.show()
    sys.exit(app.exec_())


