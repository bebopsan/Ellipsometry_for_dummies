#!/usr/bin/env python3

# Programa simple para aprender a usar Qt


from __future__ import with_statement
import sys
import matplotlib
matplotlib.use('Qt4Agg')
from PyQt4 import QtGui, QtCore
from propagator import Ui_MainWindow
from polarization_routines import plot_ellipse, getAnglesFromEllipse, getAnglesFromJones
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from numpy import pi, sqrt
import numpy as np

__version__ = '1.0'

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)
    
class MainWindow(QtGui.QMainWindow, Ui_MainWindow):
    def __init__(self, parent = None ):
        super(MainWindow, self).__init__(parent)
        self.setupUi(self)
        self.setWindowTitle("Polarization stat propagator")
        self.action_About.triggered.connect(self.about)
        self.action_Close.triggered.connect(self.fileQuit)
        
        # Initialize some global variables and generate default plots
        self.Jones_M = np.zeros((2,2),dtype='complex')
        self.in_Jones = np.matrix([[0.0],[0.0]], dtype='complex')
        self.out_Jones = np.matrix([[0.0],[0.0]], dtype='complex')
        self.Generate_in_Jones()
        self.Propagate()
        
        # Generate plots in case the propagate button is clicked
        self.propagate.clicked.connect(self.Propagate)
        
        # Synchronize sliders
        self.directionSlider.valueChanged.connect(self.directionSliderChanged)
        self.eangleSlider.valueChanged.connect(self.eangleSliderChanged)
        
        #Synchronize lineEdits
        self.directionLineEdit.editingFinished.connect(self.directionLineEditChanged)
        self.eangleLineEdit.editingFinished.connect(self.eangleLineEditChanged)
        
    def directionLineEditChanged(self):
        """ Changes value of ellipse direction if the lineEdit changes """
        if isinstance(eval(self.directionLineEdit.text()),float):
            self.Generate_in_Jones()
        
        self.Propagate()
    def eangleLineEditChanged(self):
        """ Changes value of ellipticity angle if the lineEdit changes """
        if isinstance(eval(self.eangleLineEdit.text()),float):
            self.Generate_in_Jones()
        self.Propagate()
        
    def directionSliderChanged(self):
        """ Changes value of ellipse direction if the slider changes """
        value = str(self.directionSlider.value()/100.*np.pi-np.pi/2.)
        self.directionLineEdit.setText(_translate("MainWindow", value, None))        
        self.Generate_in_Jones()
        if self.Auto_Propagate_checkBox.isChecked() == True:
                self.Propagate()
    def eangleSliderChanged(self):
        """ Changes value of ellipticity angle if the slider changes """
        value = str(self.eangleSlider.value()/100.*np.pi/2-np.pi/4.)
        self.eangleLineEdit.setText(_translate("MainWindow", value, None))
        self.Generate_in_Jones()
        if self.Auto_Propagate_checkBox.isChecked() == True:
                self.Propagate()
    
    def Build_raw_M(self):
        X,Y,Z,W = eval(self.raw_MatrixLineEdit.text())
        print(X,Y,Z,W)
        self.M_x1_LineEdit.setText(_translate("MplMainWindow", str(X), None))
        self.M_x2_LineEdit.setText(_translate("MplMainWindow", str(X), None))        
        self.M_y1_LineEdit.setText(_translate("MplMainWindow", str(Y), None))
        self.M_y2_LineEdit.setText(_translate("MplMainWindow", str(Y), None))        
        self.M_z1_LineEdit.setText(_translate("MplMainWindow", str(Z), None))
        self.M_z2_LineEdit.setText(_translate("MplMainWindow", str(Z), None))        
        self.M_w1_LineEdit.setText(_translate("MplMainWindow", str(W), None))
        self.M_w2_LineEdit.setText(_translate("MplMainWindow", str(W), None))        

    def Build_M(self):
        """" Assembles a Jones matrix from its real elements x,y,z,w """
        M = {}
        M['X1'] = eval(self.M_x1_LineEdit.text())
        M['X2'] = eval(self.M_x2_LineEdit.text())
        M['Y1'] = eval(self.M_y1_LineEdit.text())
        M['Y2'] = eval(self.M_y2_LineEdit.text())
        M['Z1'] = eval(self.M_z1_LineEdit.text())
        M['Z2'] = eval(self.M_z2_LineEdit.text())
        M['W1'] = eval(self.M_w1_LineEdit.text())
        M['W2'] = eval(self.M_w2_LineEdit.text())
        for i in M.keys(): 
            if isinstance(M[i],int):
                M[i] = float(M[i])
            else: assert isinstance(M[i],float)
            
        self.Jones_M = np.matrix([ [ M['X1'] + 1.j*M['Y1'], M['Z1'] + 1.j*M['W1'] ],\
                                 [ -M['Z2'] + 1.j*M['W2'], M['X2'] - 1.j*M['Y2'] ]],dtype = 'complex')
        
        
        #print(self.Jones_M)

    def TranslateEllipse(self):
        """ This function takes the values of ellipse parameters representation,
         translates them into the jones vector and lab representations and updates those fields.           """
        
        phi =  eval(self.directionLineEdit.text())
        theta = eval(self.eangleLineEdit.text())
        angles =  getAnglesFromEllipse(phi, theta)

        #Convert these to globa variables.
        return angles['psi'], angles['alpha'], angles['Jones']
        
    def Generate_in_Jones(self):
        self.polAngle, self.qwpAngle, self.in_Jones = self.TranslateEllipse()
        self.mpl_1.plot_ellipses(self.polAngle, self.qwpAngle)
        
    def Propagate(self):
        """ Reads values of angles and sends them to the mplWidget """
        
        if self.rawMatrix_checkBox.isChecked() == True:
            self.Build_raw_M()
        
        self.Build_M()
        
        self.out_Jones = self.Jones_M*self.in_Jones
        self.mpl_2.plot_ellipses(self.polAngle, self.qwpAngle,Jones=self.out_Jones)
        
        self.Jones_out_Y.setText(_translate("MplMainWindow", str(self.out_Jones[0,0]), None))
        self.Jones_out_X.setText(_translate("MplMainWindow", str(self.out_Jones[1,0]), None))
        

    def fileQuit(self):
        self.close()

#    def keyPressEvent(self, event):
#        key = event.key()
#        modifiers = QtGui.QApplication.keyboardModifiers()
#        if key == QtCore.Qt.Key_X:
#            #self.Generate_in_Jones()
#                
#            if modifiers == QtCore.Qt.ShiftModifier:
#                self.M_x2_LineEdit.setFocus
#            else:
#                self.M_x1_LineEdit.setFocus
            
    def about(self):
        QtGui.QMessageBox.about(self, "About",
    """Platform for propagating polarization ellipses through polarizing elements
    Santiago Echeverri Chacón, 2014
    Grupo de Óptica Aplicada

    Universidad EAFIT
    
    This program was intended for plotting the polarization ellipse of
    a state of polarization given by arbitrary direction and ellipticity.
    It also shows the polarization vector and ellipse resulting from the 
    propagation of the latter trough a polarizing or phase retarding optical 
    element. 

    It may be used and modified with no restriction; raw copies as well as
    modified versions may be distributed without limitation."""
    )

if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    frame = MainWindow()
    frame.show()
    sys.exit(app.exec_())
