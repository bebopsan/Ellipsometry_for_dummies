# -*- coding: utf-8 -*-
"""
Library with simple funcions for learning and better understanding polarization concepts.

It contains:

A function for rotating optical elements such as polarizationizers and waveplates 
defined by a Jones Matrix.

A function for plotting polarizarion ellipses of a given polarization state.

"""

__author__      = "Santiago Echeverri Chacon"
__copyright__   = "Universidad EAFIT"
__license__ = "MIT"
__version__ = "1"
__maintainer__ = "Santiago Echeverri"
__email__ = " santiag77e@gmail.com"
__status__ = "ongoing"

import numpy as np
from numpy import pi
H = np.matrix([[ 1 , 0],\
               [ 0 , 0]])

HWP = np.matrix([[ np.exp(1j*pi/2) , 0],\
                     [ 0 , np.exp(-1j*pi/2)]], dtype = 'complex')

QWP = np.matrix([[ np.exp(1.j*pi/4) , 0],\
                 [ 0 , np.exp(-1.j*pi/4)]])    

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

def plot_ellipse( J, name = '', show = True, retrieve = False):
    """ Given a Jones vector representation of a wave
    this routine plots its correspondent polarization ellipse.

    Use a vector in the form:

    J =  numpy.matrix([ [ 1 ],  \
                                    [ 0] ])
    
    """
    
    #import sympy as sy
    import matplotlib as mpl
    from mpl_toolkits.mplot3d import Axes3D
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec

    
    assert isinstance(J, np.matrix)
    assert J.shape == (2,1)
    

    "Convert elements of the vector to complex numbers"
    J_aux = np.zeros((2,1), dtype = complex)
    for i in [0,1]:
        J_aux[i] = J[i]
        
    J = J_aux

    #plt.xkcd()
    
    mpl.rcParams['legend.fontsize'] = 10

    fig = plt.figure()
    fig2 = plt.figure()
    
    gs1 = gridspec.GridSpec(1, 1)

    axes = [ ]
    axes.append(fig.add_subplot( gs1[0,0], projection='3d'))
    axes.append(fig2.add_subplot( gs1[0,0] ))

    theta = np.linspace(-4 * np.pi, 4 * np.pi, 100)

    Ex = np.abs( J[0])
    Ey = np.abs( J[1])

    phi1 = np.angle( J[0] )
    phi2 =  np.angle( J[1] )

    X = np.linspace(-2, 2, 100)
    Y = Ex * np.cos(theta + phi1)
    Z = Ey * np.cos(theta + phi2)

    lim = np.sqrt(Ex**2+Ey**2)
       
    axes[0].set_xlim(-2, 2)
    
    axes[0].set_ylim(-1, 1)
    axes[0].set_zlim(-1, 1)
     #axes[0].view_init(elev=0, azim=0)
    axes[0].set_aspect('equal')

    axes[0].set_axis_off()

    #axes[0].set_autoscale_on(False)
    #axes[0].legend()    axes[0].set_yticklabels([])
    axes[0].set_xticklabels([])
    axes[0].set_zticklabels([])
   
    #axes[1].ylim([-1,1])
    #axes[1].xlim([-1,1])
    #axes[1].zlim([-1,1])
    axes[1].axis([-lim,lim,-lim,lim])
    axes[1].set_xlabel('Normalized $x$ component of $E$')
    axes[1].set_ylabel('Normalized $y$ component of $E$')
    axes[1].set_aspect('equal')
    axes[0].plot([-4, 4],[0, 0],zs=[0,0], color='green')
    axes[0].plot([-2, -2],[0, 0],zs=[-2,2], color='green')
    axes[0].plot([-2, -2],[-2, 2],zs=[0,0],color='green')
    #axes[1].plot([-1, 1],[-1, 1], color='blue')
    #axes[1].plot([-1, 1],[-0.41421356237309503, 0.41421356237309503], color='blue')
    axes[1].plot([-2, 2],[0, 0], color='green')
    axes[1].plot([0, 0],[-2, 2], color='green')

    

    axes[0].plot(X, Y, Z, zdir ='Z', color='red')
    axes[1].plot(Y, Z, color='red')
    
    
    if show == True:
        plt.show()
    else:
        plt.close()
    if name != '':
        fig.savefig('trayectory_'+name)
        fig2.savefig('_'+name)
    if retrieve:
        return fig,fig2
    else:
        return "done"

def ellipse_gen(x0,phi, theta, give_delta = False, verbose = False,debug =False):
    """" This function is meant to be used in a minimization routine. 
    The objective is to check if a rotated polarizer followed by a QWP 
    retarder returns a Jones vector that has a polarization ellipse given 
    by the parametrization angles $\phi$ and $\theta$. 
        
    Inputs:
    x0: Array, contains the rotation angles for the polarizer and the QWP: [psi, alpha] 
        psi: Rotation angle of the polarizer
        alpha: Rotation angle of the lam/4 retarder
    phi: Angle in radians between the coordinate system and the mayor axes of the ellipse
    theta: Ellipticity angle of the polarization ellipse. It is function of the ratio between 
           the minor and mayor axes: theta = arctan(b/a).
    
    Outputs:
    result: Scalar, is the result of summing the square differences of polarization angles.
            Each angle (phi_c, theta_c) is extracted from the Jones vector after propagation 
            through optical elements, and compared to the required angle. 
    delta: Optional. Is the phase difference between y and x components of the field described 
           a Jones Vector.
       
    Author: Santiago Echeverri Chacón   
    """

    [psi,alpha] = x0
    
    # Matrices of horizantal oriented retarders to be rotated:
    
    #H = np.matrix([[ 1 , 0],\
    #               [ 0 , 0]], dtype = 'complex')
    
    HWP = np.matrix([[ np.exp(1j*pi/2) , 0],\
                     [ 0 , np.exp(-1j*pi/2)]], dtype = 'complex')
    
    QWP = np.matrix([[   np.exp(1j*pi/4) , 0],\
                     [ 0 ,  np.exp(-1j*pi/4) ]], dtype = 'complex')
    
    ## def Rotate( element, theta):   
    ##     """Applies a rotation on any polarizing element given 
    ##     its Jones matrix representation and a angle of rotation.
    ##     The angle is to be given in radians. The matrix is a 2*2 matrix.  
    ##     """
        
    ##     assert isinstance(element, np.matrix)
    ##     assert element.shape == (2,2)
    
    ##     rotator = np.matrix([[ np.cos(theta), np.sin(theta)],\
    ##                          [-np.sin(theta), np.cos(theta)]],dtype = 'complex')
    ##     return rotator.transpose()*element*rotator

    # Little test...
    #if psi == pi/2:
    #   psi = alpha
    #   alpha = pi/2
    
    
    # Definition of an imput Jones vector oriented along the axis of the polarizer,
    # this way the output vector keeps a normalized intensity.
    # It may be redundant when multiplying by the polarizer but this way I can change 
    # the definition later if I want other kinds of inputs.
    
    In = Rotate(HWP,psi/2)*np.matrix([[1],[0]],dtype = 'complex')
        
    Out = Rotate(QWP, alpha)*In
    
    # Extract the phase differences and field amplitudes from "Out" and calculate phi,theta:
    
    delta1 = np.angle( Out[0,0] )
    delta2 = np.angle( Out[1,0] )
    delta = delta2 - delta1
    
    Ex = abs(Out[0,0])
    Ey = abs(Out[1,0])
        
    phi_c = 1/2*np.arctan2(2*Ex*Ey*np.cos(delta), (Ex**2 - Ey**2))
    theta_c = -1/2*np.arcsin(2*Ex*Ey*np.sin(delta)/(Ex**2 + Ey**2))
    
    result =((phi - phi_c)**2 + (theta - theta_c)**2)

    if verbose:
        print('phi_c: ', phi_c,'phi: ', phi,'theta_c: ', theta_c,'theta: ', theta )
    if give_delta:
        
        if debug:
            return result,delta, Out, In
        else:
            return result, delta
    else:
        if debug:
             return result, Out, In
        else:
            return result
    

def getAnglesFromEllipse(phi, theta, bnds = ((-np.pi, np.pi),(-np.pi, np.pi))):
    """ Minimization procedure for finding polarization angles.
    Based on desired values of ellipse orientation and ellipticity this function uses
    the lbfgs minimization algorithm to find the combination of polarization and
    retarder angles that reproduce it.
    
    Inputs:
    phi: Angle in radians between the coordinate system and the mayor axes of the ellipse
    theta: Ellipticity angle of the polarization ellipse. It is function of the ratio between 
           the minor and mayor axes: theta = arctan(b/a).
    
    Outputs:
    dict: Dictionary that contains:
         psi: angle for input polarizer
         alpha: angle for input quarter waveplate retarder
         theta: ellipticity angle
         x: Jones vector horizontal component
         y: Jones vector vertical component
         
    Author: Santiago Echeverri Chacón
    """
    
    #import numpy as np
    from numpy import matrix
    from numpy.linalg import norm
    from scipy.optimize import fmin_l_bfgs_b
    ## if phi == np.pi/2:
    ##     psi, alpha = 0,0#0, 67.5*np.pi/180
    ## else:
        ## psi, alpha = 0.0
    psi, alpha = 0.0, 0.0
    res = fmin_l_bfgs_b(ellipse_gen,[ psi, alpha], fprime=None,approx_grad=1,\
                                    args = (phi, theta), pgtol=1e-05,bounds =bnds,epsilon=1e-10)
    psi, alpha = res[0]
    minimum, delta = ellipse_gen(res[0], phi, theta, give_delta = True, verbose = False)
    rHWP = Rotate(HWP, psi/2)
    rQWP = Rotate(QWP, alpha)
    b = (rQWP*rHWP*matrix([[1],[0]],dtype = 'complex'))
    b = b/norm(b[0])
    x = b[0,0]
    y =b[1,0]
    #print('El vector de Jones asociado a este estado es:\n {0}.'.format(b))
    #alpha =(res[0][0]*180/pi)
    #theta = (res[0][1]*180/pi)
    return {'psi':psi, 'alpha':alpha,'phi':phi,'theta':theta,'x':x,'y':y,'delta' : delta, 'Jones':b}

def ellipse_jon(x0, psi, delta, give_ellipse = False, verbose = False):
    """" This function is meant to be used in a minimization routine. 
    The objective is to check if a rotated polarizer followed by a QWP 
    retarder returns a certain Jones vector constructed with parameters
    psi and delta. 
        
    Inputs:
    x0: Array, contains the rotation angles for the polarizer and the QWP: [psi, alpha] 
        psi: Rotation angle of the polarizer
        alpha: Rotation angle of the lam/4 retarder
    phi: Angle in radians between the coordinate system and the mayor axes of the ellipse
    theta: Ellipticity angle of the polarization ellipse. It is function of the ratio between 
           the minor and mayor axes: theta = arctan(b/a).
    
    Outputs:
    result: Scalar, is the result of summing the square differences of polarization angles.
            Each angle (phi_c, theta_c) is extracted from the Jones vector after propagation 
            through optical elements, and compared to the required angle. 
    delta: Optional. Is the phase difference between y and x components of the field described 
           a Jones Vector.
       
    Author: Santiago Echeverri Chacón   
    """

    [Psi, alpha] = x0
    
    # Matrices of horizantal oriented polarizer and retarder to be rotated:
    
    #H = np.matrix([[ 1 , 0],\
    #               [ 0 , 0]], dtype = 'complex')
    HWP = np.matrix([[ np.exp(1j*pi/2) , 0],\
                     [ 0 , np.exp(-1j*pi/2)]], dtype = 'complex')
    QWP = np.matrix([[ np.exp(1.j*pi/4) , 0],\
                 [ 0 , np.exp(-1.j*pi/4)]], dtype = 'complex')

    
    ## def Rotate( element, theta):   
    ##     """Applies a rotation on any polarizing element given 
    ##     its Jones matrix representation and a angle of rotation.
    ##     The angle is to be given in radians. The matrix is a 2*2 matrix.  
    ##     """
        
    ##     assert isinstance(element, np.matrix)
    ##     assert element.shape == (2,2)
    
    ##     rotator = np.matrix([[ np.cos(theta), np.sin(theta)],\
    ##                          [-np.sin(theta), np.cos(theta)]],dtype = 'complex')
    ##     return rotator.transpose()*element*rotator

    #M = Rotate(QWP, alpha)
    #P = Rotate(H, psi)
    
    # Definition of an imput Jones vector oriented along the axis of the polarizer,
    # this way the output vector keeps a normalized intensity.
    # It may be redundant when multiplying by the polarizer but this way I can change 
    # the definition later if I want other kinds of inputs.

    In = Rotate(HWP,psi/2)*np.matrix([[1],[0]],dtype = 'complex')
        
    Out = Rotate(QWP, alpha)*In
    
    #In = np.matrix([[np.cos(psi)],[np.sin(psi)*np.exp(1j*0)]],dtype = 'complex')
        
    #Out = M*P*In
    
    # Extract the phase differences and field amplitudes from "Out" and calculate phi,theta:
    
    delta1 = np.angle( Out[0,0] )
    delta2 = np.angle( Out[1,0] )
    
    delta_test = delta2 - delta1
    Ex = abs(Out[0,0])
    Ey = abs(Out[1,0])
        
    phi_c = 1/2*np.arctan2(2*Ex*Ey*np.cos(delta), (Ex**2 - Ey**2))
    theta_c = -1/2*np.arcsin(2*Ex*Ey*np.sin(delta)/(Ex**2 + Ey**2))
    
    result =((delta_test - delta)**2+(Psi-psi)**2)

    if verbose:
        print('delta_test ', delta_test, 'psi: ', psi )
    if give_ellipse:
        return result, phi_c, theta_c
    else:
        return result

def getAnglesFromJones(psi, delta, bnds = ((-np.pi, np.pi),(-np.pi, np.pi))):
    """ Minimization procedure for finding polarization angles.
    Based on desired values of the jones vector this function uses
    the lbfgs minimization algorithm to find the combination of polarization and
    retarder angles that reproduce it.
    
    Inputs:
    psi: Angle that determines the direction of polarization at input.
    delta: Phase difference between x and y component of the wave.

    Outputs:
    dict: Dictionary that contains:
         psi: angle for input polarizer
         alpha: angle for input quarter waveplate retarder
         theta: Ellipticity angle of the polarization ellipse. It is function of the ratio between 
           the minor and mayor axes: theta = arctan(b/a).
         x: Jones vector horizontal component
         y: Jones vector vertical component
         delta: Phase difference between x and y component of the wave.
         
    Author: Santiago Echeverri Chacón
    """
    
    import numpy as np
    from scipy.optimize import fmin_l_bfgs_b
    Psi, alpha =psi, psi
    res = fmin_l_bfgs_b(ellipse_jon, [Psi,alpha], fprime=None,approx_grad=1,\
                                    args = (psi,delta), pgtol=1e-05, bounds = bnds, epsilon=1e-10)
    
    [Psi, alpha] = res[0]
    minimum, phi, theta = ellipse_jon(res[0], Psi,delta, give_ellipse = True, verbose = False)
    rHWP = Rotate(HWP, Psi/2)
    rQWP = Rotate(QWP, alpha)
    b = (rQWP*rHWP*np.matrix([[1],[0]],dtype = 'complex'))
    b = b/b[0]
    x = b[0,0]
    y =b[1,0]
    #print('El vector de Jones asociado a este estado es:\n {0}.'.format(b))
    #alpha =(res[0][0]*180/pi)
    #theta = (res[0][1]*180/pi)
    return {'psi':psi, 'alpha':alpha,'phi':phi,'theta':theta,'x':x,'y':y,'delta' : delta}

def ellipse_angle(x0,psi, alpha, give_ellipse = False):
    """" This function is meant to be used in a minimization routine.
        The objective is to check if a rotated polarizer followed by a QWP
        retarder returns a certain Jones vector constructed with parameters
        psi and delta.
        
        Inputs:
        x0: Array, contains the rotation angles for the polarizer and the QWP: [psi, alpha]
        psi: Rotation angle of the polarizer
        alpha: Rotation angle of the lam/4 retarder
        phi: Angle in radians between the coordinate system and the mayor axes of the ellipse
        theta: Ellipticity angle of the polarization ellipse. It is function of the ratio between
        the minor and mayor axes: theta = arctan(b/a).
        
        Outputs:
        result: Scalar, is the result of summing the square differences of polarization angles.
        Each angle (phi_c, theta_c) is extracted from the Jones vector after propagation
        through optical elements, and compared to the required angle.
        delta: Optional. Is the phase difference between y and x components of the field described
        a Jones Vector.
        
        Author: Santiago Echeverri Chacón
        """
    [delta] = x0
    
    # Matrices of horizantal oriented polarizer and retarder to be rotated:
    
    ## H = np.matrix([[ 1 , 0],\
    ##                [ 0 , 0]], dtype = 'complex')
    HWP = np.matrix([[ np.exp(1j*pi/2) , 0],\
                     [ 0 , np.exp(-1j*pi/2)]], dtype = 'complex')
    QWP = np.matrix([[   np.exp(1j*pi/4) , 0],\
                     [ 0 ,  np.exp(-1j*pi/4) ]], dtype = 'complex')
    
    ## def Rotate( element, theta):
    ##     """Applies a rotation on any polarizing element given
    ##         its Jones matrix representation and a angle of rotation.
    ##         The angle is to be given in radians. The matrix is a 2*2 matrix.
    ##         """
        
    ##     assert isinstance(element, np.matrix)
    ##     assert element.shape == (2,2)
        
    ##     rotator = np.matrix([[ np.cos(theta), np.sin(theta)],\
    ##                          [-np.sin(theta), np.cos(theta)]],dtype = 'complex')
    ##     return rotator.transpose()*element*rotator
    
    #M = Rotate(QWP, alpha)
    #P = Rotate(H, psi)
    
    # Definition of an imput Jones vector oriented along the axis of the polarizer,
    # this way the output vector keeps a normalized intensity.
    # It may be redundant when multiplying by the polarizer but this way I can change
    # the definition later if I want other kinds of inputs.
    In = Rotate(HWP,psi/2)*np.matrix([[1],[0]],dtype = 'complex')
        
    Out = Rotate(QWP, alpha)*In
   
    print('Out',Out)
    Out = Out/Out[0,0]
    print('Out',Out)
    
    JonesOut = np.matrix([[np.cos(psi)],[np.sin(psi)*np.exp(1j*delta)]],dtype = 'complex')
    print('JonesOut', JonesOut, delta)
    JonesOut = JonesOut / JonesOut[0,0]
    print('JonesOut', JonesOut, delta)
    
    result = (Out[0,0].real - JonesOut[0,0].real)**2 + (Out[0,0].imag - JonesOut[0,0].imag)**2+\
             (Out[1,0].real - JonesOut[1,0].real)**2 + (Out[1,0].imag - JonesOut[1,0].imag)**2

    return result

def getJonesFromAngles(psi, alpha, bnds = [(-np.pi, np.pi)]):
    """ Minimization procedure for finding phase retardance between components
        given polarization angles that produce a Jones vector.
        Based on desired values of the polarizer and quarter waveplate this function uses
        the lbfgs minimization algorithm to jones vector representation in the form:
        [[cos(psi)],[sin(psi)*exp(1j*delta)]]
        
        Inputs:
        psi: Angle that determines the direction of polarization at input.
        alpha: Angle for input quarter waveplate retarder.
        
        Outputs:
        
        dict: Dictionary that contains:
            psi: angle for input polarizer
            alpha: angle for input quarter waveplate retarder
            theta: Ellipticity angle of the polarization ellipse. It is function of the ratio between
            the minor and mayor axes: theta = arctan(b/a).
            x: Jones vector horizontal component
            y: Jones vector vertical component
            delta: Phase difference between x and y component of the wave.
        
        Author: Santiago Echeverri Chacón
        """
    
    import numpy as np
    from scipy.optimize import fmin_l_bfgs_b
    delta = np.pi/4
    res = fmin_l_bfgs_b(ellipse_angle,[delta], fprime=None,approx_grad=1,\
                        args = (psi, alpha), pgtol=1e-05, bounds = bnds, epsilon=1e-10)
    [delta] = res[0]
    
    polarizer = Rotate(H, psi)
    retarder = Rotate(QWP, alpha)
    b = (retarder*polarizer*np.matrix([[np.cos(pi/8)],[np.sin(pi/8)]],dtype = 'complex'))
    b = np.matrix([[np.cos(psi)],[np.sin(psi)*np.exp(1j*delta)]],dtype = 'complex')
    #b = b/b[0]
    x = b[0,0]
    y =b[1,0]
    print('El vector de Jones asociado a este estado es:\n {0}.'.format(b))
    #alpha =(res[0][0]*180/pi)
    #theta = (res[0][1]*180/pi)
    return {'psi':psi, 'alpha':alpha,'x':x,'y':y,'delta' : delta,'b':b}
