# -*- coding: utf-8 -*-
"""
Created on Thu Oct 16 15:18:24 2014

@author: Santiago Echeverri CHacón
"""

# Used successfully in Python2.5 with matplotlib 0.91.2 and PyQt4 (and Qt 4.3.3)
from distutils.core import setup
import py2exe

# We need to exclude matplotlib backends not being used by this executable.  You may find
# that you need different excludes to create a working executable with your chosen backend.
# We also need to include include various numerix libraries that the other functions call.

opts = {
    'py2exe': { "includes" : ["sip", "PyQt4", "matplotlib.backends",  "matplotlib.backends.backend_qt4agg",
                               "matplotlib.figure","pylab", "numpy","propagator","polarization_routines","numpy",
                               "scipy.optimize",'scipy.optimize.minpack2'],
                'excludes': ['_gtkagg', '_tkagg', '_agg2', '_cairo', '_cocoaagg',
                             '_fltkagg', '_gtk', '_gtkcairo' ],
                'dll_excludes': ['libgdk-win32-2.0-0.dll',
                                 'libgobject-2.0-0.dll']
              }
       }

# for console program use 'console = [{"script" : "scriptname.py"}]
setup(name='Propagator',
      version='1.0',console=['propagator_app.py'], options=opts)