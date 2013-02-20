# This script invokes the module for executing the algorithm. By default, if ExecuteAlgorithmWrapITK.py is available, it is used. Otherwise ExecuteAlgorithmCSwig.py is used.

import sys
try:
    from ExecuteAlgorithmWrapITK import *
except:
    try:
        from ExecuteAlgorithmCSwig import *
    except:
        print "Could not load module for executing. Please Check if ExecuteAlgorithmWrapITK.py or ExecuteAlgorithmCSwig.py is reachable by Python."
