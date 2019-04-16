

import matplotlib
matplotlib.use('Agg')
import yt
import math
import matplotlib.pyplot as plt
import numpy as np
import pdb
nar = np.array

from importlib import reload

import looper
reload(looper)
import loop_tools
reload(loop_tools)



directory = '/scratch2/dcollins/Paper19_48/B02/u05-r4-l4-128'
ds = yt.load("%s/DD%04d/data%04d"%(directory,109,109))

leaves = loop_tools.get_leaf_indices(ds) #,h5_name="NEW_PEAKS.h5" )
