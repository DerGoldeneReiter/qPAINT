# Load packages
import h5py as h5py #hdf5 handling
import numpy as np #numpy data formats and operators
import matplotlib.pyplot as plt #plotting
import os.path as ospath #platform independent paths
import importlib
from tqdm import tqdm
import javabridge
import bioformats

# Define path
dir_name='E:\\Projects\\qPAINT\\data\\20171110_Max_2x2Grid_varN'
file_name='P1_Cy3b_PPT_10mW_2_5nM_4C_Origami_1_MMStack_Pos0.ome.tif'
path=ospath.join(dir_name,file_name)



# Start javabridge for bioformats
javabridge.start_vm(class_path=bioformats.JARS)

meta=bioformats.get_omexml_metadata(path)

# Kill javabridge for bioformats
javabridge.kill_vm()


