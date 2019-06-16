## Imports 
##########

from __future__ import absolute_import, division, print_function

import pathlib

from IPython import display
from matplotlib import cm
from matplotlib import gridspec

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sklearn.metrics as metrics

import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers


## Entry Point
##############
print(tf.__version__)

tf.logging.set_verbosity(tf.logging.ERROR)
pd.options.display.max_rows = 10
pd.options.display.float_format = '{:.1f}'.format

# import California housing data 
california_housing_dataframe = pd.read_csv("../data/set1/intersections.csv", sep=",")

# randomize data
california_housing_dataframe = california_housing_dataframe.reindex(
    np.random.permutation(california_housing_dataframe.index))
print("Data: {}".format(california_housing_dataframe))