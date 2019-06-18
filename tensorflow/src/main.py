## Imports 
##########

# help backward/forward comaptibility with python 2/3
from __future__ import absolute_import, division, print_function

# tensorflow 
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers

# Math utilities
import pathlib
from IPython import display
from matplotlib import cm
from matplotlib import gridspec
import matplotlib.pyplot as plt
import numpy as np
# panda (column-oriented data analysis API)
import pandas as pd
import sklearn.metrics as metrics
# System
import sys

# GLOBALS
#########
showDistributionPlots = False
outputDataSample = False

## Entry Point
##############
print(tf.__version__)

tf.logging.set_verbosity(tf.logging.ERROR)
pd.options.display.max_rows = 10
pd.options.display.float_format = '{:.1f}'.format

# import intersection data file
intersections_dataframe = pd.read_csv("../../data/set1/intersections.csv", sep=",") # relational data table of intersection tests
# Output Intersection Data Overview
print(intersections_dataframe.describe())

# Show Distribution
if(showDistributionPlots) :
    intersections_dataframe.hist('origin.x')
    intersections_dataframe.hist('origin.y')
    intersections_dataframe.hist('origin.z')
    intersections_dataframe.hist('direction.x')
    intersections_dataframe.hist('direction.y')
    intersections_dataframe.hist('direction.z')
    intersections_dataframe.hist('hit')
    plt.show()

# Output Data
#############
if(outputDataSample) :
    print("Data: {}".format(intersections_dataframe))

# randomize data
intersections_dataframe = intersections_dataframe.reindex(
    np.random.permutation(intersections_dataframe.index))

print("#################################")
print("Intersection Tests Amount:", len(intersections_dataframe))
print("#################################")

# Define the input features
inputFeature = intersections_dataframe[["origin.x","origin.y","origin.z", "direction.x", "direction.y","direction.z"]]
#Configure a numeric feature column for hit
feature_columns = [tf.feature_column.numeric_column("origin.x"),
                   tf.feature_column.numeric_column("origin.y"),
                   tf.feature_column.numeric_column("origin.z"),
                   tf.feature_column.numeric_column("direction.x"),
                   tf.feature_column.numeric_column("direction.y"),
                   tf.feature_column.numeric_column("direction.z")]

print("Data: {}".format(inputFeature))

# Define the Target / Label
label = intersections_dataframe["hit"]

#############################
# Create Model 
#############################

# Use gradient descent as the optimizer for training the model.
my_optimizer=tf.train.GradientDescentOptimizer(learning_rate=0.0000001)
my_optimizer = tf.contrib.estimator.clip_gradients_by_norm(my_optimizer, 5.0)

# Configure the linear regression model with our feature columns and optimizer.
# Set a learning rate of 0.0000001 for Gradient Descent.
linear_regressor = tf.estimator.LinearRegressor(feature_columns=feature_columns,
                                                optimizer=my_optimizer)

