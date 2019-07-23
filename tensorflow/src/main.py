## Imports 
##########

# help backward/forward comaptibility with python 2/3
from __future__ import absolute_import, division, print_function

# tensorflow 
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.python.data import Dataset

#from IPython import display

# Math utilities
import math
import pathlib
from matplotlib import cm
from matplotlib import gridspec
from matplotlib import pyplot as plt
import numpy as np

# panda (column-oriented data analysis API)
import pandas as pd
# classifivation, regression and clusterin algorithms
from sklearn import metrics
# System
import sys

# GLOBALS
#########

showDistributionPlots = True       # Flag for showing input distributions
randomizeData = True               # Flag for randomization of input data

trainingValidationPercentage = 0.8  # [0,1] Percentage of trainingData from complete Dataset
describeSplittedData = True        # Flag Show Information about splitted sets

printDataSet = False                # Flag for printing inputDataset

# Model
PERIODS = 10
LEARNING_RATE = 0.0003
STEPS = 500
BATCH_SIZE = 5

######################
## Functions
######################

def preprocess_features(intersections_dataframe):
    """ Prepares input features from intersections data set

    Args:
        intersections_dataframe: A (panda) Dataframe expected to contain data from
            the intersections data set.
    Returns:
        A Dataframe that contains the features to be used for the model, including
        synthetic feature.
    """
    selected_features = intersections_dataframe[
        ["origin.x",
         "origin.y",
         "origin.z",
         "direction.x",
         "direction.y",
         "direction.z"]]
    preprocess_features = selected_features.copy()
    # Create a synthetic feature
    #processed_features["angAlpha"] = dir.x / (sqrt(dir.x*dir.x+dir.y*dir.y+dir.z*dir.z))
    #processed_features["angleBeta"] = dir.y / (sqrt(dir.x*dir.x+dir.y*dir.y+dir.z*dir.z))
    #processed_features["angleGamma"] = dir.z / (sqrt(dir.x*dir.x+dir.y*dir.y+dir.z*dir.z))
    return preprocess_features

def preprocess_targets(intersections_dataframe) :
    """ Prepares target features (i.e., labels) from intersections data set.
       
    Args:
        intersections_dataframe: A (panda) Dataframe expected to contain data from the intersections data set.
    Returns:
        A  DataFrame that contains the target feature.
    """
    output_targets = pd.DataFrame()
    output_targets["hit"] = intersections_dataframe["hit"]
    return output_targets

def plot_features(validation_examples,validation_targets,training_examples,training_targets) :
    plt.figure(figsize=(13, 8))

    ax = plt.subplot(1, 2, 1)
    ax.set_title("Validation Data")

    ax.set_autoscaley_on(False)
    ax.set_ylim([32, 43])   #   limit y
    ax.set_autoscalex_on(False)
    ax.set_xlim([-126, -112])   # limit x
    plt.scatter(validation_examples["origin.x"],
                validation_examples["origin.y"],
                cmap="coolwarm",
                c=validation_targets["hit"])

    ax = plt.subplot(1,2,2)
    ax.set_title("Training Data")

    ax.set_autoscaley_on(False) 
    ax.set_ylim([32, 43])   # limit y
    ax.set_autoscalex_on(False)
    ax.set_xlim([-126, -112])   # limit x
    plt.scatter(training_examples["origin.x"],
                training_examples["origin.y"],
                cmap="coolwarm",
                c=training_targets["hit"])
    _ = plt.plot()
    plt.show()

def input_fn(features, targets, batch_size=1, shuffle=True, num_epochs=None) :
    """Trains a linear regression model of multiple features.
    
    Args:
        feature: pandas DataFrame of features
        targets: pandas DataFrame of targets
        batch_size: Size of batches to be passed to the model
        shuffle: True or False. Whether to shuffle the data.
        num_epochs: Number of epochs for which data should be repeated. None = repeat indefinitely
    Returns:
        Tuple of (features, labels) for next data batch
    """

    # Convert pandas data into a dict of np arrays.
    features = {key:np.array(value) for key,value in dict(features).items() }

    # Construct a dataset, and configure batching/repeating.
    ds = Dataset.from_tensor_slices((features, targets))        # Note: wanring 2 GB LIMIT!
    ds = ds.batch(batch_size).repeat(num_epochs)

    if shuffle :
        ds = ds.shuffle(10000)

    # Return the next batch of data.
    features,labels = ds.make_one_shot_iterator().get_next()
    return features, labels

def construct_feature_columns(input_features) :
    """Construct the TensorFlow (multiple) Feature Columns.

    Args:
        input_features: The names of the numerical input features to use.
    Returns:
        A set of feature columns
    """ 
    return set([tf.feature_column.numeric_column(my_feature)
                for my_feature in input_features])

def train_model(learning_rate,
    steps,
    batch_size,
    training_examples,
    training_targets,
    validation_examples,
    validation_targets) :
    """ Trains a linear regression model of multiple features.

    In addition to training, this function also prints training progress information, 
    as well as a plot of the training and validation loss over time.

    Args:
        learning_rate: A `float`, the learning rate.
        steps: A non-zero `int`, the total number of training steps. A training step
                consists of a forward and backward pass using a single batch.
        batch_size: A non-zero `int`, the batch size.
        training_examples: A `DataFrame` containing one or more columns from 
            `intersections_dataframe` to use as input features for training.
        training_targets: A `DataFrame` containing exactly one column from
            `intersections_dataframe`  to use as target for training.
        validation_examples: A `DataFrame` containing one or more columns from 
            `intersections_dataframe` to use as input features for validation.
        validation_examples: A `DataFrame` containing exactly one column from 
            `intersections_dataframe` to use as target for validation.
    Returns:
        A `Linear Regressor object trainedon the training data`
    """
    periods = PERIODS
    steps_per_period = steps / periods

    # Create a linear regressor object
    optimizer = tf.train.GradientDescentOptimizer(learning_rate=learning_rate)
    optimizer = tf.contrib.estimator.clip_gradients_by_norm(optimizer, 5.0)
    linear_regressor = tf.estimator.LinearRegressor(
        feature_columns=construct_feature_columns(training_examples),
        optimizer=optimizer
    )

    # 1. Create input functions
    training_input_fn = lambda: input_fn(training_examples, training_targets["hit"], batch_size=batch_size)
    predict_training_input_fn = lambda: input_fn(training_examples, training_targets["hit"], num_epochs=1, shuffle=False)
    predict_validation_input_fn = lambda: input_fn(validation_examples, validation_targets["hit"], num_epochs=1, shuffle=False)
    
    # Train the model, but do so inside a loop so that we can periodically assess
    # loss metrics.
    print("Training model ...")
    print("RMSE (on training data):")
    training_rmse = []
    validation_rmse = []
    for period in range (0, periods) :
        # Train the model, starting from prior state.
        linear_regressor.train(
            input_fn=training_input_fn,
            steps=steps_per_period,
        )

        # 2. Take a break and compute predictions.
        training_predictions = linear_regressor.predict(input_fn=predict_training_input_fn)
        training_predictions = np.array([item['predictions'][0] for item in training_predictions])

        validation_predictions = linear_regressor.predict(input_fn=predict_validation_input_fn)
        validation_predictions = np.array([item['predictions'][0] for item in validation_predictions])


        # Compute training and validation loss.
        training_root_mean_squared_error = math.sqrt(
            metrics.mean_squared_error(training_predictions, training_targets))
        validation_root_mean_squared_error = math.sqrt(
            metrics.mean_squared_error(validation_predictions, validation_targets))
        # Occasinally print the current loss.
        print(" period %02d : %0.2f" % (period, training_root_mean_squared_error))
        # Add the loss metrics from this period to our list.
        training_rmse.append(training_root_mean_squared_error)
        validation_rmse.append(validation_root_mean_squared_error)
    print("Model training finished.")

    # Output a graph of loss metrics over periods.
    plt.ylabel("RMSE")
    plt.xlabel("Periods")
    plt.title("Root Mean Squared Error vs. Periods")
    plt.tight_layout()
    plt.plot(training_rmse, label="training")
    plt.plot(validation_rmse, label="validation")
    plt.legend()
    plt.show()
    
    return linear_regressor

######################
## Entry Point
######################

print(tf.__version__)

#Tensorflow log errors 
tf.logging.set_verbosity(tf.logging.ERROR)
#Panda set display options
pd.options.display.max_rows = 10
pd.options.display.float_format = '{:.1f}'.format

# import intersection data file
###############################
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
if(printDataSet) :
    print("Data: {}".format(intersections_dataframe))

if(randomizeData) :
    # randomize data
    intersections_dataframe = intersections_dataframe.reindex(
        np.random.permutation(intersections_dataframe.index))

print("#################################")
print("Intersection Tests Amount:", len(intersections_dataframe))
print("#################################")

# Define the input features
#inputFeature = intersections_dataframe[["origin.x","origin.y","origin.z", "direction.x", "direction.y","direction.z"]]
#Configure a numeric feature column for hit
#feature_columns = [tf.feature_column.numeric_column("origin.x"),
#                   tf.feature_column.numeric_column("origin.y"),
#                   tf.feature_column.numeric_column("origin.z"),
#                   tf.feature_column.numeric_column("direction.x"),
#                   tf.feature_column.numeric_column("direction.y"),
#                   tf.feature_column.numeric_column("direction.z")]

#print("Data: {}".format(inputFeature))



# Define the Target / Label
#label = intersections_dataframe["hit"]

###################################
## Split Data
###################################

intersections_dataframe = intersections_dataframe.head(500000)
print(len(intersections_dataframe))

amountTrainingData = int(len(intersections_dataframe)*trainingValidationPercentage)
amountValidationData = len(intersections_dataframe) - amountTrainingData

training_examples = preprocess_features(intersections_dataframe.head(amountTrainingData))
training_targets = preprocess_targets(intersections_dataframe.head(amountTrainingData))

validation_examples = preprocess_features(intersections_dataframe.head(amountValidationData))
validation_targets = preprocess_targets(intersections_dataframe.head(amountValidationData))

if(describeSplittedData) :
    print("##################")
    print("###Training Set###")
    print("##################")
    print(training_examples.describe())
    print(training_targets.describe())

    print("##################")
    print("###Training Set###")
    print("##################")
    print(validation_examples.describe())
    print(validation_targets.describe())

################
# Plot Features
################
################
#plot_features(validation_examples,
#              validation_targets,
#              training_examples,
#              training_targets)

sys.exit(0)

##################################
# Create Model (linear Regressor)
##################################
linear_regressor = train_model(
    # TWEAK THESE VALUES TO SEE HOW MUCH YOU CAN IMPROVE THE RMSE
    learning_rate=LEARNING_RATE,
    steps=STEPS,
    batch_size=BATCH_SIZE,
    training_examples=training_examples,
    training_targets=training_targets,
    validation_examples=validation_examples,
    validation_targets=validation_targets)

##################################
# Evaluate on Test Data
##################################
intersections_test_data = intersections_dataframe.tail(500000)

test_examples = preprocess_features(intersections_test_data)
test_targets = preprocess_targets(intersections_test_data)

predict_test_input_fn = lambda: input_fn(
    test_examples,
    test_targets["hit"],
    num_epochs=1,
    shuffle=False)

test_predictions = linear_regressor.predict(input_fn=predict_test_input_fn)
test_predictions = np.array([item['predictions'][0] for item in test_predictions])

root_mean_squared_error = math.sqrt(
    metrics.mean_squared_error(test_predictions, test_targets))

print("Final RMSE (on test data): %0.2f" % root_mean_squared_error)