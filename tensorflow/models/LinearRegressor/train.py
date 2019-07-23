""" Training the model """

___author___ = "Jonas Sorgenfrei"

########################
## Dependencies
########################

from pathlib import Path
import logging
import sys
import time

# tensorflow
import tensorflow as tf
# panda (column-oriented data analysis API)
import pandas as pd
# NumPython
import numpy as np

# Project Files
from dataHandling import preprocess_features, preprocess_targets
from dataPlotter import plot_features
from model import input_fn, train_model, serving_input_receiver_fn

########################
## Globals
########################

# Program Flags
describeData = False                 # Flag for descripion of import data
describeTrainingSets = False         # Flag for 
showDistributionPlots = False       # Flag for showing input distributions

describeSplittedData = False        # Flag Show Information about splitted sets
printDataSet = False                # Flag for printing inputDataset
plotFeatures = False

# Dataset
trainingTestPercentage = 0.8        # [0,1] Percentage of Training Test of complete 
randomizeTrainingData = True       # Flag for randomization of Training-Validation-Data
trainingValidationPercentage = 0.8  # [0,1] Percentage of training/validationData from complete trainingSet
createSyntAngles = False            # Flag for creating syntectic Angles of vectors
useFeatureCrosses = False           # Flag for crossing features

# Model Parameter
trainModel = True                  # Flag Train Model
PERIODS = 10
LEARNING_RATE = 0.0001
STEPS = 500
BATCH_SIZE = 10

########################
## Functions
########################
def construct_feature_columns(input_features) :
    """Construct the TensorFlow (multiple) Feature Columns.

    Args:
        input_features: The Dataset containing the input features.
    Returns:
        A DataFrame that contains the target feature.
    """ 
    return set([tf.feature_column.numeric_column(my_feature)
                for my_feature in input_features])

########################
## Main
########################
if __name__ == '__main__':
    print(tf.__version__)

    # Logging
    #########
    Path('model').mkdir(exist_ok=True)
    tf.compat.v1.logging.set_verbosity(logging.ERROR)
    handlers = [
        logging.FileHandler('model/train_v01.log'),
        logging.StreamHandler(sys.stdout)
    ]
    logging.getLogger('tensorflow').handlers = handlers

    #Panda set display options
    pd.options.display.max_rows = 10
    pd.options.display.float_format = '{:.1f}'.format

    # Import Data
    #############
    dataframe = pd.read_csv("../../data/set1/smallSet.csv", sep=",") # relational data table of intersection tests
    print("#################################")
    print("Imported Data; Records:", len(dataframe))
    print("#################################")

    # Print Data Overview
    if(describeData) :
        print(dataframe.describe())

    # Show Distribution
    if(showDistributionPlots) :
        dataframe.hist('origin.x')
        dataframe.hist('origin.y')
        dataframe.hist('origin.z')
        dataframe.hist('direction.x')
        dataframe.hist('direction.y')
        dataframe.hist('direction.z')
        dataframe.hist('hit')
        plt.show()
    
    # Show Data Content
    if(printDataSet) :
        print("Data: {}".format(dataframe))

    # Split Data
    # Split Training-Validaiton <--> Testing
    amountTrainingData = int(len(dataframe)*trainingTestPercentage)  # amount of Data used for training and validation
    amountTestData = len(dataframe) - amountTrainingData             # Data used for testing

    trainingValData = dataframe.head(amountTrainingData)             # Data used for training and validation

    # Split Training <--> Validation
    amountTrainingTrainginData = int(len(trainingValData)*trainingValidationPercentage)
    amountTrainingValidationData = len(trainingValData) - amountTrainingTrainginData

     # randomize data
    if(randomizeTrainingData) :
        trainingValData = trainingValData.reindex(np.random.permutation(trainingValData.index))

    print("#################################")
    print("Splitted-Data:")
    print("Training-Data Amount:", amountTrainingData)
    print("  Training Amount:", amountTrainingTrainginData)
    print("  Validation Amount:", amountTrainingValidationData)
    print("#################################")

    training_examples = preprocess_features(trainingValData.head(amountTrainingTrainginData), createSyntAngles)
    training_targets = preprocess_targets(trainingValData.head(amountTrainingTrainginData))

    validation_examples = preprocess_features(trainingValData.tail(amountTrainingValidationData), createSyntAngles)
    validation_targets = preprocess_targets(trainingValData.tail(amountTrainingValidationData))

    # Describe Splitted Data
    if(describeSplittedData) :
        print("#################################")
        print("Training Set")
        print("#################################")
        print(training_examples.describe())
        print(training_targets.describe())

        print("#################################")
        print("Validation Set")
        print("#################################")
        print(validation_examples.describe())
        print(validation_targets.describe())

    # plot features
    if(plotFeatures) :
        plot_features(validation_examples,
                      validation_targets,
                      training_examples,
                      training_targets)

    # Train model
    #############
    if(trainModel) :
        print("=============== TRAIN MODEL ===============")
        tic = time.time()

        # Create a linear regressor object
        optimizer = tf.train.GradientDescentOptimizer(learning_rate=LEARNING_RATE)
        optimizer = tf.contrib.estimator.clip_gradients_by_norm(optimizer, 5.0)
        linear_regressor = tf.estimator.LinearRegressor(
            feature_columns=construct_feature_columns(training_examples),
            model_dir='model/v001',
            optimizer=optimizer
        )

        linear_regressor = train_model(
            linear_regressor,
            steps=STEPS,
            batch_size=BATCH_SIZE,
            periods_amount=PERIODS,
            training_examples=training_examples,
            training_targets=training_targets,
            validation_examples=validation_examples,
            validation_targets=validation_targets)

        toc = time.time()
        print("=============== MODEL TRAINING FINISHED. ===============")

        print('Average time in train.py: {}s'.format((toc - tic) / 10))

        feature_columns = [tf.feature_column.numeric_column(key=key)
                   for key in training_examples.keys()]
        feature_spec = tf.feature_column.make_parse_example_spec(feature_columns)
        serving_input_receiver_fn2 = tf.estimator.export.build_parsing_serving_input_receiver_fn(feature_spec)
        linear_regressor.export_saved_model('saved_model', serving_input_receiver_fn2);