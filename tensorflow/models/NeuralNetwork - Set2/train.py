""" Training the model """
from __future__ import print_function

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
gpu_options = tf.GPUOptions(per_process_gpu_memory_fraction=0.8)
sess = tf.Session(config=tf.ConfigProto(gpu_options=gpu_options))
from keras import backend as K
import tensorflow as tf
config = tf.ConfigProto()
config.gpu_options.per_process_gpu_memory_fraction = 0.8
session = tf.Session(config=config)
K.set_session(session)
# panda (column-oriented data analysis API)
import pandas as pd
# NumPython
import numpy as np

# Matplotlib
from matplotlib import pyplot as plt

from sklearn import metrics

# Project Files
from dataHandling import preprocess_features, preprocess_targets
from dataPlotter import plot_features
from model import input_fn, train_model

########################
## Globals
########################

# Program Flags
describeData = True                 # Flag for descripion of import data
describeTrainingSets = True         # Flag for 
showDistributionPlots = False       # Flag for showing input distributions

describeSplittedData = False        # Flag Show Information about splitted sets
printDataSet = True                # Flag for printing inputDataset
plotFeatures = False

evaluateModel = True
exportModel = True

# Dataset
trainingTestPercentage = 0.7        # [0,1] Percentage of Training Test of complete 
randomizeTrainingData = True       # Flag for randomization of Training-Validation-Data
trainingValidationPercentage = 0.8  # [0,1] Percentage of training/validationData from complete trainingSet
createSyntAngles = False            # Flag for creating syntectic Angles of vectors
useFeatureCrosses = False           # Flag for crossing features

# Model Parameter
trainModel = True                  # Flag Train Model
PERIODS = 10
LEARNING_RATE = 0.002
STEPS = 900
BATCH_SIZE = 60
REGULARIZATION_STRENGTH = 0.1

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
        logging.FileHandler('model/train_set02_v01.log'),
        logging.StreamHandler(sys.stdout)
    ]
    logging.getLogger('tensorflow').handlers = handlers

    #Panda set display options
    pd.options.display.max_rows = 10
    pd.options.display.float_format = '{:.1f}'.format

    # Import Data
    #############
    dataframe = pd.read_csv("../../data/set2/completeSet.csv", sep=",") # relational data table of intersection tests
    print("#################################")
    print("Imported Data; Records:", len(dataframe))
    print("#################################")

    # Print Data Overview
    if(describeData) :
        print(dataframe.describe())

    # Show Distribution
    if(showDistributionPlots) :
        dataframe.hist('A.x')
        dataframe.hist('A.y')
        dataframe.hist('A.z')
        dataframe.hist('B.x')
        dataframe.hist('B.y')
        dataframe.hist('B.z')
        dataframe.hist('d.x')
        dataframe.hist('d.y')
        dataframe.hist('d.z')
        dataframe.hist('l')
        dataframe.hist('v')
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
    print("  Test Amount:", amountTestData)
    print("#################################")

    # Training Data
    training_examples = preprocess_features(trainingValData.head(amountTrainingTrainginData), createSyntAngles)
    training_targets = preprocess_targets(trainingValData.head(amountTrainingTrainginData))

    # Validation Data 
    validation_examples = preprocess_features(trainingValData.tail(amountTrainingValidationData), createSyntAngles)
    validation_targets = preprocess_targets(trainingValData.tail(amountTrainingValidationData))

    # Test Data 
    testData = dataframe.tail(amountTestData)
    test_examples = preprocess_features(testData, createSyntAngles)
    test_targets = preprocess_targets(testData)

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
        # Create a DNNClassifier regressor object
        hidden_units=[600, 600, 400, 400]
        #optimizer = tf.train.GradientDescentOptimizer(learning_rate=LEARNING_RATE)      # TODO: TRY DIFFERENT OPTIMIZER
        optimizer = tf.train.AdamOptimizer(learning_rate=LEARNING_RATE)
        optimizer = tf.contrib.estimator.clip_gradients_by_norm(optimizer, 5.0)
        dnn = tf.estimator.DNNClassifier(
            feature_columns=construct_feature_columns(training_examples),
            model_dir='model_set02/v003_synth',
            hidden_units=hidden_units,        # defines the structure of the NN, provides a list of ints, where each int corresponds to a hidden layer and indicates the number of nodes in it
            optimizer=optimizer
        )

        dnn = train_model(
            dnn,
            steps=STEPS,
            batch_size=BATCH_SIZE,
            periods_amount=PERIODS,
            training_examples=training_examples,
            training_targets=training_targets,
            validation_examples=validation_examples,
            validation_targets=validation_targets)

        print("=============== MODEL TRAINING FINISHED. ===============")
        
        if(evaluateModel) :
            print("=============== EVALUATE ===============")
            # Validation Data
            predict_validation_input_fn = lambda: input_fn(validation_examples, validation_targets["v"], num_epochs=1, shuffle=False)
            evaluation_metrics = dnn.evaluate(input_fn=predict_validation_input_fn)

            print("AUC on the validation set: %0.2f" % evaluation_metrics['auc'])
            print("Accuracy on the validation set: %0.2f" % evaluation_metrics['accuracy'])
            
            validation_probabilities = dnn.predict(input_fn=predict_validation_input_fn)
            # Get just the probabilities for the positive class.
            validation_probabilities = np.array([item['probabilities'][1] for item in validation_probabilities])

            false_positive_rate, true_positive_rate, thresholds = metrics.roc_curve(
                                                                    validation_targets, validation_probabilities)                               # obtain true positive and false positive rates
            plt.plot(false_positive_rate, true_positive_rate, label="our model")
            plt.plot([0, 1], [0, 1], label="random classifier")
            _ = plt.legend(loc=2)
            plt.title("ROC Curve")          # plot ROC Curve
            plt.show()

            # Test Data
            predict_test_input_fn = lambda: input_fn(test_examples, test_targets["v"], num_epochs=1, shuffle=False)
            evaluation_test_metrics = dnn.evaluate(input_fn=predict_test_input_fn)

            print("AUC on the test set: %0.2f" % evaluation_test_metrics['auc'])
            print("Accuracy on the test set: %0.2f" % evaluation_test_metrics['accuracy'])


        if(exportModel) :
            print("=============== Export ===============")
            feature_columns = [tf.feature_column.numeric_column(key=key)
                       for key in training_examples.keys()]
            feature_spec = tf.feature_column.make_parse_example_spec(feature_columns)
            serving_input_receiver_fn2 = tf.estimator.export.build_parsing_serving_input_receiver_fn(feature_spec)
            dnn.export_saved_model('saved_model_set02_3', serving_input_receiver_fn2);

        