
"""Reload and serve a saved model"""

__author__ = "Jonas Sorgenfrei"

########################
## Dependencies
########################

from pathlib import Path
import time
# panda (column-oriented data analysis API)
import pandas as pd

# Math Lib
import math
import numpy as np
# classifivation, regression and clusterin algorithms
from sklearn import metrics

# tensorflow
import tensorflow as tf
from tensorflow.contrib import predictor

# Project Files
from dataHandling import preprocess_features, preprocess_targets
# Time
import time
import datetime
########################
## Globals
########################
trainingTestPercentage = 1.0        # [0,1] Percentage of Training Test of complete

########################
## Functions
########################
def my_service():
    """Some service yielding numbers"""
    start, end = 100, 110
    for number in range(start, end):
        yield number

########################
## Main
########################
if __name__ == '__main__':
    export_dir = 'saved_model_ds2_1'
    subdirs = [x for x in Path(export_dir).iterdir()
               if x.is_dir() and 'temp' not in str(x)]
    latest = str(sorted(subdirs)[-1])
    predict_fn = predictor.from_saved_model(latest)

   
    # Import Data
    #############
    dataframe = pd.read_csv("../completeSet.csv", sep=",")    # relational data table of intersection tests

    amountTrainingData = int(len(dataframe)*trainingTestPercentage)     # amount of Data used for training and validation
    amountTestData = len(dataframe) - amountTrainingData                # Data used for testing

    testData = dataframe.tail(amountTestData)                           # Data used for training and validation

    test_examples = preprocess_features(testData, False)
    test_targets = preprocess_targets(testData)

    examples = []
    for index, row in test_examples.iterrows():
        feature = {}
        for col, value in row.iteritems():
            feature[col] = tf.train.Feature(float_list=tf.train.FloatList(value=[value]))
        example = tf.train.Example(
            features=tf.train.Features(
                feature=feature
            )
        )
        examples.append(example.SerializeToString())
    print(datetime.datetime.now())
    test_predictions = predict_fn({'inputs': examples})
    print(datetime.datetime.now())
#    test_probabilities = np.array([item['scores'][1] for item in test_predictions])
#    false_positive_rate, true_positive_rate, thresholds = metrics.roc_curve(
#                                                        test_targets, test_probabilities)                               # obtain true positive and false positive rates
#    plt.plot(false_positive_rate, true_positive_rate, label="our model")
#    plt.plot([0, 1], [0, 1], label="random classifier")
#    _ = plt.legend(loc=2)
#    plt.title("ROC Curve")          # plot ROC Curve
#    plt.show()
