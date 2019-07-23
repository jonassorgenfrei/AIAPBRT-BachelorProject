
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
# classifivation, regression and clusterin algorithms
from sklearn import metrics

# tensorflow
import tensorflow as tf
from tensorflow.contrib import predictor

# Project Files
from dataHandling import preprocess_features, preprocess_targets

########################
## Globals
########################
trainingTestPercentage = 0.8        # [0,1] Percentage of Training Test of complete 

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
    export_dir = 'saved_model'
    subdirs = [x for x in Path(export_dir).iterdir()
               if x.is_dir() and 'temp' not in str(x)]
    latest = str(sorted(subdirs)[-1])
    predict_fn = predictor.from_saved_model(latest)

    tic = time.time()

    # Import Data
    #############
    dataframe = pd.read_csv("../../data/set1/smallSet.csv", sep=",")    # relational data table of intersection tests

    amountTrainingData = int(len(dataframe)*trainingTestPercentage)     # amount of Data used for training and validation
    amountTestData = len(dataframe) - amountTrainingData                # Data used for testing

    testData = dataframe.tail(amountTestData)                           # Data used for training and validation

    test_examples = preprocess_features(testData, False)
    test_targets = preprocess_targets(testData)


    #predict_test_input_fn = lambda: input_fn( test_examples,
    #test_targets["hit"],
    #num_epochs=1,
    #shuffle=False)

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

    test_predictions = predict_fn({'inputs': examples})['outputs']
        
    root_mean_squared_error = math.sqrt(metrics.mean_squared_error(test_predictions, test_targets))
    print("Final RMSE (on test data): %0.2f" % root_mean_squared_error)
    
    acc = matrics.accuracy_score(test_targets, test_predictions)
    print("Accuracy: %0.2f" % acc)

    toc = time.time()
    print('Average time in serve.py: {}s'.format((toc - tic) / 10))