
"""Model_fn"""

__author__ = "Jonas Sorgenfrei"

########################
## Dependencies
########################

# Tensorflow
import tensorflow as tf
from tensorflow.python.data import Dataset

# Matplotlib
from matplotlib import pyplot as plt
# NumPython
import numpy as np
# Math Lib
import math
# Time
import time
import datetime

# classifivation, regression and clusterin algorithms
from sklearn import metrics

########################
## Functions
########################

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

def train_model(estimator,
    steps,
    batch_size,
    periods_amount,
    training_examples,
    training_targets,
    validation_examples,
    validation_targets) :
    """ Trains a linear regression model of multiple features.

    In addition to training, this function also prints training progress information, 
    as well as a plot of the training and validation loss over time.

    Args:
        estimator: TODO!!!!
        steps: A non-zero `int`, the total number of training steps. A training step
                consists of a forward and backward pass using a single batch.
        batch_size: A non-zero `int`, the batch size.
        periods_amount: Amount of Periods
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
    periods = periods_amount
    steps_per_period = steps / periods


    # 1. Create input functions
    training_input_fn = lambda: input_fn(training_examples, training_targets["hit"], batch_size=batch_size)
    predict_training_input_fn = lambda: input_fn(training_examples, training_targets["hit"], num_epochs=1, shuffle=False)
    predict_validation_input_fn = lambda: input_fn(validation_examples, validation_targets["hit"], num_epochs=1, shuffle=False)
    
    # Train the model, but do so inside a loop so that we can periodically assess
    # loss metrics.
    print("Train on %d samples, validate on %d samples" % (len(training_examples),len(validation_examples)))

    training_log_losses = []
    validation_log_losses = []
    
    ticMain = time.time()   # start timer
    print(datetime.datetime.now())
    for period in range (0, periods) :
        tic = time.time()   # start Timer
        # Train the model, starting from prior state.
        estimator.train(
            input_fn=training_input_fn,
            steps=steps_per_period,
        )

        # 2. Take a break and compute predictions.
        training_probabilities = estimator.predict(input_fn=predict_training_input_fn)
        training_probabilities = np.array([item['probabilities'] for item in training_probabilities])
    
        validation_probabilities = estimator.predict(input_fn=predict_validation_input_fn)
        validation_probabilities = np.array([item['probabilities'] for item in validation_probabilities])

        # Calculate Log loss using the probabilities
        # LogLoss = SUM [x,y]eD ( -y * log(ypred) - (1 - y) * log(1 - ypred))
        training_log_loss = metrics.log_loss(training_targets, training_probabilities)
        validation_log_loss = metrics.log_loss(validation_targets, validation_probabilities)

        # Add the loss metrics from this period to our list.
        training_log_losses.append(training_log_loss)
        validation_log_losses.append(validation_log_loss)
        
        toc = time.time()   # end Timer

        # Occasinally print the current loss.
        print("Epoch %d/%d" % (period+1, periods))
        print("            %ds LogLoss: %0.2f, val_LogLoss:  %0.2f" % ((toc - tic) / 10, training_log_loss, validation_log_loss))
    
    tocMain = time.time()   # end timer
    print('Average time in model.py: {}s'.format((tocMain - ticMain) / 10))

    print(datetime.datetime.now())

    if(True) :
        # Output a graph of loss metrics over periods.
        plt.ylabel("LogLoss")
        plt.xlabel("Periods")
        plt.title("LogLoss vs. Periods")
        plt.tight_layout()
        plt.plot(training_log_losses, label="training")
        plt.plot(validation_log_losses, label="validation")
        plt.legend()
        plt.show()
    
    return estimator