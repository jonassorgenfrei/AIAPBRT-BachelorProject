
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
# classifivation, regression and clusterin algorithms
from sklearn import metrics
import datetime

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

    training_rmse = []
    validation_rmse = []
    print(datetime.datetime.now())
    for period in range (0, periods) :
        tic = time.time()   # start Timer
        # Train the model, starting from prior state.
        estimator.train(
            input_fn=training_input_fn,
            steps=steps_per_period,
        )

        # 2. Take a break and compute predictions.
        training_predictions = estimator.predict(input_fn=predict_training_input_fn)
        training_predictions = np.array([item['predictions'][0] for item in training_predictions])

        validation_predictions = estimator.predict(input_fn=predict_validation_input_fn)
        validation_predictions = np.array([item['predictions'][0] for item in validation_predictions])

        # Compute training and validation loss.
        training_root_mean_squared_error = math.sqrt(metrics.mean_squared_error(training_predictions, training_targets))
        validation_root_mean_squared_error = math.sqrt(metrics.mean_squared_error(validation_predictions, validation_targets))
        
        toc = time.time()   # end Timer

        # Occasinally print the current loss.
        print("Epoch %d/%d" % (period, periods))
        print("            %ds RMSE: %0.2f, val_RMSE:  %0.2f" % ((toc - tic) / 10, training_root_mean_squared_error, validation_root_mean_squared_error))
        
        # Add the loss metrics from this period to our list.
        training_rmse.append(training_root_mean_squared_error)
        validation_rmse.append(validation_root_mean_squared_error)
    print(datetime.datetime.now())
    if(True) :
        # Output a graph of loss metrics over periods.
        plt.ylabel("RMSE")
        plt.xlabel("Periods")
        plt.title("Root Mean Squared Error vs. Periods")
        plt.tight_layout()
        plt.plot(training_rmse, label="training")
        plt.plot(validation_rmse, label="validation")
        plt.legend()
        plt.show()
    
    return estimator

def serving_input_receiver_fn() :
    """Serving input_fn that builds features from placeholders
    Returns
    -------
    tf.estimator.export.ServingInputReceiver
    """
    dirX = tf.placeholder(dtype=tf.float32, shape=[None, 1], name='direction.x')
    dirY = tf.placeholder(dtype=tf.float32, shape=[None, 1], name='direction.y')
    dirZ = tf.placeholder(dtype=tf.float32, shape=[None, 1], name='direction.z')
    originX = tf.placeholder(dtype=tf.float32, shape=[None, 1], name='origin.x')
    originY = tf.placeholder(dtype=tf.float32, shape=[None, 1], name='origin.y')
    originZ = tf.placeholder(dtype=tf.float32, shape=[None, 1], name='origin.z')

    
    receiver_tensors = {'direction.x': dirX, 'direction.y' : dirY, 'direction.z' : dirZ,
                        'origin.x': originX, 'origin.y' : originY, 'origin.z' : originZ}

    features = {
        key: tf.expand_dims(tensor, -1)
        for key, tensor in receiver_tensors.items()
    } 
   #tf.convert_to_tensor([[dirX, dirY, dirZ, originX, originY, originZ]])
    return tf.estimator.export.ServingInputReceiver(features, receiver_tensors)