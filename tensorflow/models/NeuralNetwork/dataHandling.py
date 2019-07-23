""" Restructor the Dat """

___author___ = "Jonas Sorgenfrei"

########################
## Dependencies
########################

# panda (column-oriented data analysis API)
import pandas as pd

########################
## Functions
########################

def preprocess_features(intersections_dataframe, createSyntAngles):
    """ Prepares input features from intersections data set

    Args:
        intersections_dataframe: A (panda) Dataframe expected to contain data from
            the intersections data set.
        createSyntAngles: Create synthetic angles of dir-Vectors
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
    if(createSyntAngles) :
        processed_features["angAlpha"] = dir.x / (sqrt(dir.x*dir.x+dir.y*dir.y+dir.z*dir.z))
        processed_features["angleBeta"] = dir.y / (sqrt(dir.x*dir.x+dir.y*dir.y+dir.z*dir.z))
        processed_features["angleGamma"] = dir.z / (sqrt(dir.x*dir.x+dir.y*dir.y+dir.z*dir.z))
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