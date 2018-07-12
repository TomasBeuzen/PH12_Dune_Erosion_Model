#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions supporting use of Gaussian Process

@author: TomasBeuzen, 2018
"""

import numpy as np
import pandas as pd
import scipy as sp
from sklearn.preprocessing import StandardScaler


def train_test_mda(df_X, df_y, sample_size, dist_measure, standardize):
    '''
    Inputs:
    df_X = dataframe of X data
    df_Y = dataframe of y data matching index of df_X
    sample_size = fraction of data to select using mda (0 to 1)
    dist_measure = distance measure for mda algorithm ('euclidean'/'mahalanobis')
    standardize = standardize prior to mda algorithm (True/False)
    Outputs:
    X_train = the X data training sample selected using mda
    X_test = the X data testing sample selected using mda
    y_train = the y data training sample selected using mda
    y_test = the y data testing sample selected using mda
    '''
    if sample_size <= 0 or sample_size >= 1:
        print('Error, wrong sample size: Please enter a value between 0 and 1.')
    else:
        # Pre-processing
        df = df_X.join(df_y) # join data for dissimilarity calcs
        sample_size = round(sample_size*df_X.shape[0])
        if standardize:
            scaler = StandardScaler()
            scaler.fit(df)
            d = sp.spatial.distance.cdist(scaler.transform(df), scaler.transform(df), metric=dist_measure)
        else:
            d = sp.spatial.distance.cdist(df, df, metric=dist_measure)  # pairwise dissimilarity
        # Select first training sample
        md = np.sum(d, axis=0)
        samp_idx = [np.argmax(md)]
        # Loop to select remaining samples
        for _ in range(1,sample_size):
            samp_idx.append(np.argmax(np.min(d[samp_idx],axis=0)))
        # Split back into X and y train and test sets for output
        X_train = df.iloc[samp_idx].drop(columns='runup')
        X_test = df.drop(columns='runup',index=df.index[samp_idx])
        y_train = df[['runup']].iloc[samp_idx]
        y_test = df[['runup']].drop(df.index[samp_idx])
            
    return X_train, X_test, y_train, y_test
