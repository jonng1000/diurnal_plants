# -*- coding: utf-8 -*-
"""
Created on Thu May 24 16:07:23 2018

@author: weixiong001
"""
'''
This script is to check that all day/night timepoints make sense
'''
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv

colormap_dict = {}
with open('colormap_values_PCA.csv', 'r', newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter='\t')
    next(reader)
    for row in reader:
        new_row = row[0].split(',')
        colormap_dict[new_row[0]] = new_row[1]
    