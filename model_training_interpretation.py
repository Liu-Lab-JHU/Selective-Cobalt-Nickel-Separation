#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Use the whole dataset to train a random forest model and use SHAP to interpret the model
(output the important substructure of bioacids).
"""
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
import shap


"""
Step 1: data processing
"""
df= pd.read_csv('dataset_w_WECFPs.csv')
df=df.drop(columns=['acid','SMILES','concentration'])

X = df.drop('selectivity', axis=1)
y = df['selectivity']

"""
Step 2: train a random forest regression model
"""
random_forest_model = RandomForestRegressor(random_state=42)
model = random_forest_model.fit(X, y)

"""
Step 3: use SHAP to interpret the model
"""
explainer = shap.Explainer(model)
shap_values = explainer(X)

# Visualization
shap.summary_plot(shap_values.values, X, max_display = 6)




