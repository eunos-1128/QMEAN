import pandas as pd
import numpy as np
from qmean.mlp_regressor import TrainRegressor

# Example training of multi-layer perceptron on a toy data set.
# The Boston House Price Dataset involves the prediction of a 
# house price in thousands of dollars given details of the house 
# and its neighborhood.
df = pd.read_csv('example_data/housing.csv')

df_train = df.loc[:400]
df_test = df.loc[400:]

# search for boston housing data set in the net... I'm sure 
# you'll find a description of the single features
features = ['CRIM', 'ZN', 'INDUS', 'CHAS', 'NOX', 'RM', 'AGE', 
            'DIS', 'RAD', 'TAX', 'PTRATIO', 'B', 'LSTAT']

# thats the median value of the houses we want to predict
target = 'MEDV'

# define architecture and training parameters
topology = [len(features), 20, 20, 1]
loss_function = 'mean_squared_error'
optimizer = 'adam'
epochs = 100
batch_size = 10

# train and predict 
regressor = TrainRegressor(df_train, features, target,loss_function, 
                           optimizer, topology, epochs, batch_size)
regressor_in = df_test[features].values
predictions = np.zeros(regressor_in.shape[0])
for idx in range(regressor_in.shape[0]):
    predictions[idx] = regressor.Predict(regressor_in[idx])

# estimate root mean square error
ref = df_test[target].values
rmse = np.sqrt(np.mean(np.square(predictions-df_test[target].values)))
print("testing rmse:", rmse)

