# Copyright (c) 2013-2020, SIB - Swiss Institute of Bioinformatics and
# Biozentrum - University of Basel
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import numpy as np


class Regressor:
    def __init__(self, n_input_features, mean=None, std=None):
        self._n_input_features = n_input_features
        if mean is None:
            self._mean = np.zeros((self._n_input_features, 1))
        else:
            m = np.array(mean)
            if m.size != self._n_input_features:
                raise ValueError(
                    "Number of elements in mean must be consistent with n_input_features"
                )
            self._mean = m.reshape(self._n_input_features, 1)

        if std is None:
            self._std = np.ones((self._n_input_features, 1))
        else:
            m = np.array(std)
            if m.size != self._n_input_features:
                raise ValueError(
                    "Number of elements in std must be consistent with n_input_features"
                )
            self._std = m.reshape(self._n_input_features, 1)

        # internally, the first layer is an input layer for data normalization
        self._n_layers = 1
        self._layer_sizes = [self._n_input_features]
        self._activation_functions = [0]
        self._bias = [None]
        self._weights = [None]

    def AddLayer(self, weights, activation_function, bias=None):
        if weights.shape[1] != self._layer_sizes[-1]:
            raise ValueError(
                "Number of cols in weights must be consistent with size of previous layer"
            )

        if activation_function not in [0, 1]:
            raise ValueError("Expect activation_function to be in [0,1]")

        if bias is None:
            self._bias.append(np.zeros((weights.shape[0], 1)))
        else:
            m = np.array(bias)
            if m.size != weights.shape[0]:
                raise ValueError(
                    "Number of rows in weights must be consistent with size of bias"
                )
            self._bias.append(m.reshape(m.size, 1))

        self._weights.append(np.array(weights))
        self._activation_functions.append(activation_function)
        self._n_layers += 1
        self._layer_sizes.append(self._weights[-1].shape[0])

    @staticmethod
    def Load(filepath):
        with open(filepath, "rb") as fh:
            magic_number = np.fromfile(fh, dtype=np.int32, count=1)[0]
            if magic_number != 444222:
                raise RuntimeError("Inconsistency when reading mlp regressor!")
            version = np.fromfile(fh, dtype=np.int32, count=1)[0]
            if version != 1:
                raise RuntimeError(
                    "Inconsistent file version when reading mlp regressor!"
                )
            n_layers = np.fromfile(fh, dtype=np.int32, count=1)[0]
            layer_sizes = np.fromfile(fh, dtype=np.int32, count=n_layers)
            activations = np.fromfile(fh, dtype=np.int32, count=n_layers)
            mean = np.fromfile(fh, dtype=np.float32, count=layer_sizes[0])
            std = np.fromfile(fh, dtype=np.float32, count=layer_sizes[0])
            bias = list()
            for i in range(1, n_layers):
                n = layer_sizes[i]
                v = np.fromfile(fh, dtype=np.float32, count=n)
                bias.append(v.reshape(layer_sizes[i], 1))
            weights = list()
            for i in range(1, n_layers):
                n = layer_sizes[i - 1] * layer_sizes[i]
                v = np.fromfile(fh, dtype=np.float32, count=n)
                weights.append(v.reshape(layer_sizes[i], layer_sizes[i - 1]))

        regressor = Regressor(layer_sizes[0], mean, std)
        for w, b, a in zip(weights, bias, activations[1:]):
            regressor.AddLayer(w, a, b)

        return regressor

    def Save(self, filepath):

        if self._layer_sizes[-1] != 1:
            raise RuntimeError("Expect size of 1 for last (output) layer.")

        with open(filepath, "wb") as fh:
            np.array([444222], dtype=np.int32).tofile(fh)
            np.array([1], dtype=np.int32).tofile(fh)
            np.array([self._n_layers], dtype=np.int32).tofile(fh)
            np.array([self._layer_sizes], dtype=np.int32).tofile(fh)
            np.array([self._activation_functions], dtype=np.int32).tofile(fh)
            self._mean.astype(np.float32).tofile(fh)
            self._std.astype(np.float32).tofile(fh)
            for b in self._bias[1:]:
                b.astype(np.float32).tofile(fh)
            for w in self._weights[1:]:
                w.astype(np.float32).tofile(fh)



    def Predict(self, features):
        features = np.array(features).reshape(self._layer_sizes[0], 1)
        layer = np.divide(features - self._mean, self._std)
        self._Activate(layer, 0)
        for i in range(1, self._n_layers):
            layer = np.dot(self._weights[i], layer) + self._bias[i]
            self._Activate(layer, i)

        return layer[(0, 0)]

    def _Activate(self, layer, layer_idx):
        if self._activation_functions[layer_idx] == 0:
            return
        elif self._activation_functions[layer_idx] == 1:
            layer[layer < 0] = 0
        else:
            raise RuntimeError(
                "Observed invalid activation function in loaded regressor!"
            )
