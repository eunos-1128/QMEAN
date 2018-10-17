# QMEAN:
# ----------
# 
# Copyright (C) 2018 University of Basel and the Authors.
# Authors: Gabriel Studer, Pascal Benkert, Marco Biasini
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import numpy as np

class Regressor:

    def __init__(self, filepath):

        fh = open(filepath, 'rb')

        magic_number = np.fromfile(fh, dtype=np.int32, count=1)[0]
        if(magic_number != 444222):
          raise RuntimeError("Inconsistency when reading mlp regressor!")

        version = np.fromfile(fh, dtype=np.int32, count=1)[0]
        if version != 1:
          raise RuntimeError("Inconsistent file version when reading mlp \
                              regressor!")

        self.n_layers = np.fromfile(fh, dtype=np.int32, count=1)[0]
        self.layer_sizes = np.fromfile(fh, dtype=np.int32, count=self.n_layers)

        if self.layer_sizes[-1] != 1:
            raise RuntimeError("Expect a layer size of 1 for last (output) "\
                               "layer.")

        self.activation_functions = np.fromfile(fh, dtype=np.int32, 
                                                count=self.n_layers)

        mean = np.fromfile(fh, dtype=np.float32, count=self.layer_sizes[0])
        self.mean = np.asmatrix(mean.reshape(self.layer_sizes[0], 1))
        std = np.fromfile(fh, dtype=np.float32, count=self.layer_sizes[0])
        self.one_over_std = np.asmatrix(std.reshape(self.layer_sizes[0], 1))
        for i in range(self.layer_sizes[0]):
            self.one_over_std[(i, 0)] = 1.0 / self.one_over_std[(i, 0)]


        self.bias = list()
        for i in range(1, self.n_layers):
            val = np.fromfile(fh, dtype=np.float32, count = self.layer_sizes[i])
            self.bias.append(np.asmatrix(val.reshape(self.layer_sizes[i], 1)))


        self.weights = list()
        for i in range(1, self.n_layers):
            val = np.fromfile(fh, dtype=np.float32, 
                              count=self.layer_sizes[i-1] * self.layer_sizes[i])
            self.weights.append(np.asmatrix(val).reshape(self.layer_sizes[i],
                                                         self.layer_sizes[i-1]))

        fh.close()


    def Predict(self, features):

        features = np.matrix(features).reshape(self.layer_sizes[0], 1)
        layer = np.multiply(features - self.mean, self.one_over_std)
        self._Activate(layer, 0)
        for i in range(1, self.n_layers):
            layer = self.weights[i-1] * layer + self.bias[i-1]
            self._Activate(layer, i)

        return layer[(0,0)]


    def _Activate(self, layer, layer_idx):

        if self.activation_functions[layer_idx] == 0:
            return
        elif self.activation_functions[layer_idx] == 1:
            layer[layer < 0] = 0
        else:
            raise RuntimeError("Observed invalid activation function in \
                                loaded regressor!")

