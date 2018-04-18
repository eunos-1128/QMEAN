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
        self.activation_functions = np.fromfile(fh, dtype=np.int32, 
                                                count=self.n_layers)

        self.mean = np.fromfile(fh, dtype=np.float32, count=self.layer_sizes[0])
        self.var = np.fromfile(fh, dtype=np.float32, count=self.layer_sizes[0])

        self.bias = list()
        for i in range(1, self.n_layers):
            self.bias.append(np.fromfile(fh, dtype=np.float32, 
                             count=self.layer_sizes[i]))

        self.weights = list()
        for i in range(1, self.n_layers):
            self.weights.append(list())
            for j in range(self.layer_sizes[i]):
                self.weights[-1].append(np.fromfile(fh, dtype=np.float32, 
                                        count=self.layer_sizes[i - 1]))

        fh.close()


    def Predict(self, features):

        if len(features) != self.layer_sizes[0]:
            raise ValueError("Number of input features must be consistent with \
                              nodes in input layer (%i)"%(self.layer_sizes[0]))

        # transform input features
        prev_layer = np.ndarray(self.layer_sizes[0])
        for i in range(self.layer_sizes[0]):
            prev_layer[i] = (features[i] - self.mean[i]) / self.var[i]

        self._Activate(prev_layer, 0)

        for i in range(1, self.n_layers):

            # assign the bias
            current_layer = self.bias[i-1].copy()

            # do the perceptron magic
            for j in range(len(current_layer)):              
                current_layer[j] += prev_layer.dot(self.weights[i-1][j]) 
        
            # apply the activation function
            self._Activate(current_layer, i)

            # this layer is done, lets proceed one layer
            prev_layer = current_layer

        return prev_layer[0]


    def _Activate(self, layer, layer_idx):
        if self.activation_functions[layer_idx] == 0:
            return
        elif self.activation_functions[layer_idx] == 1:
            layer[layer < 0] = 0
        else:
            raise RuntimeError("Observed invalid activation function in \
                                loaded regressor!")
