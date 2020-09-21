import keras
from keras.models import Sequential
from keras.layers import Dense

class Model:
    def __init__(self,input_dim):
        self.input_dim = input_dim
        self.adadelta = keras.optimizers.Adadelta(lr=1.0, rho=0.95, epsilon=None, decay=0.0)


    # create a 4 layer NN (2 hidden layers) for binary classification indicating whether window is open or closed
    def create_model(self, hid_layer1_units, hid_layer2_units):

        model = Sequential()
        model.add(Dense(units=hid_layer1_units, activation='relu', input_dim=self.input_dim))
        model.add(Dense(units=hid_layer2_units, activation='relu'))
        model.add(Dense(units=1, activation='sigmoid'))

        model.compile(loss='binary_crossentropy',
                    optimizer=self.adadelta,
                    metrics=['accuracy'])
        return model

    # 1 hidden layer has no neurons
    def create_model_simple1(self):
        model = Sequential()
        model.add(Dense(units=1, activation='relu', input_dim=self.input_dim))
        model.add(Dense(units=1, activation='sigmoid'))

        model.compile(loss='binary_crossentropy',
                    optimizer=self.adadelta,
                    metrics=['accuracy'])
        return model


    # both hidden layers have no neurons, basically a linear model plus a sigmoid activation function
    def create_model_simple0(self):
        model = Sequential()
        model.add(Dense(units=1, activation='sigmoid', input_dim=self.input_dim))

        model.compile(loss='binary_crossentropy',
                    optimizer= self.adadelta,
                    metrics=['accuracy'])
        return model
