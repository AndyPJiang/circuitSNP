import pandas as pd
import numpy as np
import keras
from sklearn.model_selection import train_test_split
from models import Model
import os

def split_data(seed,df):
    TEST_SPLIT = 0.125

    # only keep features 
    x_data = df.iloc[:,3:-1]
    y_data = df.iloc[:,-1]

    X_train, X_test, y_train, y_test = train_test_split(x_data, y_data, test_size=TEST_SPLIT, random_state=seed)
    return X_train, X_test, y_train, y_test



def train_all_models(path, X_train, y_train):
    # create 10 replicates of each 6 models

    # hyperparameters
    EPOCHS = 50
    INPUT_DIM = 1372
    BATCH_SIZE = 64
    VALIDATION_SPLIT = 0.125
    models = [(200,40),(50,10),(10,5),(5,3)]


    nn_models = Model(INPUT_DIM)
    early_stopping = keras.callbacks.EarlyStopping(monitor='val_loss', min_delta=0, patience=3, verbose=1, mode='auto')


    for i in range(10):
        for (hid_layer1_units, hid_layer2_units) in models:
            model = nn_models.create_model(hid_layer1_units,hid_layer2_units)

            model.fit(X_train, y_train, epochs=EPOCHS, batch_size=BATCH_SIZE,verbose=1,callbacks=[early_stopping],validation_split=VALIDATION_SPLIT, shuffle=True)
            model.save_weights(path+'model({}+{})_replicate_{}.h5'.format(hid_layer1_units, hid_layer2_units, i))
        

        # these models have layers with 0 neurons, so need to do them separately
        model = nn_models.create_model_simple0()
        model.fit(X_train, y_train, epochs=EPOCHS, batch_size=BATCH_SIZE,verbose=1,callbacks=[early_stopping],validation_split=VALIDATION_SPLIT, shuffle=True)
        model.save_weights(path+'model(0+0)_replicate_{}.h5'.format(i))
        
        model = nn_models.create_model_simple1()
        model.fit(X_train, y_train, epochs=EPOCHS, batch_size=BATCH_SIZE,verbose=1,callbacks=[early_stopping],validation_split=VALIDATION_SPLIT, shuffle=True)
        model.save_weights(path+'model(1+0)_replicate_{}.h5'.format(i))


def driver(ROOT_DIR, df, seed=42):

    X_train, X_test, y_train, y_test = split_data(seed, df)
    
    if not os.path.exists('{}models/models_seed={}'.format(ROOT_DIR,seed)):
        os.mkdir('{}models/models_seed={}'.format(ROOT_DIR,seed))
        
        print("Training models with seed = {}".format(seed))
        train_all_models('{}models/models_seed={}/'.format(ROOT_DIR,seed), X_train, y_train)
    else:
        print("Models already trained for seed = {}".format(seed))