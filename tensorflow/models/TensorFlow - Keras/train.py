import sklearn
from sklearn.metrics import confusion_matrix, precision_score
from sklearn.model_selection import train_test_split
from keras import optimizers
from keras.layers import Dense,Dropout
from keras.models import Sequential
from keras.regularizers import l2
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

data = pd.read_csv("D:\\Documents\\Development\\private\\AIAPBRT-BachelorProject\\data\\set2\\completeSet.csv", sep=",")
print(data.head())

# seperate Data
x = data.drop(columns=['v'])
y = data['v']

x_train, x_test, y_train, y_test = train_test_split(x,y,test_size=0.10, random_state=0)
print(x_train.shape,y_train.shape,x_test.shape,y_test.shape)

# define the model
model = Sequential()

#hid layer-1
model.add(Dense(600, activation='relu', input_dim=10, kernel_regularizer=l2(0.01)))
model.add(Dropout(0.3, noise_shape=None, seed=None))

#hid layer-2
model.add(Dense(600, activation='relu', kernel_regularizer=l2(0.01)))
model.add(Dropout(0.3, noise_shape=None, seed=None))

#hid layer-3
model.add(Dense(400, activation='elu', kernel_regularizer=l2(0.01)))
model.add(Dropout(0.3, noise_shape=None, seed=None))

# Output layer
model.add(Dense(1,activation='sigmoid'))	# since the labels to predict are either 0 or 1

opt = optimizers.Adam(lr=0.002, beta_1=0.9, beta_2=0.999, epsilon=None, decay=0.0, amsgrad=False)
model.compile(loss='binary_crossentropy',optimizer=opt, metrics=['accuracy'])	# using negative log los and adam
print(model.summary())

# Train the Model
model_output = model.fit(x_train, y_train, epochs=100,batch_size=50, verbose=1, validation_data=(x_test, y_test),)
print('Training Accuracy : ' , np.mean(model_output.history["acc"]))
print('Validation Accuracy : ', np.mean(model_output.history["val_acc"]))

# prediction and check precision
y_pred = model.predict(x_test)
rounded = [round(x[0]) for x in y_pred]
y_pred1 = np.array(rounded, dtype='int64')

confusion_matrix(y_test,y_pred1)
precision_score(y_test,y_pred1)

history_dict = model_output.history
print(history_dict.keys())

acc = history_dict['acc']
val_acc = history_dict['val_acc']
loss = history_dict['loss']
val_loss = history_dict['val_loss']

epochs = range(1, len(acc) + 1)

# "bo" is for "blue dot"
plt.plot(epochs, loss, 'bo', label='Training loss')
# b is for "solid blue line"
plt.plot(epochs, val_loss, 'b', label='Validation loss')
plt.title('Training and validation loss')
plt.xlabel('Epochs')
plt.ylabel('Loss')
plt.legend()

plt.show()

plt.clf()   # clear figure

plt.plot(epochs, acc, 'bo', label='Training acc')
plt.plot(epochs, val_acc, 'b', label='Validation acc')
plt.title('Training and validation accuracy')
plt.xlabel('Epochs')
plt.ylabel('Accuracy')
plt.legend()

plt.show()

model.save("Classifier_002.h5")
