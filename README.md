# Localization-Algorithm
Sensor fusion (UWB+IMU+Ultrasonic), using Kalman filter and 3 different multilateration algorithms (Least square and Recursive Least square and gradient descent)

> in this study we take readings from IMU( accelerometer and gyroscope) calculate linear acceleration to be an input to kalamn filter

> we take UWB distances reading feed it to multilateration algorithm (there is a choice of three: Least Square, Recurisve Least Square, Gradient Descent) to get first estimated position as an input for Kalman filter

> we feed also ultrasonic to our Kalman filter


> Main file is the main which we can specify all our options on it in addition it refere to the test paths

> We need to specify the path and the used multilateration algorithm for UWB distances reads

> In cas we choose gradient descent, we need to specify if we will consider boundaries as our whole working space or minizing it by shrink the boudaries using the previous step.

Speical study done for horizontal movement which can be found in folder /Ultrasonic_study
