function [X, P] = ekf(A,C,Q,H,R,x_t,Y_t,pHat)
 %dt: time second
 %x_t previous predicted value
 %Y_t: measured value
 %X state estimated value
 %P covariance estimated state
 

    xHat = A*x_t;
    pHat = A*pHat*A' + Q;
    KG = (pHat*H)/((H*pHat*H') + R);
    Y = C*Y_t; 
    X = xHat + KG*(Y-H*xHat);
    P = (eye(9)-KG*H) * pHat;    


    
    