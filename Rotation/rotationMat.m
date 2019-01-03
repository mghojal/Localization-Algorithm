function [Rx Ry Rz R] = rotationMat(angle)
phi  = angle(1);
theta  = angle(2);
psi= angle(3);

Rx = [  1   ,     0     ,    0     ;
        0   , cos(phi) ,-sin(phi);
        0   , sin(phi) , cos(phi)];
    
    
Ry = [cos(theta) , 0 , sin(theta);
      0           , 1 , 0          ;
      -sin(theta) , 0 , cos(theta) ];
  
  
Rz = [cos(psi) , -sin(psi) , 0 ;
      sin(psi) ,  cos(psi) , 0 ;
           0    ,      0     , 1 ];
       
R = [cos(angle(2))*cos(angle(3)) sin(angle(1))*sin(angle(2))*cos(angle(3))-cos(angle(1))*sin(angle(3)) cos(angle(1))*sin(angle(2))*cos(angle(3))+sin(angle(1))*sin(angle(3));
     cos(angle(2))*sin(angle(3)) sin(angle(1))*sin(angle(2))*sin(angle(3))+cos(angle(1))*cos(angle(3)) cos(angle(1))*sin(angle(2))*sin(angle(3))-sin(angle(1))*cos(angle(3));
        -sin(angle(2))           sin(angle(1))*cos(angle(2))                 cos(angle(1))*cos(angle(2))            ];

