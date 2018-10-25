function [A,Q,R,H] = preDefinedMatrix2(dt,algorithm)


% State transition matrix
    % State prediction
    A = [ 1 0 0 dt 0 0 dt^2/2 0 0 ;
          0 1 0 0 dt 0 0 dt^2/2 0 ; 
          0 0 1 0 0 dt 0 0 dt^2/2 ; 
          0 0 0 1 0 0 dt 0 0 ;
          0 0 0 0 1 0 0 dt 0 ;
          0 0 0 0 0 1 0 0 dt ;
          0 0 0 0 0 0 1 0 0 ;
          0 0 0 0 0 0 0 1 0 ;
          0 0 0 0 0 0 0 0 1];
  
    Q =       [.02^2 0 0 0 0 0 0 0 0;
               0 .02^2 0 0 0 0 0 0 0;
               0 0 .12^2 0 0 0 0 0 0;
               0 0 0 .01^2 0 0 0 0 0;
               0 0 0 0 .01^2 0 0 0 0;
               0 0 0 0 0 .1^2 0 0 0;
               0 0 0 0 0 0 .01^2 0 0;
               0 0 0 0 0 0 0 .01^2 0;
               0 0 0 0 0 0 0 0 .1^2];
    
if strcmp(algorithm,'Least_square')||strcmp(algorithm,'Trilateration')
    r1 = 0.1063;
    r2 = 0.1012;
    r3 = 0.3558;
elseif strcmp(algorithm,'Recursive_Trilateration')
    r1 = 0.2032;
    r2 = 0.0944;
    r3 = 0.6848;       
elseif strcmp(algorithm,'Gradient_Descent')
    r1 = 0.1503;
    r2 = 0.1298;
    r3 = 0.2314;                                                  
end  
    r7 = 0;
    r8 = 0;
    r9 = 1.0394;
    R = diag([r1 r2 r3 0 0 0 r7 r8 r9]);
    H = diag([1 1 1 1 1 1 1 1 1]);