function [A,C,Q,R,H] = preDefinedMatrix(dt,algorithm,path)


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



% Measurement matrix
    % The expected measurement given the predicted state
    C = [1 0 0 0 0 0;
         0 1 0 0 0 0;
         0 0 1 0 0 0;
         0 0 0 0 0 0
         0 0 0 0 0 0
         0 0 0 0 0 0
         0 0 0 1 0 0;
         0 0 0 0 1 0;
         0 0 0 0 0 1];
    


if strcmp(algorithm,'Gradient_Descent')
    if strcmp(path,'pathOption1.mat')
        r1 = 0.0848;
        r2 = 0.1555;
        r3 = 0.1139;
        r7 = 7.001;
        r8 = 11.7603;
        r9 = .0639;
    elseif strcmp(path,'pathOption2.mat')
        r1 = 0.1692;
        r2 = 0.1991;
        r3 = 0.1507;
        r7 = 2.6988;
        r8 = 2.4881;
        r9 = .0657;        
    elseif strcmp(path,'pathOption3.mat')
        r1 = 0.2572;
        r2 = 0.1602;
        r3 = 0.1612;
        r7 = 2.925;
        r8 = 1.5712;
        r9 = .0567;             
    else
        r1 = 0.3645;
        r2 = 0.3755;
        r3 = 0.3027;
        r7 = 2.5505;
        r8 = 1.7688;
        r9 = 1.4553;                                     
    end
elseif strcmp(algorithm,'Recursive_Least_square')
    if strcmp(path,'pathOption1.mat')
        r1 = 0.1046;
        r2 = 0.1947;
        r3 = 0.3307;
        r7 = 7.0004;
        r8 = 11.7610;
        r9 = .063;
    elseif strcmp(path,'pathOption2.mat')
        r1 = 0.3129;
        r2 = 0.1721;
        r3 = 1.8631;
        r7 = 2.7077;
        r8 = 2.4932;
        r9 = .0657;        
    elseif strcmp(path,'pathOption3.mat')
        r1 = 0.3176;
        r2 = 0.1269;
        r3 = 1.9883;
        r7 = 2.9103;
        r8 = 1.5768;
        r9 = .0515;             
    else
        r1 = 0.3040;
        r2 = 0.3086;
        r3 = 1.4204;
        r7 = 2.5563;
        r8 = 1.7695;
        r9 = 1.4553;                                     
    end    
elseif strcmp(algorithm,'Least_square')
    if strcmp(path,'pathOption1.mat')
        r1 = 0.1183;
        r2 = 0.1772;
        r3 = 0.3891;
        r7 = 7.0004;
        r8 = 11.7613;
        r9 = .0628;
    elseif strcmp(path,'pathOption2.mat')
        r1 = 0.2258;
        r2 = 0.2414;
        r3 = 0.4;
        r7 = 2.6970;
        r8 = 2.4861;
        r9 = .064;        
    elseif strcmp(path,'pathOption3.mat')
        r1 = 0.2727;
        r2 = 0.1729;
        r3 = 0.4;
        r7 = 2.987;
        r8 = 1.5754;
        r9 = .0253;             
    else
        r1 = 0.3589;
        r2 = 0.3198;
        r3 = 0.7341;
        r7 = 2.5520;
        r8 = 1.7692;
        r9 = 1.4552;                                     
    end    
end

    R = diag([r1 r2 r3 0 0 0 r7 r8 r9]);
    Q =       [0.1^2 0 0 0 0 0 0 0 0;
               0 0.1^2 0 0 0 0 0 0 0;
               0 0 0.01^2 0 0 0 0 0 0;
               0 0 0 0.5^2 0 0 0 0 0;
               0 0 0 0 0.4^2 0 0 0 0;
               0 0 0 0 0 0.1^2 0 0 0;
               0 0 0 0 0 0 0.5^2 0 0;
               0 0 0 0 0 0 0 0.4^2 0;
               0 0 0 0 0 0 0 0 0.1^2];

    H = diag([1 1 1 1 1 1 1 1 1]);