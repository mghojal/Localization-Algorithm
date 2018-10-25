close all
clear all
clc

addpath('Rotation');
addpath('ErrorCalculation');
%% choosen path
%options are the following:

%   1: reads of 5 anchors with different heights using update firmware with 4Hz update rate of UWB
%      raw data system and have the following path (Experiment 3)
%           ---------------------
%           |                   |
%           |                   |
%           ---------------------

%   2: reads of 4 anchors with the same heights using original firmware with 10Hz update rate of
%   UWB raw data system and have the following path (Experiment1)
%           ---------------------
%           |                   |
%           |                   |
%           ---------------------

%   3: reads of 4 anchors with the same heights using original firmware with 10Hz update rate of
%   UWB raw data system and have the following path (Extra1)
%           |
%           ---------------------
%                               |
%           ---------------------
%           |
%           ---------------------

%   4: reads of 4 anchors with different heights using original firmware with 10Hz update rate of
%   UWB raw data system and have the following path (Experiment 2)
%           ---------------------
%           |                   |
%           |                   |
%           
%   5: reads of 4 anchors with different heights using original firmware with 10Hz update rate of
%   UWB raw data system and have the following path (Extra 2)
%           ----------
%           | \    / |
%           |  \  /  |
%           |   \/   |
%           |   /\   |
%           |  /  \  |
%           | /    \ |
%           ----------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path = 'pathOption1.mat'; %Exp3
% path = 'pathOption2.mat'; %Exp1
% path = 'pathOption3.mat'; %Extra1
% path = 'pathOption4.mat'; %Exp2
% path = 'pathOption5.mat'; %Extra2

load(path);

%% calculating positions using mlat (3 algorithms)
% Algorithms are four:  1- Gradient Descent
%                       2- Recursive Least Square
%                       3- Least Square
% we need to choose algorithm:

    Algorithm = 'Gradient_Descent';
%     Algorithm = 'Recursive_Least_square';
%     Algorithm = 'Least_square';

%% Options by running the code
% yes = 1
% No = 0
without_IMU = 0;
Object_surround_by_point = 0;
Collision_detection = 0;
printout_video_of_simulation = 0;

%% Beginning of the Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Code Start  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Variable definitions
% positionEstimate: estimated postion usig an algorithm of multilateration
% filtered_angles:  the angles that will be using in rotation matrix
% dt:               Time interval the system in case 4 anchros and for only
%                   sensors in case 5 anchors
% bounds:           when using gradient decent to define bounds that
%                   estimated postion should be in
% iteration:        loop iteration of Kalman filter
% anchorsNumber:    number of anchors used in the system

positionEstimate=[];
filtered_angles = [];
earth = [0 0 9.8]';
iteration = length(acc);
anchorsNumber = size(dis,2);

%% Synchornization procedure for reads from UWB system and sensors
%defining timeStamp for distances in matrix dtDis 
if anchorsNumber == 5
    for i=2:length(dis)
        dtT(i) = tI(i,2) - tI(i-1,2);
        if dtT(i)<0
            dtT(i) = dtT(i)+60;
        end
        if dtT(i)>1
            break;
        end
    end
    dt=1/SampleRateSensors;
else
    dt=1/SampleRateSensors;
end

%% starting with process

[A,C,Q,R,H] = preDefinedMatrix(dt,Algorithm,path); 
k=1;
kIteration = length(dis);

for i=1:iteration
    if i==1
        prePos = [];
    else
%         prePos = [];
        prePos = F(i-1,1:3);
    end
    if anchorsNumber == 5
        if k<=kIteration
            positionEstimate(k,:) = mlat.do_main(Algorithm, anchorLoc, dis(k,:), prePos);
        end
    else
        positionEstimate(i,:) = mlat.do_main(Algorithm, anchorLoc, dis(i,:), prePos);
    end
    %% Initial values for i==1 for kalman filtering and angles calculation
    if (i==1) 
        if anchorsNumber==5
            x_t =  [positionEstimate(k,:) 0 0 0 0 0 0]';
        else
            x_t =  [positionEstimate(i,:) 0 0 0 0 0 0]';
        end
        pHat = [dt^4/4 0 0 dt^3/2 0 0 dt^2/2 0 0;
                0 dt^4/4 0 0 dt^3/2 0 0 dt^2/2 0;
                0 0 dt^4/4 0 0 dt^3/2 0 0 dt^2/2;
                dt^3/2 0 0 dt^2 0 0 dt 0 0;
                0 dt^3/2 0 0 dt^2 0 0 dt 0;
                0 0 dt^3/2 0 0 dt^2 0 0 dt;
                dt^2/2 0 0 dt 0 0 1 0 0;
                0 dt^2/2 0 0 dt 0 0 1 0;
                0 0 dt^2/2 0 0 dt 0 0 1];

        pre_gyro_angle = [0 0 0];
        filtered_angles = [0 0 0];
        Super_filtered_angles = [0 0 0];      
    else
        x_t = X;
        pHat = P;
        [filtered_angles, gyro_angle, acc_angles] = angles(acc(i,1:3),gyro(i,1:3),dt,pre_gyro_angle);
        pre_gyro_angle = gyro_angle;
%         if anchorsNumber==5
%             if k<=kIteration
%                 Super_filtered_angles = (filtered_angles + [atan2(positionEstimate(k,2),positionEstimate(k,3)) atan2(positionEstimate(k,3),positionEstimate(k,1)) atan2(positionEstimate(k,2),positionEstimate(k,1))])/2;
%             end
%         else
%             Super_filtered_angles = (filtered_angles + [atan2(positionEstimate(i,2),positionEstimate(i,3)) atan2(positionEstimate(i,3),positionEstimate(i,1)) atan2(positionEstimate(i,2),positionEstimate(i,1))])/2;
%         end
    end
    %% calculate rotation matrix to delete gravity effect of accelerometer values
    [Rx Ry Rz Rot] = rotationMat(filtered_angles);
    %% acceleration without gravity value
    acc2(i,1:3) = (inv(Rot)*acc(i,1:3)')' - earth';
    anorm = norm(acc2(i,1:3));
    acc(i,1:3) = -acc2(i,1:3)/anorm;
    if without_IMU == 1
        acc(i,1:3) = [0 0 0];
    end
    %% reading values position after applying Multilateration algorithms and acceleration 
    if anchorsNumber==5
        j=0;
        if k<=kIteration
            if tI(k,1) ==  acc(i,5)
                if tI(k,2) == acc(i,6)
                    Y_t = [positionEstimate(k,:) acc(i,1:3)]';          
                    k=k+1;
                    j=1;
                else
                    Y_t = [positionEstimate(k,:)+(acc(i,1:3))*dt^2, acc(i,1:3)]';            
                end
            end    
        end
    else
        Y_t = [positionEstimate(i,:) acc(i,1:3)]';
    end

    %% Kalman Filtering and updating Uncertainty for measures and Estimation
    [X, P] = ekf( A, C, Q, H, R , x_t, Y_t, pHat);

    %% Results saving
    if anchorsNumber==5
        F(i,1:3) = positionEstimate(k-1,:);
    else
        F(i,1:3) = positionEstimate(i,:);
    end   
%     diff(I*length(acc) + i,:) = F(i,1:3)-posReal(k-1,:);
    F(i,4:6) = X(1:3);  

    %% Calculating error with Error Calculation
    if anchorsNumber ==5
        if j==1
            if k==2
                ErrorMatrix(k-1,:) = [0 0];
                GG(i,:) = [0 0 0];
            else
                [eM eMKF] = errorFun(posReal,F,i,k-1);
                ErrorMatrix(k-1,:) = [eM eMKF];
                GG(i,:) = abs(F(i,1:3)-posReal(k-1,1:3));
            end
        end
    else
        [eM eMKF] = errorFun(posReal,F,i,i);
        ErrorMatrix(i,1:2) = [eM eMKF];
    end
end
%% Error Calculation for x,y,z mean, standard deviation
ErrorM = ErrorCalculation(F(:,4:6),posReal, anchorsNumber);
MeanErrorMX = mean(ErrorM(:,1));
MeanErrorMY = mean(ErrorM(:,2));
MeanErrorMZ = mean(ErrorM(:,3));
SDErrorMX = std(ErrorM(:,1));
SDErrorMY = std(ErrorM(:,2));
SDErrorMZ = std(ErrorM(:,3));

%% Analysis for Error
meanAlgorithm = mean(ErrorMatrix(9:end,1));
maxAlgorithm = max(ErrorMatrix(9:end,1));
minAlgorithm = min(ErrorMatrix(9:end,1));
fprintf('Error Analysis without Kalman Filter\n mean: %f \t max: %f \t min: %f \n\n'...
,meanAlgorithm, maxAlgorithm, minAlgorithm);
fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++\n\n');
meanKF_Algorithm = mean(ErrorMatrix(9:end,2));
maxKF_Algorithm = max(ErrorMatrix(9:end,2));
minKF_Algorithm = min(ErrorMatrix(9:end,2));
fprintf('Error Analysis with Kalman Filter\n mean: %f \t max: %f \t min: %f \n'...
,meanKF_Algorithm, maxKF_Algorithm, minKF_Algorithm');

%% Plotting section
%plotting type
plotType = 0;  % 0 => 3D
               % 1 => XY plane
               % 2 => YZ plane
               % 3 => XZ plane
               
[aa,ab] = definedFigure2(anchorLoc,1,Algorithm);
k=2;

%vertices and faces of cube
if Object_surround_by_point == 1
    ObjCentralPoint = [F(90,4:5) 1];
    [V Fa] = cubeData(F(1,:),MeanErrorMX,MeanErrorMY,MeanErrorMZ,SDErrorMX,SDErrorMY,SDErrorMZ);
else
    ObjCentralPoint = [-50 -50 -50];
    [V Fa] = cubeData([0 0 0 ObjCentralPoint],MeanErrorMX,MeanErrorMY,MeanErrorMZ,SDErrorMX,SDErrorMY,SDErrorMZ);
end
[V2 Fa2] = cubeData([0 0 0 ObjCentralPoint],0.01,0.01,0.01,0.05,0.05,ObjCentralPoint(3)/3);

S1.Vertices = V;
S1.Faces = Fa;
S1.FaceVertexCData = jet(size(V,1));
S1.FaceColor = 'interp';
S1.FaceAlpha = .1;
g = hgtransform;
S1Obj = patch(S1,'Parent',g);

S2.Vertices = V2;
S2.Faces = Fa2;
S2.FaceVertexCData = jet(size(V2,1));
S2.FaceColor = 'red';
S2.FaceAlpha = .1;
S2Obj = patch(S2);
pause(5)

 for i=1:length(F)-1  
    plot_data2(1, F(:,1:3) , F(:,4:6) , posReal, i, k, plotType, aa, ab)
    if Object_surround_by_point == 1    
        if i==1
            S1Obj = patch(S1,'Parent',g);
            gg = 0;
        else
            gg = gg + F(i,4:6)-F(i-1,4:6);
            g.Matrix = makehgtform('translate',gg);
            drawnow;
            set(S1Obj,'Vertices',get(S1Obj,'Vertices'));
            pause(0.1)
        end
    end
    if anchorsNumber == 5
        if k==length(posReal)-1
            break;
        elseif tI(k,1) ==  acc(i,5)
            if tI(k,2) == acc(i,6)
                k=k+1;
            end
        end
    else
        k = i+1;
    end
    % frame for video incase we need to save
    if printout_video_of_simulation == 1 
        frameVideo(i) = getframe(figure(1));
    end
    
    % Do collision detection
    if Collision_detection == 1 
        collisionFlag = colldetect(F(i,4:6),ObjCentralPoint,[MeanErrorMX MeanErrorMY MeanErrorMZ],[SDErrorMX SDErrorMY SDErrorMZ],[0.05,0.05,ObjCentralPoint(3)/3]);
        drawnow;
        if collisionFlag
            str = {'Collision!','Warning'};
            t = text(2,2,3,str,'FontSize',30);
            break;
        end
    end
 end
 
%% create video file
if printout_video_of_simulation == 1
    movie(figure(2),frameVideo,1);
    Vid = VideoWriter('simulationCD.avi');
    Vid.FrameRate = 30;
    Vid.Quality = 90;
    open(Vid);
    writeVideo(Vid,frameVideo);
    close(Vid);
end