close all
clear all
clc

%% choosen path
% options are the following:
% this movement focus on z axis

load('input.mat');

%% calculating positions using mlat (3 algorithms)
% Algorithms are four:  1- Gradient Descent
%                       2- Recursive Least Square
%                       3- Least Square
% we need to choose algorithm:

%     Algorithm = 'Gradient_Descent';
    Algorithm = 'Recursive_Trilateration';
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
iteration = length(acc);
anchorsNumber = size(dis,2);
dt = 0.25;
%% starting with process
[A,Q,R,H] = preDefinedMatrix2(dt,Algorithm); 
tic
for i=1:iteration

    positionEstimate(i,:) = mlat.do_main(Algorithm, anchorLoc, dis(i,:));
    %% Initial values for i==1 for kalman filtering and angles calculation
    if (i==1)  
        x_t =  [positionEstimate(i,:) 0 0 0 0 0 0]';
     
        pHat = [dt^4/4 0 0 dt^3/2 0 0 dt^2/2 0 0;
                0 dt^4/4 0 0 dt^3/2 0 0 dt^2/2 0;
                0 0 dt^4/4 0 0 dt^3/2 0 0 dt^2/2;
                dt^3/2 0 0 dt^2 0 0 dt 0 0;
                0 dt^3/2 0 0 dt^2 0 0 dt 0;
                0 0 dt^3/2 0 0 dt^2 0 0 dt;
                dt^2/2 0 0 dt 0 0 1 0 0;
                0 dt^2/2 0 0 dt 0 0 1 0;
                0 0 dt^2/2 0 0 dt 0 0 1];    
    else
        x_t = X;
        pHat = P;
    end

    %% reading values position after applying Multilateration algorithms and acceleration 
    if without_IMU == 1
        Y_t = [positionEstimate(i,:) USz(i) 0 0 0]';    
    else
        Y_t = [positionEstimate(i,:) USz(i) acc(i,1:3)]';
    end

    %% Kalman Filtering
    alpha = 0.11;
    if USz(i) > 1.0
        alpha = 1-alpha;
    end
    [C] = errorMeasuresMatrix(alpha);
    [X, P] = ekf( A, C, Q, H, R , x_t, Y_t, pHat);

    %% Results saving
    F(i,1:3) = positionEstimate(i,:);        
    F(i,4:6) = X(1:3);
    aAndv(i,:) = X(4:9);
    GG(i,:) = abs(F(i,1:3)-posReal(i,1:3));
    [eM eMKF] = errorFun(posReal,F,i,i);
    ErrorMatrix(i,1:2) = [eM eMKF];
end
toc
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
%pause(5)

 for i=1:length(F)-1  
    plot_data2(1, F(:,1:3) , F(:,4:6) , posReal, i, i, plotType, aa, ab) 
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

    % frame for video incase we need to save
    if printout_video_of_simulation == 1
        frameVideo(i) = getframe(figure(1));
    end
    
    % Do collision detection
    if Collision_detection == 1
        collisionFlag = colldetect(F(i,4:6),ObjCentralPoint,[SDErrorMX SDErrorMY SDErrorMZ],[0.05,0.05,ObjCentralPoint(3)/3]);

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
    Vid = VideoWriter('simulation.avi');
    Vid.FrameRate = 10;
    open(Vid);
    writeVideo(Vid,frameVideo);
    close(Vid);
end