clear all
close all
clc

%fileName = 'C:\Users\moham\Desktop\tests\20-06\data.xlsx';
%[acc,pos,dis,gyro,posReal] = readData(fileName);
acc = [];
gyro = [];

%load('acc.mat','acc');
%load('gyro.mat','gyro');
%acc = xlsread('test.xlsx','acc','A:C');
%gyro = xlsread('test.xlsx','gyro','A:C');
load('testing\sensor_rotation_data.mat');
% load('C:\Users\moham\Desktop\tests\03-07\test\example.mat');
drift_gyro = [-0.000304 0.002140 -0.000540];
drift_acc = [0.154200 -0.138603 0.001707];

% acc = acc/norm(acc,inf);
% gyro = gyro/norm(gyro,inf);

gyro = gyro - repmat(drift_gyro,length(gyro),1);
acc = acc - repmat(drift_acc,length(acc),1);


axis = [1 1 1];

angle = [];
dt = 0.1;
filtered_angles = [];
earth = [0 0 9.8]';
for i=1:length(acc)
    if i==1
        pre_gyro_angle = [0 0 0];
        angle = [0 0 0];
    else
        [filtered_angles, gyro_angle, acc_angles] = angles(acc(i,:),gyro(i,:),dt,pre_gyro_angle);
        angle(i,:) = filtered_angles;
        pre_gyro_angle = gyro_angle;
        all_angles(i,:) = [gyro_angle acc_angles];
    end

    [Rx Ry Rz R] = rotationMat(angle);

    R1 = Rz*Ry*Rx;
%     R = [cos(angle(i,2))*cos(angle(i,3)) sin(angle(i,1))*sin(angle(i,2))*cos(angle(i,3))-cos(angle(i,1))*sin(angle(i,3)) cos(angle(i,1))*sin(angle(i,2))*cos(angle(i,3))+sin(angle(i,1))*sin(angle(i,3));
%          cos(angle(i,2))*sin(angle(i,3)) sin(angle(i,1))*sin(angle(i,2))*sin(angle(i,3))+cos(angle(i,1))*cos(angle(i,3)) cos(angle(i,1))*sin(angle(i,2))*sin(angle(i,3))-sin(angle(i,1))*cos(angle(i,3));
%             -sin(angle(i,2))           sin(angle(i,1))*cos(angle(i,2))                 cos(angle(i,1))*cos(angle(i,2))            ];


    acc2(i,:) = (R*acc(i,:)')' - (R\earth)';
    acc2(i,:) = acc2(i,:)/norm(acc(i,:));
    acc3(i) = norm(acc(i,:));
    acc4(i) = norm(acc2(i,:));
end

plot_angles(all_angles,angle,dt);





