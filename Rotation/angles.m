function [filtered_angles, gyro_angle, acc_angles] = angles(acc,gyro,dt,pre_gyro_angle)
alpha = 0.95;
%acc angles
phi =   atan2(acc(2) , acc(3)); 
theta = atan2(acc(3) , acc(1));
psi =   atan2(acc(2) , acc(1));
% phi =   atan(acc(2)/(sqrt(acc(1)^2 + acc(3)^2)));
% theta = -atan(acc(1)/(sqrt(acc(2)^2 + acc(3)^2))); 
% psi =   atan(sqrt(acc(2)^2 + acc(1)^2)/acc(3));
acc_angles = [phi theta psi];
gyro_angle = (gyro*dt)+pre_gyro_angle;
for (k=1:3)
    if gyro_angle(k) > pi
        gyro_angle(k) = gyro_angle(k) - 2*pi;
    elseif gyro_angle(k) <= -pi
        gyro_angle(k) = gyro_angle(k) + 2*pi;
    else
        gyro_angle(k) = gyro_angle(k);
    end
end

for (k=1:3)
    if acc_angles(k) > pi
        acc_angles(k) = acc_angles(k) - 2*pi;
    elseif acc_angles(k) <= -pi
        acc_angles(k) = acc_angles(k) + 2*pi;
    else
        acc_angles(k) = acc_angles(k);
    end
end
filtered_angles = [(alpha*gyro_angle(1:2))+((1-alpha)*acc_angles(1:2)), gyro_angle(3)];

