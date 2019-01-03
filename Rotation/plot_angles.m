function plot_angles(all_angles,angle,dt)
y = 0:dt:(length(all_angles)*dt-dt);
figure('Name','angles')
ax = subplot(3,1,1);
plot (ax,y,all_angles(:,4),'r');
hold on
grid on
plot (ax,y,all_angles(:,1),'b');
plot (ax,y,angle(:,1),'k');
set(gca,'YTick',-pi:pi/2:pi,'YTicklabels',({'-\pi','-\pi/2','0','\pi/2','\pi'}));
xlabel('time (sec)');
%ylabel('angle');
legend({'acc_x^°','gyro_x^°','filtered_x^°'},'location','bestoutside');
set(gca,'FontSize',18);

ay = subplot(3,1,2);
plot (ay,y,all_angles(:,5),'r');
hold on
grid on
plot (ay,y,all_angles(:,2),'b');
plot (ay,y,angle(:,2),'k');
set(gca,'YTick',-pi:pi/2:pi,'YTicklabels',({'-\pi','-\pi/2','0','\pi/2','\pi'}));
xlabel('time (sec)');
%ylabel('angle');
legend({'acc_y^°','gyro_y^°','filtered_y^°'},'location','bestoutside');
set(gca,'FontSize',18);

az = subplot(3,1,3);
plot (az,y,all_angles(:,6),'r');
hold on
grid on
plot (az,y,all_angles(:,3),'b');
plot (az,y,angle(:,3),'k');
set(gca,'YTick',-pi:pi/2:pi,'YTicklabels',({'-\pi','-\pi/2','0','\pi/2','\pi'}));
xlabel('time (sec)');
%ylabel('angle');
legend({'acc_z^°','gyro_z^°','filtered_z^°'},'location','bestoutside');
set(gca,'FontSize',18);