function [aa,ab] = definedFigure2(anchorLoc,figureNumber,MultilaterationAlgorithmName)

figure('Name',MultilaterationAlgorithmName);
aa=subplot(2,1,1);
plot3(aa,anchorLoc(:,1),anchorLoc(:,2),anchorLoc(:,3),'ko','MarkerSize',8,'lineWidth',2,'MarkerFaceColor','k');
grid on
hold on
curve = animatedline('Color','r','Marker','o');
axis([-1 4 -1 4.9 0 3]);
ax = gca;
ax.XTick = (0:0.3:3);
ax.YTick = (0:0.3:3.9);
ax.ZTick = (0:1:3);
xlabel('x axis (m)','FontSize',18)
ylabel('y axis (m)','FontSize',18)
zlabel('z axis (m)','FontSize',18)
set(gca,'fontsize',18);
title(aa,'MLat w/o filtering')
hold on

ab=subplot(2,1,2);
plot3(ab,anchorLoc(:,1),anchorLoc(:,2),anchorLoc(:,3),'ko','MarkerSize',8,'lineWidth',2,'MarkerFaceColor','k');
grid on
hold on
curve1 = animatedline('Color','r','Marker','o');
axis([-1 4 -1 4.9 0 3]);
ax1 = gca;
ax1.XTick = (0:0.3:3);
ax1.YTick = (0:0.3:3.9);
ax1.ZTick = (0:1:3);
xlabel('x axis (m)','FontSize',18)
ylabel('y axis (m)','FontSize',18)
zlabel('z axis (m)','FontSize',18)
set(gca,'fontsize',18);
title(ab,'MLat with filtering')
hold on