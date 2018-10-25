function plot_data2(figureNumber, PMlat, PKF, PR, i, k, PT, aa, ab)

%% definingtions
    % FigureNumber = FN
    % PointsMultilateration = PMlat
    % PointsKF = PKF
    % Real Position = PR
    % Iteration I = i
    % Iteration k = k
    % Plot Type = PT
    % Subfigure 1 multilateration result = aa
    % Subfigure 2 multilateration filter result = ab
 %% Ploting funtion parts
    
LineDrawing1 = [PMlat(i,1) PMlat(i,2) PMlat(i,3);...
                PMlat(i+1,1) PMlat(i+1,2) PMlat(i+1,3)];
LineDrawing2 = [PKF(i,1) PKF(i,2) PKF(i,3);...
                PKF(i+1,1) PKF(i+1,2) PKF(i+1,3)];
LineDrawing5 = [PR(k+1,:); PR(k,:)];
            
if PT == 1
    az = 0; el = 90;
elseif PT == 2
    az = 90; el = 0;
elseif PT == 0
    az = -37.5; el = 30;
elseif PT == 3
    az = 180; el = 0;
end

view(aa,az,el);
view(ab,az,el);

line(aa,LineDrawing1(:,1),LineDrawing1(:,2),LineDrawing1(:,3),'Color','r');...
     line(aa,LineDrawing5(:,1),LineDrawing5(:,2),LineDrawing5(:,3),'Color','m');
line(ab,LineDrawing2(:,1),LineDrawing2(:,2),LineDrawing2(:,3),'Color','b');...
     line(ab,LineDrawing5(:,1),LineDrawing5(:,2),LineDrawing5(:,3),'Color','m');