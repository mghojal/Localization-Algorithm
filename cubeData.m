function [V Fa] = cubeData(F,MeanErrorMX,MeanErrorMY,MeanErrorMZ,SDErrorMX,SDErrorMY,SDErrorMZ)
Xm = F(1,4)-MeanErrorMX-3*SDErrorMX;
Xp = F(1,4)+MeanErrorMX+3*SDErrorMX;
Ym = F(1,5)-MeanErrorMY-3*SDErrorMY;
Yp = F(1,5)+MeanErrorMY+3*SDErrorMY;
Zm = F(1,6)-MeanErrorMZ-3*SDErrorMZ;
Zp = F(1,6)+MeanErrorMZ+3*SDErrorMZ;
p1 = [Xm Ym Zm];
p2 = [Xm Yp Zm];
p3 = [Xp Yp Zm];
p4 = [Xp Ym Zm];
p5 = [Xm Ym Zp];
p6 = [Xm Yp Zp];
p7 = [Xp Yp Zp];
p8 = [Xp Ym Zp];

V = [p1;p2;p3;p4;...
     p2;p3;p7;p6;...
     p3;p4;p8;p7;...
     p4;p1;p5;p8;...
     p1;p2;p6;p5;...
     p5;p6;p7;p8];
Fa = [1 2 3 4;
      5 6 7 8;
      9 10 11 12;
      13 14 15 16;
      17 18 19 20;
      21 22 23 24];