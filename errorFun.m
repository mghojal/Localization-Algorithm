% calculating mean error in 3D space
function [errorMLat errorMLatKF] = errorFun(posReal,F,i,k)
    errorMLat = sqrt((posReal(k,1)-F(i,1))^2+(posReal(k,2)-F(i,2))^2+(posReal(k,3)-F(i,3))^2);
    errorMLatKF = sqrt((posReal(k,1)-F(i,4))^2+(posReal(k,2)-F(i,5))^2+(posReal(k,3)-F(i,6))^2);
end