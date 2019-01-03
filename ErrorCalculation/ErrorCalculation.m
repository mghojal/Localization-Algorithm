function ErrorM = ErrorCalculation(posEst, posReal, AnNo)    
    j=1;
    for i=1:length(posEst)
        if i==1
            ErrorM(i,1) = abs(posEst(i,1) - posReal(j,1));
            ErrorM(i,2) = abs(posEst(i,2) - posReal(j,2));
            ErrorM(i,3) = abs(posEst(i,3) - posReal(j,3));
            j=j+1;
        elseif AnNo == 5 && rem(i-1,5)==0 || rem(i-1,5)==2
            ErrorM(i,1) = abs(posEst(i,1) - posReal(j,1));
            ErrorM(i,2) = abs(posEst(i,2) - posReal(j,2));
            ErrorM(i,3) = abs(posEst(i,3) - posReal(j,3));
            j=j+1;
        elseif  AnNo == 5
            ErrorM(i,1) = abs(posEst(i,1) - posReal(j-1,1));
            ErrorM(i,2) = abs(posEst(i,2) - posReal(j-1,2));
            ErrorM(i,3) = abs(posEst(i,3) - posReal(j-1,3));
        else
            ErrorM(i,1) = abs(posEst(i,1) - posReal(j,1));
            ErrorM(i,2) = abs(posEst(i,2) - posReal(j,2));
            ErrorM(i,3) = abs(posEst(i,3) - posReal(j,3));
            j=j+1;
        end
    end
    
end
