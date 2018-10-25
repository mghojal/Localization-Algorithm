function ErrorM = ErrorCalculation(posEst, posReal, AnNo)    
    for i=1:length(posEst)
            ErrorM(i,1) = posEst(i,1) - posReal(i,1);
            ErrorM(i,2) = posEst(i,2) - posReal(i,2);
            ErrorM(i,3) = posEst(i,3) - posReal(i,3);
    end    
end