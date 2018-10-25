function collisionFlag = colldetect(posEst,posObj,MposEst,SDposEst,SDposObj)
    A1 = posEst+MposEst+3*SDposEst;
    A2 = posEst-MposEst-3*SDposEst;
    B1 = posObj+3*SDposObj;
    B2 = posObj-3*SDposObj;
    if((A1(1) >= B2(1)) && (A2(1) <= B1(1)))
        if((A1(2) >= B2(2)) && (A2(2) <= B1(2)))
            if((A1(3) >= B2(3)) && (A2(3) <= B1(3)))
                collisionFlag = 1;
            else
                collisionFlag = 0;
            end
        else
            collisionFlag = 0;
        end
    else
        collisionFlag = 0;
    end
end