function [C] = errorMeasuresMatrix(alpha)

    % Measurement matrix
    % The expected measurement given the predicted state
    C = [1 0 0 0 0 0 0;
         0 1 0 0 0 0 0;
         0 0 alpha 1-alpha 0 0 0;
         0 0 0 0 0 0 0
         0 0 0 0 0 0 0
         0 0 0 0 0 0 0
         0 0 0 0 1 0 0;
         0 0 0 0 0 1 0;
         0 0 0 0 0 0 1];
     
end
