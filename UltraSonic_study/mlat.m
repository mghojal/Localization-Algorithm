classdef mlat
  % to run Mlat you need to use main funciton 
  % > do_main('algorithm', anchor Location, distances, varargin)
  % > 'algorithms' are four: gradient descent, recursive trilateration,
  % simple least square, normal trilateration
  % > varargin: using bounds for gradient descent and weighting matrix (error
  % variance matrix) for recursive and simple trilateration algorithms 
  % > MLAT is class contains 4 different algorithms for Multilateration and
  % those are:
  %     1- Gradient Decent using two funtions:
  %             a- gd()
  %             b- do_gdescent()
  %     2- Recursive Trilateration using three funtions:
  %             a- RecTrilateration()
  %             b- distancen()
  %             c- lsrec()
  %             d- do_RecTri()
  %     3- Simple Least Square Error algorithm using one funtion do_LLS()
  %     4- normal Trilateration algorithm using one funtion do_Tri()
  
  methods (Access = public ,Static)
      %% main function
      function position = do_main(Algorithm, anchor_Location, ranges_in, varargin)
          bounds_in = findBounds(anchor_Location);
          W = eye(length(ranges_in)); 
          if length(varargin) ~= 0
              for i = 1:2:length(varargin{1})
                switch (varargin{1}{i})
                  case 'bounds'
                    bounds_in = varargin{1}{i+1};
                  case 'weight'
                    W = varargin{1}{i+1};
                  case 'Null'
                    continue;
                end
              end
          end
          
          switch Algorithm
              case 'Gradient_Descent'
                position = mlat.do_gdescent(anchor_Location, ranges_in, 'bounds',bounds_in);
              case 'Recursive_Trilateration'
                  position = mlat.do_RecTri(anchor_Location', ranges_in, W);
              case 'Least_square'
                  position = mlat.do_LLS(anchor_Location', ranges_in);
              case 'Trilateration'
                  position = mlat.do_Tri(anchor_Location',ranges_in, W);
          end
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %% Gradinet Descent Algorithm
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Gradient Descent Algorithm %
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function result_table = gd(anchors_in, ranges_in, varargin)
      % varargin: 'bounds_in', 'n_trial', 'alpha', 'time_threshold'
      bounds_in = [];
      n_trial = 1000;
      alpha = 0.0001;
      time_threshold = 1/ n_trial;
      for i = 1:2:length(varargin{1})
        switch (varargin{1}{i})
          case 'bounds'
            bounds_in = varargin{1}{i+1};
          case 'trial'
            n_trial = varargin{1}{i+1};
          case 'alpha'
            alpha = varargin{1}{i+1};
          case 'time'
            time_threshold = varargin{1}{i+1};
        end
      end

      [n, dim] = size(anchors_in);
      bounds_temp = [anchors_in; bounds_in];
      bounds(1, :) = min(bounds_temp);
      bounds(2, :) = max(bounds_temp);
      
      ranges = nan(1, n);
      result_table = nan(n_trial, dim + 1);
      
      for i = 1:n_trial
        estimator0 = nan(1, dim);
%         estimator0(1:2) = [1.5 2.2];
%         estimator0(3) = (bounds(2, 3) - bounds(1, 3)) * rand + bounds(1, 3);
        for j = 1:dim
          estimator0(j) = (bounds(2, j) - bounds(1, j)) * rand + bounds(1, j);
        end
        estimator = estimator0;
        
        t0 = tic;
        while true
          for j = 1:n
            ranges(j) = norm(anchors_in(j, :) - estimator);
          end
          err = norm(ranges_in - ranges);
          
          delta = zeros(1, dim);
          for j = 1:n
            delta = delta + (ranges_in(j) - ranges(j)) / ranges(j) * (estimator - anchors_in(j, :));
          end
          delta = 2 * alpha * delta;
          
          estimator_next = estimator - delta;
          for j = 1:n
            ranges(j) = norm(anchors_in(j, :) - estimator_next);
          end
          err_next = norm(ranges_in - ranges);
          if err_next < err
            estimator = estimator_next;
          else
            result_table(i, 1:dim) = estimator;
            result_table(i, end) = err;
            break;
          end
          if toc - t0 > time_threshold
            break;
          end
        end
      end
    end

    function estimator = do_gdescent(anchors_in, ranges_in, varargin)
      result_table = mlat.gd(anchors_in, ranges_in, varargin);
      [~, I] = min(result_table(:, end));
      estimator = round(result_table(I, 1:end-1),2);
    end
    
    %% Simple Least Square Error algorithm    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Least Square Error Simple Trilateration %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function Nmat = do_LLS(anchors_in,ranges_in)
        A = []; b = [];
        [m n] = size(anchors_in);
        x = anchors_in(1,:); y = anchors_in(2,:); z = anchors_in(3,:);
        for i=2:n
            A(i-1,:) = [x(i)-x(1) , y(i)-y(1) , z(i)-z(1)];
            b(i-1,:) = ranges_in(1)^2 - ranges_in(i)^2 + x(i)^2 + y(i)^2 + z(i)^2 - x(1)^2 - y(1)^2 - z(1)^2;
        end
        %A = 2*A;
        Nmat = pinv(A'*A)*A'*b/2;
    end
    
    %% Recursive Trilateration
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Recursive Trilateration %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    function Nmat = RecTrilateration(P,S,W)
        % paper "An algebraic solution to the multilateration problem"
        % Author: Norrdine, Abdelmoumen  (norrdine@hotmail.de)
        % https://www.researchgate.net/publication/275027725_An_Algebraic_Solution_to_the_Multilateration_Problem
        % usage: [N1 N2] = RecTrilateration(P,S,W) 
        % P = [P1 P2 P3 P4 ..] Reference points matrix
        % S = [s1 s2 s3 s4 ..] distance matrix.
        % W : Weights Matrix (Statistics).
        % N : calculated solution
        [mp,np] = size(P);
        ns = length(S);
        if (ns~=np)
            error('number of reference points and ranges and not the same');
        end
        A=[]; b=[];
        for i1=1:np
            x = P(1,i1); y = P(2,i1); z = P(3,i1);
            s = S(i1);
            A = [A ; 1 -2*x  -2*y  -2*z]; 
            b= [b ; s^2-x^2-y^2-z^2 ];
        end
        if (np==3)
            warning off;
            Xp= pinv(A)*b;  % Gaussian elimination
            % or Xp=pinv(A)*b; P
            % the matrix  inv(A'*A)*A' or inv(A'*C*A)*A'*C or pinv(A)
            % depend only on the reference points
            % it could be computed only once
            xp = Xp(2:4,:);
            Z = null(A);
            z = Z(2:4,:);
            if rank (A)==3
                %Polynom coeff.
                a2 = z(1)^2 + z(2)^2 + z(3)^2 ;
                a1 = 2*(z(1)*xp(1) + z(2)*xp(2) + z(3)*xp(3))-Z(1);
                a0 = xp(1)^2 +  xp(2)^2+  xp(3)^2-Xp(1);
                p = [a2 a1 a0];
                t = roots(p);
                %solution
                N1 = Xp + t(1)*Z;
                N2 = Xp + t(2)*Z;
                Nmat(:,1) = N1;
                Nmat(:,2) = N2;
            end
        end
        A0 = A(1:3,:);
        if  (np>3)
            P10=P(:,1:end-1);S10=S(:,1:end-1);W0=W(1:end-1,1:end-1);
            N0mat = mlat.RecTrilateration(P10,S10,W0);
            N01 = N0mat(:,1);
            N02 = N0mat(:,2);
            %select N0
            C = W'*W;
            Xpdw = pinv(A'*C*A)*A'*C*b;
            % the matrix  inv(A'*A)*A' or inv(A'*C*A)*A'*C or pinv(A)
            % depend only on the reference points
            % it could be computed only once
            NormErrorXpdw = Xpdw(1)-norm(Xpdw(2:4))^2;
            if (norm(Xpdw(2:4)-N01(2:4))<norm(Xpdw(2:4)-N02(2:4)))
                N0 = N01;
            else
                N0 = N02;
            end
            Nmat(:,1)= N01;
            Nmat(:,2)= N02;
            W0 = W(1:3,1:3);
            C0 = W0*W0';
            P_0 = inv(A0'*C0*A0);
            %Start solution
            invP_i_1 = inv(P_0); 
            xi_1 = N0;
            % recursive Least square (Introduction to applied Math Strang pp 147)
            x0 = N0;
            [x,P] = mlat.lsrec(x0,1);
            for i=1:np-3
                An = A(i+3,:);
                Wn = W(i+3,i+3);
                yn = b(i+3);
                [xn,Pn] = mlat.lsrec(yn,An,1,x,P);
                x=xn; P=Pn;
                Nmat(:,i+2) = xn;
            end
              Nmat(:,i+3)= Xpdw;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [Sn , F] = distancen(Nn,P,S)
        % Calculate distance between Nn and referece points Pi
        % calculates distances between Nn and the reference points Pi
        % P= [P1 P2 ...]  ; measured distances : S=[s1 s2 ...]
        % Sn = []: calculated distances 
        % F : Error norm
        global P S
        for i1=1:length(S)
            Sn(i1)=norm(P(:,i1)-Nn);
        end
        F = norm(S-Sn);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [xn,Pn]=lsrec(varargin)
        %LSREC Recursive Least Squares.
        % [x,P] = LSREC(x0,W) initializes a recursive solution by returning the
        % initial solution x = x0 having a scalar weight 0 < W <= 1. If x0 is a
        % very good first estimate, use W near 1. If x0 is a poor first estimate
        % use W near 0.  If W is not given, W = 1e-12 is used. P is a matrix of size
        % length(x0)-by-length(x0) that is required for future recursive calls.
        %
        % [xn,Pn] = LSREC(yn,An,Wn,x,P) computes the recursive least squares
        % solution xn, given new equations yn = An*x, where size(An,1) >= 1 and
        % size(An,2) = length(x). Wn is the weight associated with the new data,
        % which is typically equal to 1. If Wn is a scalar it applies to all new
        % equations; if it is a vector the i-th element of Wn applies to the i-th
        % equation. x and P are the output from the most recent recursive function
        % call. xn and Pn are the updated solution vector and P matrix for future
        % recursive calls.
        %
        % This function is useful when one wants to update a least squares solution
        % repeatedly as new data becomes available, such as after each pass through
        % some iterative process.
        %
        % Reference: "Modern Control Theory," 3rd ed., William L. Brogan
        % Prentice Hall, 1991.
        %
        % See also MLDIVIDE, LSCOV, LSQNONNEG.

        % D.C. Hanselman, University of Maine, Orono, ME 04469
        % MasteringMatlab@yahoo.com
        % Mastering MATLAB 7
        % 2006-11-8

        if nargin==1                                % initialize recursive solution
           xn=varargin{1}(:);
           Pn=diag(1e12+zeros(size(xn)));

        elseif nargin==2                            % initialize recursive solution
           xn=varargin{1}(:);
           if numel(varargin{2})~=1
              error('LSREC:scalar','Scalar Weight Required.')
           else
              W=varargin{2};
              if W<=eps || W>1
                 error('LSREC:OutofBound','Weight Must be Between 0 and 1.')
              end
              Pn=diag((1/W)+zeros(size(xn)));
           end

        elseif nargin==5                                           % recursive call

           yn=varargin{1}(:); % make sure yn is a column vector
           An=varargin{2};
           Wn=varargin{3}(:);
           x=varargin{4}(:);
           P=varargin{5};
           if length(yn)~=size(An,1)
              error('LSREC:nonconform',...
                    'yn Must Have as Many Rows as An.')
           end
           if size(An,2)~=length(x)
              error('LSREC:nonconform',...
                    'An Must Have as Many Columns as x has elements.')
           end
           if size(P,1)~=size(P,2) || size(P,1)~=length(x)
              error('LSREC:nonform',...
                    'P Must be a Square Matrix of Dimension Equal to length(x).')
           end
           if length(Wn)~=1 && length(Wn)~=length(yn)
              error('LSREC:conform',...
                    'Wn Must be a Scalar or Have the Same Number of Elements as yn.')
           end
           if any(Wn<=eps) || any(Wn>1)
              error('LSREC:OutofBound','Weights Must be Between 0 and 1.')
           end
           if numel(Wn)==1   % expand scalar weight if needed
              Wn=repmat(Wn,size(yn));
           end

           K=P*An'/(An*P*An'+diag(1./Wn));
           xn=x+K*(yn-An*x);
           if nargout>1  % compute new P
              Pn=P-K*An*P;
           end
        else
           error('LSREC:rhs','Recursive Calls Require 5 Input Arguments.')
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function NFinal = do_RecTri(P,S,W)
        Nmat = mlat.RecTrilateration(P,S,W);
        N1 = Nmat(:,1:end);
        [n m]=size(N1);
        for i1=1:m
            Nn = N1(:,i1);
            Nn = Nn(2:4);
            [Sn(i1,:) F(i1)] = mlat.distancen(Nn,P,S);
            diff(i1) = N1(1,i1) - norm(Nn).^2;
        end
        [Fmin Imin] = min(F);
        [NrecDel IDel] = min(N1(1,1:3));
        NrecDel = N1(2:4,IDel);
        Nrec = N1(2:4,Imin);
        if imag(Nrec(2))==0
            NFinal = Nrec';
        else
            for jj=1:size(N1,2)
                if imag(N1(2,jj))==0
                    NFinal = N1(2:4,jj)';
                end
            end
        end
    end  
    
    %% Normal Trilateration
    %%%%%%%%%%%%%%%%%
    % Trilateration %
    %%%%%%%%%%%%%%%%%
    function Nmat = do_Tri(P,S,W)
        % paper "An algebraic solution to the multilateration problem"
        % Author: Norrdine, Abdelmoumen  (norrdine@hotmail.de)
        % https://www.researchgate.net/publication/275027725_An_Algebraic_Solution_to_the_Multilateration_Problem
        % usage: [N1 N2] = Trilateration(P,S,W) 
        % P = [P1 P2 P3 P4 ..] Reference points matrix
        % S = [s1 s2 s3 s4 ..] distance matrix.
        % W : Weights Matrix (Statistics).
        % N : calculated solution
        % THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY!! 
        [mp,np] = size(P);
        ns = length(S);
        if (ns~=np)
            error('Number of reference points and distances are different');
        end
        A=[]; b=[];
        for i1=1:np
            x = P(1,i1); y = P(2,i1); z = P(3,i1);
            s = S(i1);
            A = [A ; 1 -2*x  -2*y  -2*z]; 
            b= [b ; s^2-x^2-y^2-z^2 ];
        end

        if (np==3)
            warning off;
            Xp= A\b;  % Gaussian elimination
           % or Xp=pinv(A)*b; 
           % the matrix  inv(A'*A)*A' or inv(A'*C*A)*A'*C or pinv(A)
           % depend only on the reference points
           % it could be computed only once
            xp = Xp(2:4,:);
            Z = null(A,'r');
            z = Z(2:4,:);
            if rank (A)==3
                %Polynom coeff.
                a2 = z(1)^2 + z(2)^2 + z(3)^2 ;
                a1 = 2*(z(1)*xp(1) + z(2)*xp(2) + z(3)*xp(3))-Z(1);
                a0 = xp(1)^2 +  xp(2)^2+  xp(3)^2-Xp(1);
                p = [a2 a1 a0];
                t = roots(p);

                %Solutions
                N1 = Xp + t(1)*Z;
                N2 = Xp + t(2)*Z;
                if N1(1)<N2(1)
                    Nmat = N1;
                else
                    Nmat = N2;
                end
            end
        end
        if  (np>3)
        %Particular solution

            if W~=diag(ones(1,length(W)))
                C = W'*W;
                Xpdw =inv(A'*C*A)*A'*inv(C)*b; % Solution with Weights Matrix
            else
                Xpdw=inv(A'*A)*A'*b; % Solution without Weights Matrix
            end

            % the matrix  inv(A'*A)*A' or inv(A'*C*A)*A'*C or pinv(A)
            % depend only on the reference points
            % it could be computed only once
            N1 = Xpdw;
            N2 = N1;
            Nmat = N1(2:4)';
        end
    end
        
  end
  
end