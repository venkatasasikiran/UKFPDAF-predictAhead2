function FilterEstimateUKF = UKF_Estimate2(FilterEstimateUKF,FilterMat,loopNum)
    
    PkPlus = FilterEstimateUKF.PkMinus;
    XkPlus = FilterEstimateUKF.XkMinus;
    n = length(XkPlus);
    alpha = 10^-3;
    kappa = 0;
    beta = 2;
    lamda = alpha^2*(n+kappa) - n;    

    % Dimension the weight and sample arrays
    %

    noPoints    = 2*n+1;             % number of samples
    wPts        = zeros(1,noPoints); % sample weightings
    XPlusPts    = zeros(n,noPoints); % samples

    % Calculate the weight vector
    %
    wPts(1)             = lamda/(n+ lamda);
    for i=1:2*n
        wPts(i+1)       = 1/(2*(n+lamda));
    end
    % Calculate sigma points
    %
    [U,Eigen,~]               = svd(PkPlus);
    Cxx                       = U*sqrt(Eigen);
    c                         = sqrt(n+lamda);
    XPlusPts(:,1)             = XkPlus;
    XPlusPts(:,2:n+1)         = XkPlus*ones(1,n)+c*Cxx; 
    XPlusPts(:,n+2:2*n+1)     = XkPlus*ones(1,n)-c*Cxx;
    
    % Transform the sample points
    %
    XMinusPts = Trans_F(XPlusPts,FilterMat.F);
    m1 = size(XMinusPts,1);
    % Calculate the X_k_minus and P_k_minus
    %
    XkMinus = XMinusPts*wPts';
    PkMinus = zeros(m1);
    for i = 0: 2*n
        if i == 0
            PkMinus = PkMinus...
                      + (wPts(i+1)+(1-alpha^2+beta))...
                      *(XMinusPts(:,i+1)-XkMinus)*(XMinusPts(:,i+1) - XkMinus)';
        else
            PkMinus = PkMinus + wPts(i+1)...
                      *(XMinusPts(:,i+1)-XkMinus)*(XMinusPts(:,i+1) - XkMinus)';
        end
    end
    PkMinus = PkMinus + FilterMat.Q;
    FilterEstimateUKF.XkMinus = XkMinus;
    FilterEstimateUKF.PkMinus = PkMinus;
    FilterEstimateUKF.XkMinusVec(:,loopNum) = XkMinus;
    FilterEstimateUKF.PkVec(loopNum).Minus = PkMinus;
    FilterEstimateUKF.XMinusPts = XMinusPts;
end

function output = Trans_F(input,F)
    
    [Row,Col] = size(input);
    output = zeros(Row, Col);    
    for i = 1:Col
        output(:,i) = F*input(:,i);        
    end
end


