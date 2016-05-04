function [Model]    = LRM_Panel_GibbsSampling(Model)
global tmp Model
counter =   1;

for s = 1:Model.S
    
    if s <= Model.S1 & rem(s,Model.PrintIteration) == 0
        disp(['Burn in Gibbs sampling iteration  ' num2str(s) ' of ' num2str(Model.S1)])
    end
    if s > Model.S1 & rem(s,Model.PrintIteration) == 0
        disp(['Retained Gibbs sampling iteration ' num2str(s-Model.S1) ' of ' num2str(Model.S2)])
    end
    % DRAW THE VARIANCE OF THETA
    % Draw the V_Theta_Inv form a Wishart Distribution
    % Eq 7.27 Line 5
    Model.V_Theta     = zeros(Model.K,Model.K);                     %Initialize   
    for n = 1:Model.N
        ErrTheta    = Model.ThetaDraw(:,n)-Model.MTheta;%A draw of Theta - Mean of Theta (Start with OLS)
        Model.V_Theta = Model.V_Theta + ErrTheta*ErrTheta';
    end
    
    tmp                 =   inv(Model.V_Theta + Model.D_VTheta_Pr * Model.R);    %Page 156, Eq 7.23
    [ tmp,RequiredForce ] = ForcePSD(tmp);
    Model.V_Theta_Inv   =   wish_rnd(tmp, Model.D_VTheta_Ps);
    Model.V_Theta_Drw   =   inv(Model.V_Theta_Inv);
    
    % DRAW THE VARIANCE OF ERROR (OR PRECISION OF ERROR = h)  
    % Draw h conditional on other parameters
    Model.S2_Error_Ps = 0;
    for n = 1:Model.N
        T           = (n-1)*Model.T+1:n*Model.T;    %All observaton of the individual (data is stacked)
        xuse        = Model.X(T,:);                 %Relevant X matrix for the individual
        yuse        = Model.Y(T,1);                 %Relevant Y matrix for the individual
        thetuse     = Model.ThetaDraw(:,n);         %Relevant Theta for the individual
        e           = (yuse-xuse*thetuse);          %Errors for the individual
        Model.S2_Error_Ps = Model.S2_Error_Ps + e'*e;%Add sum of squared errors due to individual
    end
    %Page 157, Eq 7.28 Line 3
    Model.S2_Error_Ps = (Model.S2_Error_Ps + Model.DS2_Error_Pr)/Model.D_Error_Ps;                     
    Model.H_Error_Drw = gamm_rnd(1,1, .5*Model.D_Error_Ps, .5*Model.D_Error_Ps*Model.S2_Error_Ps);        
 
    % DRAW THE MEAN OF THETA
    % Draw theta0 (mean in hierarchical prior) conditional on other params
    Model.V_MTheta_Ps = inv(Model.N * Model.V_Theta_Inv + Model.V_MTheta_Inv_Pr);   %Eq 7.27 Line 2
    tmp1        = Model.N * Model.V_Theta_Inv * mean(Model.ThetaDraw')';            %Eq 7.27 Line 3
    tmp2        = Model.V_MTheta_Inv_Pr * Model.M_MTheta_Pr;                        %Eq 7.27 Line 3
    Model.M_MTheta_Ps = tmp1 + tmp2;
    Model.MTheta   = Model.V_MTheta_Ps*Model.M_MTheta_Ps + norm_rnd(Model.V_MTheta_Ps);      %Eq 7.26

      
    % DRAW THETA
    %Now draw theta-i s conditional on other parameters
    for n = 1:Model.N
        T           = (n-1)*Model.T+1:n*Model.T;%All observaton of the individual (data is stacked)
        xuse        = Model.X(T,:);                   
        yuse        = Model.Y(T,1);
        Model.V_Theta_Ps  = inv(Model.H_Error_Drw * xuse'*xuse + Model.V_Theta_Inv);                     %Eq 7.25 Line 2
        Model.M_Theta_Ps  = ...
            Model.V_Theta_Ps*(Model.H_Error_Drw*xuse'*yuse + Model.V_Theta_Inv*Model.MTheta);     %Eq 7.25 Line 3
        Model.ThetaDraw(:,n) = Model.M_Theta_Ps + norm_rnd(Model.V_Theta_Ps);
    end
    
    if s > Model.S1      %After Burnout Period, Store all Draws
        Model.H_Error(counter)      =   Model.H_Error_Drw;
        Model.M_Theta(counter,:)    =   Model.MTheta';
        Model.SV_Theta(counter,:)   =   reshape(Model.V_Theta_Drw, 1, []);
        Model.ThetaMean             =   Model.ThetaMean + Model.ThetaDraw;
        Model.ThetaStd              =   Model.ThetaStd + Model.ThetaDraw.^2;
    
        counter = counter + 1;
%         if imlike==1
%             %For Chib method, this calculates posterior chunk relating to theta0
%             logpost = -.5*kx*log(2*pi) -.5*kx*log(det(capd0))...
%                 -.5*(th0chib-th0mean)'*inv(capd0)*(th0chib-th0mean);
%             postth0=postth0+logpost;
%         end
    end
end

Model.ThetaMean = Model.ThetaMean ./ Model.S2;
Model.ThetaStd  = Model.ThetaStd  ./ Model.S2;
Model.ThetaStd  = sqrt(Model.ThetaStd - Model.ThetaMean.^2);
if isfield(Model,'Ystar')
    for n = 1:Model.N
        Model.Ystarhat(n)       = Model.Xstar(n,:)*Model.ThetaMean(:,n);
        Model.YstarhatOLS(n)    = Model.Xstar(n,:)*Model.Thetas_OLS(:,n);
    end
    C = corr([Model.Ystar,Model.Ystarhat,Model.YstarhatOLS]);
    Model.OutOfSampleCorrelation_Bayesian   = C(1,2);
    Model.OutOfSampleCorrelation_OLS        = C(1,3);
end
errors_Bayes    =   zeros(Model.N*Model.T,1);

for n = 1:Model.N
    T           = (n-1)*Model.T+1:n*Model.T;    %All observaton of the individual (data is stacked)
    xuse        = Model.X(T,:);                 %Relevant X matrix for the individual
    yuse        = Model.Y(T,1);                 %Relevant Y matrix for the individual
    thetuse     = Model.ThetaMean(:,n);         %Relevant Theta for the individual
    e           = (yuse-xuse*thetuse);          %Errors for the individual in the panel
    errors_Bayes(T,1)   = e;
end

Model.errors_Bayes  =   errors_Bayes;
