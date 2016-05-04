function [Model]    = LRM_Panel_Setup2(InputModel)
Model.Y     =   InputModel.Y;               %Dependent variable 
Model.X     =   InputModel.X;               %Independent variables
Model.N     =   InputModel.N;               %Number of Individuals
Model.T     =   InputModel.T;               %Number of Observations per individual
Model.K     =   InputModel.K;               %Number of coefficents
%Model.Obs   =   InputModel.N*InputModel.T;  %Total number of observations
Model.Obs   =   length(find(sum(isnan([Model.X Model.Y])')'==0));%Total number of observations(missing data)
Model.S1    =   InputModel.S1;              %Number of burn in replications
Model.S2    =   InputModel.S2;              %Number of retained replications
Model.S     =   Model.S1 + Model.S2;        %Total number of replications
Model.PrintIteration    =   InputModel.PrintIteration;%Print message indicating iteration reached

Model.GoodObs       =   sum(isnan([Model.X Model.Y])')'==0;

Model.V_Error_Pr    =   0; 
Model.H_Error_Pr    =   0;
Model.D_Error_Pr    =   0; 
Model.D_Error_Ps    =   0; 
Model.H2_Error_Pr   =   0;
Model.S2_Error_Pr   =   0;            
Model.DS2_Error_Pr  =   0;
Model.D_VTheta_Pr   =   0;
Model.D_VTheta_Ps   =   0; 
Model.SE_OLS        =   0;

Model.ThetaDraw     =   zeros(Model.K, Model.N);  %K by N
Model.ThetaMean     =   zeros(Model.K, Model.N);  %K by N
Model.ThetaStd      =   zeros(Model.K, Model.N);  %K by N
Model.SV_Theta      =   zeros(Model.S2,Model.K^2);


%Compute OLS estimates coefficents for each individual
errors_OLS      =   zeros(Model.N*Model.T,1)*NaN;
Thetas_OLS      =   zeros(Model.K,Model.N)*NaN;
for n = 1:InputModel.N
    T           =   (n-1)*InputModel.T+1:n*InputModel.T;  %All observatons of the individual (data is stacked)
    good        =   find(Model.GoodObs(T)==1);
    if ~isempty(good)
        xuse            =   InputModel.X(T,:);           %X matrix for the individual
        yuse            =   InputModel.Y(T,1);           %Y matrix for the individual
        good            =   find(sum(isnan([xuse yuse])')'==0);
        Thetas_OLS(:,n) =   (inv(xuse(good,:)'*xuse(good,:))*xuse(good,:)'*yuse(good));   %Coefficients for the individual estimated by OLS
        xuse            =   InputModel.X(T,:);              %X matrix for the individual
        yuse            =   InputModel.Y(T,1);              %Y matrix for the individual
        e_OLS           =   (yuse-xuse*Thetas_OLS(:,n));    %Errors for the individual estimated by OLS
        errors_OLS(T,1) =   e_OLS;
    end
end
Model.Thetas_OLS    =   Thetas_OLS;
Model.errors_OLS    =   errors_OLS;

%==========Define all the Hyperparameters=======================================
% Desc: Fills in the priors needed to calculate all the hyperparametrs
% 
% Hyperparameters are :
%   1. Theta:   The coefficients  (dK)
%   2. MTheta:  Mean of the Coefficients (dK)
%   3. VTheta:  Var  of the Coefficients (dKK)
%   4. VError:  Var  of the Errors, (or Error Precision: h) (Sc.V_Error)

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Mean of Theta = M_Theta   (NORMAL DISTRIBUTION)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%       Theta   ~ N(M_Theta,  V_Theta)          Eq 7.21
%       M_Theta ~ N(M_MTheta, V_MTheta)         Eq 7.22

Model.M_MTheta_Pr       =   zeros(Model.K,1);   %Center prior over 0 , for no effect
Model.V_MTheta_Pr       =   5*eye(Model.K);     %Allow substatial variation in prior
Model.V_MTheta_Inv_Pr   =   inv(Model.V_MTheta_Pr);



%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Var of Theta = V_Theta   (Inverse of V_Theta = WISHART DISTRIBUTION)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%       Theta      ~ N(M_Theta,  V_Theta)       Eq 7.21
%       V_ThetaInv ~ W(rho, V_ThetaInv_Pr)      Eq 7.23
% degree of freedom = rho, scale matrix R --- implying mean = rho*R
rho                 =   2;
Model.D_VTheta_Pr   =   2;
NotNaN              =   find(sum(isnan(Model.Thetas_OLS)) ==0);
Model.R             =   rho*inv(cov(Model.Thetas_OLS(:,NotNaN)'));
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Var of Error = V_Error   (Inverse of V_Error (h) = GAMMA DISTRIBUTION)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% V_Error_Inv ~ G(S_Error, D_Error)             Eq 7.24
% S_Error = Point Estimate, D_Error = D.O.F

%Model.D_Error_Pr    =   1;
%Model.H2_Error_Pr   =   5;
%Model.S2_Error_Pr   =   1/Model.H2_Error_Pr;

%Use OLS regression to get starting values of Theta
ValidObs        =   find(Model.GoodObs==1);
Model.ValidObs  =   ValidObs;
Model.Theta_OLS =   inv(Model.X(ValidObs,:)'*Model.X(ValidObs,:))*Model.X(ValidObs,:)'*Model.Y(ValidObs);   
%inv(Model.X(ValidObs,:)'*Model.X(ValidObs,:))*Model.X(ValidObs,1)'*Model.Y(ValidObs);             %OLS coefficients across all data
Model.e_OLS     =   (Model.Y(ValidObs) - Model.X(ValidObs,:)*Model.Theta_OLS);                %OLS errors
Model.SE_OLS    =   (Model.e_OLS'*Model.e_OLS)/(Model.Obs - Model.K);   %Standard Error (variance of OLS errors)
Model.MTheta    =   Model.Theta_OLS;                                    %Start Theta Mean with OLS estimate

%Model.V_Error_Pr    =   Model.SE_OLS;
%Model.H_Error_Pr    =   1/Model.SE_OLS;%starting value of h, the standard errors from the OLS regression
Model.D_Error_Pr    =   1;
Model.H2_Error_Pr   =   1/Model.SE_OLS;     %starting value of h, the standard errors from the OLS regression
Model.S2_Error_Pr   =   1/Model.H2_Error_Pr;

Model.ThetaDraw     =   repmat(Model.Theta_OLS, 1, Model.N);% Init Theta Draws for Gibbs Sampling as OLS estimates

Model.D_Error_Ps        =   Model.D_Error_Pr + Model.Obs;          %Posterior D.o.F for Error (Eq 7.28, Ln 2)
NT                      =   Model.N*Model.T;
N2                      =   ceil(NT/Model.Obs);                       %Correction for unbalanced panels
Model.D_VTheta_Ps       =   Model.D_VTheta_Pr + N2;             %Posterior D.o.F for Var of Theta (Eq 7.27, Ln 4)
%Model.D_VTheta_Ps       =   Model.D_VTheta_Pr + Model.N;           %Posterior D.o.F for Var of Theta (Eq 7.27, Ln 4)
Model.DS2_Error_Pr      =   Model.D_Error_Pr * Model.S2_Error_Pr; 

