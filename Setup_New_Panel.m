clear all
Model.NbyT =    1;                  %1 for data organized NbyT, 0 for data organized TbyN
Model.LastPeriodOutOfSample = 1;    %1 last period out of sample, 0 otherwise
Model.N    =    30;                 %Number of Individuals
Model.T    =    7;                  %Number of Observations per individual
Model.K    =    5;                  %Number of coefficents, including constant
Model.Obs  =    Model.N*Model.T;    %Total number of observations
Model.S1   =    1000;               %Number of burn in replications
Model.S2   =    2000;               %Number of retained replications
Model.S    =    Model.S1 + Model.S2;%Total number of replications
Model.PrintIteration = 500;
Model.Path = 'C:\Teaching\StandardCourses\MSF_567_Bayesian_Econometrics\AAA\Class_III\RandomCoefficientModel_Koop\Example_Simulated_Data\';
Model.DataFileName = 'SampleData.csv';
Model      =   Load_Panel_Data(Model);
Model      =   LRM_Panel_Setup(Model);
Model      =   LRM_Panel_GibbsSampling(Model);

%plot(Model.Ystar,Model.Ystarhat,'.b',Model.Ystar,Model.YstarhatOLS,'.g')
C = corr([Model.Ystar,Model.Ystarhat,Model.YstarhatOLS]);