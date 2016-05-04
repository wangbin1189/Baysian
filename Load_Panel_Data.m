function InputModel    = Load_Panel_Data(InputModel)
T   =   InputModel.T;
N   =   InputModel.N;
NT  =   N*T;
if InputModel.NbyT == 1
    Data            =   csvread([InputModel.Path InputModel.DataFileName],1,0);
    InputModel.Y    =   Data(:,1);               %Dependent variable 
    InputModel.X    =   Data(:,2:end);           %Independent variables, including a constant
end
if InputModel.NbyT == 0
    Data                =   csvread([InputModel.Path InputModel.DataFileName]);
    IndividualVector    =   nan(NT,1);
    Counter             =   0;
    for i = 1:T
        for j = 1:N
            Counter     = Counter + 1;
            IndividualVector(Counter) = j;
        end
    end
    Data2               =   [IndividualVector Data];
    Data2               =   sortrows(Data2,1);
    InputModel.Y        =   Data2(:,2);               %Dependent variable 
    InputModel.X        =   Data2(:,3:end);           %Independent variables, including a constant
    InputModel.Data2    = Data2;
    InputModel.Data     = Data;
end
InputModel              =   LRM_Panel_Setup_RemoveBadRank(InputModel);
T                       =   InputModel.T;
N                       =   InputModel.N;
NT                      =   N*T;
if InputModel.LastPeriodOutOfSample == 1
   FullSample           =   1:NT;
   OutSample            =   (NT-N+1):NT;
   InSample             =   setdiff(FullSample,OutSample);
   Y                    =   InputModel.Y;
   X                    =   InputModel.X;
   InputModel.Y         =   Y(InSample);
   InputModel.X         =   X(InSample,:);
   InputModel.Ystar     =   Y(OutSample);
   InputModel.Xstar     =   X(OutSample,:);
   InputModel.T         =   T-1;
end