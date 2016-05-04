function OutputModel = LRM_Panel_Setup_RemoveBadRank(InputModel)
OutputModel = InputModel;
InputModel.RankGood = zeros(InputModel.N,1);
for n = 1:InputModel.N
    T               =   (n-1)*InputModel.T+1:n*InputModel.T;  %All observatons of the individual (data is stacked)
    xuse            =   InputModel.X(T,:);           %X matrix for the individual
    yuse            =   InputModel.Y(T,1);           %Y matrix for the individual
    if rank(xuse) == size(xuse,2)
        InputModel.RankGood(n) = 1;
    end
end
OutputModel.N   = sum(InputModel.RankGood);
ncounter        = 0;
for n = 1:InputModel.N
    T               =   (n-1)*InputModel.T+1:n*InputModel.T;  %All observatons of the individual (data is stacked)
    xuse            =   InputModel.X(T,:);           %X matrix for the individual
    yuse            =   InputModel.Y(T,1);           %Y matrix for the individual
    if InputModel.RankGood(n) == 1
        ncounter            = ncounter+1;
        T1                  =   (ncounter-1)*InputModel.T+1:ncounter*InputModel.T;
        OutputModel.X(T1,:) =   xuse;
        OutputModel.Y(T1,1) =   yuse;
    end
end