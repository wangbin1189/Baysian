N = 3;
T = 5;
NT = N*T;
FullSample   =   1:NT;
OutSample    =   (NT-N+1):NT;
InSample     =   setdiff(FullSample,OutSample);