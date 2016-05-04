function [EWBenchMark] = CalcEQBenchMark(Returns)

[ N T ] = size(Returns);
EWBench = nanmean(Returns);

EWBenchMark =   zeros(N,T);

for n=1:N
    EWBenchMark(n,:) =EWBench;
end

