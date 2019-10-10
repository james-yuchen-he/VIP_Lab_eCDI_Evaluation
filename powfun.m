function newHist = powfun(oldHist,a)
%POWFUN Summary of this function goes here
%   Detailed explanation goes here

newHist = zeros(1,length(oldHist));

for i = 1:length(oldHist)
    if a > 0
        newHist(i) = oldHist(i)^a;
    elseif a == 0
        newHist(i) = log(oldHist(i));
    else
        newHist(i) = -1*oldHist(i)^a;
    end
end

