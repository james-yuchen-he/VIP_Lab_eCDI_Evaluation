[patients, ~] = getPatientData();
%%
adcCancerPixels = [];
adcNonCancerPixels = [];
cdiCancerPixels = [];
cdiNonCancerPixels = [];
for i = 1:(numPatients - 1)
    if patients(i).numTumor ~= 0
        cancerPixelsTmp = getCancerPixels(patients(i),'adc');
        adcCancerPixels = vertcat(adcCancerPixels, cancerPixelsTmp);
        cancerPixelsTmp = getCancerPixels(patients(i),'cdi')
        cdiCancerPixels = vertcat(cdiCancerPixels,cancerPixelsTmp);
    end
    nonCancerPixelsTmp = getNonCancerPixels(patients(i),'adc');
    adcNonCancerPixels = vertcat(adcNonCancerPixels,nonCancerPixelsTmp);
    nonCancerPixelsTmp = getNonCancerPixels(patients(i),'cdi');
    cdiNonCancerPixels = vertcat(cdiNonCancerPixels,nonCancerPixelsTmp);
end

x = [1: 2000];

paramEstCancerADC = mle(adcCancerPixels);% compute the mean and std through mle
paramEstNonCancerADC = mle(adcNonCancerPixels);
yCancerADC = normpdf(x,paramEstCancerADC(1),paramEstCancerADC(2));
yNonCancerADC = normpdf(x,paramEstNonCancerADC(1),paramEstNonCancerADC(2));


plot(x,yCancerADC)
hold on;
plot(x,yNonCancerADC);

% x = [1: 10];
% paramEstCancerCDI = mle(cdiCancerPixels);% compute the mean and std through mle
% paramEstNonCancerCDI = mle(cdiNonCancerPixels);
% yCancerCDI = normpdf(x,paramEstCancerCDI(1),paramEstCancerCDI(2));
% yNonCancerCDI = normpdf(x,paramEstNonCancerCDI(1),paramEstNonCancerCDI(2));
% 
figure(2);
% plot(x,yCancerCDI)
% hold on;
% plot(x,yNonCancerCDI);
histogram(cdiNonCancerPixels,20)
figure(3);
histogram(cdiCancerPixels,20)

[~, threshold] = min(abs(yCancerADC./yNonCancerADC-1));
[~, threshold] = min(abs(yCancerADC./yNonCancerADC-1));
truePositiveADC = patients(numPatients).adc((patients(numPatients).adc < threshold & patients(numPatients).pMask == 1) &  patients(numPatients).tumorCombined == 1);
falsePositiveADC = patients(numPatients).adc((patients(numPatients).adc < threshold & patients(numPatients).pMask == 1) &  patients(numPatients).tumorCombined == 0);
trueNegativeADC = patients(numPatients).adc((patients(numPatients).adc > threshold & patients(numPatients).pMask == 1) &  patients(numPatients).tumorCombined == 0);
falseNegativeADC = patients(numPatients).adc((patients(numPatients).adc > threshold & patients(numPatients).pMask == 1) &  patients(numPatients).tumorCombined == 1);

ntruePositiveADC = length(truePositiveADC);
nfalsePositiveADC = length(falsePositiveADC);
ntrueNegativeADC = length(trueNegativeADC);
nfalseNegativeADC = length(falseNegativeADC);

specificity = ntrueNegativeADC / (ntrueNegativeADC + nfalsePositiveADC)
sensitivity = ntruePositiveADC / (ntruePositiveADC + nfalseNegativeADC)
accuracy = (ntrueNegativeADC + ntruePositiveADC) / (nfalseNegativeADC + nfalsePositiveADC + ntrueNegativeADC + ntruePositiveADC)

