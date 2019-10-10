%NOTE: This file must be placed at the immediate parent directory of /100_original_anonymized_images_ExportedMatlab
[patients, numPatients, ~] = getPatientData();
adcCancerPixels = [];
adcNonCancerPixels = [];
cdiCancerPixels = [];
cdiNonCancerPixels = [];
for i = 1:(numPatients - 1)
    if patients(i).numTumor ~= 0
        cancerPixelsTmp = getCancerPixels(patients(i),'adc');
        adcCancerPixels = vertcat(adcCancerPixels, cancerPixelsTmp);
        cancerPixelsTmp = getCancerPixels(patients(i),'cdi');
        cdiCancerPixels = vertcat(cdiCancerPixels,cancerPixelsTmp);
    end
    nonCancerPixelsTmp = getNonCancerPixels(patients(i),'adc');
    adcNonCancerPixels = vertcat(adcNonCancerPixels,nonCancerPixelsTmp);
    nonCancerPixelsTmp = getNonCancerPixels(patients(i),'cdi');
    cdiNonCancerPixels = vertcat(cdiNonCancerPixels,nonCancerPixelsTmp);
end

x = [1: 2000];

%%
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

%%
% Do for both cdiNonCancerPixels and cdiCancerPixels
% Do for both fine and coarse alpha searches

figure(1);
start_ = -1;
end_ = 1;
a = start_ : (end_-start_)/8 : end_;

for i = 1:9
    subplot(3,3,i)
    histogram(powfun(cdiNonCancerPixels,a(i)),200)
    title("alpha = "+num2str(a(i)))
end

sgtitle('CDI Non Cancer Transformed Distributions')

%%
% Do for both cdiNonCancerPixels and cdiCancerPixels

figure(2);
histogram(powfun(cdiNonCancerPixels,0),200)
title('Log Transform')

%%
figure(3);
sgtitle("Log Transforms")

subplot(2,2,1)
histogram(powfun(cdiNonCancerPixels,0),200)
title('CDI Non Cancer')

subplot(2,2,2)
histogram(powfun(cdiNonCancerPixels,0),200)
hold on
histogram(powfun(cdiCancerPixels,0),200)
title('Not Normalized')

subplot(2,2,3)
histogram(powfun(cdiCancerPixels,0),200,'FaceColor',[0.8500, 0.3250, 0.0980])
title('CDI Cancer')

subplot(2,2,4)
histogram(powfun(cdiNonCancerPixels,0),200,'Normalization','pdf')
hold on
histogram(powfun(cdiCancerPixels,0),200,'Normalization','pdf')
title('Normalized (PDF)')

%%
[~, threshold] = min(abs(yCancerADC./yNonCancerADC-1));
[~, threshold] = min(abs(yCancerADC./yNonCancerADC-1));
truePositiveADC = patients(numPatients).adc((patients(numPatients).adc < threshold & patients(numPatients).pMask == 1) &  getCombinedCancerMask(patients(numPatients)) == 1);
falsePositiveADC = patients(numPatients).adc((patients(numPatients).adc < threshold & patients(numPatients).pMask == 1) & getCombinedCancerMask(patients(numPatients)) == 0);
trueNegativeADC = patients(numPatients).adc((patients(numPatients).adc > threshold & patients(numPatients).pMask == 1) &  getCombinedCancerMask(patients(numPatients)) == 0);
falseNegativeADC = patients(numPatients).adc((patients(numPatients).adc > threshold & patients(numPatients).pMask == 1) & getCombinedCancerMask(patients(numPatients)) == 1);

ntruePositiveADC = length(truePositiveADC);
nfalsePositiveADC = length(falsePositiveADC);
ntrueNegativeADC = length(trueNegativeADC);
nfalseNegativeADC = length(falseNegativeADC);

specificity = ntrueNegativeADC / (ntrueNegativeADC + nfalsePositiveADC)
sensitivity = ntruePositiveADC / (ntruePositiveADC + nfalseNegativeADC)
accuracy = (ntrueNegativeADC + ntruePositiveADC) / (nfalseNegativeADC + nfalsePositiveADC + ntrueNegativeADC + ntruePositiveADC)

