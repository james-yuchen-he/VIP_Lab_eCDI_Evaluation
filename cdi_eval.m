patientID = load('./patientID.mat').patientID;
patientIDWithTumor = load('./posID.mat').caseID;
numPatients = length(patientID);
numPatientsWithTumor = length(patientIDWithTumor);
%%This outer for loop organize data into an array of struct. 
%Each struct represent a patient.
for i = 1:numPatients
    path = append('./100_original_anonymized_images_ExportedMatlab/',patientID(i,:),'/');
    patients(i).adc = load(append(path,'ADC_',patientID(i,:),'.mat')).ADC;
    patients(i).cdi = load(append(path,'CDI_matlab_',patientID(i,:),'.mat')).CDI_matlab;
    patients(i).pMask = logical(load(append(path,'PMask0_',patientID(i,:),'.mat')).Mask0);%convert mask into logical array
    patients(i).tumors = struct('lesion',{},'gleasonScore',{});
    patients(i).numTumor = 0;
    patients(i).patientID = patientID(i,:);
    patients(i).tumorCombined = logical(zeros(size(patients(i).pMask)));
    lesionFilePath = append(path,'Lesion_',patientID(i,:),'.mat');
    if isfile(lesionFilePath)%if there is tumor, fill in the tumor struct
        lesionData = load(lesionFilePath).Lesion;
        lesionInfoPath = append(path,'Lesion_info_',patientID(i,:),'.mat');
        lesionInfo = load(lesionInfoPath).Lesion_info;
        patients(i).numTumor = length(lesionInfo);
        for numTumors = 1:patients(i).numTumor
           patients(i).tumors(numTumors).gleasonScore = lesionInfo(numTumors);
           patients(i).tumors(numTumors).lesion = logical(lesionData{numTumors});
           %currently two classes: GS >=6 vs  GS = 0.
           %need to add more finer classification. i.e. GS = 6 vs GS < 6
           %and GS >= 7 vs GS = 0 
        end
        patients(i).tumorCombined = getCombinedCancerMask(patients(i));
    end
end

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

