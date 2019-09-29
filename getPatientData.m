function [patients,numPatientsWithTumor] = getPatientData()
%GETPATIENTDATA Summary of this function goes here
%   Detailed explanation goes here
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
        %patients(i).tumors = struct('lesion',{},'gleasonScore',{});
        patients(i).numTumor = 0;
        patients(i).patientID = patientID(i,:);
        lesionFilePath = append(path,'Lesion_',patientID(i,:),'.mat');
        if isfile(lesionFilePath)%if there is tumor, fill in the tumor struct
            lesionData = load(lesionFilePath).Lesion;
            lesionInfoPath = append(path,'Lesion_info_',patientID(i,:),'.mat');
            lesionInfo = load(lesionInfoPath).Lesion_info;
            patients(i).numTumor = length(lesionInfo);
            for numTumors = 1:patients(i).numTumor
               patients(i).tumors(numTumors).gleasonScore = lesionInfo(numTumors);
               patients(i).tumors(numTumors).lesion = logical(lesionData{numTumors});
            end
        end
    end
    return
end

