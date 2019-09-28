function cancerPixels = getCancerPixels(patient, modality)
%Takes a patient struct as input and return an array that contains all the
%non-cancerous pixels
%@input: patient: a struct containing all the relavaent info about a patient
%@input: modality: a string indicating adc or cdi
%@output: cancerPixels: an array containing all the cancer pixels of the
%given modality
        cancerMask = getCombinedCancerMask(patient);
        if modality == 'cdi'
            cancerPixels = patient.cdi(patient.pMask == 1 & cancerMask == 1);
        else
            cancerPixels = patient.adc(patient.pMask == 1 & cancerMask == 1);
        end
end

