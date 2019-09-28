function nonCancerPixels = getNonCancerPixels(patient,modality)
%Takes a patient struct as input and return an array that contains all the
%non-cancerous pixels
%@input: patient: a struct containing all the relavent info about a patient
%@input: modality: a string indicating adc or cdi
%@output: cancerPixels: an array containing all the non-cancer pixels of the
%given modality

%BAD design of function structure. NEED to refactor!!!!
    cancerMask = getCombinedCancerMask(patient);
    if modality == 'cdi'
        nonCancerPixels = patient.cdi(patient.pMask == 1 & cancerMask == 0);
    else
        nonCancerPixels = patient.adc(patient.pMask == 1 & cancerMask == 0);
    end
    return
end


