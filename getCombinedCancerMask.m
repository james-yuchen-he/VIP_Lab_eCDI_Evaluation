function combinedCancerMask = getCombinedCancerMask(patient)
    if patient.numTumor == 0 
        combinedCancerMask = logical(zeros(size(patient.pMask)));
        return
    else
        cancerMask = patient.tumors(1).lesion;
        for i = 1:patient.numTumor
            cancerMask(cancerMask == 1 | patient.tumors(i).lesion == 1) = 1;
        end
        combinedCancerMask = cancerMask;
        return
    end
end

