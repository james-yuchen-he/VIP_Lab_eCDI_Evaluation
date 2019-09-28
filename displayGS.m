function displayGS(patients)
    for i = 1:numel(patients)
        if patients(i).numTumor ~= 0
            msg = [patients(i).patientID ':  '];
            for j = 1:patients(i).numTumor
                msg = append(msg,'   ',int2str(patients(i).tumors(j).gleasonScore));
            end
            disp(msg);
        end
    end
end

