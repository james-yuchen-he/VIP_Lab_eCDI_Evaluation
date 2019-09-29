function displayGS(patients)
%Prints the gleason score of every tumor of every patient
%input: patients: An array of Struct. Each struct contains the info of one
%patient
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

