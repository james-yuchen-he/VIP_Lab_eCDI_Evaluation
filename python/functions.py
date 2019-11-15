# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'

# %%
import os.path
import scipy.io
import numpy as np

# %% [markdown]
# ## getPatientData()

# %%
def getPatientData():
    import scipy.io
    """
    Read the data from the stray .mat files and store it in a
    list of dictionaries, with each list element containing a patient
    e.g.
        patients = 
    
        1x104 list with each entry being a dict with fields:
        
            adc       | (144x144x26) ndarray uint16
            cdi       | (144x144x26) ndarray float64
            pMask     | (144x144x26) ndarray bool
            numTumor  | int
            patientID | str
            tumors    | dict
            
        patients[2]['tumors'] = 
        
        a dict with fields:
        
            gleasonScore | list of int
            lesion       | list of (144x144x26) ndarray bool
    
    """
    patientID = scipy.io.loadmat('../patientID.mat')['patientID']
    numPatients = len(patientID)

    patientIDWithTumor = scipy.io.loadmat('../posID.mat')['caseID']
    numPatientsWithTumor = len(patientIDWithTumor)
    
    patients = [dict() for x in range(numPatients)]
    
    for i in range(numPatients):
        path = "../100_original_anonymized_images_ExportedMatlab/" + patientID[i] + '/'
        patients[i] = {
            "adc": scipy.io.loadmat(path + 'ADC_' + patientID[i] + '.mat')['ADC'],
            "cdi": scipy.io.loadmat(path + 'CDI_matlab_' + patientID[i] + '.mat')['CDI_matlab'],
            "pMask": np.array(scipy.io.loadmat(path + 'PMask0_' + patientID[i] + '.mat')['Mask0'],dtype=bool),
            "numTumor": 0,
            "patientID": patientID[i]
        }
        lesionFilePath = path + 'Lesion_' + patientID[i] + '.mat'
        if os.path.exists(lesionFilePath):
            lesionData = scipy.io.loadmat(lesionFilePath)['Lesion']
            lesionInfo = scipy.io.loadmat(path + 'Lesion_info_' + patientID[i] + '.mat')['Lesion_info']
            patients[i]["numTumor"] = len(lesionInfo[0])
            patients[i]["tumors"] = {
                "gleasonScore": [lesionInfo[0][j] for j in range(len(lesionInfo[0]))],
                "lesion": [np.array(lesionData[0][j],dtype=bool) for j in range(len(lesionInfo[0]))]
            }
    return [patients, numPatients, numPatientsWithTumor]
    

# %% [markdown]
# ## getCancerPixels()

# %%
def getCancerPixels(patient,modality,cancer):
    """
    Takes a patient dictionary as input and return an array that contains 
    all the pixels identified as either cancerous or non-cancerous
    @input: patient: a dictionary containing all relevant info about a patient
    @input modality: a string indicating adc or cdi
    @input: cancer: a bool indicating whether the pixels should be cancerous or not
    @output: an array containing all the identified pixels of the given modality
    """
    cancerMask = getCombinedCancerMask(patient)
    try:
        return patient[modality][np.logical_and(patient["pMask"]==True,cancerMask==cancer)]
    except:
        print("Unrecognized modality.")
        return 0
    

# %% [markdown]
# ## getCombinedCancerMask()

# %%
def getCombinedCancerMask(patient):
    if patient['numTumor'] == 0:
        return np.zeros(np.shape(patient['pMask']))
    else:
        cancerMask = patient["tumors"]["lesion"][0]
        for i in range(patient["numTumor"]):
            cancerMask[np.logical_or(cancerMask == 1,patient["tumors"]["lesion"][i] == 1)] = 1
        return cancerMask
    


# %%


