# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'

# %%
import matplotlib.pyplot as plt
import scipy.stats
import numpy as np
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import KFold
import math
import pdb
from functions import getCancerPixels, getCombinedCancerMask, getPatientData

# %%
def train(training_set=None):
    assert type(patients == "numpy.ndarray"), "The input parameter must be a numpy array"
    assert training_set.all() != None, "empty training set"
    assert  len(training_set) != 0, "empty training set"
    adcCancerPixels,    cdiCancerPixels    = np.array([]), np.array([])
    adcNonCancerPixels, cdiNonCancerPixels = np.array([]), np.array([])

    for patient in training_set:
        if patient["numTumor"] != 0:
            cancerPixelsTmp = getCancerPixels(patient,"adc",True)
            adcCancerPixels = np.append(adcCancerPixels, cancerPixelsTmp)
            cancerPixelsTmp = getCancerPixels(patient,"cdi",True)
            cdiCancerPixels = np.append(cdiCancerPixels, cancerPixelsTmp)
        nonCancerPixelsTmp = getCancerPixels(patient,"adc",False)
        adcNonCancerPixels = np.append(adcNonCancerPixels, nonCancerPixelsTmp)
        nonCancerPixelsTmp = getCancerPixels(patient,"cdi",False)
        cdiNonCancerPixels = np.append(cdiNonCancerPixels, nonCancerPixelsTmp)
    
    dist = scipy.stats.norm
    mu_adc_cancer,     std_adc_cancer     = dist.fit(adcCancerPixels)
    mu_adc_non_cancer, std_adc_non_cancer = dist.fit(adcNonCancerPixels)

    x_adc_cancer_lower_bound = mu_adc_cancer    -4*std_adc_cancer
    x_adc_cancer_upper_bound = mu_adc_cancer    +4*std_adc_cancer
    x_adc_non_cancer_lower_bound = mu_adc_non_cancer-4*std_adc_non_cancer
    x_adc_non_cancer_upper_bound = mu_adc_non_cancer+4*std_adc_non_cancer

    x_adc_lower_bound = x_adc_cancer_lower_bound if x_adc_cancer_lower_bound < x_adc_non_cancer_lower_bound else x_adc_non_cancer_lower_bound
    x_adc_upper_bound = x_adc_cancer_upper_bound if x_adc_cancer_upper_bound > x_adc_non_cancer_upper_bound else x_adc_non_cancer_upper_bound

    x_adc = np.linspace(x_adc_lower_bound, x_adc_upper_bound, 200)

    y_adc_cancer     = dist.pdf(x_adc,    mu_adc_cancer,    std_adc_cancer)
    y_adc_non_cancer = dist.pdf(x_adc,mu_adc_non_cancer,std_adc_non_cancer)

    mu_cdi_cancer,     std_cdi_cancer     = dist.fit(np.log(cdiCancerPixels[cdiCancerPixels>0]))
    mu_cdi_non_cancer, std_cdi_non_cancer = dist.fit(np.log(cdiNonCancerPixels[cdiNonCancerPixels>0]))

    x_cdi_cancer_lower_bound = mu_cdi_cancer    -4*std_cdi_cancer
    x_cdi_cancer_upper_bound = mu_cdi_cancer    +4*std_cdi_cancer
    x_cdi_non_cancer_lower_bound = mu_cdi_non_cancer - 4*std_cdi_non_cancer
    x_cdi_non_cancer_upper_bound = mu_cdi_non_cancer + 4*std_cdi_non_cancer

    x_cdi_lower_bound = x_cdi_cancer_lower_bound if x_cdi_cancer_lower_bound < x_cdi_non_cancer_lower_bound else x_cdi_non_cancer_lower_bound
    x_cdi_upper_bound = x_cdi_cancer_upper_bound if x_cdi_cancer_upper_bound > x_cdi_non_cancer_upper_bound else x_cdi_non_cancer_upper_bound

    x_cdi     = np.linspace(x_cdi_lower_bound, x_cdi_upper_bound, 200)

    y_cdi_cancer     = dist.pdf(x_cdi,    mu_cdi_cancer,    std_cdi_cancer)
    y_cdi_non_cancer = dist.pdf(x_cdi,mu_cdi_non_cancer,std_cdi_non_cancer)

    adc_threshold = x_adc[np.argmin(np.abs(y_adc_cancer / y_adc_non_cancer -1))]
    cdi_threshold = math.exp(x_cdi[np.argmin(np.abs(y_cdi_cancer / y_cdi_non_cancer -1))])
    print(f'ADC Decision boundary: Cancer: x < {adc_threshold} Non-Cancer: x > {adc_threshold}')
    print(f'CDI Decision boundary: Cancer: x > {cdi_threshold} Non-Cancer: x < {cdi_threshold}')

    return adc_threshold, cdi_threshold


# %%
def evaluate(test_set, adc_threshold, cdi_threshold):
    cdi_cm = np.zeros(dtype=float,shape=(2,2))
    adc_cm = np.zeros(dtype=float,shape=(2,2))
    for patient in test_set:
        adc_prediction = np.zeros(shape=patient["pMask"].shape, dtype=bool)
        adc_prediction[np.logical_and(patient["pMask"] == True, patient["adc"] > adc_threshold)] = False
        adc_prediction[np.logical_and(patient["pMask"] == True, patient["adc"] < adc_threshold)] = True
        adc_prediction[patient["pMask"] == False] = False
        adc_cm += confusion_matrix(getCombinedCancerMask(patient)[patient["pMask"] == True].flatten(), adc_prediction[patient["pMask"] == True].flatten())
        
        cdi_prediction = np.zeros(shape=patient["pMask"].shape, dtype=bool)
        cdi_prediction[np.logical_and(patient["pMask"] == True, patient["cdi"] > cdi_threshold)] = True
        cdi_prediction[np.logical_and(patient["pMask"] == True, patient["cdi"] < cdi_threshold)] = False
        cdi_prediction[patient["pMask"] == False] = False
        cdi_cm += confusion_matrix(getCombinedCancerMask(patient)[patient["pMask"] == True].flatten(), cdi_prediction[patient["pMask"] == True].flatten())
        
    return adc_cm, cdi_cm


# %%
[patients, numPatients, numPatientsWithTumor] = getPatientData()
patients = np.array([item for item in patients if item["patientID"] != "P00000249"])
adc_confusion_matrix_sum = np.zeros(shape=(2,2))
cdi_confusion_matrix_sum = np.zeros(shape=(2,2))

kf = KFold(n_splits=5, shuffle=False)

for num, (patients_train_index, patients_test_index) in enumerate(kf.split(patients)):
    print(f'Trial {num}: training sample size:{len(patients_train_index)} test sample size:{len(patients_test_index)}')
    adc_threshold, cdi_threshold = train(patients[patients_train_index])
    adc_confusion_matrix, cdi_confusion_matrix = evaluate(patients[patients_test_index],adc_threshold, cdi_threshold)
    adc_confusion_matrix_sum += adc_confusion_matrix
    cdi_confusion_matrix_sum += cdi_confusion_matrix
    print(f'ADC tn={adc_confusion_matrix[0,0]}, fp={adc_confusion_matrix[0,1]}, fn={adc_confusion_matrix[1,0]}, tp={adc_confusion_matrix[1,1]} total voxel={np.sum(adc_confusion_matrix)}')
    print(f'CDI tn={cdi_confusion_matrix[0,0]}, fp={cdi_confusion_matrix[0,1]}, fn={cdi_confusion_matrix[1,0]}, tp={cdi_confusion_matrix[1,1]} total voxel={np.sum(cdi_confusion_matrix)}')
adc_sensitivity = adc_confusion_matrix[1,1] / (adc_confusion_matrix[1,1] + adc_confusion_matrix[1,0])
adc_specificity = adc_confusion_matrix[0,0] / (adc_confusion_matrix[0,0] + adc_confusion_matrix[0,1])
adc_accuracy = (adc_confusion_matrix[0,0] + adc_confusion_matrix[1,1])/(np.sum(adc_confusion_matrix))
print("ADC:")
print(f'sensitivity = {adc_sensitivity:.6f} specificity = {adc_specificity:.6f} accuracy = {adc_accuracy:.6f}')

cdi_sensitivity = cdi_confusion_matrix[1,1] / (cdi_confusion_matrix[1,1] + cdi_confusion_matrix[1,0])
cdi_specificity = cdi_confusion_matrix[0,0] / (cdi_confusion_matrix[0,0] + cdi_confusion_matrix[0,1])
cdi_accuracy = (cdi_confusion_matrix[0,0] + cdi_confusion_matrix[1,1])/(np.sum(cdi_confusion_matrix))
print("CDI:")
print(f'sensitivity = {cdi_sensitivity:.6f} specificity = {cdi_specificity:.6f} accuracy = {cdi_accuracy:.6f}')

