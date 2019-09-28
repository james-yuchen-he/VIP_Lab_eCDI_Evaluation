patientID.mat:	list of all valid cases
posID.mat	a subset of patientID which have ROI contour

Each case folder contain the following files:

Lesion_P		- the mask for region of interest (ROI), which might be tumour
Lesion_info		- Gleason score (GS) for the ROI
			GS>=7	 : clinically significant cancer
			GS=6	 : Clinically insignificant cancer
			GS=0	 : not cancer

ADC_info_P*.mat		- ADC info, extracted from DICOM
ADC_matlab_P*.mat	- ADC calculated in Matlab
ADC_P*.mat 		- ADC calculated by ProCanVAS 
These two ADCs are slightly different (use matlab)

CDI_matlab_P*.mat	- CDI calculated in Matlab
dualCDI_P*.mat		- dual-CDI calcualted in Matlab

DWI_info_P*.mat		- DWI info, extracted from DICOM
DWI_P*.mat		- DWI images
DWI_Reg_P*.mat		- Registered DWI images

HBV_info_P*.mat		- HBV info, extracted from DICOM (HBV: High-b value image - b1600)
HBV_matlab_P*.mat 	- HBV calculated in Matlab
HBV_P*.mat		- HBV calculated by ProCanVAS
HBV_matlab and HBV should be the same but might be slightly differnt due to implementation (use matlab)

newADC_P*.mat		- ADC calculated from registered DWI in matlab
newCDI_P*.mat 		- CDI calculated from registered DWI in matlab
newHBV_P*.mat		- HBV calculated from registered DWI in matlab

optCDI_P*.mat 		- optimized CDI calcuated in Matlab

PMask0_P*.mat		- Prostate Mask
Reg_info_P*.mat		- Registration information file for registering T2 to DWI

T2_info_P*.mat 		- T2 info, extracted from DICOM

T2_P*.mat		- T2 images
T2_Reg_P*.mat		- T2 images registered to DWI
TZMask0_P*.mat		- TZ mask

	