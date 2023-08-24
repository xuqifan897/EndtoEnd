# Beam Orientation Optimization for 4pi Radiotherapy (CPU & GPU)
This project aims to be an efficient beam angle selection algorithm for 4pi Radiotherapy

## Getting Started
### Dose calculation
Step 1: Make folder: patFolder/dicomdata. Put all dicom files (CT, RTstructure, etc) in the folder

Step 2: Run `PreDosecalc_IMRT.m` to generate `beamlist.txt`, `structures.json`, and `CERR_bz2` (optional)

		Note that the following variables need to be changed:
		    'patFolder', 'patientName', 'PresriptionDose', 'beamlogfile' 
		Delete structures in `structures.json` if necessary
		
Step 3: Copy `config.json` and `run_dosecalc.sh` files from <dosecalc_template> folder

		Open run_dosecalc.sh in gedit; change name of BBOX (body contour) if necessary

Step 4: Run dose calculation on microway

		Move <patFolder> to microway
		Open terminal
		change directory to <patFolder>
			cd <patFolder>
		Run dosecalc script
			bash run_dosecalc.sh
		dose results (M-matrix) will be saved to Dose_Coefficients.h5
		structure masks will be saved to Dose_Coefficients.mask

For details on dose calculation, refer to [neph-gpu-dosecalc-paper].


### Optimization
Step 1: Run `PreOptimize_IMRT.m` to generate `params0.mat`, `StructureInfo0.mat`, `[patientName]_M.mat`, `[patientName]_M_DS.mat`
* `params0.mat:`  parameters for FISTA
* `StructureInfo0.mat:`  Structure weightings
* `[patientName]_M.mat:` full-sampled M-matrix
* `[patientName]_M_DS.mat:` down-sampled M-matrix by replacing elements below a threshold with 0 (can be adjusted with 'thresh')

Step 2: Change weightings in StructureInfo and parameters in params, save to files `params1.mat`, `StructureInfo1.mat`
* **StructureInfo.maxDose:** maximal dose (recommend to set it to prescription dose for PTV)
* **StructureInfo.maxWeights:** weightings for max dose penalization (increase if there is a hot spot in PTV/other structures)
* **StructureInfo.minDoseTarget:** target min dose (recommend to set it to prescription dose for PTV, and NaN for OARs)
* **StructureInfo.minDoseTargetWeights:** weightings for target min dose penalization (increase if there is a cold spot in PTV)
* **StructureInfo.OARWeights:** structure weightings of OARs (penalize doses to OAR)
* **StructureInfo.IdealDose:** prescription dose at PTV, 0 elsewhere
* Structures can be deleted by deleting rows of **StructureInfo**
* **params.eta:** increase to smooth the fluence map
* **params.maxIter:** maximal number of iterations (The optimization algorithm ends when the number of beams is less than desired, or when the maximal number of iterations is reached).
* **params.beamWeight:**  increase it if too many beams are selected, decrease it otherwise (The beamWeight is tuned automatically during the optimization, but may need manual tuning in some cases)
* **params.ChangeWeightsTrigger:**  The beamWeight is updated automatically every `ChangeWeightsTrigger` iterations
* The default values of parameters in **params** should work well in most cases


Step 3: Run `Main_4piIMRT` with `params1.mat` and `StructureInfo1.mat`. The optimization result is saved in patFolder/optimize
* Selected beam angles: `.csv` file, `result.gantryVarianIEC` and `result.couchVarianIEC`
* fluence map: `result.xPolish`
* dose: `result.dose`

For details on optimization algorithm, refer to [O'Connor-boo-paper] and [Lyu-boo-paper].


## Useful things to know
* Example data can be found in 'ShengNas\SharedProjectData\4pi BOO'
* It is recommended to run optimization using GPU mode whenever GPU memory fits. 
* **DoseInfo** includes the dose of all plans. To compare different plans, use 
```matlab
        plotDVH_QL(DoseInfo([Ind]), strNum, StructureInfo, numBins, 0);
```
where **Ind** specifies the indexes of the plans to compare, and **strNum** are the indexes of the structures  


## Resolving Common Issues:
### Dose Calculation GPU Memory Errors
Most issues are related to insufficient GPU memory availability for the selected quality parameters. In these events, the code will produce an error resembling the following:

> ___CUDA error at [...] code=2(cudaErrorMemoryAllocation) "cudaMalloc(...)"___  
    - or -  
> ___Request to allocate %d bytes of device memory exceeded remaining memory (%d bytes)___

To work around the error: Open the `config.yml` file that is supplied to `dosecalc-preprocess [...] --config="<config-file>" [...]` and decrease the numbers in `"max-rev-size": [<int>, <int>, <int>]` if possible (if no new errors are produced)

### Optimization Memory Errors
* Reduce 'kernel-extent' in `config.json` and rerun dose calculation; 
    * In this case, it is recommended to finalize dose calculation and fluence map optimization (polish step) using the selected 20 beams from optimization and  larger 'kernel-extent' (>2).
* Increase 'thresh' in `PreOptimize_IMRT.m` 
    * Never use thresh > 0.05; Recommend thresh < 0.01


[Neph-gpu-dosecalc-paper]: https://aapm.onlinelibrary.wiley.com/doi/10.1002/mp.13651
[O'Connor-boo-paper]: https://arxiv.org/abs/1710.05308
[Lyu-boo-paper]: https://iopscience.iop.org/article/10.1088/1361-6560/ab63b8

