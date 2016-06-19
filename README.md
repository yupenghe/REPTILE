# REPTILE
### - Regulatory Element Prediction based on TIssue-specific Local Epigenetic marks

https://github.com/yupenghe/REPTILE

REPTILE is a tool to identify the precise location of enhancers by integrating histone modification data 
and base-resolution DNA methylation profiles.

Please contact [Yupeng He](mailto:yupeng.he.bioinfo@gmail.com) for for feedbacks, questions or bugs.


## Requirement

#### R
REPTILE requires R (>= 3.2.2) and several packages that REPTILE depdends on.
R is available in [R website](https://www.r-project.org/). The depedent R packages
and their versions are 

* optparse (>= 1.3.2)

* randomForest (>= 4.6-12)

* flux(>= 0.3-0)

* foreach (>= 1.4.3)

* doMC (>= 1.3.4) 

They can be installed from CRAN using the below R code.
```R
install.packages("optparse")
install.packages("randomForest")
install.packages("flux")
install.packages("foreach")
install.packages("doMC")
```

#### python
The preprocessing script and enhancer calling script in REPTILE require python2 (>= 2.7.9) or
python3 (>= 3.5.1); They also requires numpy (>= 1.10.4) and pandas (>= 0.17.1). The 
[numpy](http://www.numpy.org/) and [pandas](http://pandas.pydata.org/) are available from their 
websties. It is recommended to install [anaconda](https://www.continuum.io/downloads) for python2.7 and 
the two modules will be included.


#### bedtools
bedtools is available in [bedtools website](http://bedtools.readthedocs.io/en/latest/).
It is recommended to use versions newer than v2.25.0. Older versions may work but they have not been tested.

  Please include the path to the bedtools execuitable in `PATH` environmental variable. If the bedtools
execuitable is at `/path/to/bedtools/`, you can add the below command to your `~/.bashrc` file (your shell config file).
```
export PATH=/path/to/bedtools/:$PATH
```


#### bigWigAverageOverBed
`bigWigAverageOverBed` execuitable can be downloaded from [UCSC genome browser utilities]( http://hgdownload.soe.ucsc.edu/admin/exe/), including binaries for
[linux](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigAverageOverBed) and
[macOS](http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/bigWigAverageOverBed). 

  Please include the path to the bigWigAverageOverBed execuitable in `PATH` environmental variable. If the 
bigWigAverageOverBed execuitable is at `/path/to/bedtools/`, you can add the below command to `~/.bashrc`
file (your shell config file).
```
export PATH=/path/to/bigWigAverageOverBed/:$PATH
```


## INSTALLATION
REPTILE can be downloaded and installed using the commands below. No root access is needed.
```
git clone https://github.com/yupenghe/REPTILE.git
cd REPTILE/
R CMD INSTALL REPTILE.tar.gz
```

## UNINSTALLATION
To uninstall REPTILE, please run the below command to remove the installed R package. You
may want to remove the REPTILE folder and related files to fully clean up.
```
R CMD UNINSTALL REPTILE
```

## Test REPTILE
Below command can be used To test whether REPTILE is correctly installed and all requirements
are met. 
```
cd test
./run_test.py
```


## Use REPTILE
#### Preprocessing
Preprocessing `./REPTILE_preprocess.py -h` to get help information.
```
./REPTILE_preprocess.py \
		data_info_file \
		training_region \
		-d dmr_file \
		tmp/training_region \
		-p 8
```

#### Training an enhancer model
```
./REPTILE_train.R \
	-i data_info_file \
	-a tmp/training_region.region_with_epimark.tsv \
	-d tmp/training_region.DMR_with_epimark.tsv \
	-l training_labels.tsv \
	-s mESC \
	-o tmp/REPTILE_model
```
`./REPTILE_train.R -h` to get help information.

#### Generate enhancer scores
```
./REPTILE_compute_score.R \
	-i data_info_file \
	-m tmp/REPTILE_model.reptile \
	-a tmp/test_region.region_with_epimark.tsv \
	-d tmp/test_region.DMR_with_epimark.tsv \
	-s E11_5_FB \
	-o results/E11_5_FB_pred
```
`./REPTILE_compute_score.R -h` to get help information.

#### Get genome-wide enhancer predictions
```
./REPTILE_call_enhancers.py \		
   	tmp/E11_5_FB_pred.R.bed \
	-d tmp/E11_5_FB_pred.DMR.bed \
   	-o results/enhancer_E11_5_FB.bed \
   	-p 0.5
```
`./REPTILE_call_enhancers.py -h` to get help information.
