# UK Biobank resources

<p align="center">
    <img style="max-width: 40%; height: auto; " src="https://upload.wikimedia.org/wikipedia/en/d/d5/UK_biobank_logo.png"><br>
    <b>A curated list of UK Biobank related resources</b><br><br>
    <a href="https://awesome.re"><img src="https://awesome.re/badge.svg"></a>
    <a href="http://creativecommons.org/publicdomain/zero/1.0/"><img src="https://img.shields.io/badge/License-CC0_1.0-lightgrey.svg"></a><br><br>
    <img src="https://img.shields.io/badge/contributions-welcome-9cf.svg?style=flat">
    <img src="https://img.shields.io/github/last-commit/s-bell/awesome-uk-biobank?color=9cf">
</p>

## Contents
- [Official quicklinks](#official-quicklinks)
- [Phenomewide association studies](#phenomewide-association-studies)
- [Derivation of variables](#derivation-of-variables)
    - [Questionnaire](#questionnaire)
    - [Primary care](#primary-care)
    - [Clinical and biochemical](#clinical-and-biochemical)
    - [Accelerometer](#accelerometer)
- [Genetics](#genetics)
- [Imaging](#imaging)
    - [Brain](#brain-mri)
    - [Cardiac](#cardiac-mri)
    - [Abdominal](#abdominal-mri)
    - [Optical coherence tomography and fundus](#optical-coherence-tomography-and-fundus)
- [Data processing](#data-processing)
- [Misc](#misc)
- [Contributing](#contributing)

## Official quicklinks
- [Researcher login](https://bbams.ndph.ox.ac.uk/ams/)
- [Data showcase](https://biobank.ndph.ox.ac.uk/showcase/)
- [DNAnexus Research Analysis Platform](https://ukbiobank.dnanexus.com/login/) (RAP)
- [Twitter feed](https://twitter.com/uk_biobank)
- [Genomics mailing list](https://www.jiscmail.ac.uk/cgi-bin/webadmin?A0=UKB-GENETICS)
- [Neuroimaging mailing list](https://www.jiscmail.ac.uk/cgi-bin/wa-jisc.exe?A0=UKB-NEUROIMAGING)

## GWAS summary data / other resources
- [Neale lab (4,203 phenotypes)](https://www.nealelab.is/uk-biobank)
- [PanUK Biobank (7,228 phenotypes)](https://pan.ukbb.broadinstitute.org/)
- [GeneAtlas (778 phenotypes)](http://geneatlas.roslin.ed.ac.uk)
- [Global Biobank Engine (~3,800 phenotypes)](https://biobankengine.stanford.edu/)
- [UK Biobank TOPMed-imputed PheWeb (1,419 phenotypes)](https://pheweb.org/UKB-TOPMed/)
- [Lee lab (~1,400 phenotypes)](https://www.leelabsg.org/resources)
- [MRC IEU (2,514 phenotypes)](https://gwas.mrcieu.ac.uk/datasets/?gwas_id__icontains=ukb-b)
- [Genebass (4,529 phenotypes)](https://app.genebass.org)
- [AstraZeneca PheWAS portal (~15,500 binary and ~1,500 continuous phenotypes)](https://azphewas.com/)
- [Oxford Brain Imaging Genetics Server (~4,000 phenotypes)](https://open.win.ox.ac.uk/ukbiobank/big40/)
- [UNC School of Medicine GWAS summary statistics for brain imaging Phenotypes](https://www.med.unc.edu/bigs2/data/gwas-summary-statistics/)
- [Sex-stratified GWAS of 530 phenotypes](https://datashare.ed.ac.uk/handle/10283/3908)
- [ExPheWas browser (1,746 phenotypes)](https://exphewas.statgen.org/)
- [GWAS summary statistics for 82 heart imaging traits](https://zenodo.org/record/5606168)
- [Gene-based tests for 1,403 phecodes and 827 continuous phenotypes](https://zenodo.org/record/4477703)
- [Rivas lab (35 blood and urine biomarkers)](https://github.com/rivas-lab/biomarkers)
- [Nightingale Health NMR measures and their ratios in ~115,000 participants](https://gwas.mrcieu.ac.uk/datasets/?gwas_id__icontains=met-d)
- [Mitochondrial PheWas](https://zenodo.org/record/4609973)
- [29 blood cell traits](https://t.ly/Nnga)
- [Telomere length](https://figshare.com/s/caa99dc0f76d62990195)

Many of these, and more, are available via the fantastic [MRC IEU OpenGWAS](https://gwas.mrcieu.ac.uk/datasets/) project.

## Perform phenome-wide association studies
| Tool | Description|
|------|------------|
| [PHESANT](https://github.com/MRCIEU/PHESANT) | PHEnome Scan ANalysis Tool (pheWAS, Mendelian randomisation (MR)-pheWAS etc) in UK Biobank |
| [createUKBphenome](https://github.com/umich-cphds/createUKBphenome) | Create a PheWAS code based phenome using ICD9 and ICD10 data from baskets of the UK Biobank |

## Derivation of variables

### Questionnaire
| Tool | Description|
|------|------------|
| [ukbb-srmed](https://github.com/PhilAppleby/ukbb-srmed) | Match terms in the UK Biobank self-reported medication data coding table with terms in both the Anatomical Therapeutic Chemical (ATC) classification system and in the British National Formulary (BNF) |
| [ukbfrailty](https://github.com/dylwil/ukbfrailty) | Stata code to generate a frailty index based on 49 items in UK Biobank |
| [ukbiobank_lca_model_predictions](https://github.com/dcmuller/ukbiobank_lca_model_predictions) | Stata programs to calculate risk of lung cancer based on a UK Biobank prediction model |

### Electronic health records
| Tool | Description|
|------|------------|
| [ukb-biomarker-phenotypes](https://github.com/spiros/ukb-biomarker-phenotypes) | Phenotyping algorithms for 31 commonly-measured biomakers in primary care   electronic health records |
| [ukbb-ehr-data](https://github.com/philipdarke/ukbb-ehr-data) | Clean and prepare UK Biobank primary care electronic health records for research |
| [ukbbhelpr](https://github.com/philipdarke/ukbbhelpr) | A collection of helper functions primarily for working with linked primary care electronic health record data in UK Biobank |
| [phemap](https://github.com/spiros/phemap) | Functions to map between ICD-10 terms and PheCodes for UK Biobank hospital electronic health records |
| [UKBCaseFinder](https://github.com/yeyixuan/UKBCaseFinder) | A tool for identifying patients in UK Biobank given the definition of disease phenotypes (ICD-10, ICD-9, OPCS or cancer histology) |
| [ukbcodes](https://github.com/SurgicalInformatics/ukbcodes) | Returns suggestions of equivalent codes based on UK Biobank mappings |
| [UKB.COVID19](https://github.com/bahlolab/UKB.COVID19) | UK Biobank COVID-19 data processing and risk factor association tests |

### Clinical and biochemical
| Tool | Description |
|---|---|
| [ukbnmr](https://github.com/sritchie73/ukbnmr) | Tools for processing Nightingale NMR biomarker data in UK Biobank |
| [CAHDbyECG](https://github.com/andywxf/Predict_CAHD_by_ECG_on_UK_Biobank) | Predict coronary atherosclerotic heart disease from electrocardiogram |
| [NeuralCVD](https://github.com/thbuerg/NeuralCVD) | Contains code to preprocess UK Biobank data and train the NeuralCVD score (neural network-based integration of polygenic and clinical information of a prediction model for 10 year risk of major adverse cardiac events) | 

### Accelerometer
| Tool | Description |
|------|-------------|
| [biobankAccelerometerAnalysis](https://github.com/OxWearables/biobankAccelerometerAnalysis) | A tool to extract meaningful health information from large accelerometer datasets |
| [BiobankActivityCSF](https://github.com/CASSON-LAB/BiobankActivityCSF) | Process the accelerometer CWA files from the UK Biobank activity monitoring study in Python |
| [snake-ukb-accel](https://github.com/tgbrooks/ukbb) | Snakemake pipeline to process and aggregate the actigraphy data in the UK Biobank |
| [axivity-ax3-tool](https://github.com/jlc-christie/axivity-ax3-tool) | Command line tool for easily extracting light and temperature data from Axivity AX3 accelerometers used in the UK Biobank |
| [UKBiobankKineticEnergyHarvesting](https://github.com/CASSON-LAB/UKBiobankKineticEnergyHarvesting) | Software for processing the UK Biobank accelerometer datasets in the paper 'Estimation of kinetic energy harvesting potential for self-powered wearable devices' | 

## Genetics
| Tool | Description|
|------|------------|
| [MoChA](https://github.com/freeseek/mocha) | MOsaic CHromosomal Alterations (MoChA) caller |
| [MT_UKB](https://github.com/GrassmannLab/MT_UKB) | Compute mtDNA abundance estimates |
| [mt_PheWAS](https://github.com/clody23/UKBiobank_mtPheWas/) | Custom scripts used for mtSNVs variant recalling in UK Biobank |
| [ddpop](https://github.com/marcustuke/ddpop) | Detect large DNA duplications/deletions in UK Biobank Affymetrix Axiom array data |
| [cnv-ukb](https://github.com/priestlab/cnv-ukb) | Contains scripts used to call CNVs using UK Biobank data |
| [ukbb-pytorch-dataloader](https://github.com/hrushikeshloya/ukbb-pytorch-dataloader) | Load genotype data from UK Biobank into PyTorch |
| [go-bgen](https://github.com/carbocation/bgen) | BGEN file format reader for parsing UK Biobank and other genotype data in Go |
| [liftover_plink_beds](https://github.com/dnanexus-rnd/liftover_plink_beds) | Converting UK Biobank genome-wide genotyping data from one reference build to another |

## Imaging phenotypes

### Brain MRI
| Tool | Description |
|------|-------------|
| [ukb-brain-imaging](https://www.fmrib.ox.ac.uk/ukbiobank/fbp) | UK Biobank brain imaging processing pipeline |
| [tvb-ukbb](https://github.com/McIntosh-Lab/tvb-ukbb) | TheVirtualBrain implementation of the UK Biobank imaging pipeline |
| [qsm-pipeline](https://git.fmrib.ox.ac.uk/cwang/uk_biobank_qsm_pipeline) | MATLAB scripts for QSM processing of the UK Biobank swMRI data |
| [XTRACT_atlases](https://github.com/SPMIC-UoN/XTRACT_atlases) | Extract 42 white matter tract atlases |
| [ukbiobank-spinalcord-csa](https://github.com/sct-pipeline/ukbiobank-spinalcord-csa) | Analysis pipeline to normalise spinal cord cross-sectional area |
| [brain-age](https://github.com/ha-ha-ha-han/UKBiobank_deep_pretrain) | Pretrained neural networks for brain age prediction in UK Biobank from MRI images |
| [slice2vol](https://github.com/alexbagur/slice2vol) | Slice-to-Volume Registration for UK Biobank data |
| [brain-sex](https://github.com/LeoBman/brain-sex-classification) | Predict participant sex from brain MRI images |
| [PSMD](http://www.psmd-marker.com/) | PSMD is a robust, fully-automated and easy-to-implement marker for cerebral small vessel disease based on diffusion tensor imaging, white matter tract skeletonization (as implemented in FSL-TBSS) and histogram analysis |
| [ePVS quantification](https://hub.docker.com/r/giorgioberardini/pvs5) | Docker image to quantify enlarged perivascular spaces (n.b. [preprocessing](https://hub.docker.com/r/giorgioberardini/pvs_preproc) should be performed first) |
| [morpho-deepsulci](https://github.com/brainvisa/morpho-deepsulci) | Deep learning methods for Morphologist sulci recognition |
| [ENIGMA-SULCI](https://hub.docker.com/r/fpizzaga/sulci) | Docker image for the sulci protocol from the ENIGMA-SULCI working group that combines Freesurfer and BrainVISA to robuslty segment and label 123 sulci across the whole brain |
| [adni_phenotypes](https://github.com/tjiagoM/adni_phenotypes) | Identifying healthy individuals with Alzheimer neuroimaging phenotypes in the UK Biobank |

### Cardiac MRI
| Tool | Description |
|------|-------------|
| [ukbb_cardiac](https://github.com/baiwenjia/ukbb_cardiac) | A toolbox used for processing and analysing cardiovascular magnetic resonance (CMR) images |
| [ukbb_lax_4ch](https://github.com/priestlab/ukbb_lax_4ch_segmentation) | Segmenting and measuring long-axis 4-chamber view cardiac MRI data sets |
| [CardioRadiomics](https://github.com/iremcetin/radiomics_cardio_risk_factors) | Radiomics signatures of cardiovascular risk factors from cardiac MRI |
| [ukb-cardiac-mri](https://github.com/HazyResearch/ukb-cardiac-mri) | Weakly supervised classification of aortic valve malformations using unlabeled cardiac MRI sequences |
| [cardiac-fractals](https://zenodo.org/record/3698268) | Artificial intelligence derived fractal structure of the heart |

### Abdominal MRI
| Tool | Description |
|------|-------------|
| [ReCOH_pipeline](https://github.com/recoh/pipeline) | Image processing and quality control for abdominal MRI in the UK Biobank |
| [ukbb-mri-sseg](https://github.com/calico/ukbb-mri-sseg/) | Semantic segmentation for Dixon abdominal MRI and other modalities by Calico |
| [ukb-kidney-segmentation](https://github.com/tarolangner/ukb_segmentation) | Kidney morphology |https://github.com/alexbagur/pancreas-segmentation-ukbb |
| [UKBiobankDXAMRIPreprocessing](https://github.com/rwindsor1/UKBiobankDXAMRIPreprocessing) | Code for downloading and preprocessing (and stitching) MRI and DXA data from the UK Biobank | 
| [AAC_scoring](https://github.com/calico/AAC_scoring) | Map DEXA images from UK Biobank dataset to abdominal aortic calcification scores  |
| [ukb_mimir](https://github.com/tarolangner/ukb_mimir) | An inference engine for UK Biobank neck-to-knee MRI |
| [ukb_segmentation](https://github.com/tarolangner/ukb_segmentation) | PyTorch implementation for 2.5D U-Net segmentation of UK Biobank neck-to-knee body MRI |
| [iliopsoas_muscle](https://github.com/recoh/iliopsoas_muscle) | Large-scale analysis of iliopsoas muscle volumes in the UK Biobank |
| [PDFF](https://github.com/marcsous/pdff) | Proton density fat fraction calculation for MRI |

### Optical coherence tomography and fundus
| Tool | Description |
|------|-------------|
| [OCT-Converter](https://github.com/marksgraham/OCT-Converter) | Tools for extracting the raw optical coherence tomography and fundus data from proprietary file formats |
| [VAMPIRE](https://vampire.computing.dundee.ac.uk/) | VAMPIRE (Vessel Assessment and Measurement Platform for Images of the REtina) is a software application for efficient, semi-automatic quantification of retinal vessel properties with large collections of fundus camera images |
| [Retina-Seg](https://github.com/vineet1992/Retina-Seg) | Automatically segment vasculature from retinal fundus photographs |

## Data processing
| Tool | Description |
|------|-------------|
| [ukbtools](https://github.com/kenhanscombe/ukbtools) | An R package to manipulate and explore UK Biobank data |
| [ukbwranglr](https://github.com/rmgpanw/ukbwranglr) | An R package for UK Biobank data wrangling |
| [ukbwranglrextra](https://github.com/rmgpanw/ukbwranglrextra) | Extra add-on functions for ukbwranglr |
| [ukbwranglr_resources](https://github.com/rmgpanw/ukbwranglr_resources) | A collection of resources to support ukbwranglr |
| [ukbREST](https://github.com/hakyimlab/ukbrest) | Efficient and streamlined data access for reproducible research of large biobanks |
| [ukbcc](https://github.com/tool-bin/ukbcc) | Tool for curation of UK Biobank data to generate cohorts |
| [ukbm](https://github.com/neurodatascience/ukbm) | Code for management of UK Biobank bulk data |
| [bulk-ukb](https://github.com/adigherman/UKB) | Download UK Biobank bulk data  |
| [UKBioPick](https://github.com/CirculatoryHealth/UKBioPick) | Extract fields of interest from a UK Biobank tab-delimited phenotypic file (ukbxxxx.tab) |
| [BiobankRead-Bash](https://github.com/saphir746/BiobankRead-Bash) | Python scripts to extract and pre-process UK Biobank data |
| [precimed-ukb](https://github.com/precimed/ukb) | Helper tools to process UK Biobank data |
| [UKBFetchbash](https://github.com/GordonMatthewson/UKBFetchbash) | Bash script to download your application's approved bulk data files from the UK Biobank dataset |
| [ukbb-slicer](https://github.com/webhash/ukbb-slicer) | Slice large UK Biobank CSV files in a memory efficient way |
| [R-ukbiobank](https://github.com/adamleejohnson/R-ukbiobank) | Contains a collection of functions that facilitate the extraction of phenotypes from UK Biobank data |
| [ukbb_parser](https://github.com/nadavbra/ukbb_parser) | A Python module for loading phenotypic and genetic data from the UK Biobank |
| [ukb_download_and_prep_template](https://github.com/OxWearables/ukb_download_and_prep_template) | Facilitates common data pre-processing steps when working with UK Biobank data |
| [UKBiobank-Table-Conversion](https://github.com/futurologist/UKBiobank-Table-Conversion) | Conversion of data tables extracted from UK Biobank into new regression-friendly data tables  |
| [ukbpheno](https://github.com/niekverw/ukbpheno) | An R package for efficiently munging the files provided by UK Biobank to generate data tables of with unified format for further analysis |
| [ukb_healthoutcomes_db](https://github.com/ccbs-stradl/ukb_healthoutcomes_db) | Store and work with UK Biobank record-level health outcomes in a SQLite database |

## Misc
| Tool | Description |
|------|-------------|
| [UK-Biobank-scraper](https://github.com/lwaw/UK-Biobank-scraper) | A script that will download meta-data from the UK Biobank website and return a json object |
| [UKBscrape](https://github.com/asoroosh/UKBscrape) | Scrapes the UK Biobank website to provide a built-in dictionary for the field IDs. Also checks which variables are not included in your application compared to the website |
| [ukbschemas](https://github.com/bjcairns/ukbschemas) | Use R to generate a database containing the UK Biobank data [schemas](http://biobank.ctsu.ox.ac.uk/crystal/schema.cgi) |
| [docker-ukbiobank-utils](https://github.com/spiros/docker-ukbiobank-utils) | Docker versions of the utilities provided by UK Biobank |
| [tofu](https://github.com/spiros/tofu) | Tofu is a Python library for generating synthetic UK Biobank data |
| [BBS](https://github.com/choishingwan/BBS) | Simulation software for generating phenotypes from UK Biobank data |

## Contributing
Please submit a pull request or create an issue to add a new resource to the list - software, packages, tutorials or similar welcomed. 

[Return to top](#uk-biobank-resources)
