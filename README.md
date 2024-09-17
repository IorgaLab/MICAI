# MICAI Antibioresistance Predictor

This project is part of [Seq2Diag](https://anr.fr/ProjetIA-20-PAMR-0010), which plans to use AI methods to make in-silico determinations of the antibioresistance phenotypes for _Escherichia coli_, _Klebsiella pneumoniae_ and _Pseudomonas aeruginosa_. \
This repository is the extension of our [web-service](http://iorgalab.org:4567/micai). \
It enables the user to :

- Run the exact same predictions than the web-service locally
- Use other subsequence sizes
- Make the predictions of all the files of a folder, and get their results in a csv file.

# Dependencies

Our tool requires :
- numpy
- pandas
- xgboost==2.0.3
- scikit-learn==1.3.0
- biopython

# Installation

First, let's clone the repository into your computer.
```
git clone https://github.com/IorgaLab/MICAI.git
cd MICAI
```

Then, let's create a Conda environment. \
If you did not install it yet, have a look on this  [link](https://conda.io/projects/conda/en/latest/user-guide/index.html). \
If you wish to not use one, you can directly install the dependencies.

```
conda env create -p env_predictor --file environment.yml
conda activate env_predictor
```
Each time you want to use the environment, you will need this last command to activate it. \
To deactivate it, you simply need to use the command "conda deactivate".\
Finally, let's decompress the data and models. It weights about 250Mo.

```
tar -xjf data.tar.bz2
rm data.tar.bz2
```

# Configuration of use

Let's get into details of the different options to choose from.

## The species (mandatory)

The prediction can be made on 3 different species :

- _Escherichia coli_ : Ec
- _Klebsiella pneumoniae_ : Kp
- _Pseudomonas aeruginosa_ : Pa

The species can be chosen using the option "--species [code]", with [code] being the 2-letters of the species. \
Please note that a wrongly selected species will not be detected, and the tool will make non-sense predictions.

## The mode (mandatory)

Our models are trained on 2 task :

- Resistance or Susceptibility : binary
- Prediction of the MIC : regression

The mode can be chose using the option "--mode [code]", with [code] being either "binary" or "regression".

Please note that due to the lack of breakpoint for the binary task, or a data distribution not large enough in our dataset for the regression task, some antibiotics might appear in one mode and not in the other. \
Furthermore, we might lack of training data in some range of the MIC predictions. Consequently, we grouped together neighbouring MIC zones : "From x to y".

## The encoding (mandatory)

Our model work on different encoding : 

- Nucleotide format (assembly file) : nucl
- Protein format (Prodigal .faa file) : prot

The encoding can be chose using the option "--encoding [code]", with [code] being either "nucl" or "prot". \
The default format is fasta. It is possible to change it using "--format [format]" option.

## The subsequence size (optional)

In our web-service, the subsequence size can not be chosen as they show similar perfomance. \
In this repertory, the same sizes will be chosen by default, but the user can choose to use an other one.

- For the nucleotide format, the default size is 24, and can be chosen : 18, 24, 33, 42, 51
- For the protein format, the default size is 8, and can be chosen : 6, 8, 11, 14, 17 

To do so, use the option "--size [code]", with [code] being one of the options above.

# Use of the single_predictor

The minimum command for the prediction of a single file is :

```
python single_predictor.py --file [file] --species [species] --mode [mode] --encoding [encoding]
```

The options --species, --mode and --encoding are described above. The options --size and --format can aditionnaly be used. \
The option --file corresponds to the path of the file you want to analyze. \
The predictions will be printed on your screen.

# Use of the group_predictor

This code is made for the parallel analysis of multiple files (using multiprocessing). It uses the same options than the single_predictor, except :

- Instead of --file, you must use --path, with the path of a folder containing all (and only) the files you want to analyze.
- You must add the "--output [output]", which is the name of the .csv file in which all the predictions will be stored.

```
python group_predictor.py --path [path] --species [species] --mode [mode] --encoding [encoding] --output [output]
```

# Example of use
```
python single_predictor.py --file example_sample_Ec.faa --species Ec --mode binary --encoding prot
python group_predictor.py --path example_folder_assembly_Kp --species Kp --mode regression --encoding nucl --size 42 --output tab_Kp.csv
```
