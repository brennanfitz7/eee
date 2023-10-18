# eee
ensemble epistasis engine

## Installation

In a terminal, type

```
git clone https://github.com/harmslab/eee.git 
cd eee
conda env create -f environment.yml
conda activate eee
python setup.py install
```

To use the *sync_structures* script/module, make sure you have the following
other packages installed and in the path. 

+ lovoalign: https://www.ime.unicamp.br/~martinez/lovoalign/software.html
+ muscle: https://www.drive5.com/muscle/
+ foldx: https://foldxsuite.crg.eu

## Usage

See the `examples` directory for example jupyter notebooks. 

### Requirements

+ numpy
+ pandas
+ matplotlib
+ tqdm
+ pytest
+ ete3
+ scikit-learn


### External software

+ muscle
+ foldx
+ lovoalign



