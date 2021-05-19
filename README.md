![RoKAI App](rokai_app_logo.png "RoKAI App")

## Introduction
RoKAI is a computation tool for inferring kinase activities in a robust manner using functional networks. It is implemented in MATLAB and the source code is available in the [src folder](src/). If you are interested in performing RoKAI through a user-friendly online interface, please visit [RoKAI Web Application](https://rokai.io).

RoKAI operates on a heterogeneous network having kinases and phosphosites as nodes and available functional associations as edges, including protein-protein interactions, kinase-substrate annotations, co-evolution and structure distance evidence between phosphosites. The key idea of RoKAI is to propagate the phosphosite quantifications on this heterogeneous network to capture the coordinated changes in the signaling, which are used to infer the kinase activities in a more robust manner. For more information, check out the [workflow of RoKAI](README_ROKAI.md).

## Getting Started
These instructions will guide you to run RoKAI locally using the source code. For information on how to use the RoKAI web application, you can check out the provided [user manual](rokai_user_manual.pdf). 

### Installation Requirements
In order to run RoKAI locally, a MATLAB environment with Statistics and Machine Learning Toolbox is required. If you do not have one already, you can install MATLAB on your computer by following [the installation steps](https://www.mathworks.com/help/install/).

### Examples
We provide a few examples on how to run RoKAI on MATLAB. Simply run the demo file:
```
demo_rokai.m
```
In the provided examples, we utilize the prospho-proteomics data within [the Atlas compiled by Ochoa et al. (2016)](http://phosfate.com/download.html) as sample data. For more information on the description and format of the data, check out the [data preprocessing section](src/data_preprocessing/).

### Input data format
To run RoKAI on your data, all you need is a single csv file having phosphosites identifiers and phosphorylation value (as log2 fold change). Overall, RoKAI desktop application accepts three types of phosphosite identifiers:

#### Ensembl protein identifiers and position
The [demo_rokai.m](demo_rokai.m) script (as well as the web application) uses this type of identifiers. In this case, you need an input csv file that has three columns:
- Protein: The Ensembl protein (ENSP) identifier
- Position: The position of the site on the protein
- Quantification: The phosphorylation of the site as log2 fold change

For more information, check our [sample data with Ensembl identifiers](data/sample_phospho_data_ensembl.csv).

#### UniprotKB protein identifiers and position
Alternatively, you can use UniprotKB protein identifiers and position information to specify the phosphosites. See the [sample data](data/sample_phospho_data_uniprotkb.csv) and [demo_rokai_uniprotkb.m](demo_rokai_uniprotkb.m) script for an example on how to run RoKAI with this type of input data.

#### Flanking sequence 
Finally, you can use the +-7 flanking sequence as phosphosite identifiers. These will be mapped to the reference phosphosites provided by [PhosphositePlus](https://www.phosphosite.org/staticDownloads). See the script [demo_rokai_flankseq.m](demo_rokai_flankseq.m) and [the sample data with this type of identifiers](data/sample_phospho_data_flanking.csv)

### For more advanced use cases
The easiest way to run RoKAI is through the web application (at http://rokai.io) or by modifing the provided demo files to load your input files instead of the sample files. For more advanced cases that require an adaptation of source code (such as running RoKAI on networks other than the ones provided), you can refer to the provided [code schematic](rokai_source_code_schematic.pdf) that explains the code structure used. 

#### Replicating the results of our study
We provide all materials (data and code) to reproduce the results of our study in [Figshare](https://doi.org/10.6084/m9.figshare.12644864). These resources might be helpful for incorporating RoKAI into your analysis pipeline readily.

## References
For more information on our study:
- Yılmaz S., Ayati M., Schlatzer D., Çiçek A. E., Chance M. R., Koyutürk M. (2021) [Robust Inference of Kinase Activity Using Functional Networks](https://doi.org/10.1038/s41467-021-21211-6). Nature Communications, 12 (1117)

Resources used for functional networks:
- Hornbeck, P. V. et al. (2015). [Phosphositeplus, 2014: mutations, ptms and recalibrations](https://academic.oup.com/nar/article/43/D1/D512/2439467). Nucleic acids research, 43(D1), D512–D520
- Minguez, P. et al. (2012). [Ptmcode: a database of known and predicted functional associations between post-translational modifications in proteins](https://academic.oup.com/nar/article/41/D1/D306/1069950). Nucleic acids research, 41(D1), D306–D311
- Szklarczyk, D. et al. (2014). [String v10: protein–protein interaction networks, integrated over the tree of life](https://academic.oup.com/nar/article/43/D1/D447/2435295). Nucleic acids research, 43(D1), D447–D452

