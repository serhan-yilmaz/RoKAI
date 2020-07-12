![RoKAI App](rokai_app_logo.png "RoKAI App")

## Introduction
RoKAI is a computation tool for inferring kinase activities in a robust manner using functional networks. It is implemented in MATLAB and the source code is available in the [src folder](src/). If you are interested in performing RoKAI through a user-friendly online interface, please visit [RoKAI Web Application](https://rokai.ngrok.io/webapps/home/session.html?app=rokai).

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
In the provided examples, we utilize the prospho-proteomics data within [the Atlas compiled by Ochoa et al. (2016)](http://phosfate.com/download.html) along with four types of functional networks. For more information on the descriptions and format of the data, check out the [data preprocessing section](src/data_preprocessing/).

## License
This project is licensed under MIT - see the [LICENSE](LICENSE) file for details.

## References
Yılmaz S., Ayati M., Schlatzer D., Çiçek E., Chance M. R., Koyutürk M. (2020) [Robust Inference of Kinase Activity Using Functional Networks](https://www.biorxiv.org/content/10.1101/2020.05.01.062802v1). Preprint in bioRxiv.

Resources used for functional networks:
Hornbeck, P. V. et al. (2015). [Phosphositeplus, 2014: mutations, ptms and recalibrations](https://academic.oup.com/nar/article/43/D1/D512/2439467). Nucleic acids research, 43(D1), D512–D520
o Minguez, P. et al. (2012). [Ptmcode: a database of known and predicted functional associations between post-translational modifications in proteins](https://academic.oup.com/nar/article/41/D1/D306/1069950). Nucleic acids research, 41(D1), D306–D311 [Link]
o Szklarczyk, D. et al. (2014). [String v10: protein–protein interaction networks, integrated over the tree of life](https://academic.oup.com/nar/article/43/D1/D447/2435295). Nucleic acids research, 43(D1), D447–D452 [Link]

