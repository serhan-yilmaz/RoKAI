![RoKAI App](rokai_app_logo.png "RoKAI App")

## Introduction
RoKAI is a computation tool for inferring kinase activities in a robust manner using functional networks. It is implemented in MATLAB and the source code is available in the [src folder](src/). If you are interested in performing RoKAI through a user-friendly online interface, please visit [RoKAI Web Application](https://rokai.ngrok.io/webapps/home/session.html?app=rokai).

RoKAI operates on a heterogeneous network having kinases and phosphosites as nodes and available functional associations as edges, including protein-protein interactions, kinase-substrate annotations, co-evolution and structure distance evidence between phosphosites. The key idea of RoKAI is to propagate the phosphosite quantifications on this heterogeneous network to capture the coordinated changes in the signaling, which are used to infer the kinase activities in a more robust manner.

## Getting Started
This section is under construction. In the meanwhile, you can check out the [user manual](rokai_user_manual.pdf). 
For the data used in [the demo](demo_rokai.m), check out the [data preprocessing section](src/data_preprocessing/).
