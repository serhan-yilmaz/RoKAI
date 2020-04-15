## Data Preprocessing
Here, we prepare the functional networks and phosphorylation data utilized in [RoKAI demo](../../demo_rokai.m).
### 1 - Loading known phosphosites and kinase-substrate annotations
We obtain a list of known phosphosites and known kinase-substrate annotations from [PhosphositePlus (PSP)](https://www.phosphosite.org/staticDownloads). We use a [MATLAB script](load_psp_kinase_substrates.m) to parse PSP data and save as a compressed data (.mat) file for use in RoKAI.
### 2 - Preparing kinase-kinase interaction network
For the kinase-kinase interactions, we use the protein-protein interaction (PPI) network obtained from [STRING](https://string-db.org/cgi/download.pl). The interactions are weighted by the 'combined_score' provided by STRING. See [STRING FAQ](http://version10.string-db.org/help/faq/) and [MATLAB script used to parse PPI network](load_string_ppi_network.m) for more information. 
### 3 - Preparing phosphosite level interaction network
For phosphosite-phosphosite interaction network, we use the co-evolution and structure distance evidence provided by [PTMcode](https://ptmcode.embl.de/data.cgi). For co-evolution network, we only use sites with relative residue conservation score (rRCS) greater than 90%. See [PTMcode help page](https://ptmcode.embl.de/help.cgi) and the [MATLAB script used to parse PTMcode data](load_ptmcode_networks.m) for more information.
### 4 - Combining (mapping) network data 

### 5 - Preparing sample phospho-proteomics data

