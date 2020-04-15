## Data Preprocessing
Here, we prepare the functional networks and phosphorylation data utilized in [RoKAI demo](../../demo_rokai.m).
### 1 - Loading known phosphosites and kinase-substrate annotations
We obtain a list of known phosphosites and known kinase-substrate annotations from [PhosphositePlus (PSP)](https://www.phosphosite.org/staticDownloads). We use a [MATLAB script](load_psp_kinase_substrates.m) to parse PSP data and save as a compressed data (.mat) file for use in RoKAI.
### 2 - Loading protein-protein interaction network
For the kinase-kinase interactions, we use the protein-protein interaction (PPI) network obtained from [STRING](https://string-db.org/cgi/download.pl). The interactions are weighted by the 'combined_score' provided by STRING. See [STRING FAQ](http://version10.string-db.org/help/faq/) and [MATLAB script used to parse PPI network](load_string_ppi_network.m) for more information. 
### 3 - Preparing phosphosite level interaction network
For phosphosite-phosphosite interaction network, we use the co-evolution and structure distance evidence provided by [PTMcode](https://ptmcode.embl.de/data.cgi). For co-evolution network, we only use sites with relative residue conservation score (rRCS) greater than 90%. See [PTMcode help page](https://ptmcode.embl.de/help.cgi) and the [MATLAB script used to parse PTMcode data](load_ptmcode_networks.m) for more information.
### 4 - Combining (mapping) network data 
To obtain the network data used in RoKAI, we use a [MATLAB script](combine_functional_networks.m) to combine the kinase-substrate annotations from PhosphositePlus protein-protein interactions from STRING, and phosphosite level interactions from PTMcode. For mapping these datasets together, we employ several steps:
##### Mapping Ensembl protein identifiers to UniprotKb
In order to create the kinase-kinase interaction network, we map the string protein identifiers (ENSP) to psp kinase identifiers (UniprotKB). We obtain the ensembl protein to uniprotkb id mappings using [Uniprot mapping tool](https://www.uniprot.org/uploadlists/). The [results of the query (run on 2020-04-04)](../../data/string_proteins_uniprotkb.tab) are provided with the repository.. 

### 5 - Preparing sample phospho-proteomics data













