## Data Preprocessing
Here, we prepare the functional networks and phosphorylation data utilized in [RoKAI demo](../../demo_rokai.m).
### Loading known phosphosites and kinase-substrate annotations
We obtain a list of known phosphosites and known kinase-substrate annotations from [PhosphositePlus (PSP)](https://www.phosphosite.org/staticDownloads). We use a [MATLAB script](load_psp_kinase_substrates.m) to parse PSP data and save as a MATLAB data (.mat) file. 
### Preparing kinase-kinase interaction network
For the kinase-kinase interactions, we use the protein-protein interaction (PPI) network obtained from [STRING](https://string-db.org/cgi/download.pl). The interactions are weighted by the 'combined_score' provided by STRING. See [STRING FAQ](http://version10.string-db.org/help/faq/) and [MATLAB script used to parse PPI network](load_string_ppi_network.m) for more information. 
### Preparing phosphosite level interaction network

### Combining 

