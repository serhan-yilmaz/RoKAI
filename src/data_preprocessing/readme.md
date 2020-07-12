## Data Preprocessing
Here, we prepare the functional networks and the sample phospho-proteomics data utilized in [RoKAI demo](../../demo_rokai.m).
### 1 - Loading known phosphosites and kinase-substrate annotations
We obtain a list of known phosphosites and known kinase-substrate annotations from [PhosphositePlus (PSP)](https://www.phosphosite.org/staticDownloads). We use a [MATLAB script](load_psp_kinase_substrates.m) to parse PSP data and save as a compressed data (.mat) file for use in RoKAI.
### 2 - Loading protein-protein interaction network
For the kinase-kinase interactions, we use the protein-protein interaction (PPI) network obtained from [STRING](https://string-db.org/cgi/download.pl). The interactions are weighted by the 'combined_score' provided by STRING. See [STRING FAQ](http://version10.string-db.org/help/faq/) and [MATLAB script used to parse PPI network](load_string_ppi_network.m) for more information. 
### 3 - Preparing phosphosite level interaction network
For phosphosite-phosphosite interaction network, we use the co-evolution and structure distance evidence provided by [PTMcode](https://ptmcode.embl.de/data.cgi). For co-evolution network, we only use sites with relative residue conservation score (rRCS) greater than 90%. See [PTMcode help page](https://ptmcode.embl.de/help.cgi) and the [MATLAB script used to parse PTMcode data](load_ptmcode_networks.m) for more information.
### 4 - Combining (mapping) network data 
To obtain the network data used in RoKAI, we use a [MATLAB script](combine_functional_networks.m) to combine the kinase-substrate annotations from PhosphositePlus protein-protein interactions from STRING, and phosphosite level interactions from PTMcode. For mapping these datasets together, we employ several steps:
#### Mapping Ensembl protein identifiers to UniprotKb
In order to create the kinase-kinase interaction network, we map the Ensembl protein identifiers (ENSP) used by STRING to PSP kinase identifiers (UniprotKB). We obtain the Ensembl protein to UniprotKB id mappings using [Uniprot mapping tool](https://www.uniprot.org/uploadlists/). The [results of the query](../../data/string_proteins_uniprotkb.tab) are provided with the repository.
#### Mapping PTMcode Gene Identifiers to UniprotKb
In order to map the phosphosite level interaction network of PTMcode to the known phosphosites in PSP, we first map the gene identifiers used in PTMcode to PSP protein identifiers (UniprotKB). We obtain the gene symbol and ENSG to UniprotKB id mappings using [Uniprot mapping tool](https://www.uniprot.org/uploadlists/). The [results of the query](../../data/ptmcode_genes_uniprotkb.tab) are provided with the repository.
### 5 - Preparing sample phospho-proteomics data
For use in the demo, we utilize the prospho-proteomics data within [the Atlas compiled by Ochoa et al. (2016)](http://phosfate.com/download.html). Here, we include only a subset of the data, a single experiment of Pan et al. (2009) comparing EGF treated samples with the untreated. We provide [this sample data](../../data/sample_phospho_data_ensembl.csv) with the repository. In this data, the ENSP IDs are used as protein identifiers. To map these to UniprotKb IDs used by PSP, we use [Uniprot mapping tool](https://www.uniprot.org/uploadlists/). The [results of the query](../../data/phospho_data_proteins_uniprotkb.tab) are provided with the repository. See the [MATLAB script used to parse this data](load_sample_phospho_data.m) for more information.
### References
- Ochoa, David, et al. "An atlas of human kinase regulation." Molecular systems biology 12.12 (2016).
- Pan, Cuiping, et al. "Global effects of kinase inhibitors on signaling networks revealed by quantitative phosphoproteomics." Molecular & Cellular Proteomics 8.12 (2009): 2796-2808. 




