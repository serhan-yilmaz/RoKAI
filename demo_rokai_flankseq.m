%% Running RoKAI with Flanking Sequence Identifiers
% Here, we will run RoKAI using input data with flanking sequence
% identifiers for phosphosites.
% This time our sample input has only two columns:
% - SiteFlanking: The +-7 flanking sequence of the site
% - Quantification: The phosphorylation of the site as log2 fold change
filePath = [folder, 'sample_phospho_data_flanking.csv'];
ds = datastore(filePath);
ds.TextscanFormats = {'%q', '%f'};
T = readall(ds);

% Load network data prepared for Uniprotkb
load([folder, 'rokai_network_data_uniprotkb.mat']);

% Set the site identifier to flanking sequence
NetworkData.Site.Identifier = NetworkData.Site.Flanking;

% We can RoKAI by specifying the identifier type to be 'flanking'
[KinaseTable, SiteTable] = rokai(T, NetworkData, 'Identifier', 'flanking');

% Add the source files to search path
addpath('src/rokai');
%% Running RoKAI with default options
% The first output of RoKAI is a table containing the inferred activity of
% kinases. Each row is a kinase and there are total of four columns:
% (1) KinaseID: The UniProt accession id of the kinase
% (2) KinaseName: The name/description of the kinase
% (3) Gene: The gene symbol of the kinase
% (4) Activity: The inferred activity of the kinase based on the
%  phosphorylation of sites in their functional neighborhood.
% Note that, the table is sorted by kinase activites in descending order.
[KinaseTable] = rokai(T, NetworkData, 'Identifier', 'flanking');

% The second output of RoKAI is a table containing the refined
% phosphorylation profile after network propagation.
% In addition to site identifiers, it contains the following columns:
% - Raw_Q: The raw phosphorylation (log2 fold change) of the site. 
% These are equal to the input quantifications without any modification.
% - RoKAI_Q: The refined phosphorylation of the site after network
%  propagation. These are an aggregate of the phosphorylation of the sites
%  themselves as well as the sites in their functional neighborhood.
% Note that, the table is sorted by the absolute values of the refined
% profile obtained by RoKAI. Thus, the sites on the top exhibit higher 
% difference with the control sample.
[KinaseTable, SiteTable] = rokai(T, NetworkData, 'Identifier', 'flanking');
%% Running RoKAI with customized options
options = struct();

% The inference method to be used. The available options are:
% - Mean: The mean phosphorylation of kinase substrates
% - Zscore: The mean phosphorylation normalized for statistical significance
% - Linear: A linear model where the phosphorylation of a site is modeled 
% as the summation of activities of kinases that phosphorylate that site
% - GSEA: Gene set enrichment analysis adapted for kinase activity inference
options.InferenceMethod = 'linear';

% Whether to consider protein-proteins interactions between kinases
options.IncludePPI = false;

% Whether to consider structure distance in phosphosite interactions
options.IncludeStructureDistance = true;

% Whether to consider coevolution evidence in phosphosite interactions
options.IncludeCoevolution = false;

% Whether to keep missing sites without quantifications in the network
options.IncludeMissingSites = true;

% Set the input data identifier
options.Identifier = 'Flanking';

[KinaseTableb, SiteTableb] = rokai(T, NetworkData, options);








