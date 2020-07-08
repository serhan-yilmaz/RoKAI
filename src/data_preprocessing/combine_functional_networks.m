%% Load networks
% Kinase substrate annotations obtained from PhosphositePlus. 
% See 'load_psp_kinase_substrates.m' for more information. 
dataFolder = '../../data/processed/';
PSP = load([dataFolder, 'psp_sites_and_kinase_substrate_network.mat']);

% Protein-Protein Interaction (PPI) network obtained from STRING.
% See 'load_string_ppi_network.m' for more information. 
STRING = load([dataFolder, 'string_ppi_network.mat']);

% Phosphosite-level interaction network based on co-evolution and/or 
% structure distance evidence obtained from PTMcode. 
% See 'load_ptmcode_networks.m' for more information. 
PTMcode = load([dataFolder, 'ptmcode_networks.mat']);
%% Mapping Ensembl protein identifiers to UniprotKb
% In order to create the kinase-kinase interaction network, we map the
% string protein identifiers (ENSP) to psp kinase identifiers (UniprotKB).
% We obtain the ensembl protein to uniprotkb id mappings using Uniprot
% mapping tool: https://www.uniprot.org/uploadlists/
% The 'string_proteins_uniprotkb.tab' contains the results of the query
% run on 2020-04-04, using the ids given in 'STRING.Proteins' as input. 
dataFolder = '../../data/raw/';
filePath = [dataFolder, 'string_proteins_uniprotkb.tab'];
ds = tabularTextDatastore(filePath, 'FileExtensions', '.tab');
ds.TextscanFormats = repmat({'%q'}, 1, length(ds.VariableNames));
Mapping = ds.readall();

% Filter for kinase identifiers
validRows = ismember(Mapping.Entry, PSP.Kinase.ID);
M = Mapping(validRows, :);

% Map to primary identifier
M.PrimaryEnsp = regexprep(M.Ensembl_Protein, ',[^.]*', '');

[~, kin_indices] = ismember(M.Entry, PSP.Kinase.ID);
[~, ppi_protein_indices] = ismember(M.PrimaryEnsp, STRING.Proteins);
Mkin2ppi = sparse(kin_indices, ppi_protein_indices, ...
        true, height(PSP.Kinase), length(STRING.Proteins));

% Create kinase-kinase interaction network, a subset of STRING ppi network
% Edges are weighted according to the combined score provided by STRING.
Kinase_PPI = Mkin2ppi * STRING.PPI * Mkin2ppi';
%% Mapping PTMcode Gene Identifiers to UniprotKb
% In order to create the phosphosite level interaction network, we map the
% PTMcode gene identifiers to psp protein identifiers (UniprotKB).
% We obtain the gene name and ensg to uniprotkb id mappings using 
% Uniprot mapping tool: https://www.uniprot.org/uploadlists/
% The 'ptmcode_genes_uniprotkb.tab' contains the results of the query
% run on 2020-04-05, using the ids given in 'PTMcode.Gene' as input. 
dataFolder = '../../data/raw/';
filePath = [dataFolder, 'ptmcode_genes_uniprotkb.tab'];
ds = tabularTextDatastore(filePath, 'FileExtensions', '.tab');
ds.TextscanFormats = repmat({'%q'}, 1, length(ds.VariableNames));
Mapping = ds.readall();

% Map to primary identifiers 
Mapping.Gene = regexprep(Mapping.Gene, ',[^.]*', '');

[b, indices] = ismember(Mapping.Entry, PSP.Protein);
M = Mapping(b, :);
Mmapping2pspprotein = sparse((1:height(M))', indices(b), ...
        1, height(M), length(PSP.Protein));

[~, indices] = ismember(M.Gene, regexprep(PTMcode.Gene, ' ', ''));
Mmapping2ptmcodegene = sparse((1:height(M))', indices, ...
        1, height(M), length(PTMcode.Gene));

% Mapping of PSP proteins to PTM Gene Identifiers
Mpspprotein2ptmcodegene = Mmapping2pspprotein' * Mmapping2ptmcodegene;

% Map the ptmcode sites to psp sites according to protein/gene ids
[~, pspProteinIndices] = ismember(PSP.Site.Protein, PSP.Protein);
Mpsp2pspprotein = sparse((1:height(PSP.Site))', ...
    pspProteinIndices, 1, height(PSP.Site), length(PSP.Protein));

Mptmcode2ptmcodegene = sparse((1:height(PTMcode.Site))', ...
    PTMcode.Site.GeneIndex, 1, height(PTMcode.Site), length(PTMcode.Gene));

% Mapping of PSP sites to PTMcode sites according to protein/gene ids
Mprotein_psp2ptmcode = Mpsp2pspprotein ...
                     * Mpspprotein2ptmcodegene ...
                     * Mptmcode2ptmcodegene';

% Map the ptmcode sites to psp sites according to residue
residues = [PSP.Site.Residue; PTMcode.Site.Residue];
[~, pspResidueIndices] = ismember(PSP.Site.Residue, residues);
[~, ptmcodeResidueIndices] = ismember(PTMcode.Site.Residue, residues);

Mpsp2residue = sparse((1:height(PSP.Site))', ...
    pspResidueIndices, 1, height(PSP.Site), length(residues));
Mptmcode2residue = sparse((1:height(PTMcode.Site))', ...
    ptmcodeResidueIndices, 1, height(PTMcode.Site), length(residues));

% Mapping of P  SP sites to PTMcode sites according to residue
Mresidue_psp2ptmcode = Mpsp2residue * Mptmcode2residue';

% Combined mapping of PSP sites to PTMcode sites
Mpsp2ptmcode = double(Mprotein_psp2ptmcode & Mresidue_psp2ptmcode);

% Map the co-evolution network to PSP sites
Wcoev = logical(Mpsp2ptmcode * PTMcode.Wcoev * Mpsp2ptmcode');

% Map the structure distance network to PSP sites
Wsd = logical(Mpsp2ptmcode * PTMcode.Wsd * Mpsp2ptmcode');

% Cleanup
clear Mmapping2pspprotein Mmapping2ptmcodegene Mpspprotein2ptmcodegene
clear Mpsp2pspprotein Mptmcode2ptmcodegene
clear Mpsp2residue Mptmcode2residue
clear Mprotein_psp2ptmcode Mresidue_psp2ptmcode

% Prepare output 
Kinase = table();
Kinase.KinaseID = PSP.Kinase.ID;
Kinase.KinaseName = PSP.Kinase.Name;
Kinase.Gene = PSP.Kinase.Gene;

position = regexprep(PSP.Site.Residue, '[^0-9]', '');

Site = table();
Site.Identifier = cellstr(join([PSP.Site.Protein, position], '_'));
Site.Protein = PSP.Site.Protein;
Site.Position = PSP.Site.Residue;
Site.Flanking = PSP.Site.Flanking;

NetworkData = struct();
NetworkData.Site = Site;
NetworkData.Kinase = Kinase;
NetworkData.Wkin2site = PSP.KS;
NetworkData.Wkin2kin = Kinase_PPI;
NetworkData.Wsite2site_coev = Wcoev;
NetworkData.Wsite2site_sd = Wsd;
NetworkData.nSite = height(Site);
NetworkData.nKinase = height(Kinase);

% Save the results
outputFolder =  '../../data/';
save([outputFolder, 'rokai_network_data_uniprotkb.mat'], 'NetworkData');







