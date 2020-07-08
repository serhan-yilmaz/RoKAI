% We obtain co-evolution and structure distance associations from PTMcode. 
% The provided phosphosite associations are organized in two files:
% (1) Within-Proteins: 'PTMcode2_associations_within_proteins.txt'
% (2) Between-Proteins: 'PTMcode2_associations_between_proteins.txt'
% Due to their excessive file size, they are not included in the
% repository. But, they are available at:   
% https://ptmcode.embl.de/data.cgi
%% Load PTMcode Within-Protein Associations
folder = '../../data/big/';
filename = [folder, 'PTMcode2_associations_within_proteins.txt'];
ds = tabularTextDatastore(filename, 'delimiter', '\t');
ds.VariableNames{1} = 'Protein';
T1 = ds.readall();

% Filter for humans
validRows = strcmpi(T1.Species, 'homo sapiens');
T1 = T1(validRows, :);

% Filter for phosphorylation
validRows = strcmpi(T1.PTM1, 'phosphorylation');
validRows = validRows & strcmpi(T1.PTM2, 'phosphorylation');
T1 = T1(validRows, :);

% Remove redundant fields
T1.PTM1 = [];
T1.PTM2 = [];
T1.Species = [];
T1.same_residue_competition__evidence = [];
%% Load Between-Protein Associations
folder = '../../data/big/';
filename = [folder, 'PTMcode2_associations_between_proteins.txt'];
ds = tabularTextDatastore(filename, 'delimiter', '\t');
ds.VariableNames{1} = 'Protein1';
T2 = ds.readall();

% Filter for humans
validRows = strcmpi(T2.Species, 'homo sapiens');
T2 = T2(validRows, :);

% Filter for phosphorylation
validRows = strcmpi(T2.PTM1, 'phosphorylation');
validRows = validRows & strcmpi(T2.PTM2, 'phosphorylation');
T2 = T2(validRows, :);

% Remove redundant fields
T2.PTM1 = [];
T2.PTM2 = [];
T2.Species = [];
%% Combine Within-Protein and Between-Protein Associations
[Proteins] = unique([T1.Protein; T2.Protein1; T2.Protein2]);

[~, T1.ProteinIndex1] = ismember(T1.Protein, Proteins);
T1.ProteinIndex2 = T1.ProteinIndex1;
T1.Protein = [];

[~, T2.ProteinIndex1] = ismember(T2.Protein1, Proteins);
[~, T2.ProteinIndex2] = ismember(T2.Protein1, Proteins);
T2.Protein1 = [];
T2.Protein2 = [];

Tc = [T1; T2];
clear T1 T2
%% Create networks
Tc.SiteID1 = cellstr(join([Proteins(Tc.ProteinIndex1), Tc.Residue1], '_'));
Tc.SiteID2 = cellstr(join([Proteins(Tc.ProteinIndex2), Tc.Residue2], '_'));

[SiteIDs] = unique([Tc.SiteID1; Tc.SiteID2]);
[~, Tc.SiteIndex1] = ismember(Tc.SiteID1, SiteIDs);
[~, Tc.SiteIndex2] = ismember(Tc.SiteID2, SiteIDs);
Tc.SiteID1 = [];
Tc.SiteID2 = [];

% Get associations with co-evolution evidence
% Filter for rows with rRCS1 & rRCS2 >= 0.9
validRows = Tc.Coevolution_evidence == true;
validRows = validRows & (Tc.rRCS1 >= 90) & (Tc.rRCS2 >= 90);
Tcoev = Tc(validRows, :);

% Create Co-Evolution Network
Wcoev = logical(sparse(Tcoev.SiteIndex1, Tcoev.SiteIndex2, ...
    1, length(SiteIDs), length(SiteIDs)));
% Make sure the network is undirected
Wcoev = Wcoev | Wcoev';
% Remove self-edges (if any)
Wcoev = logical(Wcoev - diag(diag(Wcoev)));

% Get associations with structure distance evidence
validRows = Tc.Structure_distance_evidence == true;
Tsd = Tc(validRows, :);

% Create Structure Distance Network
Wsd = logical(sparse(Tsd.SiteIndex1, Tsd.SiteIndex2, ...
    1, length(SiteIDs), length(SiteIDs)));
% Make sure the network is undirected
Wsd = Wsd | Wsd';
% Remove self-edges (if any)
Wsd = logical(Wsd - diag(diag(Wsd)));

% Filter for sites in the networks
validSites = sum(Wcoev | Wsd, 2) > 0;
Wcoev = Wcoev(validSites, validSites);
Wsd = Wsd(validSites, validSites);

% Create Site Information Table
genes = regexprep(SiteIDs(validSites), '_[^\char]*', '');
residues = regexprep(SiteIDs(validSites), '[^\char]*_', '');
[Gene, ~, geneIndices] = unique(genes);

Site = table();
Site.GeneIndex = geneIndices;
Site.Residue = residues;

% Save the results
outputFolder =  '../../data/processed/';
save([outputFolder, 'ptmcode_networks.mat'], ...
    'Wcoev', 'Wsd', 'Site', 'Gene');
%%


















