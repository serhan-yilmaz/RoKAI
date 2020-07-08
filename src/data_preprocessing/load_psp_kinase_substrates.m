%% Load Kinase-Substrate network
% We obtain the kinase-substrate annotations from PhosphositePlus at:
% https://www.phosphosite.org/staticDownloads
dataFolder = '../../data/raw/';
filename = [dataFolder, 'Kinase_Substrate_Dataset'];
ds = datastore(filename);
ds.TextscanFormats = repmat({'%q'}, 1, length(ds.VariableNames));
KStable = ds.readall();

% Filter for KIN_ORGANISM = 'human' and SUB_ORGANISM = 'human'
validRows = strcmpi(KStable.KIN_ORGANISM, 'human');
validRows = validRows & strcmpi(KStable.SUB_ORGANISM, 'human');
KStable = KStable(validRows, :);

% Remove isoform identifiers
KStable.KIN_ACC_ID = regexprep(KStable.KIN_ACC_ID, '-[^.]', '');
KStable.SUB_ACC_ID = regexprep(KStable.SUB_ACC_ID, '-[^.]', '');
% Extract Kinases
[kins, ib, ic] = unique(KStable.KIN_ACC_ID);
KStable.KinaseIndex = ic;

Kinase = table();
Kinase.ID = kins;
Kinase.Gene = KStable.GENE(ib);
Kinase.Name = KStable.KINASE(ib);

% Create unique substrate identifier
KStable.SubstrateID = cellstr(join([KStable.SUB_ACC_ID, KStable.SUB_MOD_RSD], '_'));

% Extract Substrates
[subst, ib, ic] = unique(KStable.SubstrateID);
KStable.SubstrateIndex = ic;

Substrate = table();
Substrate.ID = subst;
Substrate.Gene = KStable.SUB_GENE(ib);
Substrate.Name = KStable.SUBSTRATE(ib);
Substrate.Protein = KStable.SUB_ACC_ID(ib);
Substrate.Residue = KStable.SUB_MOD_RSD(ib);
Substrate.Flanking = KStable.SITE____7_AA(ib);

% Create Kinase-Substrate network
KS = logical(sparse(KStable.KinaseIndex, KStable.SubstrateIndex, ...
        1, height(Kinase), height(Substrate)));
%% Load Known Phosphorylation Sites
% We obtain a list of known phosphosites from PhosphositePlus.
% Due to its excessive file size, 'Phosphorylation_site_dataset'
% is not included in this repository, but is available at: 
% https://www.phosphosite.org/staticDownloads
dataFolder = '../../data/big/';
filename = [dataFolder, 'Phosphorylation_site_dataset'];
ds = datastore(filename);
ds.TextscanFormats = repmat({'%q'}, 1, length(ds.VariableNames));
Sitetable = ds.readall();

% Filter for ORGANISM = 'human'
validRows = strcmpi(Sitetable.ORGANISM, 'human');
Sitetable = Sitetable(validRows, :);

% Create a unique identifier
Sitetable.Residue = regexprep(Sitetable.MOD_RSD, '-[^.]*', '');
Sitetable.ID = join([Sitetable.ACC_ID Sitetable.Residue], '_');

% Combine the site table with substrate table
siteIds = unique([Sitetable.ID; Substrate.ID]);
Site = table();
Site.Protein = cell(length(siteIds), 1);
Site.Residue = cell(length(siteIds), 1);
Site.Flanking = cell(length(siteIds), 1);

% Map the sitetable to combined 'Site' table
[b, indices] = ismember(siteIds, Sitetable.ID);
Site.Protein(b) = Sitetable.ACC_ID(indices(b));
Site.Residue(b) = Sitetable.Residue(indices(b));
Site.Flanking(b) = Sitetable.SITE____7_AA(indices(b));

% Map the substrates to combined 'Site' table
[b, indices] = ismember(siteIds, Substrate.ID);
Site.Protein(b) = Substrate.Protein(indices(b));
Site.Residue(b) = Substrate.Residue(indices(b));
Site.Flanking(b) = Substrate.Flanking(indices(b));

% Map the kinase-substate network to combined site table
[~, indices] = ismember(Substrate.ID, siteIds);
Msubstrate2site = sparse((1:height(Substrate))', ... 
    indices, 1, height(Substrate), height(Site));
KS = logical(KS * Msubstrate2site);

% Create a list of unique proteins
Protein = unique(Site.Protein);

% % Remove redundant fields
clear KStable Sitetable

% Save the results
outputFolder =  '../../data/processed/';
save([outputFolder, 'psp_sites_and_kinase_substrate_network.mat'], ...
    'Protein', 'Site', 'Kinase', 'KS');
%%









