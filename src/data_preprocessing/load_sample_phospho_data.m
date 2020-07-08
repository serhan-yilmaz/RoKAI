% For use in the demo, we utilize the phosphorylation data within
% the Atlas compiled by Ochoa et al. (2016), obtained from: 
% http://phosfate.com/download.html
% Here, we include only a subset of the data, a single experiment of 
% Pan et al. (2009) comparing EGF treated samples with the untreated.
% This data is provided with the repository in 'sample_phospho_data.csv'.
% 
% References:
% Ochoa, David, et al. "An atlas of human kinase regulation." 
%  Molecular systems biology 12.12 (2016).
% Pan, Cuiping, et al. "Global effects of kinase inhibitors on  
%  signaling networks revealed by quantitative phosphoproteomics." 
%  Molecular & Cellular Proteomics 8.12 (2009): 2796-2808.
dataPath = '../../data/';
filename = [dataPath, 'sample_phospho_data_ensembl.csv'];
ds = tabularTextDatastore(filename, 'delimiter', ',');
T = ds.readall();

% Obtain a list of unique protein identifiers
Proteins = unique(T.Protein);
%% Mapping of Ensembl Protein Identifiers to UniprotKb
% We map the ensembl protein identifiers (ENSP) to uniprotkb in order to
% use the phosphorylation data together with PhosphositePlus
% kinase-substrate annotations and other functional networks.
% We obtain the ensembl protein to uniprotkb id mappings using Uniprot
% mapping tool: https://www.uniprot.org/uploadlists/
% The 'phospho_data_proteins_uniprotkb.tab' contains the results of the 
% query run on 2020-04-13, using the ids given in 'Proteins' as input. 
dataPath = '../../data/raw/';
filepath = [dataPath, 'phospho_data_proteins_uniprotkb.tab'];
ds = tabularTextDatastore(filepath, ...
    'FileExtensions', '.tab');
ds.TextscanFormats = repmat({'%q'}, 1, length(ds.VariableNames));
Mapping = ds.readall();

% Map to primary identifier
Mapping.PrimaryEnsp = regexprep(Mapping.Ensembl_Protein, ',[^.]*', '');

% Map to uniprotkb protein identifiers
[b, ib] = ismember(T.Protein, Mapping.Ensembl_Protein);

% Filter out rows that are not mapped
T = T(b, :);

% Assign the uniprotkb Id
T.Protein = Mapping.Entry(ib(b), :);

% Save the results
outputPath = '../../data/';
% save([outputPath, 'sample_phospho_data.mat'], 'T');
writetable(T, [outputPath, 'sample_phospho_data_uniprotkb.csv']);

