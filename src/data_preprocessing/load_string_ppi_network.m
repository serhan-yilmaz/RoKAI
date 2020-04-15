%% Load PPI network
% We obtain the protein-protein interaction network from STRING. 
% Due to its excessive file size, '9606.protein.links.v11.0.txt'
% is not included in this repository, but is available at: 
% https://string-db.org/cgi/download.pl
dataFolder = '../../data/big/';
filename = [dataFolder, '9606.protein.links.v11.0.txt'];
ds = datastore(filename);
STRING_PPI = ds.readall();

Proteins = unique([STRING_PPI.protein1; STRING_PPI.protein2]);
[~, proteinIndices1] = ismember(STRING_PPI.protein1, Proteins);
[~, proteinIndices2] = ismember(STRING_PPI.protein2, Proteins);
Proteins = regexprep(Proteins, '9606.', '');

PPI = sparse(proteinIndices1, proteinIndices2, ...
    STRING_PPI.combined_score, length(Proteins), length(Proteins));
clear STRING_PPI proteinIndices1 proteinIndices2

% Save the results
outputFolder =  '../../data/processed/';
save([outputFolder, 'string_ppi_network.mat'], 'Proteins', 'PPI');
%%













