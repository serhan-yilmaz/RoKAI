function [ X ] = rokai_input_mapping( T, NetworkData, identifier )
    if(nargin < 3)
        identifier = 'protein-position';
    end
    switch(identifier)
        case 'protein-position'
            T.Properties.VariableNames{1} = 'Protein';
            T.Properties.VariableNames{2} = 'Position';
            T.Properties.VariableNames{3} = 'Quantification';
            position = regexprep(T.Position, '[^0-9]', '');
            T.Identifier = cellstr(join([T.Protein, position], '_'));
        case 'flanking'
            T.Properties.VariableNames{1} = 'SiteFlanking';
            T.Properties.VariableNames{2} = 'Quantification';
            T.Identifier = T.SiteFlanking;
        otherwise
            error('Invalid identifier type.');
    end
    
    
    [b, ib] = ismember(upper(T.Identifier), upper(NetworkData.Site.Identifier));
    X = nan(NetworkData.nSite, 1);
    X(ib(b)) = T.Quantification(b);
end

