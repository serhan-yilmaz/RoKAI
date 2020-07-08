function [A] = rokai_inference( X, Wkin2site, method )
    if(nargin < 3)
       method = 'mean'; 
    end
    nKinase = size(Wkin2site, 1);
    nCondition = size(X, 2);
    A = nan(nKinase, nCondition);
    for iCondition = 1:nCondition
        validSites = ~isnan(X(:, iCondition));
        V = X(validSites, iCondition);
        kinase_substrates = Wkin2site(:, validSites);
        switch(method)
            case 'mean'
                kinaseScores = (kinase_substrates * V) ...
                                    ./ (kinase_substrates * ones(size(V)));
            case 'zscore'
                S = std(V);
                kinaseScores = (kinase_substrates * V) ...
                                ./ (S*sqrt(kinase_substrates * ones(size(V))));
            case 'linear'
                k = 0.1;
                kinaseScores = linear_kinase_activity_inference(V, kinase_substrates, k);
            case 'gsea'
                nPerm = 10000;
                kinaseScores = gsea_kinase_activity_inference(V, kinase_substrates, nPerm);
            otherwise
                error('Invalid method: %s', method);
        end
        A(:, iCondition) = kinaseScores;
    end
end

