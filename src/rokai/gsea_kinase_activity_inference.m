function [kinaseBeta] = gsea_kinase_activity_inference(V, substrates, varargin)
    p = inputParser;
    p.CaseSensitive = false;
    addParameter(p, 'nPerm', 10000, @(x) true);
    addParameter(p, 'Cache', [], @(x) true);
    parse(p, varargin{:});
    param = p.Results;

    nKinase = size(substrates, 1);
    kinaseBeta = nan(nKinase, 1);
    validKinases = sum(substrates, 2) > 0;
    kinaseIndices = find(validKinases);
    [Vs, si] = sort(V, 'descend');
    subs = substrates(:, si);
    
    if(isempty(param.Cache))
        nPerm = param.nPerm;
        GSEA_perm = gsea_inference_permutation(Vs, subs, nPerm);
    else
        nPerm = size(param.Cache, 1);
        GSEA_perm = param.Cache;
    end
    
    subs = subs';
    for iValidKinase = 1:length(kinaseIndices)
        iKinase = kinaseIndices(iValidKinase);
        subst = subs(:, iKinase);
        nSubs = nnz(subst);
        [score] = gsea(Vs, subst);
        pval = (1 + nnz(GSEA_perm(:, nSubs) >= score)) / (1 + nPerm);
        kinaseBeta(iKinase) = -log10(pval);
    end
end

