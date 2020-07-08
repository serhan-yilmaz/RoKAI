function [kinaseBeta] = gsea_kinase_activity_inference(V, substrates, nPerm)
    nKinase = size(substrates, 1);
    kinaseBeta = nan(nKinase, 1);
    validKinases = sum(substrates, 2) > 0;
    kinaseIndices = find(validKinases);
    [Vs, si] = sort(V, 'descend');
    subs = substrates(:, si)';
    
    GSEA_perm = gsea_inference_permutation(Vs, subs, nPerm);
    
%     nPerm = size(perm, 1);
    for iValidKinase = 1:length(kinaseIndices)
        iKinase = kinaseIndices(iValidKinase);
        subst = subs(:, iKinase);
        nSubs = nnz(subst);
        [score] = gsea(Vs, subst);
        pval = (1 + nnz(GSEA_perm(:, nSubs) >= score)) / (1 + nPerm);
        kinaseBeta(iKinase) = -log10(pval);
    end
end

function [G] = gsea_inference_permutation(V, Wkin2site, nPerm)
    n = length(V);
    nMax = max(sum(Wkin2site, 2));    
    Vs = sort(V, 'descend');
    G = nan(nMax, nPerm);
    for iPermutation = 1:nPerm
        perms = randperm(n, nMax);
        I = false(n, 1);
        for iN = 1:nMax
            I(perms(iN)) = true;
            score = gsea(Vs, I);
            G(iN, iPermutation) = score;
        end        
    end
end

