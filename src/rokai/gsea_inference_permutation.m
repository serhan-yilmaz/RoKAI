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
    G = G';
end