function [kinaseBeta] = linear_kinase_activity_inference(V, substrates, k)
    if(nargin < 3)
       k = 0.1; 
    end
    validKinases = sum(substrates, 2) > 0;
    X = full(substrates(validKinases, :)');
    Y = V;
    [n, m] = size(X);    
%     k = 0.1;    
    pseudo = sqrt(k) * eye(m+1);
    Xplus = [ones(n, 1) X; pseudo];
    Yplus = [Y; zeros(m + 1, 1)];
    beta = Xplus\Yplus;
    
%     kinaseBeta = linsolve(full(kinase_substrates_b), Vb);
    kinaseBeta = nan(size(substrates, 1), 1);
    kinaseBeta(validKinases) = beta(2:end);
%     kinaseBeta(invalidKinases) = NaN;
end

