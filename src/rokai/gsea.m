function [score] = gsea(Xs, Is)
    Xp = abs(Xs(Is));
    X1 = zeros(size(Xs));
    X1(Is) = Xp / sum(Xp);    
    X0 = ~Is / (length(Is) - nnz(Is));
    scores = cumsum(X1 - X0);
%     Xp = abs(Xs(Is));
%     Xc = ~Is / ((length(Is) - nnz(Is))*-1);
%     Xc(Is) = Xp / sum(Xp);
%     scores = cumsum(Xc);
    score = max(abs(scores));
end



















