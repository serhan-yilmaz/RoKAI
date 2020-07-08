function [Kinase] = kinasetable(Kinase, A)
    Kinase.Activity = A;
    valids = ~isnan(A);
    Kinase = Kinase(valids, :);
    [~, si] = sort(abs(Kinase.Activity), 'descend');
    Kinase = Kinase(si, :);
end

