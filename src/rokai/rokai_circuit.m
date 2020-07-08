function [ V, R, I, P, Pd ] = rokai_circuit( C, W, varargin )
    p = inputParser;
    p.CaseSensitive = false;
    validVector = @(x) validateattributes(x, {'numeric'}, ...
        {'vector','nonempty','real'});
    validSquareMatrix = @(x) validateattributes(x, {'logical', 'numeric'}, ...
        {'2d','square','nonempty','real','nonnan','finite'});
    validNonnegScalar = @(x) validateattributes(x, {'numeric'}, ...
        {'scalar', 'nonempty', 'real', 'nonnan', 'nonnegative', '<=', 1});
    addRequired(p, 'C', validVector);
    addRequired(p, 'W', validSquareMatrix);
    addParameter(p, 'Tau', 1e8, validNonnegScalar);
    addParameter(p, 'DampingFactor', 0.5, validNonnegScalar);
    parse(p, C, W, varargin{:});
    param = p.Results;
    df = param.DampingFactor;
    W = W - diag(diag(W));
    tau = param.Tau;
    
    n = length(C);
    r = (1-df) / df;
    Rd = ~isnan(C) + 1/(r*tau) * isnan(C);
    Rs = spdiags(Rd + sum(W, 2)/r, 0, n, n);
    C(isnan(C)) = 0;
    R = Rs - W / r;
    V = R \ C;
    
    if(nargout >= 3)
        indices = find(R);
        [i1, i2] = ind2sub(size(R), indices);
        valids = i1 ~= i2;
        i1 = i1(valids);
        i2 = i2(valids);
        indices = indices(valids);
        Vv = (V(i1) - V(i2));
        Iv = -R(indices) .* Vv;
        I = sparse(i1, i2, Iv, size(R,1), size(R,2)); 
        Pv = Iv .* Vv;
        P = sparse(i1, i2, Pv, size(R,1), size(R,2));
        Pd = full(Rd .* V .* V);
    end
end

