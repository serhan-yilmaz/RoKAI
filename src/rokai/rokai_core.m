function [Xs] = rokai_core(X, Wkin2site, Wkin2kin, Wsite2site, varargin)
    p = inputParser;
    p.CaseSensitive = false;
    validPhosphorylation = @(x) validateattributes(x, {'numeric'}, ...
        {'2d', 'nonempty', 'real'});   
    validNetworkMatrix = ...
        @(x) validateattributes(x, {'logical', 'numeric'}, ...
        {'2d', 'nonempty','real','nonnan','finite'});
    validSquareNetworkMatrix = ...
        @(x) validateattributes(x, {'logical', 'numeric'}, ...
        {'2d','square', 'real','nonnan','finite'});
    addRequired(p, 'X', validPhosphorylation);
    addRequired(p, 'Wkin2site', validNetworkMatrix);
    addRequired(p, 'Wkin2kin', validSquareNetworkMatrix);
    addRequired(p, 'Wsite2site', validSquareNetworkMatrix);
    addParameter(p, 'KeepMissingSites', true, @islogical);
    parse(p, X, Wkin2site, Wkin2kin, Wsite2site, varargin{:});
    param = p.Results;
    [nKinase, nSite] = size(Wkin2site);
    
    if(isempty(Wkin2kin))
        Wkin2kin = sparse(nKinase, nKinase);
    end
    
    if(isempty(Wsite2site))
        Wsite2site = sparse(nSite, nSite);
    end
    
    if(size(X, 1) ~= nSite)
       error(['The number of rows in phosphorylation data X must match ', ...
           'the number of columns in kinase-phosphosite network.']);
    end
    
    if(size(Wkin2kin, 1) ~= nKinase)
       error(['The number of kinases in the kinase-phosphosite network ', ...
           'Wkin2site must match the kinase-kinase network Wkin2kin.']);
    end
        
    if(size(Wsite2site, 1) ~= nSite)
       error(['The number of phosphosites in the kinase-phosphosite network ', ...
           'Wkin2site must match the phosphosite-phosphosite network Wsite2site.']);
    end
    
    Xs = nan(size(X));
    nCondition = size(X, 2);
    for iCondition = 1:nCondition
        Xv = X(:, iCondition);
        if(param.KeepMissingSites)
            validsites = true(size(Xv));
        else
            validsites = ~isnan(Xv);
        end
        Wk2s = Wkin2site(:, validsites);
        Ws2s = Wsite2site(validsites, validsites);
        
        [n, m] = size(Wk2s);
        site_indices = 1:m;
        %kinase_indices = m + (1:n);
        network = [Ws2s Wk2s'; Wk2s Wkin2kin];
        
        I = [Xv(validsites); nan(n, 1)];
        Q = rokai_circuit(I, network);
        Xs(validsites, iCondition) = Q(site_indices);
        Xs(isnan(Xv), iCondition) = NaN;
    end
end

