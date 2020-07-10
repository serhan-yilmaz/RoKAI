function [Wkin2site, Wkin2kin, Wsite2site] = rokai_networks( NetworkData, varargin )
    p = inputParser;
    p.CaseSensitive = false;
    addRequired(p, 'NetworkData', @isstruct);
    addParameter(p, 'IncludePPI', true, @islogical);
    addParameter(p, 'IncludeCoevolution', true, @islogical);
    addParameter(p, 'IncludeStructureDistance', true, @islogical);
    parse(p, NetworkData, varargin{:});
    param = p.Results;

    Wkin2site = NetworkData.Wkin2site;
    nSite = size(Wkin2site, 2);
    if(param.IncludePPI)
        Wkin2kin = NetworkData.Wkin2kin*1e-3;
    else
        Wkin2kin = [];
    end
    Wsite2site = sparse(nSite, nSite);
    if(param.IncludeCoevolution)
        Wsite2site = Wsite2site + NetworkData.Wsite2site_coev;
    end
    
    if(param.IncludeStructureDistance)
        Wsite2site = Wsite2site + NetworkData.Wsite2site_sd;
    end
end

