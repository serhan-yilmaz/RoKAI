function [KT, ST] = rokai(T, NetworkData, varargin)
    p = inputParser;
    p.CaseSensitive = false;
    validInferenceMethod = @(x) ...
        any(validatestring(x, {'mean', 'zscore', 'linear', 'gsea'}));
    validIdentifier = @(x) ...
        any(validatestring(x, {'protein-position', 'flanking'}));
    addRequired(p, 'T', @istable);
    addRequired(p, 'NetworkData', @isstruct);
    addParameter(p, 'Identifier', 'protein-position', validIdentifier);
    addParameter(p, 'InferenceMethod', 'mean', validInferenceMethod);
    addParameter(p, 'IncludeMissingSites', true, @islogical);
    addParameter(p, 'IncludePPI', true, @islogical);
    addParameter(p, 'IncludeCoevolution', true, @islogical);
    addParameter(p, 'IncludeStructureDistance', true, @islogical);
    parse(p, T, NetworkData, varargin{:});
    param = p.Results;
    
    X = rokai_input_mapping(T, NetworkData, param.Identifier);
      
    [Wkin2site, Wkin2kin, Wsite2site] = rokai_networks(NetworkData, ...
        'IncludePPI', param.IncludePPI, ...
        'IncludeCoevolution', param.IncludeCoevolution, ...
        'IncludeStructureDistance', param.IncludeStructureDistance);

    [Xs] = rokai_core(X, Wkin2site, Wkin2kin, Wsite2site, ...
        'KeepMissingSites', param.IncludeMissingSites);
    ST = sitetable(NetworkData.Site, X, Xs, Wkin2site);

    A = rokai_inference(Xs, Wkin2site, param.InferenceMethod);
    KT = kinasetable(NetworkData.Kinase, A);
end
























