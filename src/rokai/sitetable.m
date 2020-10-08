function [Site] = sitetable(Site, Xraw, Xrefined, Wkin2site)
    Site.RoKAI_Q = Xrefined;    
    Site.Raw_Q = Xraw;
    
    valids = ~isnan(Site.Raw_Q);
    if(nargin >= 4)
        valids = valids & (sum(Wkin2site, 1) > 0)';
    end
    
    Site = Site(valids, :);
    [~, si] = sort(abs(Site.RoKAI_Q), 'descend');
    Site = Site(si, :);
end

