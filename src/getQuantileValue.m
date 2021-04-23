function q = getQuantileValue(ref,val)
    if (isempty(ref))
        q = NaN;
        return;
    end

    quantiles = atob(0,1,min(length(ref),100));
    
    qq = nanquantile(ref(:),quantiles);
    
    if (val < qq(1))
        q = quantiles(1);
    elseif (val > qq(end))
        q = quantiles(end);
    else
        q = interp1(qq,quantiles,val,'linear');
    end
end