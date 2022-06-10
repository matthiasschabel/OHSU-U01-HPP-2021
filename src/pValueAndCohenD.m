function [p,d,c] = pValueAndCohenD(ds1,ds2,robust)
    if (nargin<3)
        robust = true;
    end
    
    if (robust)
        center = @(x)median(flatten(x));
        spread = @(x)iqr(flatten(x))/1.349;
    else
        center = @(x)mean(flatten(x));
        spread = @(x)std(flatten(x));
    end
    
    ds1 = ds1(isfinite(ds1));
    ds2 = ds2(isfinite(ds2));

    [~,p] = kstest2(ds1,ds2);
    
    c(1) = center(ds1);
    c(2) = center(ds2);
   
    pooledStdDev = spread([ds1 ds2]);
    
    d = abs(diff(c))/pooledStdDev;
end