% fit is output from nlfit
function roc = generateROC(fit,msk)
    if (nargin < 2)
        msk = true(size(fit(2).x));
    end
    
    [ux,fux] = unique(fit(1).x);

    predloalphaint = [];
    predhialphaint = [];

    for idx=1:size(fit(1).predloalpha,2) 
        predloalphaint(:,idx)=interp1(fit(1).x(fux),fit(1).predloalpha(fux,idx),fit(2).x);
        predhialphaint(:,idx)=interp1(fit(1).x(fux),fit(1).predhialpha(fux,idx),fit(2).x); 
    end

    ttmp = [predloalphaint fliplr(predhialphaint)];

    x = [];
    y = [];

    for idx=1:size(ttmp,2) 
        x(idx) = idx;
        y(idx) = sum(fit(2).y(msk)<=ttmp(msk,idx));
    end

    x = x/max(x);
    y = y/sum(msk);

    roc = analyzeROC(x,y);
    
    % use max Youden J statistic to find optimal percentile cutoff
    %   the computed J index is 1-FPR location of maximum point relative to diagonal
    %   we are using confidence interval percentiles, so values below/above those need to be halved
    roc.predPercentiles = (1+[-fliplr(fit(1).options.prediction_alpha) fit(1).options.prediction_alpha])/2;
    
    roc.optimalPercentile = roc.predPercentiles(roc.Jindex); 
end
