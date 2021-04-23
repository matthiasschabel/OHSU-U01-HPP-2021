function roc = analyzeROC(FPR,TPR)
    roc.FPR = FPR;
    roc.TPR = TPR;
    roc.Sensitivity = TPR;
    roc.Specificity = 1-FPR;
    roc.AUC = trapz(FPR,TPR);
    
    % Youden's J
    roc.J = roc.Sensitivity - (1-roc.Specificity);
    
    f=find(roc.J==max(roc.J),1,'first');
    
    roc.Jmax = roc.J(f);
    roc.Jindex = f;
end