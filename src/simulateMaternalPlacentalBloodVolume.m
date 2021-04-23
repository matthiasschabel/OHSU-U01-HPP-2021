function vmpb = simulateMaternalPlacentalBloodVolume(gwfit,Hbfit,SpO2fit,T2sfit,numberOfRealizations,noiseFactor)
    if (nargin<5)
        numberOfRealizations = 10000;
    end
    
    if (nargin<6)
        noiseFactor = 1;
    end

    Hbfit(1).noise = correlatedNoise(numberOfRealizations,Hbfit(1).cov)*noiseFactor;
    Hbfit(2).noise = correlatedNoise(numberOfRealizations,Hbfit(2).cov)*noiseFactor;
    SpO2fit(1).noise = correlatedNoise(numberOfRealizations,SpO2fit(1).cov)*noiseFactor;
    SpO2fit(2).noise = correlatedNoise(numberOfRealizations,SpO2fit(2).cov)*noiseFactor;
    T2sfit(1).noise = correlatedNoise(numberOfRealizations,T2sfit(1).cov)*noiseFactor;
    T2sfit(2).noise = correlatedNoise(numberOfRealizations,T2sfit(2).cov)*noiseFactor;

    for idx=1:numberOfRealizations
        Hbfitpar1(idx,:) = Hbfit(1).beta+Hbfit(1).noise(idx,:);
        Hbfitpar2(idx,:) = Hbfit(2).beta+Hbfit(2).noise(idx,:);
        SpO2fitpar1(idx,:) = SpO2fit(1).beta+SpO2fit(1).noise(idx,:);
        SpO2fitpar2(idx,:) = SpO2fit(2).beta+SpO2fit(2).noise(idx,:);
        T2sfitpar1(idx,:) = T2sfit(1).beta+T2sfit(1).noise(idx,:);
        T2sfitpar2(idx,:) = T2sfit(2).beta+T2sfit(2).noise(idx,:);

        Hbfitval1 = Hbfit(1).model(Hbfitpar1(idx,:),gwfit);
        Hbfitval2 = Hbfit(2).model(Hbfitpar2(idx,:),gwfit);
        SpO2fitval1 = SpO2fit(1).model(SpO2fitpar1(idx,:),gwfit);
        SpO2fitval2 = SpO2fit(2).model(SpO2fitpar2(idx,:),gwfit);

        % 0.6206 factor converts [Hb] in g/dl to mmol/l
        dR2s = 20.2e-3*0.6206*(Hbfitval1.*(100-SpO2fitval1)/100-Hbfitval2.*(100-SpO2fitval2)/100);

        T2sfitval1 = T2sfit(1).model(T2sfitpar1(idx,:),gwfit);
        T2sfitval2 = T2sfit(2).model(T2sfitpar2(idx,:),gwfit);

        vmpb(idx,:) = (1./T2sfitval1-1./T2sfitval2)./dR2s;
    end
end