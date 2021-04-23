% analyze human placenta T2* mapping data
exportFiles = true;
exportFiles = false;

rootDirectory = '/Volumes/Data/U01';

figureDirectory = [rootDirectory '/Figures'];
% figureDirectory = [rootDirectory '/Figures_quality_masked'];
% figureDirectory = [rootDirectory '/Figures_all_data'];

% get adjudication groups/severity from REDCap
%   -2 is invalid/unrecorded value, -1 is abnormal, 0 is uncomplicated, 1 is adverse
category = adjudicationGroup;

% set up various data masks

% quality control mask
qualityMask = true(size(category));
% qualityMask = qualityMetricROI >= 0.5;
% qualityMask = uncertaintyMetricT2s<=0.25;
% qualityMask = qualityMetricROI>=0.95;
% qualityMask = (qualityMetricROI >= quantile(qualityMetricROI,0.1) & ...
%                uncertaintyMetricT2s <= quantile(uncertaintyMetricT2s,0.9));

iugrMask = ((studyIDs>=600 & studyIDs<700) | siteNumber==3);                    
% iugrMask = ((studyIDs>=600 & studyIDs<700) | siteNumber==3) & category == 1;

sgaStudyIDs = [201 272 282 320 345 349 351 378 603 604 605 607 608 610 611 614 407 490 527 801 802 808];

sgaMask = false(size(studyIDs));

for idx=1:length(sgaStudyIDs) 
    sgaMask(studyIDs==sgaStudyIDs(idx))=true; 
end

% preMask = true(size(category));
% preMask = qualityMask;
% preMask = qualityMask & siteNumber ~= 3;
preMask = qualityMask & ~iugrMask & siteNumber ~= 3;    % exclude all BCNatal and IUGR studies

controlMask = category==0 & preMask;
adverseMask = category==1 & preMask;
abnormalMask = category==-1 & preMask;

ohsuMask = siteNumber==1 & preMask;
utahMask = siteNumber==2 & preMask;

% the following analysis excludes IUGR-specific studies
siteMasks = {controlMask & ohsuMask,controlMask & utahMask,controlMask & (ohsuMask | utahMask)};
categoryMasks = {controlMask,adverseMask,abnormalMask};
categoryMasksIncludingIUGR = {controlMask,adverseMask | iugrMask,abnormalMask};
adverseSiteMasks = {adverseMask & ohsuMask,adverseMask & utahMask};
severityMasks = {adverseMask & adjudicationSeverity==0,adverseMask & adjudicationSeverity==1};

fetalSexMasks = {fetalSex==1 & controlMask,fetalSex==2 & controlMask};
maternalAgeMasks = {controlMask & maternalAge<35,controlMask & maternalAge>=35};

bmiMasks = {controlMask & prepregnancyBMI<18.5,...
            controlMask & prepregnancyBMI>=18.5&prepregnancyBMI<25,...
            controlMask & prepregnancyBMI>=25&prepregnancyBMI<30,...
            controlMask & prepregnancyBMI>=30};
        
lowBirthWeightMask = preMask & birthWeightPercentile<=10;
veryLowBirthWeightMask = preMask & birthWeightPercentile<=5;

birthWeightPercentileMasks = {controlMask & birthWeightPercentile>10,...
                              lowBirthWeightMask,...
                              veryLowBirthWeightMask};

% site-specific masks ( could add fourth mask for IUGR cases)
ohsuCategoryMasks = {ohsuMask & controlMask,ohsuMask & adverseMask,ohsuMask & abnormalMask};
utahCategoryMasks = {utahMask & controlMask,utahMask & adverseMask,utahMask & abnormalMask};
    
% regression models
line = @(p,x)(p(1)+p(2)*x);
sigmoid = @(p,x)(p(1)./(1+exp(p(2)*(x-p(3))))+p(4));
sigmoid_derivative = @(p,x)(-(exp(p(2).*(x+p(3))).*p(1).*p(2))./(exp(p(2).*p(3))+exp(p(2).*x)).^2);

sigmoidGuessT2s = [-60 -.25 30 80];

uniqueStudyIDs = unique(studyIDs);

for gidx=1:length(uniqueStudyIDs)
    currentStudyID = uniqueStudyIDs(gidx);
    currentNumberOfStudies = length(find(studyIDs==currentStudyID));
    
    numberOfStudies(studyIDs==currentStudyID) = currentNumberOfStudies;
    
    repeatGestationalDays = gestationalDay(studyIDs==currentStudyID);
    repeatT2sValues = medianT2s(studyIDs==currentStudyID);
    
    medianGD(studyIDs==currentStudyID,1:currentNumberOfStudies-1) = repmat((repeatGestationalDays(1:end-1)+repeatGestationalDays(2:end))/2,[currentNumberOfStudies 1]);
    deltaGD(studyIDs==currentStudyID,1:currentNumberOfStudies-1) = repmat(diff(repeatGestationalDays),[currentNumberOfStudies 1]);
    deltaT2s(studyIDs==currentStudyID,1:currentNumberOfStudies-1) = repmat(diff(repeatT2sValues),[currentNumberOfStudies 1]);
end

medianGD(medianGD==0|deltaGD==0|deltaT2s==0) = NaN;
deltaGD(medianGD==0|deltaGD==0|deltaT2s==0) = NaN;
deltaT2s(medianGD==0|deltaGD==0|deltaT2s==0) = NaN;

% fit SLN shift parameter vs. gestational week by site
figure;
[ph,siteFit.SLNshift] = compareModelFits(gestationalDay/7,shifted_lognormal_distfit_T2s(:,1)',siteMasks,line,[0 1]);
title('SLN shift by site');
xlabel('Gestational weeks');
ylabel('T2* shift (ms)');
legend([ph{1}(2) ph{2}(2) ph{3}(2)],{'OHSU','Utah','Both'})
if (exportFiles) 
    exportfig(gcf,[figureDirectory '/SLN_shift_by_site.eps'],'Color','rgb');
else
    pause(.1);
    drawnow;
end

% fit SLN shift parameter vs. gestational week by control/adverse/abnormal
figure;
[ph,categoryFit.SLNshift] = compareModelFits(gestationalDay/7,shifted_lognormal_distfit_T2s(:,1)',categoryMasks,line,[0 1]);
title('SLN shift by category');
xlabel('Gestational weeks');
ylabel('T2* shift (ms)');
legend([ph{1}(2) ph{2}(2) ph{3}(2)],{'Control','Adverse','Abnormal'})
if (exportFiles) 
    exportfig(gcf,[figureDirectory '/SLN_shift_by_category.eps'],'Color','rgb');
else
    pause(.1);
    drawnow;
end

% fit SLN mu parameter vs. gestational week by site
figure;
[ph,siteFit.SLNmu] = compareModelFits(gestationalDay/7,shifted_lognormal_distfit_T2s(:,2)',siteMasks,sigmoid,[-1 -.25 30 4]);
title('SLN \mu by site');
xlabel('Gestational weeks');
ylabel('T2* \mu (ms)');
legend([ph{1}(2) ph{2}(2) ph{3}(2)],{'OHSU','Utah','Both'})
if (exportFiles) 
    exportfig(gcf,[figureDirectory '/SLN_mu_by_site.eps'],'Color','rgb');
else
    pause(.1);
    drawnow;
end

% fit SLN mu parameter vs. gestational week by control/adverse/abnormal
figure;
[ph,categoryFit.SLNmu] = compareModelFits(gestationalDay/7,shifted_lognormal_distfit_T2s(:,2)',categoryMasks,sigmoid,[-1 -.25 30 4]);
title('SLN \mu by category');
xlabel('Gestational weeks');
ylabel('T2* \mu (ms)');
legend([ph{1}(2) ph{2}(2) ph{3}(2)],{'Control','Adverse','Abnormal'})
if (exportFiles) 
    exportfig(gcf,[figureDirectory '/SLN_mu_by_category.eps'],'Color','rgb');
else
    pause(.1);
    drawnow;
end

% fit SLN sigma parameter vs. gestational week by site
figure;
[ph,siteFit.SLNsigma] = compareModelFits(gestationalDay/7,shifted_lognormal_distfit_T2s(:,3)',siteMasks,sigmoid,[-.2 -.4 30 .5]);
title('SLN \sigma by site');
xlabel('Gestational weeks');
ylabel('T2* \sigma (ms)');
legend([ph{1}(2) ph{2}(2) ph{3}(2)],{'OHSU','Utah','Both'})
if (exportFiles) 
    exportfig(gcf,[figureDirectory '/SLN_sigma_by_site.eps'],'Color','rgb');
else
    pause(.1);
    drawnow;
end

% fit SLN sigma parameter vs. gestational week by control/adverse/abnormal
figure;
[ph,categoryFit.SLNsigma] = compareModelFits(gestationalDay/7,shifted_lognormal_distfit_T2s(:,3)',categoryMasks,sigmoid,[-.2 -.4 30 .5]);
title('SLN \sigma by category');
xlabel('Gestational weeks');
ylabel('T2* \sigma (ms)');
legend([ph{1}(2) ph{2}(2) ph{3}(2)],{'Control','Adverse','Abnormal'})
if (exportFiles) 
    exportfig(gcf,[figureDirectory '/SLN_sigma_by_category.eps'],'Color','rgb');
else
    pause(.1);
    drawnow;
end

% fit sigmoid directly to measured T2* vs. gestational week data by site
figure;
[ph,siteFit.T2s] = compareModelFits(gestationalDay/7,medianT2s,siteMasks,sigmoid,sigmoidGuessT2s);
set(gca,'xlim',[10 40],'ylim',[0 120]);
title('T2* by site');
xlabel('Gestational weeks');
ylabel('T2* (ms)');
legend([ph{1}(2) ph{2}(2) ph{3}(2)],{'OHSU','Utah','Both'})
if (exportFiles) 
    exportfig(gcf,[figureDirectory '/median_T2s_by_site.eps'],'Color','rgb');
else
    pause(.1);
    drawnow;
end

% fit sigmoid directly to measured T2* vs. gestational week data by control/adverse/abnormal
figure;
[ph,categoryFit.T2s] = compareModelFits(gestationalDay/7,medianT2s,categoryMasks,sigmoid,sigmoidGuessT2s);
set(gca,'xlim',[10 40],'ylim',[0 120]);
title('T2* by category');
xlabel('Gestational weeks');
ylabel('T2* (ms)');
legend([ph{1}(2) ph{2}(2) ph{3}(2)],{'Control','Adverse','Abnormal'})
if (exportFiles) 
    exportfig(gcf,[figureDirectory '/median_T2s_by_category.eps'],'Color','rgb');
else
    pause(.1);
    drawnow;
end

% fit sigmoid directly to measured adverse T2* vs. gestational week data by site
figure;
[ph,adverseFit.T2s] = compareModelFits(gestationalDay/7,medianT2s,adverseSiteMasks,sigmoid,sigmoidGuessT2s);
set(gca,'xlim',[10 40],'ylim',[0 120]);
title('Adverse T2* by site');
xlabel('Gestational weeks');
ylabel('T2* (ms)');
legend([ph{1}(2) ph{2}(2)],{'OHSU','Utah'});
if (exportFiles) 
    exportfig(gcf,[figureDirectory '/median_T2s_adverse_by_site.eps'],'Color','rgb');
else
    pause(.1);
    drawnow;
end

% fit sigmoid directly to measured T2* vs. gestational week data by fetal sex
figure;
[ph] = compareModelFits(gestationalDay/7,medianT2s,fetalSexMasks,sigmoid,sigmoidGuessT2s);
set(gca,'xlim',[10 40],'ylim',[0 120]);
title('T2* by fetal sex');
xlabel('Gestational weeks');
ylabel('T2* (ms)');
legend([ph{1}(2) ph{2}(2)],{'Male','Female'})
if (exportFiles) 
    exportfig(gcf,[figureDirectory '/median_T2s_by_fetal_sex.eps'],'Color','rgb');
else
    pause(.1);
    drawnow;
end

% fit sigmoid directly to measured T2* vs. gestational week data by control/adverse/abnormal, including IUGR recruits from OHSU and BCNatal
figure;
[ph] = compareModelFits(gestationalDay/7,medianT2s,categoryMasksIncludingIUGR,sigmoid,sigmoidGuessT2s);
set(gca,'xlim',[10 40],'ylim',[0 120]);
title('T2* by category including IUGR cases');
xlabel('Gestational weeks');
ylabel('T2* (ms)');
legend([ph{1}(2) ph{2}(2) ph{3}(2)],{'Control','Adverse','Abnormal'})
if (exportFiles) 
    exportfig(gcf,[figureDirectory '/median_T2s_by_category_including_IUGR_and_BCN.eps'],'Color','rgb');
else
    pause(.1);
    drawnow;
end

% fit sigmoid directly to measured T2* vs. gestational week data by control/adverse/abnormal, underweighting studyIDs with multiple scans
figure;
[ph] = compareModelFits(gestationalDay/7,medianT2s,categoryMasks,sigmoid,sigmoidGuessT2s,1./sqrt(numberOfStudies));
set(gca,'xlim',[10 40],'ylim',[0 120]);
title('T2* by category with underweighted repeat studies');
xlabel('Gestational weeks');
ylabel('T2* (ms)');
legend([ph{1}(2) ph{2}(2) ph{3}(2)],{'Control','Adverse','Abnormal'})
if (exportFiles) 
    exportfig(gcf,[figureDirectory '/median_T2s_by_category_underweight_repeat_scans.eps'],'Color','rgb');
else
    pause(.1);
    drawnow;
end

% fit sigmoid directly to measured T2* vs. gestational week data by adverse severity
figure;
[ph] = compareModelFits(gestationalDay/7,medianT2s,severityMasks,sigmoid,sigmoidGuessT2s);
set(gca,'xlim',[10 40],'ylim',[0 120]);
title('T2* by adversity');
xlabel('Gestational weeks');
ylabel('T2* (ms)');
legend([ph{1}(2) ph{2}(2)],{'Mild','Severe'})
if (exportFiles) 
    exportfig(gcf,[figureDirectory '/median_T2s_by_severity.eps'],'Color','rgb');
else
    pause(.1);
    drawnow;
end

% fit sigmoid directly to measured T2* vs. gestational week data by prepregnancy BMI
figure;
[ph] = compareModelFits(gestationalDay/7,medianT2s,bmiMasks,sigmoid,sigmoidGuessT2s);
set(gca,'xlim',[10 40],'ylim',[0 120]);
title('T2* by BMI');
xlabel('Gestational weeks');
ylabel('T2* (ms)');
legend([ph{1}(2) ph{2}(2) ph{3}(2) ph{4}(2)],{'Underweight','Normal','Overweight','Obese'})
if (exportFiles) 
    exportfig(gcf,[figureDirectory '/median_T2s_by_BMI.eps'],'Color','rgb');
else
    pause(.1);
    drawnow;
end

% bmiq = nanquantile(prepregnancyBMI,[.1 .25 .5 .75 .9],2);
% 
% figure;
% [ph] = compareModelFits(gestationalDay/7,medianT2s,{category==0&prepregnancyBMI>=bmiq(2)&prepregnancyBMI<=bmiq(4),category==0&prepregnancyBMI<bmiq(2),category==0&prepregnancyBMI>bmiq(4)},sigmoid,qFit(51).beta);
% set(gca,'xlim',[10 40],'ylim',[0 120]);
% title('T2* by BMI');
% xlabel('Gestational weeks');
% ylabel('T2* (ms)');
% legend([ph{1}(2) ph{2}(2) ph{3}(2)],{'25%-75%','<25%','>75%'})
% if (exportFigures) 
%     exportfig(gcf,[figureDirectory '/median_T2s_by_BMI.eps','Color','rgb');
% end
%     
% % fit sigmoid directly to measured T2* vs. gestational week data by maternal age
% maq = nanquantile(maternalAge,[.1 .25 .5 .75 .9],2);
% 
% figure;
% [ph] = compareModelFits(gestationalDay/7,medianT2s,{category==0&maternalAge<=maq(3),category==0&maternalAge>maq(3),category==0&maternalAge>maq(4)},sigmoid,sigmoidGuessT2s);
% set(gca,'xlim',[10 40],'ylim',[0 120]);
% title('T2* by maternal age');
% xlabel('Gestational weeks');
% ylabel('T2* (ms)');
% legend([ph{1}(2) ph{2}(2) ph{3}(2)],{'<=50%','>50%','>75%'})
% if (exportFigures) 
%     exportfig(gcf,[figureDirectory '/median_T2s_by_maternal_age.eps','Color','rgb');
% end

% normal vs. advanced maternal age
figure;
[ph] = compareModelFits(gestationalDay/7,medianT2s,maternalAgeMasks,sigmoid,sigmoidGuessT2s);
set(gca,'xlim',[10 40],'ylim',[0 120]);
title('T2* by maternal age');
xlabel('Gestational weeks');
ylabel('T2* (ms)');
legend([ph{1}(2) ph{2}(2)],{'Age<35','Age>=35 (AMA)'})
if (exportFiles) 
    exportfig(gcf,[figureDirectory '/median_T2s_by_AMA.eps'],'Color','rgb');
else
    pause(.1);
    drawnow;
end

% fit sigmoid directly to measured T2* vs. gestational week data by birthweight percentile
figure;
[ph] = compareModelFits(gestationalDay/7,medianT2s,birthWeightPercentileMasks,sigmoid,sigmoidGuessT2s);
set(gca,'xlim',[10 40],'ylim',[0 120]);
title('T2* by birthweight percentile');
xlabel('Gestational weeks');
ylabel('T2* (ms)');
legend([ph{1}(2) ph{2}(2) ph{3}(2)],{'>10%','<=10%','<=5%'})
if (exportFiles) 
    exportfig(gcf,[figureDirectory '/median_T2s_by_birthweightPercentile.eps'],'Color','rgb');
else
    pause(.1);
    drawnow;
end

% linear fit to maternal hemoglobin vs. gestational week data by site
figure;
[ph,siteFit.maternalHb] = compareModelFits(gestationalDay/7,hemoglobin,siteMasks,line,[12 0]);
set(gca,'xlim',[10 40],'ylim',[8 16]);
title('Maternal Hb by site');
xlabel('Gestational weeks');
ylabel('Hb (mg/dl)');
legend([ph{1}(2) ph{2}(2) ph{3}(2)],{'OHSU','Utah','Both'})
if (exportFiles) 
    exportfig(gcf,[figureDirectory '/Hb_by_site.eps'],'Color','rgb');
else
    pause(.1);
    drawnow;
end

% linear fit to maternal hemoglobin vs. gestational week data by control/adverse/abnormal
figure;
[ph,categoryFit.maternalHb] = compareModelFits(gestationalDay/7,hemoglobin,categoryMasks,line,[12 0]);
set(gca,'xlim',[10 40],'ylim',[8 16]);
title('Maternal Hb by category');
xlabel('Gestational weeks');
ylabel('Hb (mg/dl)');
legend([ph{1}(2) ph{2}(2) ph{3}(2)],{'Control','Adverse','Abnormal'})
if (exportFiles) 
    exportfig(gcf,[figureDirectory '/Hb_by_category.eps'],'Color','rgb');
else
    pause(.1);
    drawnow;
end

% linear fit to estimated placental volume vs. gestational week data by site
figure;
[ph,siteFit.placentalVolume] = compareModelFits(gestationalDay/7,placentalVolume,siteMasks,line,[0 20]);
title('Placental volume by site');
xlabel('Gestational weeks');
ylabel('Placental volume (cm^3)');
legend([ph{1}(2) ph{2}(2) ph{3}(2)],{'OHSU','Utah','Both'})
if (exportFiles) 
    exportfig(gcf,[figureDirectory '/placental_volume_by_site.eps'],'Color','rgb');
else
    pause(.1);
    drawnow;
end

% linear fit to estimated placental volume vs. gestational week data by control/adverse/abnormal
figure;
[ph,categoryFit.placentalVolume] = compareModelFits(gestationalDay/7,placentalVolume,categoryMasks,line,[0 20]);
title('Placental volume by category');
xlabel('Gestational weeks');
ylabel('Placental volume (cm^3)');
legend([ph{1}(2) ph{2}(2) ph{3}(2)],{'Control','Adverse','Abnormal'})
if (exportFiles) 
    exportfig(gcf,[figureDirectory '/placental_volume_by_category.eps'],'Color','rgb');
else
    pause(.1);
    drawnow;
end

% linear fit to pre-MRI SpO2 vs. gestational week data by site
figure;
[ph,siteFit.SpO2] = compareModelFits(gestationalDay/7,SpO2Pre,siteMasks,line,[100 0]);
set(gca,'xlim',[10 40],'ylim',[75 105]);
title('Maternal SpO_2 by site');
xlabel('Gestational weeks');
ylabel('Maternal SpO_2 (%)');
legend([ph{1}(2) ph{2}(2) ph{3}(2)],{'OHSU','Utah','Both'})
if (exportFiles) 
    exportfig(gcf,[figureDirectory '/maternal_SpO2_by_site.eps'],'Color','rgb');
else
    pause(.1);
    drawnow;
end

% linear fit to pre-MRI SpO2 vs. gestational week data by control/adverse/abnormal
figure;
[ph,categoryFit.SpO2] = compareModelFits(gestationalDay/7,SpO2Pre,categoryMasks,line,[100 0]);
set(gca,'xlim',[10 40],'ylim',[75 105]);
title('Maternal SpO_2 by category');
xlabel('Gestational weeks');
ylabel('Maternal SpO_2 (%)');
legend([ph{1}(2) ph{2}(2) ph{3}(2)],{'Control','Adverse','Abnormal'})
if (exportFiles) 
    exportfig(gcf,[figureDirectory '/maternal_SpO2_by_category.eps'],'Color','rgb');
else
    pause(.1);
    drawnow;
end

% linear fit to median T1 vs. gestational week data by control/adverse/abnormal
figure;
[ph,categoryFit.T1] = compareModelFits(gestationalDay/7,medianT1,categoryMasks,line,[3000 -5]);
title('T1 by category (OHSU only)');
xlabel('Gestational weeks');
ylabel('T1 (ms)');
legend([ph{1}(2) ph{2}(2) ph{3}(2)],{'Control','Adverse','Abnormal'})
if (exportFiles) 
    exportfig(gcf,[figureDirectory '/median_T1_by_category.eps'],'Color','rgb');
else
    pause(.1);
    drawnow;
end

% linear fit to median T1 vs. median T2s by control/adverse/abnormal
figure;
[ph] = compareModelFits(medianT2s,medianT1,categoryMasks,line,[3000 .05]);
title('T1 vs. T2* (OHSU only)');
xlabel('T2* (ms)');
ylabel('T1 (ms)');
legend([ph{1}(2) ph{2}(2) ph{3}(2)],{'Control','Adverse','Abnormal'})
if (exportFiles) 
    exportfig(gcf,[figureDirectory '/median_T1_vs_median_T2s_by_category.eps'],'Color','rgb');
else
    pause(.1);
    drawnow;
end

% % birthweight percentile vs. T2s percentile (relative to control fit)
% tmp = medianT2s-(siteFit.T2s(1).model(siteFit.T2s(1).beta,gestationalDay/7));
% qq = getQuantileValue(tmp,tmp)*100;
% 
% figure;
% compareModelFits(qq,birthWeightPercentile,siteMasks,line,[0 0.01]);set(gca,'xlim',[0 100],'ylim',[0 100]);
% 
% figure;
% compareModelFits(qq,birthWeightPercentile,categoryMasks,line,[0 0.01]);set(gca,'xlim',[0 100],'ylim',[0 100]);

% estimate maternal placental blood volume fraction (v_mp) vs. gestational week data for controls from difference between OHSU and Utah
% figure;
% 
% [~,Hbfit] = compareModelFits(gestationalDay/7,hemoglobin,siteMasks,line,[12 0]);
% [~,SpO2fit] = compareModelFits(gestationalDay/7,SpO2Pre,siteMasks,line,[100 0]);
% [~,T2sfit] = compareModelFits(gestationalDay/7,medianT2s,siteMasks,sigmoid,qFit(51).beta);
% 
% gwfit = 11:1/7:36;
% % 0.6206 factor converts [Hb] in g/dl to mmol/l
% dR2s = 20.2e-3*0.6206*(Hbfit(1).model(Hbfit(1).p,gwfit).*(100-SpO2fit(1).model(SpO2fit(1).p,gwfit))/100-Hbfit(2).model(Hbfit(2).p,gwfit).*(100-SpO2fit(2).model(SpO2fit(2).p,gwfit))/100);
% 
% [~,uidx] = unique(T2sfit(1).x);
% [~,vidx] = unique(T2sfit(2).x);
% 
% plot(gwfit,(1./interp1(T2sfit(1).x(uidx),T2sfit(1).fit(uidx),gwfit,'spline')-1./interp1(T2sfit(2).x(vidx),T2sfit(2).fit(vidx),gwfit,'spline'))./dR2s,'r')

% Monte Carlo to estimate uncertainty in v_mb
gwfit = 11:1/7:36;

figure;

% OHSU vs. Utah controls
[~,Hbfit] = compareModelFits(gestationalDay/7,hemoglobin,siteMasks,line,[12 0]);
[~,SpO2fit] = compareModelFits(gestationalDay/7,SpO2Pre,siteMasks,line,[100 0]);
[~,T2sfit] = compareModelFits(gestationalDay/7,medianT2s,siteMasks,sigmoid,sigmoidGuessT2s);

% % normal vs. adverse
% [~,Hbfit] = compareModelFits(gestationalDay/7,hemoglobin,categoryMasks,line,[12 0]);
% [~,SpO2fit] = compareModelFits(gestationalDay/7,SpO2Pre,categoryMasks,line,[100 0]);
% [~,T2sfit] = compareModelFits(gestationalDay/7,medianT2s,categoryMasks,sigmoid,sigmoidGuessT2s);

% % OHSU only normal vs. adverse
% [~,Hbfit] = compareModelFits(gestationalDay/7,hemoglobin,ohsuCategoryMasks,line,[12 0]);
% [~,SpO2fit] = compareModelFits(gestationalDay/7,SpO2Pre,ohsuCategoryMasks,line,[100 0]);
% [~,T2sfit] = compareModelFits(gestationalDay/7,medianT2s,ohsuCategoryMasks,sigmoid,sigmoidGuessT2s);

% % Utah only normal vs. adverse
% [~,Hbfit] = compareModelFits(gestationalDay/7,hemoglobin,utahCategoryMasks,line,[12 0]);
% [~,SpO2fit] = compareModelFits(gestationalDay/7,SpO2Pre,utahCategoryMasks,line,[100 0]);
% [~,T2sfit] = compareModelFits(gestationalDay/7,medianT2s,utahCategoryMasks,sigmoid,sigmoidGuessT2s);

vmpb = simulateMaternalPlacentalBloodVolume(gwfit,Hbfit,SpO2fit,T2sfit,10000,1);    % estimate uncertainty using Monte Carlo and measured uncertainties

ph = plot(gwfit,nanquantile(vmpb,[.05 .25 .5 .75 .95]));
set(gca,'ylim',[0 1]);
title('Maternal placental blood volume from normals OHSU vs. Utah');
xlabel('Gestational weeks');
ylabel('v_{mpb}');
legend(ph(1:5),{'5th percentile','25th percentile','50th percentile','75th percentile','95th percentile'});
if (exportFiles) 
    exportfig(gcf,[figureDirectory '/v_mp_estimate_Monte_Carlo.eps'],'Color','rgb');
else
    pause(.1);
    drawnow;
end

% generate statistical correction curve to Utah data using site-specific T2*, SpO2, and hemoglobin fits
gw_int = 10:.1:38;

vmpb_int = simulateMaternalPlacentalBloodVolume(gw_int,Hbfit,SpO2fit,T2sfit,1,0);       % simulate for observed site-specific T2*/[Hb]/SpO2 differences

T2s_int1 = siteFit.T2s(1).model(siteFit.T2s(1).beta,gw_int);
T2s_int2 = siteFit.T2s(2).model(siteFit.T2s(2).beta,gw_int);
Hb_int1 = siteFit.maternalHb(1).model(siteFit.maternalHb(1).beta,gw_int);
Hb_int2 = siteFit.maternalHb(2).model(siteFit.maternalHb(2).beta,gw_int);
SpO2_int1 = siteFit.SpO2(1).model(siteFit.SpO2(1).beta,gw_int);
SpO2_int2 = siteFit.SpO2(2).model(siteFit.SpO2(2).beta,gw_int);

dR2s_int = 20.2e-3*0.6206*vmpb_int.*(Hb_int2.*(100-SpO2_int2)/100 - Hb_int1.*(100-SpO2_int1)/100);

dR2s = interp1(gw_int,dR2s_int,gestationalDay/7);

medianT2s_uncorrected = medianT2s;
medianT2s_corrected = medianT2s;

% Only Utah measurements are corrected
correctedMask = false(size(category));
correctedMask(siteNumber==2) = true;

medianT2s_corrected(correctedMask) = 1./(1./medianT2s(correctedMask)-dR2s(correctedMask));

% stem plot showing change in T2* values when correction is applied
figure;
ph = plot(gestationalDay/7,medianT2s_uncorrected,'b.',[gestationalDay(correctedMask)/7;gestationalDay(correctedMask)/7],[medianT2s_uncorrected(correctedMask);medianT2s_corrected(correctedMask)],'r-');
set(gca,'ylim',[0 120]);
title('T2* correction');
xlabel('Gestational weeks');
ylabel('T2* (ms)');
legend([ph(1) ph(2)],{'uncorrected','corrected'});
if (exportFiles) 
    exportfig(gcf,[figureDirectory '/median_T2s_correction.eps'],'Color','rgb');
else
    pause(.1);
    drawnow;
end

% apply T2s correction for maternal SpO2 and [Hb] and fit
figure;
[~,siteFit.T2s_uncorr] = compareModelFits(gestationalDay/7,medianT2s_uncorrected,siteMasks(1:2),sigmoid,sigmoidGuessT2s);
[ph,siteFit.T2s_corr] = compareModelFits(gestationalDay/7,medianT2s_corrected,siteMasks(1:2),sigmoid,sigmoidGuessT2s);
set(gca,'ylim',[0 120]);
title('T2* by site, corrected');
xlabel('Gestational weeks');
ylabel('T2* (ms)');
legend([ph{1}(2) ph{2}(2)],{'OHSU','Utah, corrected'})
if (exportFiles) 
    exportfig(gcf,[figureDirectory '/median_T2s_by_site_corrected.eps'],'Color','rgb');
else
    pause(.1);
    drawnow;
end

figure;
[~,categoryFit.T2s_uncorr] = compareModelFits(gestationalDay/7,medianT2s_uncorrected,categoryMasks,sigmoid,sigmoidGuessT2s);
[ph,categoryFit.T2s_corr] = compareModelFits(gestationalDay/7,medianT2s_corrected,categoryMasks,sigmoid,sigmoidGuessT2s);
set(gca,'ylim',[0 120]);
title('T2* by category, corrected');
xlabel('Gestational weeks');
ylabel('T2* (ms)');
legend([ph{1}(2) ph{2}(2) ph{3}(2)],{'Control, corrected','Adverse, corrected','Abnormal, corrected'})
if (exportFiles) 
    exportfig(gcf,[figureDirectory '/median_T2s_by_category_corrected.eps'],'Color','rgb');
else
    pause(.1);
    drawnow;
end
% 
% correctedCategoryMasks = {categoryMasks{1} & correctedMask,categoryMasks{2} & correctedMask,categoryMasks{3} & correctedMask};
% 
% figure;
% [~,categoryFit.T2s_uncorr] = compareModelFits(gestationalDay/7,medianT2s_uncorrected,correctedCategoryMasks,sigmoid,sigmoidGuessT2s);
% [ph,categoryFit.T2s_corr] = compareModelFits(gestationalDay/7,medianT2s_corrected,correctedCategoryMasks,sigmoid,sigmoidGuessT2s);
% set(gca,'ylim',[0 120]);
% title('T2* by category, corrected');
% xlabel('Gestational weeks');
% ylabel('T2* (ms)');
% legend([ph{1}(2) ph{2}(2) ph{3}(2)],{'Control, corrected','Adverse, corrected','Abnormal, corrected'})
% if (exportFigures) 
%     exportfig(gcf,[figureDirectory '/median_T2s_by_category_corrected.eps'],'Color','rgb');
% end

figure;
ph = plot(siteFit.T2s_uncorr(1).x,siteFit.T2s_uncorr(1).fit,'b--',siteFit.T2s_corr(1).x,siteFit.T2s_corr(1).fit,'b',...
          siteFit.T2s_uncorr(2).x,siteFit.T2s_uncorr(2).fit,'r--',siteFit.T2s_corr(2).x,siteFit.T2s_corr(2).fit,'r');
set(gca,'ylim',[0 120]);
title('T2* by site, uncorrected vs. corrected');
xlabel('Gestational weeks');
ylabel('T2* (ms)');
legend(ph,{'OHSU, uncorrected','OHSU, corrected',...
           'Utah, uncorrected','Utah, corrected'})
if (exportFiles) 
    exportfig(gcf,[figureDirectory '/median_T2s_by_site_corrected_vs_uncorrected.eps'],'Color','rgb');
else
    pause(.1);
    drawnow;
end

figure;
ph = plot(categoryFit.T2s_uncorr(1).x,categoryFit.T2s_uncorr(1).fit,'b--',categoryFit.T2s_corr(1).x,categoryFit.T2s_corr(1).fit,'b',...
          categoryFit.T2s_uncorr(2).x,categoryFit.T2s_uncorr(2).fit,'r--',categoryFit.T2s_corr(2).x,categoryFit.T2s_corr(2).fit,'r',...
          categoryFit.T2s_uncorr(3).x,categoryFit.T2s_uncorr(3).fit,'g--',categoryFit.T2s_corr(3).x,categoryFit.T2s_corr(3).fit,'g');
set(gca,'ylim',[0 120]);
title('T2* by category, uncorrected vs. corrected');
xlabel('Gestational weeks');
ylabel('T2* (ms)');
legend(ph,{'Control, uncorrected','Control, corrected',...
           'Adverse, uncorrected','Adverse, corrected',...
           'Abnormal, uncorrected','Abnormal, corrected'})
if (exportFiles) 
    exportfig(gcf,[figureDirectory '/median_T2s_by_category_corrected_vs_uncorrected.eps'],'Color','rgb');
else
    pause(.1);
    drawnow;
end

% % generate correction curve for calibration between OHSU and Utah data using actual measured Hb and SpO2 values
% vmpb0 = simulateMaternalPlacentalBloodVolume(gwfit,Hbfit,SpO2fit,T2sfit,1,0);       % simulate for observed site-specific T2*/[Hb]/SpO2 differences
% 
% % interpolate to gestational week
% for gidx=1:length(medianT2s)
%     gw = gestationalDay(gidx)/7;    
%     v = interp1(gwfit,vmpb0,gw);
%     
%     dR2s(gidx) = v*20.2e-3*0.6206*(hemoglobin(gidx)*(100-SpO2Pre(gidx))/100 - Hbfit(1).f(Hbfit(1).beta,gw)*(100-SpO2fit(1).f(SpO2fit(1).beta,gw))/100);
%     dR2sref(gidx) = v*20.2e-3*0.6206*(Hbfit(2).f(Hbfit(2).beta,gw)*(100-SpO2fit(2).f(SpO2fit(2).beta,gw))/100 - Hbfit(1).f(Hbfit(1).beta,gw)*(100-SpO2fit(1).f(SpO2fit(1).beta,gw))/100);
%     
%     T2s1(gidx) = T2sfit(1).f(T2sfit(1).beta,gw);
%     T2s2(gidx) = T2sfit(2).f(T2sfit(2).beta,gw);
% end
% 
% medianT2s_corrected = 1./(1./medianT2s-dR2s);
% medianT2s_uncorrected = 1./(1./medianT2s-0*dR2s);   % would be the same as medianT2s except for availability of SpO2 and/or hemoglobin measurements
% 
% % apply T2s correction for maternal SpO2 and [Hb] and fit
% [~,siteFit.T2s_corr] = compareModelFits(gestationalDay/7,medianT2s_corrected,siteMasks(1:2),sigmoid,sigmoidGuessT2s);
% [~,categoryFit.T2s_corr] = compareModelFits(gestationalDay/7,medianT2s_corrected,categoryMasks,sigmoid,sigmoidGuessT2s);
% 
% [~,siteFit.T2s_uncorr] = compareModelFits(gestationalDay/7,medianT2s_uncorrected,siteMasks(1:2),sigmoid,sigmoidGuessT2s);
% [~,categoryFit.T2s_uncorr] = compareModelFits(gestationalDay/7,medianT2s_uncorrected,categoryMasks,sigmoid,sigmoidGuessT2s);
% 
% figure;
% ph = plot(siteFit.T2s_uncorr(1).x,siteFit.T2s_uncorr(1).fit,'b--',siteFit.T2s_corr(1).x,siteFit.T2s_corr(1).fit,'b',...
%           siteFit.T2s_uncorr(2).x,siteFit.T2s_uncorr(2).fit,'r--',siteFit.T2s_corr(2).x,siteFit.T2s_corr(2).fit,'r');
% set(gca,'ylim',[0 120]);
% title('T2* by site, uncorrected vs. corrected ');
% xlabel('Gestational weeks');
% ylabel('T2* (ms)');
% legend(ph,{'OHSU, uncorrected','OHSU, corrected',...
%            'Utah, uncorrected','Utah, corrected'})
% exportfig(gcf,[figureDirectory '/median_T2s_by_site_corrected_vs_uncorrected.eps'],'Color','rgb');
% 
% figure;
% ph = plot(categoryFit.T2s_uncorr(1).x,categoryFit.T2s_uncorr(1).fit,'b--',categoryFit.T2s_corr(1).x,categoryFit.T2s_corr(1).fit,'b',...
%           categoryFit.T2s_uncorr(2).x,categoryFit.T2s_uncorr(2).fit,'r--',categoryFit.T2s_corr(2).x,categoryFit.T2s_corr(2).fit,'r',...
%           categoryFit.T2s_uncorr(3).x,categoryFit.T2s_uncorr(3).fit,'g--',categoryFit.T2s_corr(3).x,categoryFit.T2s_corr(3).fit,'g');
% set(gca,'ylim',[0 120]);
% title('T2* by category, uncorrected vs. corrected');
% xlabel('Gestational weeks');
% ylabel('T2* (ms)');
% legend(ph,{'Control, uncorrected','Control, corrected',...
%            'Adverse, uncorrected','Adverse, corrected',...
%            'Abnormal, uncorrected','Abnormal, corrected'})
% exportfig(gcf,[figureDirectory '/median_T2s_by_category_corrected_vs_uncorrected.eps'],'Color','rgb');

% fit change in T2s with gestation within individual subjects
siteDeltaMasks = { repmat(siteMasks{1},[2 1])', repmat(siteMasks{2},[2 1])', repmat(siteMasks{3},[2 1])' };
categoryDeltaMasks = { repmat(categoryMasks{1},[2 1])', repmat(categoryMasks{2},[2 1])', repmat(categoryMasks{3},[2 1])' };

% figure;
% [ph,categoryFit.deltaT2sdGD] = compareModelFits(medianGD(:),deltaT2s(:)./deltaGD(:),{categoryDeltaMasks{1}(:),categoryDeltaMasks{2}(:),categoryDeltaMasks{3}(:)},line,[-1 0]);
% set(gca,'xlim',[100 220],'ylim',[-1.5 .5]);
% title('Intrapatient T2* gestational rate of change - normal vs. adverse');
% xlabel('Gestational days');
% ylabel('dT2*/dGD (ms/day)');
% legend([ph{1}(2) ph{2}(2) ph{3}(2)],{'Normal','Adverse','Abnormal'})
% if (exportFigures) 
%     exportfig(gcf,[figureDirectory '/dT2sdGD_by_category.eps'],'Color','rgb');
% end
% 
% figure;
% [ph,siteFit.deltaT2sdGD] = compareModelFits(medianGD(:),deltaT2s(:)./deltaGD(:),{siteDeltaMasks{1}(:),siteDeltaMasks{2}(:), siteDeltaMasks{3}(:)},line,[-1 0]);
% set(gca,'xlim',[100 220],'ylim',[-1.5 .5]);
% title('Intrapatient T2* gestational rate of change - normal pregnancies');
% xlabel('Gestational days');
% ylabel('dT2*/dGD (ms/day)');
% legend([ph{1}(2) ph{2}(2) ph{3}(2)],{'OHSU','Utah','Both'})
% if (exportFigures) 
%     exportfig(gcf,[figureDirectory '/dT2sdGD_by_site.eps'],'Color','rgb');
% end

% % fit with line
% figure;
% [ph,siteFit.deltaT2sdGW] = compareModelFits(medianGD(:)/7,deltaT2s(:)./(deltaGD(:)/7),{siteDeltaMasks{1}(:),siteDeltaMasks{2}(:), siteDeltaMasks{3}(:)},line,[-1 0]);
% set(gca,'xlim',[15 35],'ylim',[-10 5]);
% title('Intrapatient T2* gestational rate of change - normal pregnancies');
% xlabel('Gestational weeks');
% ylabel('dT2*/dGW (ms/week)');
% legend([ph{1}(2) ph{2}(2) ph{3}(2)],{'OHSU','Utah','Both'})
% if (exportFigures) 
%     exportfig(gcf,[figureDirectory '/dT2sdGW_by_site.eps'],'Color','rgb');
% end

% fit with derivative of sigmoid
figure;
[ph,siteFit.deltaT2sdGW] = compareModelFits(medianGD(:)/7,deltaT2s(:)./(deltaGD(:)/7),{siteDeltaMasks{1}(:),siteDeltaMasks{2}(:),siteDeltaMasks{3}(:)},sigmoid_derivative,sigmoidGuessT2s(1:3));
set(gca,'xlim',[15 35],'ylim',[-10 5]);
title('Intrapatient T2* gestational rate of change - normal pregnancies');
xlabel('Gestational weeks');
ylabel('dT2*/dGW (ms/week)');
legend([ph{1}(2) ph{2}(2) ph{3}(2)],{'OHSU','Utah','Both'})
if (exportFiles) 
    exportfig(gcf,[figureDirectory '/dT2sdGW_by_site.eps'],'Color','rgb');
else
    pause(.1);
    drawnow;
end

% % fit with line
% figure;
% [ph,categoryFit.deltaT2sdGW] = compareModelFits(medianGD(:)/7,deltaT2s(:)./(deltaGD(:)/7),{categoryDeltaMasks{1}(:),categoryDeltaMasks{2}(:),categoryDeltaMasks{3}(:)},line,[-1 0]);
% set(gca,'xlim',[15 35],'ylim',[-10 5]);
% title('Intrapatient T2* gestational rate of change - normal vs. adverse');
% xlabel('Gestational weeks');
% ylabel('dT2*/dGW (ms/week)');
% legend([ph{1}(2) ph{2}(2) ph{3}(2)],{'Normal','Adverse','Abnormal'})
% if (exportFigures) 
%     exportfig(gcf,[figureDirectory '/dT2sdGW_by_category.eps'],'Color','rgb');
% end

% fit with derivative of sigmoid
figure;
[ph,categoryFit.deltaT2sdGW] = compareModelFits(medianGD(:)/7,deltaT2s(:)./(deltaGD(:)/7),{categoryDeltaMasks{1}(:),categoryDeltaMasks{2}(:),categoryDeltaMasks{3}(:)},sigmoid_derivative,sigmoidGuessT2s(1:3));
set(gca,'xlim',[15 35],'ylim',[-10 5]);
title('Intrapatient T2* gestational rate of change - normal vs. adverse');
xlabel('Gestational weeks');
ylabel('dT2*/dGW (ms/week)');
legend([ph{1}(2) ph{2}(2) ph{3}(2)],{'Normal','Adverse','Abnormal'})
if (exportFiles) 
    exportfig(gcf,[figureDirectory '/dT2sdGW_by_category.eps'],'Color','rgb');
else
    pause(.1);
    drawnow;
end

gestationalWeek = floor(gestationalDay/7);
gestationalTP = floor(gestationalDay/7/10);

% convert T2* z-scores to T2* percentiles
T2sPercentile = NaNs(size(medianT2s));

for i=1:length(studyIDs) 
    try
        refT2s = medianT2s(gestationalWeek==gestationalWeek(i)&category==0);
        T2sPercentile(i) = getQuantileValue(refT2s,medianT2s(i))*100;
    catch e
    end
    
    if (category(i)~=1)
        continue;
    end
    
    fprintf('%3d %8i Category %2d %5.1f %3.0f\n',studyIDs(i),studyDates(i),category(i),medianT2s(i),T2sPercentile(i)); 
end

% organize T2* by studyID
timepointData = struct();

[timepointData.uniqueStudyIDs,timepointData.uidx,timepointData.vidx] = unique(studyIDs);

timepointData.siteNumber = NaNs([1 length(timepointData.uniqueStudyIDs)]);
timepointData.category = NaNs([1 length(timepointData.uniqueStudyIDs)]);

timepointData.gestationalDay = NaNs([length(timepointData.uniqueStudyIDs) 3]);
timepointData.medianT2s = NaNs([length(timepointData.uniqueStudyIDs) 3]);

for idx=1:length(studyIDs) 
    sidx = timepointData.vidx(idx);
    tidx = gestationalTP(idx);
    
    timepointData.siteNumber(sidx) = siteNumber(idx);
    timepointData.category(sidx) = category(idx);
    
    timepointData.gestationalDay(sidx,tidx) = gestationalDay(idx);
    timepointData.medianT2s(sidx,tidx) = medianT2s(idx);
end

% plot(timepointData.gestationalDay(timepointData.siteNumber==1&timepointData.category==0,:)'/7,...
%      timepointData.medianT2s(timepointData.siteNumber==1&timepointData.category==0,:)','b-',...
%      timepointData.gestationalDay(timepointData.siteNumber==2&timepointData.category==0,:)'/7,...
%      timepointData.medianT2s(timepointData.siteNumber==2&timepointData.category==0,:)','r-')

medianT2s_for_ROC = medianT2s;
% medianT2s_for_ROC = medianT2s_corrected;

% ROC for predicting adverse pregnancy
rocDescription = 'Adverse only';
rocMasks = {controlMask,adverseMask};
ohsuRocMasks = {controlMask & ohsuMask,adverseMask & ohsuMask};
utahRocMasks = {controlMask & utahMask,adverseMask & utahMask};
% ohsuRocMasks = {controlMask,adverseMask & ohsuMask};
% utahRocMasks = {controlMask,adverseMask & utahMask};

% % ROC for predicting adverse pregnancy with low birthweight
% rocDescription = 'Adverse and low birth weight';
% rocMasks = {controlMask,adverseMask & lowBirthWeightMask};
% ohsuRocMasks = {controlMask & ohsuMask,adverseMask & ohsuMask & lowBirthWeightMask};
% utahRocMasks = {controlMask & utahMask,adverseMask & utahMask & lowBirthWeightMask};

% % ROC for predicting adverse or abnormal pregnancy
% rocDescription = 'Adverse or abnormal';
% rocMasks = {controlMask,adverseMask | abnormalMask};
% ohsuRocMasks = {controlMask & ohsuMask,(adverseMask | abnormalMask) & ohsuMask};
% utahRocMasks = {controlMask & utahMask,(adverseMask | abnormalMask) & utahMask};

% % ROC for predicting low birth weight (below 10th percentile)
% rocDescription = 'Low birthweight including IUGR';
% rocMasks = {controlMask & birthWeightPercentile>=10,preMask & birthWeightPercentile<10};
% ohsuRocMasks = {ohsuMask & birthWeightPercentile>=10,ohsuMask & birthWeightPercentile<10};
% utahRocMasks = {utahMask & birthWeightPercentile>=10,utahMask & birthWeightPercentile<10};

for gidx=1:3
    figure;
    switch gidx
        % generate ROC curves using prediction intervals from sigmoid fit to T2* control data
        % OHSU only ROC analysis
        case 1,     [~,out] = compareModelFits(gestationalDay/7,medianT2s_for_ROC,ohsuRocMasks,sigmoid,sigmoidGuessT2s);
                    siteName = 'OHSU';
        % Utah only ROC analysis
        case 2,     [~,out] = compareModelFits(gestationalDay/7,medianT2s_for_ROC,utahRocMasks,sigmoid,sigmoidGuessT2s);
                    siteName = 'Utah';
        % both OHSU and Utah
        case 3,     [~,out] = compareModelFits(gestationalDay/7,medianT2s_for_ROC,rocMasks,sigmoid,sigmoidGuessT2s);
                    siteName = 'Both';
        otherwise,  error();
    end
        
    for gtpidx=0:3
        if (gtpidx == 0)
            roc = generateROC(out,true(size(out(2).x)));
        else            
            roc = generateROC(out,floor(out(2).x/10)==gtpidx);
        end
        
        if (gtpidx==0)
            suffix = '_all';
        else
            suffix = ['_GTP' num2str(gtpidx)];
        end
            
        figure;
        plot(roc.FPR,roc.TPR,[0 1],[0 1],'r',roc.FPR(roc.Jindex),roc.TPR(roc.Jindex),'r*');
        axis square;
        set(gca,'xlim',[0 1],'ylim',[0 1]);
        title([siteName ' ' rocDescription ' ' suffix(2:end) ...
              ' AUC : ' num2str(roc.AUC) ...
              ' J : ' num2str(roc.Jmax) ...
              ' Opt. pct. : ' num2str(roc.optimalPercentile)]);
        
        pause(.1);

        if (exportFiles)
            switch gidx
                case 1, exportfig(gcf,[figureDirectory '/ROC' suffix '_OHSU.eps'],'Color','rgb');
                case 2, exportfig(gcf,[figureDirectory '/ROC' suffix '_Utah.eps'],'Color','rgb');
                case 3, exportfig(gcf,[figureDirectory '/ROC' suffix '.eps'],'Color','rgb');
                otherwise error();
            end
        else
            pause(.1);
            drawnow;
        end
    end
end

gwbins = 10.5:39.5;
fgw = floor(gestationalDay/7);

% plot histograms of study gestational ages by site
h1 = hist(fgw(ohsuMask),gwbins);
h2 = hist(fgw(utahMask),gwbins);
h3 = hist(fgw(ohsuMask | utahMask),gwbins);

figure;
ph = plot(gwbins,h1,'b.-',gwbins,h2,'r.-',gwbins,h3,'k.-');
set(gca,'xlim',[10 40]);
title('Gestational weeks at MRI by site');
xlabel('Gestational weeks at MRI');
ylabel('Number of studies');
legend(ph,{'OHSU','Utah','Both'});
if (exportFiles) 
    exportfig(gcf,[figureDirectory '/MRI_enrollment_by_site.eps'],'Color','rgb');
else
    pause(.1);
    drawnow;
end

% plot histograms of study gestational ages by category
h1 = hist(fgw(controlMask),gwbins);
h2 = hist(fgw(adverseMask),gwbins);
h3 = hist(fgw(abnormalMask),gwbins);

figure;
ph = plot(gwbins,h1,'b.-',gwbins,h2,'r.-',gwbins,h3,'g.-');
set(gca,'xlim',[10 40]);
title('Gestational weeks at MRI by category');
xlabel('Gestational weeks at MRI');
ylabel('Number of studies');
legend(ph,{'Normal','Adverse','Abnormal'});
if (exportFiles) 
    exportfig(gcf,[figureDirectory '/MRI_enrollment_by_category.eps'],'Color','rgb');
else
    pause(.1);
    drawnow;
end

% plot z-score for T2* data
zscore = [];

for idx=1:2
    if (idx==1)
        ff = categoryFit.T2s;
        suffix = '';
    else
        ff = categoryFit.T2s_corr;
        suffix = '_corrected';
    end
    
    [~,uidx] = unique(ff(1).x);
    
    zscore{1} = (ff(1).y-ff(1).fit)./(0.25*(ff(1).predhi-ff(1).predlo));
    zscore{2} = (ff(2).y-ff(1).f(ff(1).beta,ff(2).x))./interp1(ff(1).x(uidx),0.25*(ff(1).predhi(uidx)-ff(1).predlo(uidx)),ff(2).x);
    zscore{3} = (ff(3).y-ff(1).f(ff(1).beta,ff(3).x))./interp1(ff(1).x(uidx),0.25*(ff(1).predhi(uidx)-ff(1).predlo(uidx)),ff(3).x);
    
    for jdx=1:3
        if (idx==1)
            categoryFit.T2s(jdx).zscore = zscore{jdx};
        else
            categoryFit.T2s_corr(jdx).zscore = zscore{jdx};
        end
    end
    
    figure;
    ph = plot(ff(1).x,zscore{1},'b.',ff(2).x,zscore{2},'r.',ff(3).x,zscore{3},'g.');
    set(gca,'ylim',[-6 6]);
    title(['T2* z-score by category ' suffix(2:end)]);
    xlabel('Gestational weeks at MRI');
    ylabel('z-score');
    legend(ph,{'Normal','Adverse','Abnormal'});
    if (exportFiles) 
        exportfig(gcf,[figureDirectory '/T2s_zscore_by_category' suffix '.eps'],'Color','rgb');
    else
        pause(.1);
        drawnow;
    end
    
    % z-score histograms by category
    figure;
    zbins=-6:.5:6;
    h1=normalize(hist(zscore{1},zbins));
    h2=normalize(hist(zscore{2},zbins));
    h3=normalize(hist(zscore{3},zbins));
    ph = plot(zbins,h1,'b',zbins,h2,'r',zbins,h3,'g');
    title(['T2* z-score by category ' suffix(2:end)]);
    xlabel('z-score');
    ylabel('relative frequency');
    legend(ph,{'Normal','Adverse','Abnormal'});
    if (exportFiles) 
        exportfig(gcf,[figureDirectory '/zscore_histogram_by_category' suffix '.eps'],'Color','rgb');
    else
        pause(.1);
        drawnow;
    end
end

% combined z-score
zscore = NaN(size(medianT2s));

for idx=1:3
    zscore(categoryMasks{idx}) = categoryFit.T2s(idx).zscore;
end

% compute T2* percentile from z-score
normalDistribution = makedist('Normal',0,1);

T2sPercentile = cdf(normalDistribution,zscore)*100;

% histograms of T2* percentiles (should be uniform in ideal case)
T2sBins = 2.5:5:97.5;
T2sHistogram = [];

for idx=1:3 
    T2sHistogram(idx,:) = normalize(hist(T2sPercentile(categoryMasks{idx}),T2sBins)); 
    T2sHistogram(3+idx,:) = normalize(hist(T2sPercentile(ohsuCategoryMasks{idx}),T2sBins)); 
    T2sHistogram(6+idx,:) = normalize(hist(T2sPercentile(utahCategoryMasks{idx}),T2sBins)); 
end

figure;
ph = bar(T2sBins,T2sHistogram(1:3,:),'stacked');
ph(1).FaceColor=[0 0 1];
ph(2).FaceColor=[1 0 0];
ph(3).FaceColor=[0 1 0];
title('Relative distribution of T2* percentiles by category');
xlabel('Percentile');
ylabel('Relative frequency');
legend(ph,{'Normal','Adverse','Abnormal'});
if (exportFiles)
    exportfig(gcf,[figureDirectory '/T2s_percentile_histogram_by_category.eps'],'Color','rgb');
else
    pause(.1);
    drawnow;
end

% histograms of birthweight percentiles
birthweightBins = 5:10:95;

for sidx=1:3
    switch sidx
        case 1,     masks = categoryMasks; 
                    titleSuffix = '';
                    fileSuffix = '';
        case 2,     masks = ohsuCategoryMasks; 
                    titleSuffix = ' (OHSU)';
                    fileSuffix = '_OHSU';
        case 3,     masks = utahCategoryMasks; 
                    titleSuffix = ' (Utah)';
                    fileSuffix = '_Utah';
    end
    
    birthWeightHistogram = [];

    for idx=1:3
        birthWeightHistogram(idx,:) = normalize(hist(birthWeightPercentile(masks{idx}&scanIndex==1),birthweightBins));
    end

    figure;
    ph = bar(birthweightBins,birthWeightHistogram,'stacked');
    ph(1).FaceColor=[0 0 1];
    ph(2).FaceColor=[1 0 0];
    ph(3).FaceColor=[0 1 0];
    title(['Relative distribution of birthweight percentiles by category' titleSuffix]);
    xlabel('Percentile');
    ylabel('Relative frequency');
    legend(ph,{'Normal','Adverse','Abnormal'});
    if (exportFiles)
        exportfig(gcf,[figureDirectory '/birthweight_percentile_histogram_by_category' fileSuffix '.eps'],'Color','rgb');
    else
        pause(.1);
        drawnow;
    end
end

% generate summary output for statistical analysis
shift = transpose(shifted_lognormal_distfit_T2s(:,1));
mu = transpose(shifted_lognormal_distfit_T2s(:,2));
sigma = transpose(shifted_lognormal_distfit_T2s(:,3));

outputTable = table;

outputTable.studyID = studyIDs;
outputTable.studyDate = studyDates;
outputTable.siteNumber = siteNumber;
outputTable.gestationalDay = gestationalDay;
outputTable.gestationalWeek = gestationalWeek;
outputTable.gestationalTP = gestationalTP;
outputTable.patientSize = patientSize;
outputTable.patientWeight = patientWeight;
outputTable.maternalAge = maternalAge;
outputTable.gestationalAgeAtDelivery = gad;
outputTable.placentalWeightAtDelivery = placentalWeight_g;
outputTable.prepregnancyBMI = prepregnancyBMI;
outputTable.deliveryBMI = deliveryBMI;
outputTable.birthWeight = birthWeight;
outputTable.birthWeightPercentile = birthWeightPercentile;
outputTable.fetalSex = fetalSex;
outputTable.hemoglobin = hemoglobin;
outputTable.hemoglobinFingerstick = hemoglobinFingerstick;
outputTable.SpO2Pre = SpO2Pre;
outputTable.heartRatePre = heartRatePre;
outputTable.SpO2Post = SpO2Post;
outputTable.heartRatePost = heartRatePost;
outputTable.systolicBloodPressure = systolicBloodPressure;
outputTable.diastolicBloodPressure = diastolicBloodPressure;
outputTable.adjudicationGroup = adjudicationGroup;
outputTable.adjudicationSeverity = adjudicationSeverity;
outputTable.category = category;
outputTable.scanIndex = scanIndex;
outputTable.medianT2s = medianT2s;
outputTable.medianR2s = medianR2s;
outputTable.SLNFitShift = shift;
outputTable.SLNFitMu = mu;
outputTable.SLNFitSigma = sigma;
outputTable.SLNFitMedian = exp(mu);
outputTable.SLNFitMean = exp(mu+sigma.^2/2);
outputTable.SLNFitVariance = (exp(sigma.^2)-1).*exp(2*mu+sigma.^2);
outputTable.T2sPercentile = T2sPercentile;
outputTable.uncertaintyMetricT2s = uncertaintyMetricT2s;
outputTable.qualityMetricROI = qualityMetricROI;
outputTable.medianT1 = medianT1;
outputTable.zscoreT2s = zscore;
outputTable.normalMask = categoryMasks{1};
outputTable.adverseMask = categoryMasks{2};
outputTable.abnormalMask = categoryMasks{3};
outputTable.ohsuMask = ohsuMask;
outputTable.utahMask = utahMask;

outputStruct = table2struct(outputTable);

outputFields = fieldnames(outputStruct);

if (exportFiles)
    save([rootDirectory '/U01_summaryData.mat'],'outputTable','outputStruct','outputFields');
end

% need to transpose all fields to write to xlsx
f = fieldnames(outputStruct);

for i=1:length(f)
    outputStruct.(f{i}) = transpose(outputStruct.(f{i}));
end

if (exportFiles)
    writetable(struct2table(outputStruct),[rootDirectory '/U01_summaryData.xlsx']);
end

% demographics
demographics=struct();

for idx=1:length(redcapRecords) 
    demographics(idx).preterm = str2num(redcapRecords(idx).gestational_age_at_delivery)/7<37; 
    demographics(idx).gravidity = str2num(redcapRecords(idx).gravidity);
    demographics(idx).pregnancy_complications = cell2mat(redcapRecords(idx).pregnancy_complications);
    demographics(idx).racial_background = cell2mat(redcapRecords(idx).racial_background);
    demographics(idx).hispanic_origin = cell2mat(redcapRecords(idx).hispanic_origin);
    demographics(idx).pregnancy_complications = cell2mat(redcapRecords(idx).pregnancy_complications);
    demographics(idx).stillbirth_pregnancy_loss = cell2mat(redcapRecords(idx).stillbirth_pregnancy_loss);
    demographics(idx).hypertensive_diagnosis = cell2mat(redcapRecords(idx).hypertensive_diagnosis);
    demographics(idx).other_outcomes = cell2mat(redcapRecords(idx).other_outcomes);
    demographics(idx).sga = birthWeightPercentile(idx)<10;
end

fn=fieldnames(demographics); 

for idx=1:length(demographics) 
    for jdx=1:length(fn) 
        cf=demographics(idx).(fn{jdx});
        
        if (~ischar(cf) & isempty(cf)) 
            demographics(idx).(fn{jdx}) = 0; 
        end
    end
end

msk = categoryMasks{2} & scanIndex==1;
% msk = ohsuCategoryMasks{2} & scanIndex==1;
% msk = utahCategoryMasks{2} & scanIndex==1;

spl = int32([demographics(msk).stillbirth_pregnancy_loss])-48;
oth = int32([demographics(msk).other_outcomes])-48;
htn = int32([demographics(msk).hypertensive_diagnosis])-48;
ptm = [demographics(msk).preterm];

demographicData = struct('PIH',sum(sum(htn([2 3 4 5 7],:)>0)),...
                         'GHTN',sum(htn(7,:),2),...
                         'PEWO',sum(htn(4,:),2),...
                         'PEW',sum(sum(htn([2 3 5],:))),...
                         'SGA',sum([demographics(msk).sga]),...
                         'SBFL',sum(spl(:)),...
                         'PA',sum(oth(2,:)),...
                         'PIHSGA',sum(sum(htn([2 3 4 5 7],:))>0&[demographics(msk).sga]),...
                         'PB',sum(ptm));

valpct = @(x)[num2str(x) ' (' num2str(100*x/sum(msk)) ')'];
mnsd = @(x)[nanmean(x) nanstd(x)];

disp(['PIH : ' valpct(demographicData.PIH)]);
disp(['Gestational HTN : ' valpct(demographicData.GHTN)]);
disp(['Pre-e w/o sev : ' valpct(demographicData.PEWO)]);
disp(['Pre-e w/sev : ' valpct(demographicData.PEW)]);
disp(['SGA : ' valpct(demographicData.SGA)]);
disp(['Stillbirth or fetal loss : ' valpct(demographicData.SBFL)]);
disp(['Placental abruption : ' valpct(demographicData.PA)]);
disp(['Both PIH and SGA : ' valpct(demographicData.PIHSGA)]);
disp(['Preterm birth <37 weeks : ' valpct(demographicData.PB)]);

% compute p-values using chi-square test from prop_test.m

% racial background
racial_backgrounds = {'White','African Descent','Native American','Asian Indian','Other Asian','Native Hawaiian','Pacific Islander','Other','Refuse','Hispanic'}

for idx=1:3
    switch idx
        case 1, disp('Uncomplicated Normal');
        case 2, disp('Primary Adverse');
        case 3, disp('Secondary Abnormal');
    end
    
    for jdx=1:3
        switch jdx
            case 1,     msk = categoryMasks{idx} & scanIndex==1;
            case 2,     msk = ohsuCategoryMasks{idx} & scanIndex==1;
            case 3,     msk = utahCategoryMasks{idx} & scanIndex==1;
        end
    
        rb = int32([demographics(msk).racial_background])-48;
        ho = int32([demographics(msk).hispanic_origin])-48;
    
        disp([racial_backgrounds; num2cell(sum(rb')); num2cell(100*sum(rb')/sum(rb(:)))]');
        
%         disp(num2str(mnsd(maternalAge(msk))));
%         disp(num2str(mnsd(prepregnancyBMI(msk))));
                
%         gr = int32([demographics(msk).gravidity]);
%         
%         for kdx=0:4
%             if (kdx<4)
%                 disp(num2str([kdx sum(gr==kdx)]));
%             else
%                 disp(num2str([kdx sum(gr>=kdx)]));
%             end
%         end
    end
end

% % push postprocessed values to REDCap
% % for idx=794:794
% for idx=1:length(studyIDs)
%     switch siteNumber(idx)
%         case 1, rcstudyID = ['PIPO-' sprintf('%03i',studyIDs(idx))];
%         case 2, rcstudyID = ['PIPU' sprintf('%03i',studyIDs(idx))];
%         case 3, rcstudyID = ['BCNatal-' sprintf('%03i',studyIDs(idx))];
%     end
%         
%     rcstudyDate = num2str(studyDates(idx));
%     rcplacentalVolume = num2str(placentalVolume(idx));
%     rcT2s = num2str(medianT2s(idx));
%     
%     % MacOS
%     [status,systemOutput] = system(['/usr/local/bin/python3 U01_redcap_update.py ' rcstudyID ' ' rcstudyDate ' ' rcplacentalVolume ' NaN NaN NaN ' rcT2s]);        
% %     disp(['/usr/local/bin/python3 U01_redcap_update.py ' rcstudyID ' ' rcstudyDate ' ' rcplacentalVolume ' NaN NaN NaN ' rcT2s]);        
% 
%     if (status ~= 0)
%         disp(['Problem running U01_redcap_update.py for ' rcstudyID ' ' rcstudyDate]);
%         disp(systemOutput);
%     else
%         disp(['U01_redcap_update.py for ' rcstudyID ' ' rcstudyDate ' complete']);
%     end
% end
