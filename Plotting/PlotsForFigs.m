% PenCorrelationLength.m 2_25 updates

clear; clc; close all;
load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Pipeline\7_8_23_Analysis\CombinedCorrmat.mat"); %Load combined results (Corrmat)
load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Pipeline\7_8_23_Analysis\CombinedKFmat.mat"); %Load combined results (KFmat)

%% Table 1
KFtodel = ~logical(CombinedKFmat(:,6));
CombinedKFmat(KFtodel,:) = [];

wts = 1./(CombinedKFmat(:,3).^2);
kvec = abs(CombinedKFmat(:,2));
fvec = CombinedKFmat(:,1);

%variance-weighted k_vaso
% r2 = pvcomb(:,2).^2;
kvarnum = dot(wts,kvec);
kvardenom = sum(wts);
kvarmean = kvarnum/kvardenom

r2wtd2ndmoment = dot(wts,kvec.^2)/sum(wts);
r2wtduncert = sqrt(r2wtd2ndmoment - kvarmean^2)

%% FIG 3: f vs k scatter and marginal histograms
cd('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Pipeline\7_8_23_Analysis\Fig3')
KFtodel = ~logical(CombinedKFmat(:,6));
CombinedKFmat(KFtodel,:) = []; %Delete PAs with insignificant vasomotor coherence

%F vs K scatter
ispos = CombinedKFmat(:,2) > 0;
posplot = CombinedKFmat(ispos,:);
negplot = CombinedKFmat(~ispos,:);
f1 = figure;
scatter(abs(posplot(:,2)),posplot(:,1),'filled','MarkerFaceAlpha',1,'MarkerFaceColor','r','MarkerEdgeColor','none'); hold on;
scatter(abs(negplot(:,2)),negplot(:,1),'filled','MarkerFaceAlpha',1,'MarkerFaceColor','b','MarkerEdgeColor','none');
xlim([0 3.5]); ylim([0 0.18]);
xlabel('Magnitude K (rad/mm)','Interpreter','latex');
ylabel('Vessel peak frequency (Hz)','Interpreter','latex');
title({'Penetrating arteriole f vs k, 7 animals, 14 trials, 140 penetrating arterioles','Blue = deep precedes shallow Red = shallow precedes deep'},'Interpreter','latex','FontSize',13)
ax = gca;
ax.TickLabelInterpreter = 'latex';
f1.Position = [4608 45 759 665];
print(gcf, '-depsc2', 'PAFvsK_RedBlue');

%Calculate median points for velocity fit and plot
wts = 1./(CombinedKFmat(:,3).^2);
kvec = abs(CombinedKFmat(:,2));
fvec = CombinedKFmat(:,1);
addpath('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\ExtractPGCode')
[wtdmedf,wtdmedk] = weightedMedFit(kvec,fvec,wts,10);
hold on;
scatter(wtdmedk,wtdmedf,'filled','k');
% [xData, yData] = prepareCurveData( wtdmedk, wtdmedf );
% % Set up fittype and options.
% ft = fittype( 'poly1' );
% opts = fitoptions( 'Method', 'LinearLeastSquares' );
% opts.Lower = [-Inf 0];
% opts.Upper = [Inf 0];
% % Fit model to data.
% [fitresult, gof] = fit( xData, yData, ft, opts );
% ci = confint(fitresult);
% df = length(wtdmedf) - 1;
% tcum = tinv(0.95,df);
% se_fit = ((ci(2,1) - ci(1,1))/2)/tcum;
[b,SEb_origin,R2Origin] = FitSlope_NoIntercept(wtdmedk,wtdmedf)
% se_v = (2*pi*fitresult.p1)*(se_fit)/fitresult.p1;
%Bootstrap SE for speed
[b_bootstrap] = fvsk_bootstrap(fvec,kvec,wts);
figure; histogram(b_bootstrap);
SE_c = 2*pi*std(b_bootstrap)

xplot = 0:0.01:1;
% yplot = xplot*fitresult.p1;
yplot = xplot*b;
plot(xplot,yplot,'k')
print(gcf, '-depsc2', 'PAFvsK_RedBlue_WtdMedFit_8_20_23');

%Plot weights
f2 = figure
scatter(kvec,wts,'filled','MarkerFaceAlpha',1,'MarkerFaceColor','k','MarkerEdgeColor','none');
xlim([0 3.5]);
xlabel('Magnitude K (rad/mm)','Interpreter','latex');
ylabel('Weight, $1/(dk)^2$','Interpreter','latex');
title({'Penetrating arteriole f vs k, 7 animals, 14 trials, 141 penetrating arterioles','Blue = deep precedes shallow Red = shallow precedes deep'},'Interpreter','latex','FontSize',13)
ax = gca;
ax.TickLabelInterpreter = 'latex';
% f2.Position = [4608 45 759 665];
print(gcf, '-depsc2', 'PA_Weights_v_K');

%Plot marginal histograms
figure
histogram(abs(posplot(:,2)),'BinWidth',0.25,'FaceColor','r','FaceAlpha',1); hold on
histogram(abs(negplot(:,2)),'BinWidth',0.25,'FaceColor','b','FaceAlpha',1)
xlim([0 3.5]);
xlabel('Magnitude K (rad/mm)','Interpreter','latex');
print(gcf, '-depsc2', 'K marginal Hist');

figure
h1 = histogram(posplot(:,1),'BinWidth',0.01,'FaceColor','r','FaceAlpha',1); hold on
h2 = histogram(negplot(:,1),'BinWidth',0.01,'FaceColor','b','FaceAlpha',1)
set(h1,'FaceAlpha',1); set(h1,'FaceAlpha',1);
set(h2,'FaceAlpha',1); set(h1,'FaceAlpha',1);
h1.Orientation = 'horizontal';
h2.Orientation = 'horizontal';
ylabel('Frequency (Hz)','Interpreter','latex');
ylim([0 0.18])
xlabel('Count','Interpreter','latex')
ax = gca;
ax.TickLabelInterpreter = 'latex';
print(gcf, '-depsc2', 'f marginal Hist');

%Plot cumulative distributions
figure
[e1,x1] = ecdf(kvec);
plot(x1,e1,'k');
xlabel('magnitude k','Interpreter','latex');
ylabel('cumulative probability','Interpreter','latex');
xlim([0 3.5]);
print(gcf, '-depsc2', 'PA_K_CDF');
figure
[e2,x2] = ecdf(fvec);
hold on; plot(e2,x2,'k');
ylim([0 0.18]);
ylabel('Frequency (Hz)','Interpreter','latex');
xlabel('Cumulative probability','Interpreter','latex');
print(gcf, '-depsc2', 'PA_fvaso_CDF');


%% Figure 4

%Plot Combined fluorescence image
%Average time series (Ca and FWHM) with GCaMP axis flipped
%Normalized cross correlation summary (average correlation curve across vessels)
%Normalized maximum correlation for each vessel (fxn of average vessel diameter)
%Also coherence (averaged over locations and vessels) fxn of freq, look for fall off at 0.3Hz like Thomas's data show

%Comb Fluorescence image (Re-make it after correcting 
%ROI 10)
clear; close all;
% data_folder = 'Y:\Data backup\20230613 JD221211F1';
data_folder = 'C:\Users\duckw\Desktop\FileTransfer';
animal = 'JD221211F1'; %For file naming

cd(data_folder);
files = dir('*001.tif');

file = 2;
namestr = extractBefore(files(file).name,'_00001.tif');

%Load images
cd('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P')
outvars = pyrunfile("TiffDimDetector.py 'C:\Users\duckw\Desktop\FileTransfer\roi10_00001.tif'",'test');
ims = str2double(outvars.char)

%Correct for noise floor before making fused image
cd('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P')
tmpimsD = pyrunfile("TiffImReader.py 'C:\Users\duckw\Desktop\FileTransfer\roi10_00001.tif'",'test', r1 = int32(1), r2 = int32(round(ims+1)), r3 = int32(2)); 
D_data = double(tmpimsD); 
% D_manual_offset = min(D_data(:)) %Calcualate an offset manually. Scanimage can apply an offset to center noise around zero (zero voltage = 0 photons detected)
%Track frames with saturation in cy5.5 channel & exclude (usually due to
%large dye particle)
figure; imagesc(squeeze(D_data(100,:,:))); daspect([1,4,1]); title('Select rectangle for noise quantification');
[noisecol,noiserow] = ginput(2);
D_noise = D_data(:,round(noiserow(1)):round(noiserow(2)),round(noisecol(1)):round(noisecol(2)));
figure; histogram(D_noise); title('Choose Bounsd for gaussian fit');%Fit gaussian to noise
D_noise = D_noise(:);
[gausscol, ~] = ginput(2); 
D_noise_tofit = D_noise(and(D_noise > gausscol(1),D_noise < gausscol(2)));
pd = fitdist(D_noise_tofit,'Normal');
noise_mean = pd.mu %Noise level to be applied as offset manually.
clearvars D_data tmpimsD D_noise D_noise_tofit 

tmpimsCa = pyrunfile("TiffImReader.py 'C:\Users\duckw\Desktop\FileTransfer\roi10_00001.tif'",'test', r1 = int32(0), r2 = int32(round(ims)), r3 = int32(2)); 
Ca_data = double(tmpimsCa); 
% D_manual_offset = min(D_data(:)) %Calcualate an offset manually. Scanimage can apply an offset to center noise around zero (zero voltage = 0 photons detected)
%Track frames with saturation in cy5.5 channel & exclude (usually due to
%large dye particle)
figure; imagesc(squeeze(Ca_data(100,:,:))); daspect([1,4,1]); title('Select rectangle for noise quantification');
[noisecol,noiserow] = ginput(2);
Ca_noise = Ca_data(:,round(noiserow(1)):round(noiserow(2)),round(noisecol(1)):round(noisecol(2)));
figure; histogram(Ca_noise); title('Choose Bounsd for gaussian fit');%Fit gaussian to noise
Ca_noise = Ca_noise(:);
[gausscol, ~] = ginput(2); 
Ca_noise_tofit = Ca_noise(and(Ca_noise > gausscol(1),Ca_noise < gausscol(2)));
pdCa = fitdist(Ca_noise_tofit,'Normal');
noise_mean_Ca = pdCa.mu %Noise level to be applied as offset manually.
clearvars Ca_data tmpimsCa Ca_noise Ca_noise_tofit 

%ROI10 Acquisitin rate = 31.97Hz
%Load images 140 - 290 where diameter is relatively constant -> 150 images,
%4.7 seconds.
%Make vessel masks.
outvarsD = pyrunfile("TiffImReader.py 'C:\Users\duckw\Desktop\FileTransfer\roi10_00001.tif'",'test', r1 = int16(141), r2 = int16(291), r3 = int16(2)); 
Dmask_tmp = double(outvarsD) - noise_mean; Dmask = squeeze(mean(Dmask_tmp,1));
outvarsCa = pyrunfile("TiffImReader.py 'C:\Users\duckw\Desktop\FileTransfer\roi10_00001.tif'",'test', r1 = int16(140), r2 = int16(290), r3 = int16(2));
Camask_tmp = double(outvarsCa) - noise_mean_Ca; Camask = squeeze(mean(Camask_tmp,1));
% figure; subplot(2,1,1); imagesc(Dmask); daspect([1,4,1]); axis off; subplot(2,1,2); imagesc(Camask); daspect([1,4,1]); axis off;
fusemask = imfuse(Dmask,Camask); figure; imshow(fusemask); daspect([1,4,1]); %300 x 480 (original size) 
fusemask_resize = imresize(fusemask,[size(fusemask,1),size(fusemask,2)*4],'nearest'); %Try not median filtering first
figure; imshow(fusemask_resize); daspect([1,1,1]); %1200 x 480 (original size) 

cd(['Z:\Thomas_1P\JD_GC_Diam2P\6_13_23\7_24Analysis\',namestr]);
savefig('Ca_Lumen_imfust_noisecorrected.fig');
%Save as eps
cd('Z:\Thomas_1P\JD_GC_Diam2P\6_13_23\Plots7_24\');
print(gcf, '-depsc2', 'JD221211F1ROI10_CombChannels_8_20Version');

%Plot with lines
cd(['Z:\Thomas_1P\JD_GC_Diam2P\6_13_23\7_24Analysis\',namestr]);
% load('seg_inds.mat'); %Need to load and order inds using code above
load('resultsmat_sf3.mat');
load('seg_inds.mat');
mask = logical(rgb2gray(imread('Mask.tif')));
indstodel = resultsmat.indstodel;
todel = resultsmat.todel;
inds(logical(indstodel)) = [];

numpts = sum(~todel);
colormat = zeros(numpts,3); colormat(:,1) = 0; colormat(:,3) = 1; colormat(:,2) = flip((1/numpts):(1/numpts):1);
colormat(:,4) = 0.5;
coloriter = 1;
%Plot locations on vessel image
h1 = openfig('Ca_Lumen_imfust_noisecorrected.fig','reuse'); 
ax1 = gca;
for i = 1:length(inds)/10
    if todel(i) == 0
        hold on;
        linenum = i*10-5;
        ind1 = inds{1,linenum}(1);
        ind2 = inds{1,linenum}(end);
        [row1,col1] = ind2sub(size(mask),ind1);
        [row2,col2] = ind2sub(size(mask),ind2);
        line([col1,col2],[row1,row2],'Color',colormat(coloriter,:),'LineWidth',1.5)
        coloriter = coloriter + 1;
    end
end
cd('Z:\Thomas_1P\JD_GC_Diam2P\6_13_23\Plots7_24\');
print(gcf, '-depsc2', 'JD221211F1ROI10_CombChannels_SingleLines_8_20Version');

%Plot all lines
cd(['Z:\Thomas_1P\JD_GC_Diam2P\6_13_23\7_24Analysis\',namestr]);
h1 = openfig('Ca_Lumen_imfust_noisecorrected.fig','reuse'); 
ax1 = gca;
coloriter = 1;
for i = 1:length(inds)/10 %Update this - plotting less than 10 lines per group
    if todel(i) == 0
        hold on;
        for j = 1:10
            linenum = i*10-9;
            linenum2 = j-1;
            ind1 = inds{1,linenum+linenum2}(1);
            ind2 = inds{1,linenum+linenum2}(end);
            [row1,col1] = ind2sub(size(mask),ind1);
            [row2,col2] = ind2sub(size(mask),ind2);
            line([col1,col2],[row1,row2],'Color',colormat(coloriter,:))
        end
        coloriter = coloriter + 1;
    end
end
%Plot scale bar
line([41 80],[443 443],'Color','white');
cd('Z:\Thomas_1P\JD_GC_Diam2P\6_13_23\Plots7_24\');
print(gcf, '-depsc2', 'JD221211F1ROI10_CombChannels_AllLines_8_20Version');


%Average time series (Ca and FWHM) with GCaMP axis flipped
time = resultsmat.time;
diam = resultsmat.diam;
wave = resultsmat.wave;
figure
yyaxis('left')
plot(time,mean(diam,1),'k');
ylabel('Diam (um)','Interpreter','latex','Color','k');
set(gca,'ycolor','k'); ylim([5 20]);
yyaxis right
plot(time,mean(wave),'g');
ylabel('Ca2+ dF/F','Interpreter','latex','Color','g');
set(gca,'ycolor',[0.4660 0.6740 0.1880]); ylim([-1 1]);
set(gca,'YDir','reverse')
% legend({'Diam FWHM','Calcium dF/F'})
xlabel('Time (s)','Interpreter','latex');
cd('Z:\Thomas_1P\JD_GC_Diam2P\6_13_23\Plots7_24\');
print(gcf, '-depsc2', 'ROI10AverageDiam_GCaMPdF_F_8_20Version');

lumwave = resultsmat.lumwave;
figure
yyaxis('left')
plot(time,mean(lumwave,1),'k');
ylabel('Sum lumen intensity (arb)','Interpreter','latex','Color','k');
set(gca,'ycolor','k') 
yyaxis right
plot(time,mean(wave),'g');
ylabel('Ca2+ dF/F','Interpreter','latex','Color','g');
set(gca,'ycolor',[0.4660 0.6740 0.1880]) 
set(gca,'YDir','reverse')
% legend({'Diam FWHM','Calcium dF/F'})
xlabel('Time (s)','Interpreter','latex');
cd('Z:\Thomas_1P\JD_GC_Diam2P\6_13_23\Plots7_24\');
print(gcf, '-depsc2', 'AverageIntensity_GCaMPdF_F');

%%
%Fig4D
%Calcium-SumIntensity lag vs average vessel diameter
clear; clc; close all;
cd('Z:\Thomas_1P\JD_GC_Diam2P\6_13_23\7_24Analysis\');
files = dir('*/*resultsmat.mat');

allwts = [];
allLags = [];
allCorr = [];
counter = 1;
for i = 1:length(files)
    clearvars D_Ca_lagvec D_Ca_corrvec
    %Load results mat
    load([files(i).folder,'/',files(i).name])
    lowp = designfilt('lowpassiir','PassbandFrequency',0.3,...
        'StopbandFrequency',0.35,'StopbandAttenuation',10,'SampleRate',resultsmat.rate);
    % lowp = designfilt('lowpassiir','PassbandFrequency',1,...
    %     'StopbandFrequency',1.05,'StopbandAttenuation',10,'SampleRate',resultsmat.rate);

    wave = resultsmat.wave;
    diam = resultsmat.diam;
%     if i == 9
%         wave = wave(:,750:end);
%         diam = diam(:,750:end);
%     end
    diammean = resultsmat.diammean;
    wavemean = resultsmat.wavemean;
    diammeanmat = repmat(diammean,[1,size(diam,2)]);
    wavemeanmat = repmat(wavemean,[1,size(wave,2)]);
    diam_zmean = diam - diammeanmat;
    wave_zmean = wave - wavemeanmat;

    wave_lp = transpose(filtfilt(lowp,wave_zmean'));
    diam_lp = transpose(filtfilt(lowp,diam_zmean'));
    
    numpts = size(diam,1);
    corrmat_todel_trial = zeros(numpts,1);
    for j = 1:numpts %Calculate cross correlations from time series
        tmpD = diam_lp(j,:);
        tmpca = wave_lp(j,:);
        [corrtmp,lagstmp] = xcorr(tmpD,-tmpca,'normalized'); %Max lag positive when intensity lags ca2+

        [maxcorr,maxlag] = max(corrtmp);
        D_Ca_lagvec(j) = lagstmp(maxlag)/resultsmat.rate;
        D_Ca_corrvec(j) = maxcorr;

        corrmat(counter,1:length(lagstmp)) = corrtmp; %For plotting the average curve.
        if maxcorr < 0.3
            corrmat_todel(counter) = 1;
            corrmat_todel_trial(j) = 1;
        end
        counter = counter + 1;
    end

    %For each vessel, calculate the |Corr|-weighted lag and
    %weighted-uncertainty.

    todelD = abs(D_Ca_corrvec) < 0.3;
    D_Ca_corrvec(todelD) = [];
    D_Ca_lagvec(todelD) = [];
    diammean(todelD) = [];

    nD = numel(D_Ca_lagvec);
    avgdiamD(i) = mean(diammean); %um
    stddiamD(i) = std(diammean);
    sediamD(i) = std(diammean)/sqrt(nD); %Standard error of mean diameter across vessel segment
    
    D_CaAvgLag(i) = mean(D_Ca_lagvec);
    wts = abs(D_Ca_corrvec);
    D_CaAvgLag_Wtd(i) = sum(wts.*D_Ca_lagvec)/sum(wts);
    wtd2ndmom = sum(wts .* D_Ca_lagvec.^2)/sum(wts);
    wtdmean = D_CaAvgLag_Wtd(i);
    D_CaLagSE_Wtd(i) = sqrt(wtd2ndmom - wtdmean^2)*(1/sqrt(nD));
    D_CaLagSE(i) = std(D_Ca_lagvec)/sqrt(nD);   

    D_CaAvgCorr(i) = mean(D_Ca_corrvec);
    D_CaCorrSE(i) = std(D_Ca_corrvec)/sqrt(numel(D_Ca_corrvec));

    allwts = [allwts,wts];
    allLags = [allLags,D_Ca_lagvec];
    allCorr = [allCorr,D_Ca_corrvec];

    numlines(i) = numpts - sum(corrmat_todel_trial);
    lagmat(i,1:length(lagstmp)) = lagstmp/resultsmat.rate;

    i
end

%% Temporary, visualize each line raw data and cross correlation
figure;
maxval = max(diam(:));
minval = min(diam(:));
camaxval = max(wave(:));
caminval = min(wave(:));
lowp = designfilt('lowpassiir','PassbandFrequency',0.3,...
    'StopbandFrequency',0.35,'StopbandAttenuation',10,'SampleRate',resultsmat.rate);
diammean = resultsmat.diammean;
wavemean = resultsmat.wavemean;
diammeanmat = repmat(diammean,[1,size(diam,2)]);
wavemeanmat = repmat(wavemean,[1,size(wave,2)]);
diam_zmean = diam - diammeanmat;
wave_zmean = wave - wavemeanmat;

wave_lp = transpose(filtfilt(lowp,wave_zmean'));
diam_lp = transpose(filtfilt(lowp,diam_zmean'));

% wave_lp = wave_lp(:,794:end);
% diam_lp = diam_lp(:,794:end);

for i = 1:size(diam,1)
%     plotdiam = diam(i,794:end);
%     plotwave = wave(i,794:end);
%     plottime = resultsmat.time(794:end);
    plotdiam = diam(i,1:end);
    plotwave = wave(i,1:end);
    plottime = resultsmat.time(1:end);

    subplot(3,1,1);
    plot(plottime,plotdiam); ylim([minval,maxval]);
    subplot(3,1,2);
    plot(plottime,plotwave); ylim([caminval,camaxval]);

    tmpD = diam_lp(i,:);
    tmpca = wave_lp(i,:);
    [corrtmp,lagstmp] = xcorr(tmpD,-tmpca,'normalized'); %Max lag positive when intensity lags ca2+
    [maxcorr,maxlag] = max(corrtmp);

    subplot(3,1,3);
    plot(lagstmp/resultsmat.rate,corrtmp);
    xlim([-20 20]); ylim([-0.5 1]);
    xline(lagstmp(maxlag)/resultsmat.rate);

    title(['line ',num2str(i),' lag = ',num2str(round(lagstmp(maxlag)/resultsmat.rate,2)),' s']);

    D_Ca_lagvec(i) = lagstmp(maxlag)/resultsmat.rate;
    D_Ca_corrvec(i) = maxcorr;

%     pause();


end
mean(D_Ca_lagvec)
mean(D_Ca_corrvec)

%%

totwtdLag = sum(allwts.*allLags)/sum(allwts);
totAvgCorr = mean(allCorr);
totAvgCorrSE = std(allCorr)/sqrt(numel(allCorr));

%Plot weighted results
f1 = figure;
errorbar(avgdiamD,D_CaAvgLag_Wtd,D_CaLagSE_Wtd,D_CaLagSE_Wtd,sediamD,sediamD,'o')
title({'Calcium to Diameter (FWHM) lags, 11 Vessels',['Weighted average = ',num2str(round(totwtdLag,2)),' seconds']},'Interpreter','latex');
xlabel('Vessel Diameter ($\mu$m)','Interpreter','latex');
ylabel('Calcium-to-Diameter Lag (s)','Interpreter','latex');
ylim([0 5]);
yline(totwtdLag)
xlim([5 20])
xticks([0 5 10 15 20])
f1.Position = [4551 422 677 309];
% cd('Z:\Thomas_1P\JD_GC_Diam2P\6_13_23\Plots7_24\');
% print(gcf, '-depsc2', 'PeakCorrelationLag_Errorbars');

%Plot correlation averages
f2=figure;
errorbar(avgdiamD,D_CaAvgCorr,D_CaCorrSE,D_CaCorrSE,sediamD,sediamD,'o')
title({'Calcium to Diameter (FWHM) peak correlation, 11 Vessels',['Weighted average = ',num2str(round(totAvgCorr,2))]},'Interpreter','latex');
xlabel('Vessel Diameter ($\mu$m)','Interpreter','latex');
ylabel('Calcium-to-Diameter norm. peak correlation (s)','Interpreter','latex');
xlim([5 20])
ylim([0 1]);
yline(totAvgCorr)
xticks([0 5 10 15 20])
f2.Position = [4551 422 677 309];
% cd('Z:\Thomas_1P\JD_GC_Diam2P\6_13_23\Plots7_24\');
% print(gcf, '-depsc2', 'PeakCorrelation_Errorbars');

%Interpolate and plot Correlation curve
corrmat(logical(corrmat_todel),:) = [];
[~,isminlag] = min(max(lagmat,[],2));
lagreference = lagmat(isminlag,:);
lagreference = lagreference(lagreference ~= 0);

intpcorr = [];
for i = 1:length(numlines)
    numlags = sum(lagmat(i,:) ~= 0);
    if i==1
        v = corrmat(1:numlines(1),1:numlags);
    else
        v = corrmat((sum(numlines(1:(i-1))) + 1):sum(numlines(1:i)),1:numlags);
    end
    intpcorr = [intpcorr,interp1(lagmat(i,1:numlags),v',lagreference)];

end

figure
plot(lagreference,mean(intpcorr,2),'b'); 
xlim([-20 20]); ylim([-0.2 1]);
yticks([0 0.5 1])
xline(0); yline(0);
[max_intpcorr,max_intpcorr_loc] = max(mean(intpcorr,2));
xline(lagreference(max_intpcorr_loc))
hold on;
%Calc SE after atanh transformation
sigma12 = squeeze(std(atanh(intpcorr),0,2));%sigma12 = sqrt(dim-1)*squeeze(std(atanhCxyk,1,1));
Cu = atanh(mean(intpcorr,2)) + sigma12/sqrt(size(intpcorr,2)); %Transformed variable confidence limits
Cl = atanh(mean(intpcorr,2)) - sigma12/sqrt(size(intpcorr,2)); %Transformed variable confidence limits

Cerr(1,:) = tanh(Cl);
Cerr(2,:) = tanh(Cu);

plot(lagreference,Cerr(1,:),'b')
plot(lagreference,Cerr(2,:),'b')
xlabel('Lag (s)','Interpreter','latex');
ylabel('Normalized Cross Correlation','Interpreter','latex');
title([num2str(size(intpcorr,2)),' locations across 11 vessels $\pm$ SEM'],'Interpreter','latex');
cd('Z:\Thomas_1P\JD_GC_Diam2P\6_13_23\Plots7_24\');
print(gcf, '-depsc2', 'AverageCorrVSLag_Lowpassfilt');

%% Figure 4new
clear; clc; close all;
cd('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Pipeline\7_8_23_Analysis\Fig4\9_4_23')
%F vs K pial vessels
load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\AllSegments\1_30_23_Results\CombinedResults\9_4_23\pvcomb_vesselfv_tapha_01_750um_8869ves.mat");
kvec = abs(pvcomb(:,1));
fvec = pvcomb(:,3);
r2vec = pvcomb(:,2).^2;

minlength = 0.75; %mm
f1 = figure;
alpha = 0.01;
scatter(kvec,fvec,'filled','MarkerFaceColor','k','MarkerEdgeColor','none','MarkerFaceAlpha',1);
xlim([0 2]); ylim([0 0.18]);
xlabel('Phase gradient magnitude (rad/mm)','Interpreter','latex');
ylabel('Vessel peak vasomotor frequency (Hz)','Interpreter','latex');
str = sprintf('Pial arteriole f vs k, 24 animals, %.0f vessels',length(pvcomb));
str1 = sprintf('%.2f mm minimum vessel length, T-test',minlength);
str2 = sprintf('%.2f',alpha);
title({str,[str1,' $\alpha$ = ',str2]},'Interpreter','latex')
ax = gca;
ax.TickLabelInterpreter = 'latex';
f1.Position = [4608 45 759 665];
print(gcf, '-depsc2', 'fvsk_75length_01Tfilt');

[wtdmedf,wtdmedk] = weightedMedFit(kvec,fvec,r2vec,10);
hold on;
scatter(wtdmedk,wtdmedf,'filled','r');

% [xData, yData] = prepareCurveData( wtdmedk, wtdmedf );
% % Set up fittype and options.
% ft = fittype( 'poly1' );
% opts = fitoptions( 'Method', 'LinearLeastSquares' );
% opts.Lower = [-Inf 0];
% opts.Upper = [Inf 0];
% % Fit model to data.
% [fitresult, gof] = fit( xData, yData, ft, opts );
% ci = confint(fitresult);
% df = length(wtdmedf) - 1;
% tcum = tinv(0.95,df);
% se_fit = ((ci(2,1) - ci(1,1))/2)/tcum;
% se_v = (2*pi*fitresult.p1)*(se_fit)/fitresult.p1;
[b,SEb_origin,R2Origin] = FitSlope_NoIntercept(wtdmedk,wtdmedf);
[b_bootstrap] = fvsk_bootstrap(fvec,kvec,r2vec);
figure; histogram(b_bootstrap);
SE_c = 2*pi*std(b_bootstrap)
c = 2*pi*b

xplot = 0:0.01:1;
% yplot = xplot*fitresult.p1;
yplot = xplot*b;


plot(xplot,yplot,'r')
print(gcf, '-depsc2', 'fvsk_75length_01Tfilt_WtdMedFit_9_4_23');




kvarnum = dot(r2vec,kvec);
kvardenom = sum(r2vec);
kvarmean = kvarnum/kvardenom

r2wtd2ndmoment = dot(r2vec,kvec.^2)/sum(r2vec);
r2wtduncert = sqrt(r2wtd2ndmoment - kvarmean^2)

%Marginal Histograms
figure
h1 = histogram(fvec,'BinWidth',0.01);
set(h1,'FaceAlpha',1);
h1.Orientation = 'horizontal';
ylabel('Frequency (Hz)','Interpreter','latex');
ylim([0 0.18])
xlabel('Count','Interpreter','latex')
ax = gca;
ax.TickLabelInterpreter = 'latex';
% title({'Vasomotion frequency distribution 7 animals 138 penetrating arterioles','Freq. at max power of shallow vessel segment spectrum'},'Interpreter','latex');
savefig('VasofreqMarginalHist.fig');
print(gcf, '-depsc2', 'VasofreqMarginalHist');


figure
h1 = histogram(abs(kvec),'BinWidth',0.05);
xlabel('Magnitude phase gradient (rad/mm)','Interpreter','latex');
ylabel('Count','Interpreter','latex')
xlim([0 2])
ax = gca;
ax.TickLabelInterpreter = 'latex';
set(h1,'FaceAlpha',1);
xline(2*pi/25); xline(2*pi/9); %Example wavelengths in plot.
% title('Phase gradient distribution 7 animals 138 penetrating arterioles','Interpreter','latex');
savefig('PhasegradMarginalHist.fig');
print(gcf, '-depsc2', 'PhaseGradMarginalHist');

[f,x] = ecdf(pvcomb(:,3));
figure
plot(f,x)
ylim([0 0.18])
ylabel('Frequency (Hz)', 'Interpreter','latex')
xlabel('Probability','Interpreter','latex')
title('Vasomotion Frequency Cumulative Distribution, Alpha = 0.01','Interpreter','latex')
ax = gca;
ax.TickLabelInterpreter = 'latex'; %This is necesssary!!! who knows why.
cd('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Pipeline\7_8_23_Analysis\Fig4\9_4_23')
print(gcf, '-depsc2', 'MarginalCumulativeDist_Alpha0_01');
savefig('VasoFreq_MarginalCumulativeDist_Alpha0_01.fig')

figure
ecdf(abs(pvcomb(:,1)))
xlim([0 2])
xlabel('Phase gradient magnitude (rad/mm)','Interpreter','latex');
ylabel('Probability','Interpreter','latex')
% str = sprintf('JD221024F2, %.0f vessels',length(pvcomb));
% str1 = sprintf('%.2f mm minimum vessel length, T-test',minlength);
% str2 = sprintf('%.2f',alpha);
% title({['Phase gradient Cumulative Distribution ',str],[str1,' $\alpha$ = ',str2]},'Interpreter','latex');
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 11
cd('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Pipeline\7_8_23_Analysis\Fig4\9_4_23')
print(gcf, '-depsc2', 'PhaseGrad_MarginalCumulativeDist_Alpha0_01');
savefig('PhaseGrad_MarginalCumulativeDist_Alpha0_01.fig')
%% Supp Fig 3
clear; clc; close all;
load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Pipeline\7_8_23_Analysis\CombinedCorrmat_WT10_OnePA18Measurement.mat"); %Load combined results (Corrmat)
load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Pipeline\7_8_23_Analysis\CombinedKFmat.mat"); %Load combined results (KFmat)

% Plot correlation figure
CohTodel = ~logical(CombinedCorrMat(:,12)); %Delete combinations where either vessel have not significant coherence
CombinedCorrMat(CohTodel,:) = [];
KFtodel = ~logical(CombinedKFmat(:,6));
CombinedKFmat(KFtodel,:) = [];
probinward = sum(sign(CombinedKFmat(:,2)) > 0)/size(CombinedKFmat,1);
proboutward = sum(sign(CombinedKFmat(:,2)) < 0)/size(CombinedKFmat,1);
prob_samesign = probinward^2 + proboutward^2 %Observed probability of two random vessels traveling in the same direction
prob_opsign = 2*probinward*proboutward %Observed probability of two random vessels traveling in the opposite direction

f1 = figure;
histogram(CombinedCorrMat(:,5),'BinWidth',0.1);
xlabel('Distance (mm)','Interpreter','latex');
ylabel('Count','Interpreter','latex');
title('Pair Distance Distribution','Interpreter','latex');
[f,d] = ecdf(CombinedCorrMat(:,5));
figure; plot(d,f);

numSDs = 0;
numbins = 10; %10 percent of pairs in each bin
tol = (f(2)-f(1))/2;
binvec = zeros(numbins,1);
for i=1:numbins
    binvectmp = find(abs(f-(i/numbins)) < tol);
    [~,minbinvec] = min(abs(f(binvectmp)-(i/numbins)));
    binvec(i) = d(binvectmp(minbinvec));
end
binvec = [0;binvec];
expecvec = zeros(length(binvec),3);
% CombCorrMat_tmp = CorrMat; %or Combined Corr Mat
CombCorrMat_tmp = CombinedCorrMat;
todel1 = abs(CombCorrMat_tmp(:,6)) < numSDs * CombCorrMat_tmp(:,7);
todel2 = abs(CombCorrMat_tmp(:,8)) < numSDs * CombCorrMat_tmp(:,9);
todel = or(todel1,todel2);
CombCorrMat_tmp(todel,:) = [];
for i=2:size(binvec,1)
    inbin1 = CombCorrMat_tmp(:,5) <= binvec(i);
    inbin2 = CombCorrMat_tmp(:,5) > binvec(i-1);
    inbin = and(inbin1,inbin2); %indices for pairs within distance bin
    x_dist(i) = median(CombCorrMat_tmp(inbin,5)); %Use median of binned distances to plot
    
    spinvec = sign(CombCorrMat_tmp(:,3).*CombCorrMat_tmp(:,4));
    spinvec = spinvec(inbin);

    findmax = [sum(spinvec<0),sum(spinvec>0)];
    [foundmax,maxind] = max(findmax);
    if maxind == 1
        probsuccess = prob_opsign;
    elseif maxind == 2
        probsuccess = prob_samesign;
    end

%     prob(i) = 1 - binocdf(foundmax,numel(spinvec),probsuccess); %1-probability to see more than (foundmax) successes
    prob(i) = binopdf(foundmax,numel(spinvec),probsuccess); %Probability density at observed value only (don't use this)
    possiblevals = 1:numel(spinvec);
    pdfvals = binopdf(possiblevals,numel(spinvec),probsuccess);
    [~,pdfmaxind] = max(pdfvals);
    bpdf(i) = possiblevals(pdfmaxind); %Value with highest binopdf
    if maxind == 1 %There are more opposite-signed pairs in the bin
        mostlikelyavg(i) = ((numel(spinvec)-bpdf(i))-bpdf(i))/numel(spinvec);
        biolims = binoinv([0.025 0.975],numel(spinvec),probsuccess);
        limsavg(i,:) = ((numel(spinvec)-biolims)-biolims)./numel(spinvec);
    elseif maxind == 2 %There are more same-signed pairs in the bin
        mostlikelyavg(i) = (-(numel(spinvec)-bpdf(i))+bpdf(i))/numel(spinvec);
        biolims = binoinv([0.025 0.975],numel(spinvec),probsuccess);
        limsavg(i,:) = (-(numel(spinvec)-biolims)+biolims)./numel(spinvec);
    end

    meanspinvec = mean(spinvec); %This is what we plot on y axis
    stdspinvec = std(spinvec);

    expecvec(i,1) = meanspinvec;
    expecvec(i,2) = stdspinvec;
    expecvec(i,3) = stdspinvec/sqrt(numel(spinvec));
    numinbin(i) = numel(spinvec);
end
expecvec(1,:) = []; %First entry has nothing
prob(1) = [];
limsavg(1,:) = [];
mostlikelyavg(1) = [];
x_dist(1) = [];
numinbin(1) = [];
%Make scatter from expecvec
% x_dist = (binvec(1:end-1) + binvec(2:end))/2; %Get mean distance of every bin for plotting
%Set significance level as constant for constant number in each bin
minnum = min(numinbin);
% maxnum = max(numinbin);
%Here we're choosing opposite signs as success and calculating 2.5% and
%97.5% bounds for <sig_i*sig_j>. Could also use probsucess = prob_samesign
%and get the same error bounds. Use min(numbin) to be conservative, and
%show that even in the best case we don't observe directionality different
%than chance.
probsuccess = prob_opsign; 
biolimsconst = binoinv([0.025 0.975],minnum,probsuccess);
limsavgconst = ((minnum-biolimsconst)-biolimsconst)./minnum; %(minnum-biolimsconst) is number of SAME SIGNS (+1), then subtract number of OPPOSITE SIGNS (-1) and divide by total to get <sig_i*sig_j>.

y_expec = expecvec(:,1); %<sig1 * sig2>
% err = expecvec(:,3); %SDM
f2 = figure;
% subplot(2,1,1);
scatter(x_dist,y_expec,'filled','blue');
hold on
% s1 = scatter(x_dist,mostlikelyavg,'_r')
% s1.LineWidth = 1.5;
% s2 = scatter(x_dist,limsavg,'_r')
% s2 = scatter(x_dist,limsavg,'_r')
yline(limsavgconst(1),'r'); yline(limsavgconst(2),'r');
% s2.SizeData = 0.75
xlim([0 1.8]); ylim([-1 1]); yline(0);
xlabel('Distance (10Percent of all pairs in each bin)','Interpreter','latex','FontSize',12);
ylabel('$\biggl \langle \sigma_i \sigma_j \biggr \rangle$','Interpreter','latex','FontSize',12)
grid on
title('Correlation vs euclidean distance, 7 animals, 692 pairs','Interpreter','latex','FontSize',15);
cd('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Pipeline\7_8_23_Analysis\SuppFig3');
print(gcf, '-depsc2', 'CorrVDist_7_11_23');

%% Calc p-value for observed number of ups and downs.
n = length(CombinedKFmat);
nin = sum(sign(CombinedKFmat(:,2)) > 0);
nout = sum(sign(CombinedKFmat(:,2)) < 0);
testp = 0.5;
testsuccess = [0:1:min(nin,nout),max(nin,nout):1:n];
y = binopdf(testsuccess,n,testp);
pval = sum(y);

%Same for stim
clear;
load("X:\DataAnalysis\VesCorrPhase\Rui_2P\Pipeline\8_4_StimAnalysis\20221208PA5_Corr_KF_struct_stim.mat");
KM1 = Corr_KF_struct(1).KFmat; 
load("X:\DataAnalysis\VesCorrPhase\Rui_2P\Pipeline\8_4_StimAnalysis\20221212PA7_Corr_KF_struct_stim.mat");
KM2 = Corr_KF_struct(1).KFmat; 
CombinedKFmat = [KM1;KM2];
n = length(CombinedKFmat);
nin = sum(sign(CombinedKFmat(:,2)) > 0);
nout = sum(sign(CombinedKFmat(:,2)) < 0);
testp = 0.5;
testsuccess = [0:1:min(nin,nout),max(nin,nout):1:n];
y = binopdf(testsuccess,n,testp);
pval = sum(y);
%Plot this binomial distribution pdf
x = 0:1:n;
y = binopdf(x,n,testp);
figure; bar(x,y,1); xlabel('Number of sucesses'); ylabel('Probability');

%Calculate p value for vessel pairs
spinvec = sign(CombCorrMat_tmp(:,3).*CombCorrMat_tmp(:,4));
n = size(CombinedCorrMat,1);
x = 0:1:n; %Should be 692 after coherence filtering.
y = binopdf(x,n,prob_samesign);
figure; bar(x,y,1); xlabel('Number of same signs'); ylabel('Probability'); title('Expected Binom dist of #Same signs given observed ups & downs')
num_samesign = sum(spinvec > 0)
findobs_same = find(x==num_samesign)
pval = 2*sum(y(findobs_same:end))


%% Supplemental figure 4
%Calculate coherence and phase for calcium to diameter 2P data
%Coherence (averaged over locations and vessels) fxn of freq, do
%jacknife confidence interval with transformation.
%Re-writeso we are not interpolating complex values 
%7/24/23 version2 - Interpolate raw data so that we have a common frequency grid and no
%interpolation is required. All trials are 600s long.

clear; clc; close all;
% cd('Z:\Thomas_1P\JD_GC_Diam2P\6_13_23');
% cd('Z:\Thomas_1P\JD_GC_Diam2P\6_13_23\7_24Analysis\');
% files = dir('*/*resultsmat.mat');
cd('Z:\Thomas_1P\JD_GC_Diam2P');
files1 = load('JD230306F3_files.mat'); files2 = load('JD221223F2_files.mat'); files3 = load('JD221211F1_files.mat');
files1 = files1.files; files2 = files2.files; files3 = files3.files;
%Exclude ROI7 for JD221211F1 and ROI3 for JD221223F2



%Get rate for all trials, then interpolate raw data to common time points
addpath(genpath('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\chronux_2_12'))
for i = 1:length(files)
     %Load results mat
    load([files(i).folder,'/',files(i).name])
    allrates(i) = resultsmat.rate;
    timetmp = resultsmat.time;
    timetmp = (timetmp - timetmp(1)) + (1/(2*resultsmat.rate));
    timestruct(i).time = timetmp;
    endtimes(i) = timetmp(end);
end
%Find smallest rate will be common to all trials after interpolation. 
[smallestrate,issmallestrate] = min(allrates);
timeref = timestruct(issmallestrate).time;
rateref = smallestrate;
timeref = (timeref - timeref(1)) + (1/(2*rateref)); %Time should start at middle of first frame Interpolate all timeseries to this.

%Find shortest time series 
[shortest_triallength,shortest_trial] = min(endtimes);
timeref(timeref > shortest_triallength) = []; %Delete any reference time points that are longer than shortest trial.
%Now interpolate all raw data to the reference time.
diamstruct = struct();
wavestruct = struct();
for i = 1:length(files)
    load([files(i).folder,'/',files(i).name]);
    diamstruct(i).diam = resultsmat.diam';
    wavestruct(i).wave = resultsmat.wave';
    diamstruct(i).intp_diam = interp1(timestruct(i).time,resultsmat.diam',timeref);
    wavestruct(i).intp_wave = interp1(timestruct(i).time,resultsmat.wave',timeref);
end


counter = 1;
addpath(genpath('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\chronux_2_12'))
for i = 1:length(files) %Iterate over trials
% for i = 1:1
    %Load results mat
    load([files(i).folder,'/',files(i).name])

    diam = diamstruct(i).intp_diam;
    diammean = resultsmat.diammean;

    wave = wavestruct(i).intp_wave;
    wavemean = resultsmat.wavemean;
    rate = smallestrate; %We interpolated everything to this rate.

    params.Fs = rate;
    params.pad = 0;
    params.fpass = [0 1]; %Hz, default is [0 Fs/2]
    params.err   = [2 .05];
    params.trialave = 0;
%     T = resultsmat.time(end);
    T = timeref(end);
    BW = 0.04; %500s trial -> TBW =~ 10 Use larger BW for coherence
    params.tapers = [round(T*BW),round(2*T*BW-1)]; %Time-BW product and number of tapers

    [tapers1,pad1,Fs,fpass1,err1,trialave1]=getparams(params);
    N1 = size(wave,1);
    nfft1=max(2^(nextpow2(N1)+pad1),N1);
    [f1,findx1]=getfgrid(Fs,nfft1,fpass1);
    tapers1=dpsschk(tapers1,N1,Fs);
    %Calculate coherence
    usedtapers(i) = params.tapers(2);
    numlines(i) = size(diam,2);
    fmat(1:length(f1),i) = f1; %This is common for all calculations
    Timevec(i) = T;
    sizef(i) = length(f1);

    for j = 1:size(diam,2) %Iterate over diameter traces
        wave0 = -wave(:,j); %Use (-) wave -> Increase means decrease in Calcium
        J0 = mtfftc(wave0,tapers1,nfft1,Fs);
        J0 = J0(findx1,:);

        wave1 = diam(:,j) - diammean(j);
        J1 = mtfftc(wave1,tapers1,nfft1,Fs);
        J1 = J1(findx1,:); %Freq x (tapers) each iteration

        %Save Fourier transforms for jacknife calcs
        if j == 1 && i == 1 
            Jwavetot=J0;
            Jdiamtot=J1;
        else
            Jwavetot = [Jwavetot,J0];
            Jdiamtot = [Jdiamtot,J1];
        end

        S12_trials(:,counter) = mean(conj(J0).*J1,2);
        S1_trials(:,counter) = mean(conj(J0).*J0,2);
        S2_trials(:,counter) = mean(conj(J1).*J1,2);

        counter = counter + 1;

        %Keep J1,J2 for error calcs. Calculate mean coherence after loop.
    end
%     numlines_tapers(i) = size(Jwaveplottmp,2);
    i
end

S12 = mean(conj(Jwavetot).*Jdiamtot,2); %Average over trials and tapers (dim 2)
S1 = mean(conj(Jwavetot).*Jwavetot,2); %Average over trials and tapers (dim 2)
S2 = mean(conj(Jdiamtot).*Jdiamtot,2); %Average over trials and tapers (dim 2)
C12 = S12./sqrt(S1.*S2);
C = abs(C12); %Magnitude coherence estimate, averaged over trials and tapers
phi = angle(C12); %Phase of coherence estimate

figure
plot(f1,C); ylim([0 1]); xlabel('Frequency (Hz)','Interpreter','latex');
ylabel('Magnitude Coherence','Interpreter','latex');
str = sprintf('Average over %.0f lines across %.0f vessels',sum(numlines),length(files));
title(str,'Interpreter','latex');



%% Do error estimates with jackknife

%Calc conf intervals 
%Adapted from Chronux
dim = size(Jdiamtot,2); %Tapers * channels
dof = 2*dim;
tcrit=tinv(0.05,dof-1);

% eJ1ktmp_all = Jwavetot.*conj(Jwavetot); %Freq x (segs * tapers) Real, magnitude squared
% eJ2ktmp_all = Jdiamtot.*conj(Jdiamtot); %Freq x (segs * tapers) Real, magnitude squared
% eJ12ktmp_all = conj(Jwavetot).*Jdiamtot; %Freq x (segs * tapers) Complex
% Coh_all = eJ12ktmp_all./sqrt(eJ1ktmp_all.*eJ2ktmp_all); %All complex coherences
% Coh_mag_all = abs(Coh_all); %Magnitude coherence, tapers * segs
% Coh_phase_all = angle(Coh_all); %Angle of coherence, tapers * segs

% Interpolate these to common f
% intp_eJ12ktmp = [];
% intp_eJ1ktmp = [];
% intp_eJ2ktmp = [];
% intp_phase = [];
% for i = 1:length(numlines) %Interpolate FTs
%     numfreqs = sum(fmat(:,i) ~= 0) + 1; %+1 for f=0
%     if i==1
%         veJ12ktmp = eJ12ktmp_all(1:numfreqs,1:numlines_tapers(1));
%         veJ1ktmp = eJ1ktmp_all(1:numfreqs,1:numlines_tapers(1));
%         veJ2ktmp = eJ2ktmp_all(1:numfreqs,1:numlines_tapers(1));
% 
%         v_phasetmp = Coh_phase_all(1:numfreqs,1:numlines_tapers(1));
%     else
%         veJ12ktmp = eJ12ktmp_all(1:numfreqs,(sum(numlines_tapers(1:(i-1))) + 1):sum(numlines_tapers(1:i)));
%         veJ1ktmp = eJ1ktmp_all(1:numfreqs,(sum(numlines_tapers(1:(i-1))) + 1):sum(numlines_tapers(1:i)));
%         veJ2ktmp = eJ2ktmp_all(1:numfreqs,(sum(numlines_tapers(1:(i-1))) + 1):sum(numlines_tapers(1:i)));
% 
%         v_phasetmp = Coh_phase_all(1:numfreqs,(sum(numlines_tapers(1:(i-1))) + 1):sum(numlines_tapers(1:i)));
%     end
%     intp_eJ12ktmp = [intp_eJ12ktmp,interp1(fmat(1:numfreqs,i),veJ12ktmp,freference)]; %Interpolate complex coherence
%     intp_eJ1ktmp = [intp_eJ1ktmp,interp1(fmat(1:numfreqs,i),veJ1ktmp,freference)]; %Interpolate magnitude J1 squared
%     intp_eJ2ktmp = [intp_eJ2ktmp,interp1(fmat(1:numfreqs,i),veJ2ktmp,freference)]; %Interpolate magnitude J2 squared
% 
%     intp_phase = [intp_phase,interp1(fmat(1:numfreqs,i),v_phasetmp,freference)];
% end

tic
parfor k = 1:dim %Takes about 8 mins, 8 pools
    indxk = setdiff(1:dim,k);
    J1k = Jwavetot(:,indxk);
    J2k = Jdiamtot(:,indxk);
% 
%     eJ1k=squeeze(sum(J1k.*conj(J1k),2));
%     eJ2k=squeeze(sum(J2k.*conj(J2k),2));
%     eJ12k=squeeze(sum(conj(J1k).*J2k,2));
% 
%     Cxyk=eJ12k./sqrt(eJ1k.*eJ2k);
%     absCxyk=abs(Cxyk);
%     logCxyk(k,:) = log((absCxyk.^2)./(1-absCxyk.^2)); %This ranges from -inf to inf
%     phasefactorxyk(k,:)=Cxyk./absCxyk;

    %Same for Spectra!
    eJ1k=sum(J1k.*conj(J1k),2);
    S1k(k,:)=eJ1k/(dim-1); % 1-drop (mean) spectrum
    eJ2k=sum(J2k.*conj(J2k),2);
    S2k(k,:)=eJ2k/(dim-1); % 1-drop (mean) spectrum
end
toc

%%
logC = log(C.^2./(1-C.^2)); %Transfromed magnitude coherence

%Calculate coherence and phase error bars
sigma12 = sqrt(dim-1)*squeeze(std(logCxyk,1,1)); %Standard deviation of transformed variable
mean12 = mean(logCxyk,1); %Mean of transformed variable
% Cu = logC + tcrit(ones(size(f1,1),1),:).*sigma12'; %Transformed variable confidence limits
% Cl = logC - tcrit(ones(size(f1,1),1),:).*sigma12'; %Transformed variable confidence limits
Cu = logC + sigma12'; %Transformed variable + standard error
Cl = logC - sigma12'; %Transformed variable + standard error

Cerr(1,:) = max((1+exp(-Cl)).^(-1/2),0); %Update this, to transform back from z
Cerr(2,:) = (1+exp(-Cu)).^(-1/2);
phistd = sqrt( (2*dim-2)*(1-abs(squeeze(mean(phasefactorxyk)))) ); %Standard error of phase

%Calculate spectra error bars
Swave=squeeze(mean(conj(Jwavetot).*Jwavetot,2));
Sdiam=squeeze(mean(conj(Jdiamtot).*Jdiamtot,2));

sigma=sqrt(dim-1)*squeeze(std(log(S1k),1,1));
conf=sigma; conf=squeeze(conf);
Swaveerr(1,:)=Swave'.*exp(-conf); Swaveerr(2,:)=Swave'.*exp(conf);

sigma=sqrt(dim-1)*squeeze(std(log(S2k),1,1));
conf=sigma; conf=squeeze(conf);
Sdiamerr(1,:)=Sdiam'.*exp(-conf); Sdiamerr(2,:)=Sdiam'.*exp(conf);

if sum(usedtapers == usedtapers(1)) == length(usedtapers)
    dof = 2*usedtapers(1);
end
df = 1./((dof/2)-1);
confC = sqrt(1 - 0.05.^df)

%% Plot results
figure
plot(f1,C,'b');
xlim([0 1]); ylim([0 1]);
yline(confC);
hold on;
plot(f1,Cerr','Color',[0,0,1,0.2])
xlabel('Freq (Hz)','Interpreter','latex');
ylabel('Magnitude Coherence','Interpreter','latex');
title([num2str(size(intpcoh,2)),' locations across 11 vessels $\pm$ jackknife 95Pct CI'],'Interpreter','latex');
cd('Z:\Thomas_1P\JD_GC_Diam2P\6_13_23\Plots7_24');
print(gcf, '-depsc2', 'AverageCoherence_DiameterCa_04BW');

%Shift phi by pi
% phitmp = phi;
% phitmp(phitmp > 0) = phitmp(phitmp > 0) - pi;
% phitmp(phitmp < 0) = phitmp(phitmp < 0) + pi;

figure
plot(f1,phi,'Color',[0,0,1]);
hold on
plot(f1,phi + tcrit*phistd','Color',[0,0,1,0.2])
plot(f1,phi - tcrit*phistd','Color',[0,0,1,0.2])
xlim([0 1]); ylim([-pi pi]); yline(0);
xlabel('Frequency (Hz)','Interpreter','latex');
ylabel('Phase (rad)')
title([num2str(sum(numlines)),' locations across 11 vessels $\pm$ jackknife SEM'],'Interpreter','latex');
cd('Z:\Thomas_1P\JD_GC_Diam2P\6_13_23\Plots7_24');
print(gcf, '-depsc2', 'AveragePhase_DiameterCa_04BW');


figure
yyaxis left
plot(f1,Cerr','Color',[0,0,1,0.2],'LineStyle','-')
hold on;
plot(f1,C,'b','LineStyle','-')
ylim([0 1]); yline(confC,'-.b','$|C|_{0.95}$','Interpreter','latex','LabelHorizontalAlignment','left');
xlabel('Freq (Hz)','Interpreter','latex');
ylabel('Magnitude Coherence','Interpreter','latex');
yyaxis right
plot(f1,phi + tcrit*phistd','Color',[0.6350 0.0780 0.1840 0.2],'LineStyle','-')
plot(f1,phi - tcrit*phistd','Color',[0.6350 0.0780 0.1840 0.2],'LineStyle','-')
hold on
plot(f1,phi,'Color',[0.6350 0.0780 0.1840],'LineStyle','-');
ylim([-pi pi]); yline(0,'Color','k','Alpha',0.3); ylabel('Phase (rad)','Interpreter','latex');
ax = gca;
ax.YAxis(1).Color = 'b';
ax.YAxis(2).Color = [0.6350 0.0780 0.1840];
set(ax.YAxis(2),'TickValues',-pi:pi/2:pi) 
ax.YAxis(2).TickLabels = {'-\pi/2','-\pi','0','\pi/2','\pi'}
cd('Z:\Thomas_1P\JD_GC_Diam2P\6_13_23\Plots7_24');
print(gcf, '-depsc2', 'AverageCoherencePhase_DiameterCa_04BW');

%Same, but plot only phase for significant coherence
sig_coh = C > confC;
sig_f = f1(sig_coh);
figure
yyaxis left
plot(f1,Cerr','Color',[0,0,1,0.2],'LineStyle','-')
hold on;
plot(f1,C,'b','LineStyle','-')
ylim([0 1]); yline(confC,'-.b','$|C|_{0.95}$','Interpreter','latex','LabelHorizontalAlignment','left');
xlabel('Freq (Hz)','Interpreter','latex');
ylabel('Magnitude Coherence','Interpreter','latex');
yyaxis right
plot(sig_f,phi(1:length(sig_f)) + tcrit*phistd(1:length(sig_f))','Color',[0.6350 0.0780 0.1840 0.2],'LineStyle','-')
plot(sig_f,phi(1:length(sig_f)) - tcrit*phistd(1:length(sig_f))','Color',[0.6350 0.0780 0.1840 0.2],'LineStyle','-')
hold on
plot(sig_f,phi(1:length(sig_f)),'Color',[0.6350 0.0780 0.1840],'LineStyle','-');
ylim([-pi pi]); yline(0,'Color','k','Alpha',0.3); ylabel('Phase (rad)','Interpreter','latex');
ax = gca;
ax.YAxis(1).Color = 'b';
ax.YAxis(2).Color = [0.6350 0.0780 0.1840];
set(ax.YAxis(2),'TickValues',-pi:pi/2:pi) 
ax.YAxis(2).TickLabels = {'-\pi/2','-\pi','0','\pi/2','\pi'}
cd('Z:\Thomas_1P\JD_GC_Diam2P\6_13_23\Plots7_24');
print(gcf, '-depsc2', 'AverageCoherencePhase_DiameterCa_04BW_SigCoherence');

%Plot Spectra
figure
plot(f1,log10(Swave),'b'); %Calcium spectrum
hold on
plot(f1,log10(Swaveerr),'Color',[0,0,1,0.2]);
xlabel('Freq (Hz)','Interpreter','latex');
ylabel('log10 Power','Interpreter','latex');
title('Average Calcium Spectrum','Interpreter','latex');
cd('Z:\Thomas_1P\JD_GC_Diam2P\6_13_23\Plots7_24');
print(gcf, '-depsc2', 'AverageCalciumSpectrum_04BW');

figure
plot(f1,log10(Sdiam),'r'); %Diameter spectrum
hold on;
plot(f1,log10(Sdiamerr),'Color',[1,0,0,0.2]);
xlabel('Freq (Hz)','Interpreter','latex');
ylabel('log10 Power','Interpreter','latex');
title('Average Diameter Spectrum','Interpreter','latex');
cd('Z:\Thomas_1P\JD_GC_Diam2P\6_13_23\Plots7_24');
print(gcf, '-depsc2', 'AverageDiameterSpectrum_04BW');

%Normalize and plot both on one axis
%Normalize by power from 0 to 1 Hz:
caint = sum(Swave);
diamint = sum(Sdiam);
Swave_norm = Swave/caint;
Sdiam_norm = Sdiam/diamint;
Swaveerr_norm = Swaveerr/caint;
Sdiamerr_norm = Sdiamerr/diamint;

figure
plot(f1,log10(Swave_norm),'b'); hold on;
plot(f1,log10(Swaveerr_norm),'Color',[0,0,1,0.2]);
xlabel('Freq (Hz)','Interpreter','latex');
ylabel('log10 Normalized Power','Interpreter','latex');
plot(f1,log10(Sdiam_norm),'r'); %Diameter spectrum
hold on;
plot(f1,log10(Sdiamerr_norm),'Color',[1,0,0,0.2]);
title({'Average Diameter Spectrum (red) and GCaMP Spectrum (blue)','normalized by power from 0 to 1 Hz'},'Interpreter','latex');
cd('Z:\Thomas_1P\JD_GC_Diam2P\6_13_23\Plots7_24');
print(gcf, '-depsc2', 'AverageDiameter_GCaMPSpectrum_04BW_Normalized01Hz');






%% Calculate spectra to double check

repdiammean = repmat(diammean,[1,size(diam,2)]);
diam_mean = diam - repdiammean;
params.trialave = 1;
[S,f] = mtspectrumc(diam_mean',params);


figure
plot(f,log10(S))

hold on
plot(f1,log10(mean(Jdiamplottmp.*conj(Jdiamplottmp),2)))



%% Supplemental figure 6

%Example of neuronal strips and arteriole phase gradients. (Show both?) 
%Overlay on quasi-raw image (average or standard deviation?)
%Need a different example vessel, this one is too short. Use function for
%plotting.

clear; close all; clc;
analyzed_folder = '\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\AllSegments\1_30_23_Results\TB201021M3\11.10.02_TB201021M3_PwrG80R75_leftpuff_0.12Hz_rightVis_0.13Hz_dualCH';
data_folder = '\\dk-server.dk.ucsd.edu\tbroggini\Prime\RAW\31.Mar.2021_11.10.02_TB201021M3_PwrG80R75_leftpuff_0.12Hz_rightVis_0.13Hz_dualCH';
cd(data_folder);
toplot = dir('*beforepuff_toplot.mat'); load(toplot.name);

Mask1 = imread('Mask1.tif');
f1 = figure; imshow(Mask1);
%Draw rectangle for roi
[col,row] = ginput(2);
pos = [col(1),row(1),diff(col),diff(row)];
hold on; rectangle('Position',pos,'EdgeColor','r');
cd('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\AllSegments\ManuscriptPlots\Supp6');
print(gcf, '-depsc2', 'ExampleMask1_11_10_02_TB201021M3');


Mask1Crop = imcrop(Mask1,pos);
% f2 = figure;
% imshow(Mask1Crop)

cd(analyzed_folder);
results = dir('*results.mat');
for i = 1:2
    if contains(results(i).name,'neural')
        isneu = i;
    else
        isves = i;
    end
end
neuresults = load(results(isneu).name); neuresults = neuresults.link_results_struct_neu;
vesresults = load(results(isves).name); vesresults = vesresults.link_results_struct;   
 
pv = vesresults.phasevec;
pvn = neuresults.n_phasevec;

%Select a vessel of interest and look at statistics.
figure(f1);
[vesx,vesy] = ginput(1);
vestestind = sub2ind(size(Mask1),round(vesy),round(vesx));
%Find ind in skeleton labels
isves = ismember(toplot.mask_ind,vestestind);
ves_skel_lab = toplot.skel_label(find(isves)); %Skeleton label for this vessel.

for i = 1:length(vesresults)
    skel_labtmp = vesresults(i).skel_lab;
    ismemtmp = ismember(skel_labtmp,ves_skel_lab);
    if ~isempty(find(ismemtmp,1))
        foundskel_lab = [i,find(ismemtmp)]; %fist is seg second is link number
        disp('found skeleton seg');
    end
end

ktmp = pv(foundskel_lab(2),1);
Rtmp = pv(foundskel_lab(2),2);
vesR2 = Rtmp^2
phitmp = vesresults(1).phi_mat(:,foundskel_lab(2));
phitmp = phitmp(~isnan(phitmp));
dtmp = vesresults(1).dist_mat(:,foundskel_lab(2));
dtmp = dtmp(~isnan(dtmp));
pix_mm = dtmp(end)/vesresults(1).link_lengths_mm(foundskel_lab(2));

figure
scatter(dtmp/pix_mm,phitmp,'filled');
xlabel('Dist (mm)','Interpreter','latex');
ylabel('Phase (rad)','Interpreter','latex');
cd('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\AllSegments\ManuscriptPlots\Supp6');
print(gcf, '-depsc2', 'ExamplePhaseProgression_11_10_02_TB201021M3');

kntmp = pvn(foundskel_lab(2),1);
Rntmp = pvn(foundskel_lab(2),2);
neuR2 = Rntmp^2
phintmp = neuresults(1).n_phi_mat(:,foundskel_lab(2));
phintmp = phintmp(~isnan(phintmp));
dntmp = neuresults(1).n_dist_mat(:,foundskel_lab(2));
dntmp = dntmp(~isnan(dntmp));

hold on
scatter(dntmp/pix_mm,phintmp,'filled');

%Plot neural envelope on top of cropped image.
figure(f2);
Mask1tmp = Mask1;
for i = 1:27
    Mask1tmp(neuresults(1).inds{i,foundskel_lab(2)}) = (i/27)*150;
end
figure
Mask1tmpCrop = imcrop(Mask1tmp,pos);
imshow(Mask1tmpCrop,[0 256]);
cd('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\AllSegments\ManuscriptPlots\Supp6');
print(gcf, '-depsc2', 'NeuEnvelopeExample_11_10_02_TB201021M3');

%Plot vessel segments
% Mask1tmpves = Mask1;
Mask1tmpves = repmat(Mask1,[1,1,3]);
%Get all skeleton labels and inds.
for i = 1:length(vesresults)
    skels(i) = vesresults(i).skel_lab(foundskel_lab(2));
    if skels(i) ~= 0
        skelinds = ismember(toplot.skel_label,skels(i));
        maskinds{i} = toplot.mask_ind(skelinds);
    end
end
skels = skels(~isnan(skels));


vesmask = zeros(size(Mask1tmpves,1),size(Mask1tmpves,2));
% alpha = 1;
numcolors = 3;
% linescmap = lines(numcolors);
% linescmap = prism(numcolors);
linescmap = (1/255*[24,69,59;123,189,0;255,255,255]);
for i = 1:length(maskinds)
    colornum = mod(i,numcolors)+1; %Need to order maskinds in distance from starting point
    colortmp = round(linescmap(colornum,:)*256);
    for j = 1:length(maskinds{i})
        [tmpvescol,tmpvesrow] = ind2sub([size(Mask1,1),size(Mask1,2)],maskinds{i}(j));
        Mask1tmpves(tmpvescol,tmpvesrow,:) = colortmp;
%         alphamask(tmpvescol,tmpvesrow) = alpha;
        vesmask(tmpvescol,tmpvesrow) = 1;
    end
end
ves_alpha = double(Mask1).*double(vesmask);
alphavals = ves_alpha(ves_alpha ~= 0);
ves_alpha = ves_alpha - min(alphavals);
ves_alpha = ves_alpha./max(ves_alpha(:))*0.7;

Mask1Crop = imcrop(Mask1,pos);
Mask1tmpvesCrop = imcrop(Mask1tmpves,pos);
ves_alphaCrop = imcrop(ves_alpha,pos);
vesmask_Crop = imcrop(vesmask,pos);

figure
imshow(Mask1Crop);
hold on;
h = imshow(Mask1tmpvesCrop);








for i = 1:length(maskinds)
    Mask1tmpves(maskinds{i}) = (i/length(maskinds))*150;
end
figure
Mask1tmpvesCrop = imcrop(Mask1tmpves,pos);

imshow(Mask1tmpvesCrop,[0 256]);
cd('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\AllSegments\ManuscriptPlots\Supp6');
print(gcf, '-depsc2', 'VesSegsExample_11_10_02_TB201021M3');


%% Supplemental Fig 5 - all of Thomas's Magnitude and Phase at vasomotion frequency

%Need to get all data used - Thomas's data used in traveling wave analysis
%and data used for Xinyue's predictive analysis
clear; clc;

cd('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\AllSegments\1_30_23_Results');
allfiles = dir(pwd);
todel = zeros(length(allfiles),1);
for i = 1:length(allfiles)
    if and(~contains(allfiles(i).name,'JD'),~contains(allfiles(i).name,'TB')) || contains(allfiles(i).name,'JD221024F2') || contains(allfiles(i).name,'JD221024F3')  %Delete under these conditions. Last animal is open window.
        todel(i) = 1;
    end
end
allfiles(logical(todel)) = [];

%Iterate over all trials analyzed, get folder locations where toplot file
%is.
WdfSummary = struct(); %struct to keep track of all results over animals and trials.
counter = 1;
for i = 1:length(allfiles)
% for i = 6:6
    cd([allfiles(i).folder,'\',allfiles(i).name]);
    resultfiles = dir('*\*results.mat');
    todel = zeros(length(resultfiles),1);
    for j = 1:length(resultfiles)
        if contains(resultfiles(j).name,'neural')
            todel(j) = 1;
        end
    end
    resultfiles(logical(todel)) = [];

    for j = 1:length(resultfiles) %iterate over result files within this animal
%     for j = 25:length(resultfiles) %iterate over result files within this animal
        load([resultfiles(j).folder,'\',resultfiles(j).name])
        WdfSummary(counter).results = resultfiles(j).name;
        WdfSummary(counter).resultsfolder = resultfiles(j).folder;
        WdfSummary(counter).datafolder = link_results_struct(1).folder;
        WdfSummary(counter).ParamFile = link_results_struct(1).file;

        cd(link_results_struct(1).folder) %Change directory to data location, look for toplot file with mag/phase info.
        toplot_files = dir('*toplot.mat'); %mag/phase could be saved in either of these file names
        mat_files = dir('*beforepuff.mat');
        
        %First look in toplot files
        toplot_todel = zeros(length(toplot_files),1);
        for k = 1:length(toplot_files)
            if contains(toplot_files(k).name,'_puff_') || contains(toplot_files(k).name,'neurons')
                toplot_todel(k) = 1;
            end
        end
        toplot_files(logical(toplot_todel)) = [];
        if ~isempty(toplot_files) && length(toplot_files) == 1 %Get magnitude and phase plotting info
            load([toplot_files(1).folder,'\',toplot_files(1).name])
            WdfSummary(counter).HasToplot = 1;
            %Check if toplot struct has necessary fields
            if isfield(toplot,'U') %This is phase and magnitude info
                WdfSummary(counter).ToplotHasMagPhase = 1;
                WdfSummary(counter).ToplotFile = toplot_files(1).folder;
                WdfSummary(counter).ToplotName = toplot_files(1).name;
            else
                WdfSummary(counter).ToplotHasMagPhase = 0;
            end
        else
            WdfSummary(counter).HasToplot = 0;
            WdfSummary(counter).ToplotHasMagPhase = 0;
        end
        %If toplot file doesn't exist or if toplot doesnt have mag/phase
        %info, check the .mat file.
        if and(isempty(toplot_files),~isempty(mat_files)) || and(~isempty(mat_files),WdfSummary(counter).ToplotHasMagPhase == 0)
            load([mat_files(1).folder,'\',mat_files(1).name])
            if isfield(toplot,'ophase') && isfield(toplot,'U') %This is phase and magnitude info
                WdfSummary(counter).MatHasMagPhase = 1;
            else
                WdfSummary(counter).MatHasMagPhase = 0;
            end
        end
        %Record if there is a beforepuff mat file or not
        if ~isempty(mat_files) && length(mat_files) == 1 %Get magnitude and phase plotting info
            WdfSummary(counter).HasMat = 1;
            WdfSummary(counter).MatFile = mat_files(1).folder;
            WdfSummary(counter).ToplotName = mat_files(1).name;
        else
            WdfSummary(counter).HasMat = 0;
            WdfSummary(counter).MatHasMagPhase = 0;
        end



    clearvars toplot
    [i,j]
    counter = counter + 1;
    end
end

%% Reanalyze necessary files
clear; clc; close all
load("X:\DataAnalysis\VesCorrPhase\AllSegments\WdfSummary_Mag_Phase.mat");
for i = 1:length(WdfSummary)
    if isempty(WdfSummary(i).MatHasMagPhase)
        WdfSummary(i).MatHasMagPhase = 0;
    end
end

toplottoanalyze = ~vertcat(WdfSummary.ToplotHasMagPhase);
mattoanalyze = ~vertcat(WdfSummary.MatHasMagPhase);
toanalyze = and(toplottoanalyze,mattoanalyze);
num_filestoreanalyze = sum(toanalyze)

toanalyze = find(toanalyze);
for iter = 1:length(toanalyze)

    file = toanalyze(iter);
    cd(WdfSummary(file).datafolder)
    load(WdfSummary(file).ToplotName)
    %Analyze trial for dominant mode, magnitude, phase
    wavefiles = dir('*beforepuff_wave.h5');
    if length(wavefiles) == 1
        wave = h5read(wavefiles(1).name,'/wave');
    end

    [U,S,V]=svd(wave,0);
    sig_modes=length(V)/2;
    Un=single(U(:,1:sig_modes));
    Sn=single(S(1:sig_modes,1:sig_modes));
    Vn=single(V(:,1:sig_modes));
    clear wave

    filename = toplot.fname;
    ext=extractBetween(filename,'_0.','Hz');
    fpeak2=str2double(erase(['0.',ext{1}],'_'));
    if size(ext,1)>2
        fpeak3=str2double(erase(['0.',ext{2}],'_'));
        fpeak4=str2double(erase(['0.',ext{3}],'_'));
        fpeak5=fpeak2*2;
        toplot.f_peak(4)=fpeak4;
        toplot.f_peak(5)=fpeak5;
    elseif size(ext,1)>1
        fpeak4=fpeak2*2;
        fpeak3=str2double(erase(['0.',ext{2}],'_'));
        toplot.f_peak(4)=fpeak4;
    else
        fpeak3=fpeak2*2;
    end
    toplot.f_peak(2)=fpeak2;
    toplot.f_peak(3)=fpeak3;
    padding_ratio = 2; %padded time base is 4 times higher then 2^nextpow2(size(wave,2));
    toplot.Delta_f = 0.02; % Hz. denoted as half-bandwidth
    if toplot.f_peak(2)<0.05
        toplot.Delta_f = 0.01
    elseif size(toplot.f_peak,2)>3 && abs(toplot.f_peak(2)-toplot.f_peak(3))<2*toplot.Delta_f
        toplot.Delta_f=(abs(toplot.f_peak(2)-toplot.f_peak(3))-0.002)/2
    end
    num_pixel = size(Un,1);
    num_frame= size(Vn,1);
    padding_ratio = 2; %padded time base is 4 times higher then 2^nextpow2(size(wave,2));
    num_frame_pad = (2 ^ ceil(log2(num_frame))) * padding_ratio;
    toplot.pad = num_frame_pad;
    p = round(num_frame * toplot.Delta_f / toplot.rate); % (num_frame / acq_frequency = acquisition time, i.e. T in the paper)
    toplot.Delta_f = p * toplot.rate / num_frame;
    disp(['Bandwidth = ', num2str(toplot.Delta_f), ' Hz'])
    toplot.num_tapers = 2 * p - 1;
    [slep,~] = dpss(num_frame, p, toplot.num_tapers);


    addpath(genpath('C:\chronux_2_12'))
    ntapers = [(toplot.num_tapers+1)/2,toplot.num_tapers];
    Fs = toplot.rate;
    nfft = toplot.pad;
    tapers = dpsschk(ntapers,size(Vn,1),Fs);
    [f,toplot.findx] = getfgrid(Fs,nfft,[0,Fs/2]);
    %FFT
    taperedFFT = complex(zeros(length(f),ntapers(2),size(Vn,2)));
    for i = 1:size(Vn,2)
        J = mtfftc(Vn(:,i),tapers,nfft,Fs);
        J = J(toplot.findx,:,:);
        taperedFFT(:,:,i) = J;
    end
    scores = Un*Sn;
    clear Un Sn i J

    S_tot = zeros(size(scores,1),size(taperedFFT,1),'single');
    for k = 1:ntapers(2)
        z = scores*squeeze(taperedFFT(:,k,:))';
        S_tot = S_tot + conj(z).*z;
    end
    S_tot = S_tot./ntapers(2);
    clear Fs ntapers nfft tapers findx
    [num_pixel] = size(toplot.mask_ind,1);
    num_frame_pad=size(S_tot,2);
    toplot.mpowr=[];
    toplot.mpowr=mean(S_tot,1);
    toplot.f = [];
    toplot.f = f;
    figure('units','normalized','outerposition',[0 0 1 1]);
    plot(f, log10(mean(S_tot, 1)), 'ko');

    [~, toplot.pwr] = max(S_tot(:, toplot.wind(1):toplot.wind(2)), [], 2); %toplot.pwr is frequency at max of mpowr within toplot.wind
    toplot.pwr = toplot.f(toplot.pwr + toplot.wind(1) - 1);

    disp(' Done');

    ntapers = [(toplot.num_tapers+1)/2,toplot.num_tapers];
    Fs = toplot.rate;
    nfft = toplot.pad;

    tapers = dpsschk(ntapers,size(V,1),Fs);
    [f,findx] = getfgrid(Fs,nfft,[0,Fs/2]);

    %FFT
    clear taperedFFT S_tot wave z
    taperedFFT = complex(zeros(length(f),ntapers(2),size(V,2)));
    for i = 1:size(V,2)
        J = mtfftc(V(:,i),tapers,nfft,Fs);
        J = J(findx,:,:);
        taperedFFT(:,:,i) = J;
    end

    scoresC = U*S;
    clear U S i J tapers
    if isunix
        rmpath(genpath(''));
    else
        rmpath(genpath('C:\chronux_2_12'));
    end



    % Calculate Coherence versus frequency
    plot_num_f_points = fix(4*toplot.rate/2 / toplot.Delta_f); % 4X oversampling
    toplot.coherence = zeros(plot_num_f_points,1);
    toplot.f_vector = linspace(0,toplot.rate/2,plot_num_f_points);
    tic
    interp_FFT = interp1(f,reshape(taperedFFT,size(taperedFFT,1),[]),toplot.f_vector);
    interp_FFT = reshape(interp_FFT,[],size(taperedFFT,2),size(taperedFFT,3));
    parfor i = 1:length(toplot.f_vector)
        m = scoresC*squeeze(interp_FFT(i,:,:))';
        s = svd(m,0);
        coherence(i) = squeeze(s(1))^2./sum(s.^2); %"Global coherence" in Observed brain dynamics p210
    end
    toplot.coherence = coherence;
    toc
    disp('Done')

    clear i m s interp_FFT coherence
    % Make U and save Prameters
    tic
    close all
    f_global = toplot.f_peak;
    interp_FFT = interp1(f,reshape(taperedFFT,size(taperedFFT,1),[]),f_global);
    interp_FFT = reshape(interp_FFT,[],size(taperedFFT,2),size(taperedFFT,3));
    toplot.U = zeros(length(toplot.skel_label),toplot.num_tapers,length(toplot.f_peak));
    %toplot.U = zeros(size(toplot.mask_ind,2),toplot.num_tapers,length(toplot.f_peak));
    for i = 1:length(toplot.f_peak)
        m = scoresC*squeeze(interp_FFT(i,:,:))';
        [u,~,~] = svd(m,0);
        toplot.U(:,:,i) = u(toplot.skel_label, :);
        %globalU(skel_label_unique,:,i) = u;
    end
    toplot.U=permute(toplot.U, [3,1,2]);
    clear i m s interp_FFT f_global
    toc

    %Save toplot
    if exist('SensorData','var')
        if ~contains(WdfSummary(file).ToplotName,'toplot')
            save([erase(WdfSummary(file).ToplotName,'.mat'),'_toplot.mat'], 'toplot' , 'recprams' ,'Sensorprams','SensorData', 'aprams', 'Vn','scores');
        elseif contains(WdfSummary(file).ToplotName,'toplot')
            save(WdfSummary(file).ToplotName, 'toplot' , 'recprams' ,'Sensorprams','SensorData', 'aprams', 'Vn','scores');
        end

    else
        if ~contains(WdfSummary(file).ToplotName,'toplot')
            save([erase(WdfSummary(file).ToplotName,'.mat'),'_toplot.mat'], 'toplot' , 'recprams' , 'Vn','scores');
        elseif contains(WdfSummary(file).ToplotName,'toplot')
            save(WdfSummary(file).ToplotName, 'toplot' , 'recprams' , 'Vn','scores');
        end
    end
    disp('Done Saving')
    file
    close all

end


%% Plot all resuls
%Get all files analyzed and ready to plot.
%Center phases using mean angle (tphase), plot on non-linear cyclic colormap
%(tanh) to accentuate phase differences.
clear; clc; close all
load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\AllSegments\WdfSummary_Mag_Phase_8_30_23.mat");
for i = 1:length(WdfSummary)
    if isempty(WdfSummary(i).MatHasMagPhase)
        WdfSummary(i).MatHasMagPhase = 0;
    end
end

vertnr = 7;
horznr = 4;
totplots = vertnr*horznr;
addpath('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\subtightplot\')
fignum = 1;
plotnum = 1;
filecounter = 1;
%%
% for file = 1:length(WdfSummary)
for file = 1:length(WdfSummary)+1
    clearvars -except file WdfSummary vertnr horznr totplots fignum plotnum filecounter

    if mod(file,totplots) == 1
        fig=figure('units','inches','outerposition',[0 0 8.5 11]); hold on;
    end
    if file == 1 %Plot colorbar as first plot
        ha = subtightplot(vertnr,horznr,plotnum,[.01 .01],[.01 .01],[.01 .01]);
        colormap hsv;
        caxis([-pi pi]);
        c = colorbar('north','Ticks',[-pi -pi/2 0 pi/2 pi],'TickLabels',{'-\pi','-\pi/2','0','\pi/2','\pi'},'FontSize',10,'TickLength',0.03);
        ha.Visible = 'off';
        figposition = ha.Position;
        cheight = c.Position(4);
        c.Position = [figposition(1) figposition(2)+(figposition(4)/2) figposition(3) cheight];
        plotnum = plotnum + 1;
        pause(1);
    else

        %Load toplot file and necessary info for plotting
        cd(WdfSummary(filecounter).datafolder)
        if ~contains(WdfSummary(filecounter).ToplotName,'toplot')
            load([erase(WdfSummary(filecounter).ToplotName,'.mat'),'_toplot.mat']);
        elseif contains(WdfSummary(filecounter).ToplotName,'toplot')
            load(WdfSummary(filecounter).ToplotName);
        end
        

        ha = subtightplot(vertnr,horznr,plotnum,[.01 .01],[.01 .01],[.01 .01]);
        %     subtightplot(vertnr,horznr,plotnum,[0 0],[0 0],[0 0]);

        phase=squeeze(angle(toplot.U(1, :, 1))); %tmp_taper is last index
        tphase=angle(sum(toplot.U(1, :, 1))); %Angle of summed vector
        ophase = mod(squeeze(phase-tphase)-pi,2*pi)-pi; %Must be between -pi and pi (inclusive)

        map = zeros(toplot.mask_size);
        map(toplot.mask_ind) = ophase;
        rescaled_pwr = nan(toplot.mask_size);
        tmp_pixel_value = map(toplot.mask_ind);
        % tmp_pixel_value = min(tmp_maxphase, max(tmp_minphase , tmp_pixel_value) ); %Should not do this, is cutting off data.
        rescaled_pwr(toplot.mask_ind) = tmp_pixel_value;


        cmap = colormap('hsv');
        % int_to_cmap = linspace(0,1,size(cmap,1));
        int_to_cmap = linspace(-pi,pi, size(cmap,1)); %256 values between -pi and pi
        non_nan_ind = find(map);
        num_nonnan = numel(non_nan_ind);
        rbg_pwr_list = zeros(num_nonnan, 3);
        rbg_pwr_list(:, 1) = interp1(int_to_cmap, cmap(:, 1), rescaled_pwr(non_nan_ind)); %Interpolate cmap values (assigned to values between -pi to pi) to arbitrary phase values from ophase.
        rbg_pwr_list(:, 2) = interp1(int_to_cmap, cmap(:, 2), rescaled_pwr(non_nan_ind));
        rbg_pwr_list(:, 3) = interp1(int_to_cmap, cmap(:, 3), rescaled_pwr(non_nan_ind));


        im_size = size(toplot.mask);
        % rgb_pwr = zeros(3, prod(im_size));
        rgb_pwr = ones(3, prod(im_size)) * 0.5; %This sets the gray level of the background
        rgb_pwr(:, non_nan_ind) = rbg_pwr_list.'; %Plot phase values
        rgb_pwr = reshape(rgb_pwr, 3, im_size(1), im_size(2));
        toplot.rgb_phase = permute(rgb_pwr, [3, 2, 1]); %Put RGB data in dim 3
        clear i j pixcolor int_to_cmap
        ftmp = imagesc(flip(toplot.rgb_phase,2)); %Flip orientation
        axis image;
        box off;
        % colormap jet;
        colormap hsv
        caxis([-pi pi]);
        daspect([1,1,1]);
        axis off
        hg = hggroup;
        set(ftmp,'Parent',hg)
        str = sprintf('%.0f',filecounter);
        set(hg,'Displayname',str);
        %     legend(hg,'FontSize',12,'TextColor','white','Position',[0.85 0.5 0.1 0.1])
        leg = legend(hg,'FontSize',12,'TextColor','k','location','northeast');
        legend('boxoff');
        %     if mod(file,totplots) == 1
        %         legposition = leg.Position;
        %         figposition = ha.Position;
        %         newlegposition = [figposition(1)+figposition(3)-0.12 figposition(2)+figposition(4)-0.022 0.0750 0.0260];
        %         set(leg,'position',newlegposition)
        %     else
        legposition = leg.Position;
        figposition = ha.Position;
        newlegposition = [figposition(1)+figposition(3)-0.12 figposition(2)+figposition(4)-0.022 0.0750 0.0260];
        set(leg,'position',newlegposition)
        %     end


        if mod(plotnum,totplots) == 0 || filecounter == length(WdfSummary) 
            cd('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\AllSegments\ManuscriptPlots\Supp5');
            str = sprintf('Phase %.0f',fignum);
            print(gcf, '-depsc2', str);
            print(gcf,str,'-dpdf','-r0','-vector')
            fignum = fignum + 1;
            close all
            pause(1);
            plotnum = 1; %Reset this for next figure
        else
            plotnum = plotnum + 1;
        end
        filecounter = filecounter + 1;
    end
    file
end

% x = (1:length(cmap))/length(cmap);
% x = (x - mean(x))*2; %Center around 0
% % y = atanh(x);
% y = tanh(x);
% nonlincmap = interp1(y,cmap,x);
% 
% int_to_cmap = linspace(-pi,pi, size(nonlincmap,1));
% non_nan_ind = find(map);
% num_nonnan = numel(non_nan_ind);
% rbg_pwr_list = zeros(num_nonnan, 3);
% rbg_pwr_list(:, 1) = interp1(int_to_cmap, nonlincmap(:, 1), rescaled_pwr(non_nan_ind));
% rbg_pwr_list(:, 2) = interp1(int_to_cmap, nonlincmap(:, 2), rescaled_pwr(non_nan_ind));
% rbg_pwr_list(:, 3) = interp1(int_to_cmap, nonlincmap(:, 3), rescaled_pwr(non_nan_ind));
% 
% im_size = size(toplot.mask);
% % rgb_pwr = zeros(3, prod(im_size));
% rgb_pwr = ones(3, prod(im_size)) * 0.5;
% rgb_pwr(:, non_nan_ind) = rbg_pwr_list.';
% rgb_pwr = reshape(rgb_pwr, 3, im_size(1), im_size(2));
% toplot.rgb_phase = permute(rgb_pwr, [3, 2, 1]);
% clear i j pixcolor cmap int_to_cmap
% figure; imagesc(flip(toplot.rgb_phase,2));
% axis image;
% box off;
% % colormap jet;
% colormap(nonlincmap)
% caxis([-pi pi]);
% daspect([1,1,1]);
% axis off



% end
%% Sup fig 7, Xinyue's v2n summaries
close all
load("C:\Users\duckw\Downloads\rest_v2n_bar_summary\rest_v2n_bar_summary\h1.mat")
openfig("C:\Users\duckw\Dropbox\Thomas_SpaceTime_DropBox\Meeting Slides\Xinyue\Figure_save\fig\rest_v2n_bar_summary.fig")
legend off
%Calc mean and sd
v2nalldat = h1{2,1};
v2nmean = mean(v2nalldat(:),'omitnan');
v2nsd = std(v2nalldat(:),[],'omitnan');
hold on
yline(v2nmean);
yline(v2nmean - v2nsd);
yline(v2nmean + v2nsd);
cd('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\AllSegments\ManuscriptPlots\Supp7');
print(gcf, '-depsc2', 'beforepuffv2n_testmean_pm_SD');
%%
close all
load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\XinyueBarCharts\h1_n2v.mat")
openfig("C:\Users\duckw\Dropbox\Thomas_SpaceTime_DropBox\Meeting Slides\Xinyue\Figure_save\fig\stim_n2v_bar_summary.fig")
legend off
%Calc mean and sd
v2nalldat = h1_n2v{2,1};
v2nmean = mean(v2nalldat(:),'omitnan');
v2nsd = std(v2nalldat(:),[],'omitnan');
hold on
yline(v2nmean);
yline(v2nmean - v2nsd);
yline(v2nmean + v2nsd);
cd('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\AllSegments\ManuscriptPlots\Supp7');
print(gcf, '-depsc2', 'stimn2v_testmean_pm_SD');
%%
clear; clc;
close all
load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\XinyueBarCharts\h1_v2n.mat")
openfig("C:\Users\duckw\Dropbox\Thomas_SpaceTime_DropBox\Meeting Slides\Xinyue\Figure_save\fig\stim_v2n_bar_summary.fig")
legend off
%Calc mean and sd
v2nalldat = h1_v2n{2,1};
v2nmean = mean(v2nalldat(:),'omitnan');
v2nsd = std(v2nalldat(:),[],'omitnan');
hold on
yline(v2nmean);
yline(v2nmean - v2nsd);
yline(v2nmean + v2nsd);
cd('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\AllSegments\ManuscriptPlots\Supp7');
print(gcf, '-depsc2', 'stimv2n_testmean_pm_SD');



%%
%Raw v1*v2 vs distance
figure
scatter(CombinedCorrMat(:,5),CombinedCorrMat(:,4).*CombinedCorrMat(:,3),'filled','MarkerFaceAlpha',0.4); ylim([-30 30]); yline(0);
xlabel('Distance (mm)','Interpreter','latex'); ylabel('v1 * v2','Interpreter','latex');
title({'7 animals, 14 trials, 736 penetrating arteriole pairs', 'pair-wise velocity product vs distance'},'Interpreter','latex');
cd('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Plots_3_6_23')
savefig('V_Product_14Experiments.fig');

%v1*v2 filtered by phase uncertainty
numSDs = 2;
todel = zeros(length(CombinedCorrMat),1);
for i=1:length(CombinedCorrMat)
    if abs(CombinedCorrMat(i,6)) < numSDs * CombinedCorrMat(i,7) || abs(CombinedCorrMat(i,8)) < numSDs * CombinedCorrMat(i,9)
        todel(i) = 1;
    end
end
distvec = CombinedCorrMat(:,5);
vel_vec = CombinedCorrMat(:,4).*CombinedCorrMat(:,3);
distvec(logical(todel)) = [];
vel_vec(logical(todel)) = [];
figure
scatter(distvec,vel_vec,'filled','MarkerFaceAlpha',0.4); ylim([-30 30]); yline(0);
xlabel('Distance (mm)','Interpreter','latex'); ylabel('v1 * v2','Interpreter','latex');
title({'7 animals, 14 trials, 125 penetrating arteriole pairs with phase $>$ 2SD from 0', 'pair-wise velocity product vs distance'},'Interpreter','latex');
ylim([-10 10]);
cd('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Plots_3_6_23')
savefig('V_Product_14Experiments_2SDFilt.fig');

%Plot number of ups and downs
x = categorical(cellstr(['Bottom precedes top';'Top precedes bottom']));
kvec1 = CombinedKFmat(:,2);
dkvec = CombinedKFmat(:,3);
todel1 = abs(kvec1) < dkvec;
kvec1(todel1) = [];
kvec2 = CombinedKFmat(:,2);
todel2 = abs(kvec2) < 2*dkvec;
kvec2(todel2) = [];
vals = [sum(CombinedKFmat(:,2)<0),sum(CombinedKFmat(:,2)>0);sum(kvec1<0),sum(kvec1>0);sum(kvec2<0),sum(kvec2>0)];
figure;
b = bar(x,vals);
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center','VerticalAlignment','bottom');
xtips2 = b(2).XEndPoints;
ytips2 = b(2).YEndPoints;
labels2 = string(b(2).YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center','VerticalAlignment','bottom');
xtips3 = b(3).XEndPoints;
ytips3 = b(3).YEndPoints;
labels3 = string(b(3).YData);
text(xtips3,ytips3,labels3,'HorizontalAlignment','center','VerticalAlignment','bottom');
legend({'All Data','Filtered 1SD','Filtered 2SD'});
ylabel('Number of vessels','Interpreter','latex','FontSize',14);
title({'Penetrating vessel traveling direction','146 vessels measured'},'Interpreter','latex');
savefig('TravelingDirectionBar_14Experiments.fig')

%Calculate added velocity metric
AddV_num = (CombinedCorrMat(:,3) + CombinedCorrMat(:,4));
AddV = (CombinedCorrMat(:,3) + CombinedCorrMat(:,4))./sqrt(CombinedCorrMat(:,3).^2 + CombinedCorrMat(:,4).^2);
figure
scatter(CombinedCorrMat(:,5),AddV,'filled','MarkerFaceAlpha',0.4);
xlabel('Distance (mm)','Interpreter','latex'); ylabel('$(v1 + v2)/sqrt(v1^2 + v2^2)$','Interpreter','latex')
title({'6 animals, 12 trials, 625 penetrating arteriole pairs', 'pair-wise velocity sum, normalized'},'Interpreter','latex');
cd('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Plots_2_28_23')
savefig('V_Sum_12Experiments.fig');

AddV(logical(todel)) = [];
figure
scatter(distvec,AddV,'filled','MarkerFaceAlpha',0.4);
xlabel('Distance (mm)','Interpreter','latex'); ylabel('$(v1 + v2)/sqrt(v1^2 + v2^2)$','Interpreter','latex')
title({'6 animals, 12 trials, 625 penetrating arteriole pairs with phase $>$ 1SD from 0', 'pair-wise velocity sum, normalized'},'Interpreter','latex');
cd('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Plots_2_28_23')
savefig('V_Sum_12Experiments_1SDFilt.fig');

%Plot pair distance histogram
figure
histogram(CombinedCorrMat(:,5),'BinWidth',0.2); xlim([0 2]);
xlabel('Pair-wise Distance (mm)','Interpreter','latex'); ylabel('Pairs','Interpreter','latex');
title('Pair-wise distance distribution, 7 animals 736 pairs','Interpreter','latex');
cd('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Plots_3_6_23')
savefig('PairDistanceDistribution_14Experiments.fig');

figure
histogram(distvec,'BinWidth',0.2); xlim([0 2]);
xlabel('Pair-wise Distance (mm)','Interpreter','latex'); ylabel('Pairs','Interpreter','latex');
title({'Pair-wise distance distribution, 7 animals 736 pairs','Red,orange,purple = Phase filtered at 1, 1.5, 2 SD'},'Interpreter','latex');
cd('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Plots_2_28_23')
savefig('PairDistanceDistribution_14Experiments_1_2_SDfilt.fig');

%For each distance bin, calculate <sig1 sig2> (try 0.2mm bins to start)
%Use combined Corr mat, calculate <> values for each distance bin
binsz = 0.2; %mm
maxdist = max(CombinedCorrMat(:,5))
expecvec = zeros(9,3);
for i=1:9
    inbin1 = CombinedCorrMat(:,5) < i*binsz;
    inbin2 = CombinedCorrMat(:,5) > (i-1)*binsz;
    inbin = and(inbin1,inbin2);
    
%     spinvec = CombinedCorrMat(:,3);
    spinvec = sign(CombinedCorrMat(:,3).*CombinedCorrMat(:,4));
    spinvec = spinvec(inbin);
    meanspinvec = mean(spinvec);
    stdspinvec = std(spinvec);

    expecvec(i,1) = meanspinvec;
    expecvec(i,2) = stdspinvec;
    expecvec(i,3) = stdspinvec/sqrt(numel(spinvec));
    
    expecvec(i,4) = numel(spinvec);
end
%Make scatter from expecvec
x_dist = 0:binsz:binsz*(size(expecvec,1)-1);
x_dist = x_dist + binsz/2
y_expec = expecvec(:,1); %<sig1 * sig2>
err = expecvec(:,3); %SDM
figure
errorbarplot = errorbar(x_dist,y_expec,err,"o",'MarkerFaceColor','b','MarkerEdgeColor','none');
errobarplot.Color = [0,0,1];
errorbarplot.LineWidth = 1.5
xlim([0 1.8]); ylim([-1.2 1.2]);
xlabel('Distance (mm)','Interpreter','latex','FontSize',15);
ylabel('$\bigl \langle \sigma_i \sigma_j \bigr \rangle$','Interpreter','latex','FontSize',15)
% title({'Correlation function, 5 trials, 174 penetrating arteriole pairs','$\Bigl \langle \sigma_i \sigma_j \bigr \rangle$ = Corr(distance) = mean(sgn $\phi_i *$ sgn $\phi_j$) for each distance bin'},'Interpreter','latex','FontSize',13)
grid on
title({'Correlation function, 7 animals, 14 trials, 736 penetrating arteriole pairs','$\bigl \langle \sigma_i \sigma_j \bigr \rangle$ = Corr(distance) = mean(sgn $v_i *$ sgn $v_j$) for each distance bin $\pm$ SDM'},'Interpreter','latex','FontSize',13)
ax = gca;
ax.TickLabelInterpreter = 'latex';
cd('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Plots_3_6_23')
savefig('CorrVsDist_14Experiments.fig');

%Correlation function only for significant phase differences
numSDs = 1.65;
binsz = 0.2; %mm
maxdist = max(CombinedCorrMat(:,5))
expecvec = zeros(9,3);
CombCorrMat_tmp = CombinedCorrMat;
todel1 = abs(CombCorrMat_tmp(:,6)) < numSDs * CombCorrMat_tmp(:,7);
todel2 = abs(CombCorrMat_tmp(:,8)) < numSDs * CombCorrMat_tmp(:,9);
todel = or(todel1,todel2);
CombCorrMat_tmp(todel,:) = [];
for i=1:9
    inbin1 = CombCorrMat_tmp(:,5) < i*binsz;
    inbin2 = CombCorrMat_tmp(:,5) > (i-1)*binsz;
    inbin = and(inbin1,inbin2);
    
    spinvec = sign(CombinedCorrMat(:,3).*CombinedCorrMat(:,4));
    spinvec = spinvec(inbin);
    meanspinvec = mean(spinvec);
    stdspinvec = std(spinvec);

    expecvec(i,1) = meanspinvec;
    expecvec(i,2) = stdspinvec;
    expecvec(i,3) = stdspinvec/sqrt(numel(spinvec));
    
end
%Make scatter from expecvec
%For 1/12/23 analysis, only have one point greater than 1.6mm so delete
%this entry
% expecvec([8,9],:) = [];
x_dist = 0:binsz:binsz*(size(expecvec,1)-1);
x_dist = x_dist + binsz/2

y_expec = expecvec(:,1); %<sig1 * sig2>
err = expecvec(:,3); %SDM
figure
errorbarplot = errorbar(x_dist,y_expec,err,"o",'MarkerFaceColor','b','MarkerEdgeColor','b');
errobarplot.Color = [0,0,1];
errorbarplot.LineWidth = 1.5
xlim([0 1.8]); ylim([-1.2 1.2]);
xlabel('Distance (mm)','Interpreter','latex','FontSize',15);
ylabel('$\Biggl \langle \sigma_i \sigma_j \biggr \rangle$','Interpreter','latex','FontSize',15)
grid on
title({'Correlation function, 7 animals, 14 trials, 125 penetrating arteriole pairs phase $>$ 2SD from 0','$\bigl \langle \sigma_i \sigma_j \bigr \rangle$ = Corr(distance) = mean(sgn $v_i *$ sgn $v_j$) for each distance bin $\pm$ SDM'},'Interpreter','latex','FontSize',13)
ax = gca;
ax.TickLabelInterpreter = 'latex';
cd('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Plots_3_6_23')
savefig('CorrVsDist_14Experiments_2SDFilt.fig');

%Correlation function only for significant phase differences
%Plot scatter with percent chance of observed directions
numSDs = 1.65;
binsz = 0.2; %mm
maxdist = max(CombinedCorrMat(:,5))
expecvec = zeros(9,3);
CombCorrMat_tmp = CombinedCorrMat;
todel1 = abs(CombCorrMat_tmp(:,6)) < numSDs * CombCorrMat_tmp(:,7);
todel2 = abs(CombCorrMat_tmp(:,8)) < numSDs * CombCorrMat_tmp(:,9);
todel = or(todel1,todel2);
CombCorrMat_tmp(todel,:) = [];
for i=1:9
    inbin1 = CombCorrMat_tmp(:,5) < i*binsz;
    inbin2 = CombCorrMat_tmp(:,5) > (i-1)*binsz;
    inbin = and(inbin1,inbin2); %indices for pairs within distance bin
    
    spinvec = sign(CombinedCorrMat(:,3).*CombinedCorrMat(:,4));
    spinvec = spinvec(inbin);

    findmax = [sum(spinvec<0),sum(spinvec>0)];
    [foundmax,maxind] = max(findmax);
    if maxind == 1
        probsuccess = 0.4696; %No filtering
%         probsuccess = 0.4758; %Filter 1SD
%         probsuccess = 0.4870; %Filter 2SD
    elseif maxind == 2
        probsuccess = 0.5304; %No filtering
%         probsuccess = 0.5242; %Filter 1SD
%         probsuccess = 0.5130; %Filter 2SD
    end

%     prob(i) = 1 - binocdf(foundmax,numel(spinvec),probsuccess); %1-probability to see more than (foundmax) successes
    prob(i) = binopdf(foundmax,numel(spinvec),probsuccess); %Probability density at observed value only
    possiblevals = 1:numel(spinvec);
    pdfvals = binopdf(possiblevals,numel(spinvec),probsuccess);
    [~,pdfmaxind] = max(pdfvals);
    bpdf(i) = possiblevals(pdfmaxind); %Value with highest binopdf
    if maxind == 1
        mostlikelyavg(i) = ((numel(spinvec)-bpdf(i))-bpdf(i))/numel(spinvec);
        biolims = binoinv([0.025 0.975],numel(spinvec),probsuccess);
        limsavg(i,:) = ((numel(spinvec)-biolims)-biolims)./numel(spinvec);
    elseif maxind == 2
        mostlikelyavg(i) = (-(numel(spinvec)-bpdf(i))+bpdf(i))/numel(spinvec);
        biolims = binoinv([0.025 0.975],numel(spinvec),probsuccess);
        limsavg(i,:) = (-(numel(spinvec)-biolims)+biolims)./numel(spinvec);
    end

    meanspinvec = mean(spinvec);
    stdspinvec = std(spinvec);

    expecvec(i,1) = meanspinvec;
    expecvec(i,2) = stdspinvec;
    expecvec(i,3) = stdspinvec/sqrt(numel(spinvec));
    
end
%Make scatter from expecvec
%For 1/12/23 analysis, only have one point greater than 1.6mm so delete
%this entry
% expecvec([8,9],:) = [];
x_dist = 0:binsz:binsz*(size(expecvec,1)-1);
x_dist = x_dist + binsz/2
y_expec = expecvec(:,1); %<sig1 * sig2>
% err = expecvec(:,3); %SDM
figure
scatter(x_dist,y_expec,'filled','blue');
hold on
s1 = scatter(x_dist,mostlikelyavg,'_r')
s1.LineWidth = 1.5;
s2 = scatter(x_dist,limsavg,'_r')
% s2.SizeData = 0.75
xlim([0 1.8]); ylim([-1 1]); yline(0);
xlabel('Distance (binned at 0.2 mm)','Interpreter','latex','FontSize',15);
ylabel('$\biggl \langle \sigma_i \sigma_j \biggr \rangle$','Interpreter','latex','FontSize',15)
grid on
% prob = round(prob,3);
% b = num2str(prob'); c = cellstr(b);
% dy = 0.2; dx = 0.07; %Text displacement
% t = text(x_dist - dx,y_expec - dy,c,'FontWeight','bold');
% for i=1:9
%     t(i).Color = 'red';
% end
% title({'Correlation function, 7 animals, 14 trials, 125 penetrating arteriole pairs','Red = Probability of observed directions occuring given P(-), P(+), and num pairs in each bin','Phase difference filtered for those 2SD from 0'},'Interpreter','latex','FontSize',11')
title({'Correlation function, 7 animals, 14 trials, 146 penetrating arteriole pairs','95 percent of binomial distribution falls within red bounds'},'Interpreter','latex','FontSize',11')
cd('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Plots_3_6_23\PenVesselPlots3_6_23')
print(gcf, '-depsc2', 'CorrVsDist_95PctBinomialDistribution');
save('CorrVsDist_95PctBinomialDistribution.fig')

ax = gca;
ax.TickLabelInterpreter = 'latex';
cd('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Plots_3_6_23')
savefig('CorrVsDist_BinomialProb_14Experiments_2SDFilt.fig');

%Marginal histogram vasomotion frequencies
figure
h1 = histogram(CombinedKFmat(:,1),'BinWidth',0.01);
set(h1,'FaceAlpha',1);
xlabel('Frequency (Hz)','Interpreter','latex');
ylabel('Count','Interpreter','latex')
title({'Vasomotion frequency distribution 146 penetrating arterioles','Freq. at max power of shallow vessel segment spectrum'},'Interpreter','latex');
savefig('VasofreqMarginalHist.fig');
figure
h1 = histogram(abs(CombinedKFmat(:,2)),'BinWidth',0.1);
xlabel('Magnitude phase gradient (rad/mm)','Interpreter','latex');
ylabel('Count','Interpreter','latex')
xlim([0 4])
title('Phase gradient distribution 146 penetrating arterioles','Interpreter','latex');
savefig('PhasegradMarginalHist.fig');
print(gcf, '-depsc2', 'PhaseGradMarginalHist');
savefig('PhaseGradMarginalHist.fig')

%Plot F vs K
figure
scatter(abs(CombinedKFmat(:,2)),CombinedKFmat(:,1),'filled','MarkerFaceAlpha',1,'MarkerEdgeColor','none');
xlim([0 4]); ylim([0 0.2]);
xlabel('Magnitude K (rad/mm)','Interpreter','latex');
ylabel('Vessel peak frequency (Hz)','Interpreter','latex');
title('Penetrating arteriole f vs k, 7 animals, 14 trials, 146 penetrating arterioles','Interpreter','latex','FontSize',13)
ax = gca;
ax.TickLabelInterpreter = 'latex';
savefig('FvsK_14Experiments.fig');
print(gcf, '-depsc2', 'FvsK_14Experiments');

kvec = abs(CombinedKFmat(:,2));
dkvec = abs(CombinedKFmat(:,3));
weights = 1./(dkvec.^2);
fvec = CombinedKFmat(:,1);

%Only fit for k < 0.3 rad/mm
kvec2 = kvec;
fvec2 = fvec;
weights2 = weights;
kvec2(kvec>0.3) = []; %As done 7/25 (page 56)
fvec2(kvec>0.3) = [];
weights2(kvec>0.3) = [];
figure
scatter(kvec2,fvec2)

kvec = CombinedKFmat(:,2);
dkvec = CombinedKFmat(:,3);
todel = abs(kvec) < dkvec;
% todel = ~logical(CombinedKFmat(:,6));
fvec = CombinedKFmat(:,1);
kvec1 = kvec(~todel);
fvec1 = fvec(~todel);
kvec2 = kvec(todel);
fvec2 = fvec(todel);
figure
scatter(abs(kvec1),fvec1,'filled','MarkerFaceAlpha',0.4,'MarkerFaceColor','red');
hold on
scatter(abs(kvec2),fvec2,'filled','MarkerFAceAlpha',0.4,'MarkerFaceColor','red')
xlim([0 4]); ylim([0 0.2]);
xlabel('Magnitude K (rad/mm)','Interpreter','latex');
ylabel('Vessel peak frequency (Hz)','Interpreter','latex');
title({'Penetrating arteriole f vs k, 7 animals, 14 trials', 'blue = 94 penetrating arterioles k $>$ 1SD from 0'},'Interpreter','latex','FontSize',13)
savefig('FvsK_14Experiments_1SDFilt.fig');

figure
errorbar(abs(kvec),fvec,CombinedKFmat(:,3),'horizontal','.','Color',[0,0,1,1],'MarkerSize',10)
xlim([-0.5 8]); ylim([0 0.2]);
xlabel('Magnitude K (rad/mm)','Interpreter','latex');
ylabel('Vessel peak frequency (Hz)','Interpreter','latex');
title('Penetrating arteriole f vs k, 146 vessels +/- SD','Interpreter','latex')
savefig('FvsK_14Experiments_K1SDBars.fig');
print(gcf, '-depsc2', 'FvsK_14Experiments_K1SDBars');


%To do: filter by significant coherence and power ratio at 2 depths
%Plot %same direction within each distance bin
%Calculate theoretical errors at each depth given number of pairs (binomial
%distribution)

%Plot f vs k filtered by significant coherence, shaded by power metric
kvec = CombinedKFmat(:,2);
fvec = CombinedKFmat(:,1);
Ctodel = ~logical(CombinedKFmat(:,6));
kvec(Ctodel) = []; fvec(Ctodel) = [];
figure
scatter(abs(kvec),fvec,'filled','MarkerFaceAlpha',0.4,'MarkerEdgeColor','none');
xlim([0 4]); ylim([0 0.2]);
xlabel('Magnitude K (rad/mm)','Interpreter','latex');
ylabel('Vessel peak frequency (Hz)','Interpreter','latex');
title({'Penetrating arteriole f vs k, 7 animals, 14 trials, 138 penetrating arterioles','Only vessels with significant coherence at vasomotor freq between depths'},'Interpreter','latex','FontSize',13)
savefig('FvsK_14Experiments_CoherenceFilt.fig')

kvec = CombinedKFmat(:,2);
fvec = CombinedKFmat(:,1);
Ctodel = ~logical(CombinedKFmat(:,6));
kvec(Ctodel) = []; fvec(Ctodel) = [];
figure
scatter(abs(kvec),fvec,'filled','MarkerFaceAlpha',1,'MarkerEdgeColor','none');
xlim([0 4]); ylim([0 0.2]);
xlabel('Magnitude K (rad/mm)','Interpreter','latex');
ylabel('Vessel peak frequency (Hz)','Interpreter','latex');
title({'Penetrating arteriole f vs k, 7 animals, 14 trials, 138 penetrating arterioles','Only vessels with significant coherence at vasomotor freq between depths'},'Interpreter','latex','FontSize',13)
ax = gca;
ax.TickLabelInterpreter = 'latex';
print(gcf, '-depsc2', 'FvsK_14Experiments_CoherenceFilt');

%Shade points by diamdeep/diamshallow 
kvec = CombinedKFmat(:,2);
fvec = CombinedKFmat(:,1);
Ctodel = ~logical(CombinedKFmat(:,6));
AlphaData = (CombinedKFmat(:,5)./CombinedKFmat(:,4)).^(1/2);
kvec(Ctodel) = []; fvec(Ctodel) = []; AlphaData(Ctodel) = [];
AlphaData(AlphaData>1) = 1;
figure
s = scatter(abs(kvec),fvec,'filled');
xlim([0 4]); ylim([0 0.2]);
xlabel('Magnitude K (rad/mm)','Interpreter','latex');
ylabel('Vessel peak frequency (Hz)','Interpreter','latex');
s.MarkerFaceColor = 'b';
s.MarkerFaceAlpha = 'flat';
s.AlphaData = AlphaData;
title({'Penetrating arteriole f vs k, 7 animals, 14 trials, 138 penetrating arterioles','Vessels with significant coherence between depths at vasomotor freq','Transparency = $\sqrt{Power_d/Power_s}$ $\in$ [0,1] $\propto$ $DiamAmp_d / DiamAmp_s$'},'Interpreter','latex','FontSize',13)
savefig('FvsK_14Experiments_CFilt_AmpRatioTransparency.fig')

%Color points by diamdeep/diamshallow  
kvec = CombinedKFmat(:,2);
fvec = CombinedKFmat(:,1);
Ctodel = ~logical(CombinedKFmat(:,6));
AlphaData = (CombinedKFmat(:,5)./CombinedKFmat(:,4)).^(1/2);
kvec(Ctodel) = []; fvec(Ctodel) = []; AlphaData(Ctodel) = [];
AlphaData(AlphaData>1) = 1;
figure
s = scatter(abs(kvec),fvec,'filled');
xlim([0 4]); ylim([0 0.2]);
xlabel('Magnitude K (rad/mm)','Interpreter','latex');
ylabel('Vessel peak frequency (Hz)','Interpreter','latex');
s.MarkerFaceColor = 'r';
s.MarkerFaceAlpha = 'flat';
s.AlphaData = AlphaData;
title({'Penetrating arteriole f vs k, 7 animals, 14 trials, 138 penetrating arterioles','Vessels with significant coherence between depths at vasomotor freq','Transparency = $\sqrt{Power_d/Power_s}$ $\in$ [0,1] $\propto$ $DiamAmp_d / DiamAmp_s$'},'Interpreter','latex','FontSize',13)
savefig('FvsK_14Experiments_CFilt_AmpRatioTransparency.fig')

kvec = CombinedKFmat(:,2);
fvec = CombinedKFmat(:,1);
Ctodel = ~logical(CombinedKFmat(:,6));
AlphaData = (CombinedKFmat(:,5)./CombinedKFmat(:,4)).^(1/2);
kvec(Ctodel) = []; fvec(Ctodel) = []; AlphaData(Ctodel) = [];
AlphaData(AlphaData>1) = 1;
cmap = flipud(colormap('gray'));
AlphaData = round(AlphaData*256); %Now these are cmap indices for the points
c = cmap(AlphaData,:);
figure
s = scatter(abs(kvec),fvec,36,c,'filled');
ax = gca;
ax.Colormap = cmap;           
cb = colorbar;
cb.Label.String = '(Power_d / Power_s)^{1/2}';
cb.TickLabelInterpreter = 'latex';
xlim([0 4]); ylim([0 0.2]);
xlabel('Magnitude K (rad/mm)','Interpreter','latex');
ylabel('Vessel peak frequency (Hz)','Interpreter','latex');
s.MarkerFaceAlpha = 1;
title({'Penetrating arteriole f vs k, 7 animals, 14 trials, 138 penetrating arterioles','Vessels with significant coherence between depths at vasomotor freq','Color = $\sqrt{Power_d/Power_s}$ $\in$ [0,1] $\propto$ $DiamAmp_d / DiamAmp_s$'},'Interpreter','latex','FontSize',13)
ax = gca;
ax.TickLabelInterpreter = 'latex';
savefig('FvsK_14Experiments_CFilt_AmpRatioColor.fig');
% print(gcf, '-depsc2', 'FvsK_14Experiments_CFilt_AmpRatioColor');

%Make marginal histograms for f and k
kvec = CombinedKFmat(:,2);
fvec = CombinedKFmat(:,1);
Ctodel = ~logical(CombinedKFmat(:,6));
kvec(Ctodel) = []; fvec(Ctodel) = [];
figure
h1 = histogram(fvec,'BinWidth',0.01);
set(h1,'FaceAlpha',1);
h1.Orientation = 'horizontal';
ylabel('Frequency (Hz)','Interpreter','latex');
ylim([0 0.2])
xlabel('Count','Interpreter','latex')
ax = gca;
ax.TickLabelInterpreter = 'latex';
title({'Vasomotion frequency distribution 7 animals 138 penetrating arterioles','Freq. at max power of shallow vessel segment spectrum'},'Interpreter','latex');
savefig('VasofreqMarginalHist.fig');
print(gcf, '-depsc2', 'VasofreqMarginalHist');



figure
h1 = histogram(abs(kvec),'BinWidth',0.1);
xlabel('Magnitude phase gradient (rad/mm)','Interpreter','latex');
ylabel('Count','Interpreter','latex')
xlim([0 4])
ax = gca;
ax.TickLabelInterpreter = 'latex';
set(h1,'FaceAlpha',1);
title('Phase gradient distribution 7 animals 138 penetrating arterioles','Interpreter','latex');
savefig('PhasegradMarginalHist.fig');
print(gcf, '-depsc2', 'PhaseGradMarginalHist');

%Plot histogram of power ratio
kvec = CombinedKFmat(:,2);
fvec = CombinedKFmat(:,1);
Ctodel = ~logical(CombinedKFmat(:,6));
AlphaData = (CombinedKFmat(:,5)./CombinedKFmat(:,4)).^(1/2);
kvec(Ctodel) = []; fvec(Ctodel) = []; AlphaData(Ctodel) = [];
% AlphaData(AlphaData>1) = 1;
figure
h1 = histogram(AlphaData,'BinWidth',0.1);
xlabel('$\sqrt{Power_d/Power_s}$','Interpreter','latex');
ylabel('Count','Interpreter','latex')
xlim([0 1.4])
ax = gca;
ax.TickLabelInterpreter = 'latex';
set(h1,'FaceAlpha',1);
title({'Square root of Diameter power ratio, 7 animals 138 penetrating arterioles','$\sqrt{Power_d/Power_s}$ $\propto$ $DiamAmp_d / DiamAmp_s$','For plots color map, $\sqrt{Power_d/Power_s}$ $>$ 1 set to 1'},'Interpreter','latex');
savefig('PowerRatioHist.fig');
print(gcf, '-depsc2', 'PowerRatioHist');

%Plot penetrating f vs k on surface data
load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\AllSegments\1_30_23_Results\CombinedResults\pvcomb_vesselfv_tapha_01.mat",'pvcomb');
minlength = 0.75;
lengthvec = pvcomb(:,5)./pvcomb(:,end);
length_todel = lengthvec < minlength;
t_todel = zeros(size(pvcomb,1),1);
for j=1:length(t_todel) %default t values calculated for alpha = 0.05
    if abs(pvcomb(j,9))<pvcomb(j,10)% || abs(pvncomb(j,9))<pvncomb(j,10) %If t<tcrit can't reject H0
        t_todel(j) = 1;
    end
end
t_todel = logical(t_todel);
todel = or(length_todel,t_todel);
pvcomb(todel,:) = [];

figure
scatter(abs(pvcomb(:,1)),pvcomb(:,3),'filled','MarkerFaceColor','b','MarkerEdgeColor','none','MarkerFaceAlpha',0.1);
xlim([0 2]); ylim([0 0.18]);
xlabel('Phase gradient magnitude (rad/mm)','Interpreter','latex');
ylabel('Vessel peak vasomotor frequency (Hz)','Interpreter','latex');
kvec = CombinedKFmat(:,2);
fvec = CombinedKFmat(:,1);
Ctodel = ~logical(CombinedKFmat(:,6));
AlphaData = (CombinedKFmat(:,5)./CombinedKFmat(:,4)).^(1/2);
kvec(Ctodel) = []; fvec(Ctodel) = []; AlphaData(Ctodel) = [];
AlphaData(AlphaData>1) = 1;
hold on;
s = scatter(abs(kvec),fvec,'filled');
xlim([0 4]); ylim([0 0.18]);
xlabel('Magnitude K (rad/mm)','Interpreter','latex');
ylabel('Vessel peak frequency (Hz)','Interpreter','latex');
s.MarkerFaceColor = [0.4940, 0.1840, 0.5560]; %Purple
% s.MarkerFaceColor = [0.85 0.3250 0.0980]; %Orange
% s.MarkerFaceColor = [0.4660 0.6740 0.188]; %Green
s.MarkerFaceAlpha = 'flat';
s.AlphaData = AlphaData;
ax = gca; ax.TickLabelInterpreter = 'latex';
title({'Blue: 12 animals 6978 surface vessels','Purple: 7 animals 138 penetrating vessels, shaded by $\sqrt{Power_d/Power_s}$'},'Interpreter','latex');
savefig('SurfAndPenetrating_FvsK.fig')

figure
scatter(abs(pvcomb(:,1)),pvcomb(:,3),'filled','MarkerFaceColor','b','MarkerEdgeColor','none','MarkerFaceAlpha',1);
xlim([0 2]); ylim([0 0.18]);
xlabel('Phase gradient magnitude (rad/mm)','Interpreter','latex');
ylabel('Vessel peak vasomotor frequency (Hz)','Interpreter','latex');

kvec = CombinedKFmat(:,2);
fvec = CombinedKFmat(:,1);
Ctodel = ~logical(CombinedKFmat(:,6));
AlphaData = (CombinedKFmat(:,5)./CombinedKFmat(:,4)).^(1/2);
kvec(Ctodel) = []; fvec(Ctodel) = []; AlphaData(Ctodel) = [];
AlphaData(AlphaData>1) = 1;
cmap = flipud(colormap('gray'));
AlphaData = round(AlphaData*256); %Now these are cmap indices for the points
c = cmap(AlphaData,:);
hold on;
s = scatter(abs(kvec),fvec,30,c,'filled');
ax = gca;
ax.Colormap = cmap;           
cb = colorbar;
cb.Label.String = '(Power_d / Power_s)^{1/2}';
cb.TickLabelInterpreter = 'latex';
xlim([0 4]); ylim([0 0.2]);
xlabel('Magnitude K (rad/mm)','Interpreter','latex');
ylabel('Vessel peak frequency (Hz)','Interpreter','latex');
s.MarkerFaceAlpha = 1;
title({'Blue: 12 animals 6978 surface vessels, Gray: 7 animals 138 penetrating arterioles','Color = $\sqrt{Power_d/Power_s}$ $\in$ [0,1] $\propto$ $DiamAmp_d / DiamAmp_s$'},'Interpreter','latex','FontSize',13)
ax = gca;
ax.TickLabelInterpreter = 'latex';
% savefig('FvsK_14Experiments_CFilt_AmpRatioColor.fig');
print(gcf, '-depsc2', 'SurfAndPenetrating_FvsK');

%Plot penetrating and surface vessel marginal histograms
kvec = CombinedKFmat(:,2);
fvec = CombinedKFmat(:,1);
Ctodel = ~logical(CombinedKFmat(:,6));
kvec(Ctodel) = []; fvec(Ctodel) = [];
load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\AllSegments\1_30_23_Results\CombinedResults\pvcomb_vesselfv_tapha_01.mat",'pvcomb');
minlength = 0.75;
lengthvec = pvcomb(:,5)./pvcomb(:,end);
length_todel = lengthvec < minlength;
t_todel = zeros(size(pvcomb,1),1);
for j=1:length(t_todel) %default t values calculated for alpha = 0.05
    if abs(pvcomb(j,9))<pvcomb(j,10)% || abs(pvncomb(j,9))<pvncomb(j,10) %If t<tcrit can't reject H0
        t_todel(j) = 1;
    end
end
t_todel = logical(t_todel);
todel = or(length_todel,t_todel);
pvcomb(todel,:) = [];

figure
h1 = histogram(pvcomb(:,3),'BinWidth',0.01,'Normalization','probability');
h1.FaceColor = 'b';
set(h1,'FaceAlpha',0.5);
hold on;
h2 = histogram(fvec,'BinWidth',0.01,'Normalization','probability');
h2.FaceColor = 'r';
set(h2,'FaceAlpha',0.5);
% h1.Orientation = 'horizontal';
xlabel('Frequency (Hz)','Interpreter','latex');
xlim([0 0.2])
% ylabel('Count','Interpreter','latex')
ylabel('Probability','Interpreter','latex')
ax = gca;
ax.TickLabelInterpreter = 'latex';
% title('Surface (blue) and penetrating (red) peak vasomotor frequency','Interpreter','latex');
title({'Surface (blue n=6978) and penetrating (red n=138) peak vasomotor frequency','Probability distributions'},'Interpreter','latex');
savefig('SurfAndPenetrating_VasofreqProbHist.fig')
print(gcf, '-depsc2', 'SurfAndPenetrating_VasofreqProbHist');

figure
h1 = histogram(abs(pvcomb(:,1)),'BinWidth',0.1);
h1.FaceColor = 'b';
set(h1,'FaceAlpha',0.5);
hold on;
h2 = histogram(abs(kvec),'BinWidth',0.1);
h2.FaceColor = 'r';
set(h2,'FaceAlpha',0.5);
% h1.Orientation = 'horizontal';
xlabel('Phase gradient magnitude (rad/mm)','Interpreter','latex');
xlim([0 2])
ylim([0 2100])
% ylabel('Count','Interpreter','latex')
ylabel('Probability','Interpreter','latex')
ax = gca;
ax.TickLabelInterpreter = 'latex';
% title('Surface (blue) and penetrating (red) peak vasomotor frequency','Interpreter','latex');
title({'Surface (blue n=6978) and penetrating (red n=138) phase gradient'},'Interpreter','latex');
savefig('SurfAndPenetrating_PhaseGradHist.fig')
print(gcf, '-depsc2', 'SurfAndPenetrating_PhaseGradHist');

%Plot cumulative distributions
figure
[e1,x1] = ecdf(pvcomb(:,3));
plot(x1,e1,'b');
[e2,x2] = ecdf(fvec);
hold on; plot(e2,x2,'r');
ylim([0 0.18]);
ylabel('Frequency (Hz)','Interpreter','latex');
xlabel('Cumulative probability','Interpreter','latex');
title('Surface (blue n=6978) and penetrating (red n=138) vasomotor frequency cdf','Interpreter','latex');
ax = gca;
ax.TickLabelInterpreter = 'latex';
savefig('SurfAndPenetrating_VasofreqCDF.fig')
print(gcf, '-depsc2', 'SurfAndPenetrating_VasofreqCDF');

figure
[e1,x1] = ecdf(abs(pvcomb(:,1)));
plot(x1,e1,'b');
[e2,x2] = ecdf(abs(kvec));
hold on; plot(x2,e2,'r');
xlim([0 4]);
xlabel('Phase Gradient Magnitude (rad/mm)','Interpreter','latex');
ylabel('Cumulative probability','Interpreter','latex');
title('Surface (blue n=6978) and penetrating (red n=138) phase gradient cdf','Interpreter','latex');
ax = gca;
ax.TickLabelInterpreter = 'latex';
savefig('SurfAndPenetrating_PhaseGradCDF.fig')
print(gcf, '-depsc2', 'SurfAndPenetrating_PhaseGradCDF');


%% Old plotting

PM1 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20221212 WT_PA_7\Pen_K_mat.mat");
PM2 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20221129 WT_PA_6\Pen_K_mat.mat");
% PM3 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20221120 WT_PA_5\Pen_K_mat.mat");
PM4 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20221128 WT_PA_5\Pen_K_mat.mat");
PM5 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20221208 WT_PA_5\Pen_K_mat.mat");

PM1 = PM1.phi_mat; PM2 = PM2.phi_mat; PM3 = PM3.phi_mat; PM4 = PM4.phi_mat; PM5 = PM5.phi_mat;
CombinedPhiMat = [PM1;PM2;PM3;PM4;PM5];






%Histogram of velocities
binedges = -10:0.5:10;
figure
histogram(CombinedPhiMat(:,6),'BinEdges',binedges); xline(0,'r','LineWidth',2);
xlabel('Traveling wave velocity (mm/s)','Interpreter','latex'); ylabel('Count');
title({'3 animals, 5 trials, 42 penetrating arterioles','Traveling wave velocity, (-) is deep preceeds shallow'},'Interpreter','latex');

%Histogram of velocities filtered by phase uncertainty
todel = zeros(length(CombinedPhiMat),1);
for i=1:length(CombinedPhiMat)
    if abs(CombinedPhiMat(i,1)) < 2*CombinedPhiMat(i,3)
        todel(i) = 1;
    end
end
Vvec1 = CombinedPhiMat(:,6);
Vvec1(logical(todel)) = []; %27/42 are 1SD away from 0 in phase

figure
histogram(Vvec1,'BinEdges',binedges); xline(0,'r','LineWidth',2);
xlabel('Traveling wave velocity (mm/s)','Interpreter','latex'); ylabel('Count');
title({'3 animals, 5 trials, 27 penetrating arterioles with phase $>$ 1SD from 0','Traveling wave velocity, (-) is deep preceeds shallow'},'Interpreter','latex');

%Bar graph of ups and down, stacked by animal 
X = categorical({'Down (+)','Up (-)'});
%Make Y with dimensions [animal1#down,animal2#down,... ;
%animal1#up,animal2#up,...]
Y = zeros(2,5);
Y(1,1) = sum(PM1(:,1)>0)
Y(1,2) = sum(PM2(:,1)>0)
Y(1,3) = sum(PM3(:,1)>0)
Y(1,4) = sum(PM4(:,1)>0)
Y(1,5) = sum(PM5(:,1)>0)
Y(2,1) = sum(PM1(:,1)<0)
Y(2,2) = sum(PM2(:,1)<0)
Y(2,3) = sum(PM3(:,1)<0)
Y(2,4) = sum(PM4(:,1)<0)
Y(2,5) = sum(PM5(:,1)<0)

figure
bar(X,Y,'stacked'); ylabel('Count','Interpreter','latex');
title('Penetrating arteriole traveling wave direction, 5 animals 42 arterioles','Interpreter','latex');
lgd = legend('12/12 PA7','11/29 PA6','11/20 PA5','11/28 PA5','12/8 PA5');
title(lgd,'Trial & animal')

%Bar graph of ups and down after filtering, stacked by animal 
Y = zeros(2,5);
PM1tmp = PM1; PM1tmp(abs(PM1(:,1)) < PM1(:,3),:) = [];
PM2tmp = PM2; PM2tmp(abs(PM2(:,1)) < PM2(:,3),:) = [];
PM3tmp = PM3; PM3tmp(abs(PM3(:,1)) < PM3(:,3),:) = [];
PM4tmp = PM4; PM4tmp(abs(PM4(:,1)) < PM4(:,3),:) = [];
PM5tmp = PM5; PM5tmp(abs(PM5(:,1)) < PM5(:,3),:) = [];
Y(1,1) = sum(PM1tmp(:,1)>0)
Y(1,2) = sum(PM2tmp(:,1)>0)
Y(1,3) = sum(PM3tmp(:,1)>0)
Y(1,4) = sum(PM4tmp(:,1)>0)
Y(1,5) = sum(PM5tmp(:,1)>0)
Y(2,1) = sum(PM1tmp(:,1)<0)
Y(2,2) = sum(PM2tmp(:,1)<0)
Y(2,3) = sum(PM3tmp(:,1)<0)
Y(2,4) = sum(PM4tmp(:,1)<0)
Y(2,5) = sum(PM5tmp(:,1)<0)

figure
bar(X,Y,'stacked'); ylabel('Count','Interpreter','latex');
title({'Penetrating arteriole traveling wave direction', '5 animals 27 penetrating arterioles with phase $>$ 1SD from 0'},'Interpreter','latex');
lgd = legend('12/12 PA7','11/29 PA6','11/20 PA5','11/28 PA5','12/8 PA5');
title(lgd,'Trial & animal')

%For each distance bin, calculate <sig1 sig2> (try 0.2mm bins to start)
%Use combined Corr mat, calculate <> values for each distance bin
binsz = 0.2; %mm
maxdist = max(CombinedCorrMat(:,5))
expecvec = zeros(9,3);
for i=1:9
    inbin1 = CombinedCorrMat(:,5) < i*binsz;
    inbin2 = CombinedCorrMat(:,5) > (i-1)*binsz;
    inbin = and(inbin1,inbin2);
    
%     spinvec = CombinedCorrMat(:,3);
    spinvec = sign(CombinedCorrMat(:,3).*CombinedCorrMat(:,4));
    spinvec = spinvec(inbin);
    meanspinvec = mean(spinvec);
    stdspinvec = std(spinvec);

    expecvec(i,1) = meanspinvec;
    expecvec(i,2) = stdspinvec;
    expecvec(i,3) = stdspinvec/sqrt(numel(spinvec));
    
    expecvec(i,4) = numel(spinvec);
end
%Make scatter from expecvec
%For 1/12/23 analysis, only have one point greater than 1.6mm so delete
%this entry
% expecvec(end,:) = [];
x_dist = 0:binsz:binsz*(size(expecvec,1)-1);
x_dist = x_dist + binsz/2

y_expec = expecvec(:,1); %<sig1 * sig2>
err = expecvec(:,3); %SDM
figure
errorbarplot = errorbar(x_dist,y_expec,err,"o",'MarkerFaceColor','b','MarkerEdgeColor','none');
errobarplot.Color = [0,0,1];
errorbarplot.LineWidth = 1.5
xlim([0 1.8]); ylim([-1.2 1.2]);
xlabel('Distance (mm)','Interpreter','latex','FontSize',15);
ylabel('$\bigl \langle \sigma_i \sigma_j \bigr \rangle$','Interpreter','latex','FontSize',15)
% title({'Correlation function, 5 trials, 174 penetrating arteriole pairs','$\Bigl \langle \sigma_i \sigma_j \bigr \rangle$ = Corr(distance) = mean(sgn $\phi_i *$ sgn $\phi_j$) for each distance bin'},'Interpreter','latex','FontSize',13)
grid on
title({'Correlation function, 3 animals, 6 experiment days, 219 penetrating arteriole pairs','$\bigl \langle \sigma_i \sigma_j \bigr \rangle$ = Corr(distance) = mean(sgn $\phi_i *$ sgn $\phi_j$) for each distance bin $\pm$ SDM'},'Interpreter','latex','FontSize',13)
ax = gca;
ax.TickLabelInterpreter = 'latex';
cd('X:\DataAnalysis\VesCorrPhase\Rui_2P\Plots_1_13_23')
print(gcf, '-depsc2', 'CorrvsDistance_11Experiments');

%Correlation function only for significant phase differences
binsz = 0.2; %mm
maxdist = max(CombinedCorrMat(:,5))
expecvec = zeros(9,3);
CombCorrMat_tmp = CombinedCorrMat;
todel1 = abs(CombCorrMat_tmp(:,6)) < CombCorrMat_tmp(:,7);
todel2 = abs(CombCorrMat_tmp(:,8)) < CombCorrMat_tmp(:,9);
todel = or(todel1,todel2);
CombCorrMat_tmp(todel,:) = [];
for i=1:9
    inbin1 = CombCorrMat_tmp(:,5) < i*binsz;
    inbin2 = CombCorrMat_tmp(:,5) > (i-1)*binsz;
    inbin = and(inbin1,inbin2);
    
    spinvec = CombCorrMat_tmp(:,3);
    spinvec = spinvec(inbin);
    meanspinvec = mean(spinvec);
    stdspinvec = std(spinvec);

    expecvec(i,1) = meanspinvec;
    expecvec(i,2) = stdspinvec;
    expecvec(i,3) = stdspinvec/sqrt(numel(spinvec));
    
end
%Make scatter from expecvec
%For 1/12/23 analysis, only have one point greater than 1.6mm so delete
%this entry
expecvec([8,9],:) = [];
x_dist = 0:binsz:binsz*(size(expecvec,1)-1);
x_dist = x_dist + binsz/2

y_expec = expecvec(:,1); %<sig1 * sig2>
err = expecvec(:,3); %SDM
figure
errorbarplot = errorbar(x_dist,y_expec,err,"o",'MarkerFaceColor','b','MarkerEdgeColor','b');
errobarplot.Color = [0,0,1];
errorbarplot.LineWidth = 1.5
xlim([0 1.6]); ylim([-1.2 1.2]);
xlabel('Distance (mm)','Interpreter','latex','FontSize',15);
ylabel('$\Biggl \langle \sigma_i \sigma_j \biggr \rangle$','Interpreter','latex','FontSize',15)
title({'Correlation function, 5 trials, 70 penetrating arteriole pairs phase $>$ 1SD from 0','$\Biggl \langle \sigma_i \sigma_j \biggr \rangle$ = Corr(distance) = mean(sgn $\phi_i *$ sgn $\phi_j$) for each distance bin $\pm$ SDM'},'Interpreter','latex','FontSize',13)
grid on

%% Plot f vs k


figure
scatter(abs(CombinedKFmat(:,2)),abs(CombinedKFmat(:,1)),'filled','MarkerFaceAlpha',0.4);
xlim([0,2]); ylim([0 0.2])



%% Plots for david 1/13/23

%Raw data image - PA5 11/28/22 PA7
data_folder = '\\dk-server.dk.ucsd.edu\rul087\Data backup\20221129 WT_PA_6';

animal = '20221129PA6'; %For file naming

cd(data_folder);
files = dir('*001.tif');
files([files.bytes] < 9e8) = []; %Only get the large files
files([files.bytes] > 1e10) = []; %Delete large files (not image sequence)
todel = logical(zeros(length(files),1));
for i=1:length(files)
    if contains(files(i).name,'test')
        todel(i) = true;
    end
end
files(todel) = [];
file = 2; 
namestr = extractBefore(files(file).name,'_7.6x');
depth_str = extractAfter(namestr,'_');
depth1 = extractBefore(depth_str,'_')
depth2 = extractAfter(depth_str,'_')
PA = extractBefore(namestr,'_')

pix_um = 300/100; 
rate = 7.25; %Hz

%Load Images - takes  mins
tic
loadims = 1;
if loadims == 1
    im_cell = cell(11000,1);
    filesize = files(file).bytes;
    im_counter = 1;
    size_counter = 0;
    im_length = 100; %Double check this for each file
    tic
    for im_counter = 1:im_length
        tmp_im = imread([files(file).folder,'\',files(file).name],im_counter);
        im_cell{im_counter,1} = tmp_im;
        if mod(im_counter,400) == 0
            disp(im_counter)
        end
    end
    toc
else
end
toc
%
im_cell = im_cell(~cellfun('isempty',im_cell));
%Average across 1 frame
tokeep_tmp1 = zeros(1,size(im_cell,1));
tokeep_tmp2 = zeros(1,size(im_cell,1));
tokeep_tmp1(1:2:end) = 1;
tokeep_tmp2(2:2:end) = 1;

im_cell1 = im_cell(logical(tokeep_tmp1));
im_cell2 = im_cell(logical(tokeep_tmp2));

im_size = [size(im_cell1{1,1},1),size(im_cell1{1,1},2)];

avg_ves1 = zeros(im_size(1),im_size(2),length(im_cell1));
avg_ves2 = zeros(im_size(1),im_size(2),length(im_cell2));
for i=1:size(avg_ves1,3)
    avg_ves1(:,:,i) = im_cell1{i};
end
for i=1:size(avg_ves1,3)
    avg_ves2(:,:,i) = im_cell2{i};
end
clearvars -except avg_ves1 avg_ves2 files file data_folder im_size loadims animal PA depth1 depth2 pix_um rate
current_cell = 2;

avg_ves_cell = cell(2,1);
avg_ves_cell{1,1} = avg_ves1;
avg_ves_cell{2,1} = avg_ves2;
% clearvars avg_ves2 avg_ves1

avg_ves_tmp = avg_ves_cell{current_cell,1};
maskedImage_t = zeros(size(avg_ves_tmp));


%Define average centroid of vessel of interest
avg_im = mean(avg_ves_tmp,3);
figure; imagesc(avg_im); daspect([1,1,1]); title('Drag ROI');
roi = drawrectangle;
rect = round(roi.Position);
rect_mask = createMask(roi);
figure; imagesc(avg_im); daspect([1,1,1]); title('Click Vessel Center');
[windx,windy] = ginput(1);
% BW = createMask(roi); %Mask for vessel region


addpath(genpath('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P'))


%Plot them now
imnumber = 20
tmp_im = avg_ves1(:,:,imnumber);
tmp_im = medfilt2(tmp_im,[5,5]);
% tmp_im_crop = imcrop(tmp_im,rect);
tmp_im = rescale(tmp_im); %rescale tmp_im to [0 1];

%% Fit with TIRS
% current_cell =2; %Change this to change depth
thresholdFWHM = 0.40; %0.3 to 0.5 from paper (0.35 used in paper)
thresholdInv = 0.20; %0.1 to 0.3 from paper (0.2 used in paper)
medFiltRange = 5;

connectivity = 8;
%Crop image
if current_cell == 1
    imgSeq = avg_ves1;
elseif current_cell == 2
    imgSeq = avg_ves2;
end
CircVec = NaN(size(imgSeq,3),1);
radonRVec = NaN(size(imgSeq,3),1);

tic;
% for i=1:size(imgSeq,3)
for i = 1:imnumber
    % for i=1:50
    tmp_im = imgSeq(:,:,i);
    tmp_im = medfilt2(tmp_im,[medFiltRange,medFiltRange]);
    tmp_im = imcrop(tmp_im,rect);
    tmp_im = tmp_im - mean(tmp_im(:));
    if i==1
        tmp_cropped = zeros(size(tmp_im,1),size(tmp_im,2),size(imgSeq,3));
        thetaRange = 0:179;
        nAngles = length(thetaRange);
        nRadii = 2*ceil(norm(size(tmp_im) - ...
            floor((size(tmp_im) - 1)/2) - 1)) + 3;
        [nRows, nCols, nFrames] = size(tmp_cropped);
        areaPixels = zeros(nFrames, 1);
        imgNorm = zeros(nRadii, nAngles, nFrames);
        idxEdges = zeros(2, nAngles, nFrames);
        imgInv = zeros(size(tmp_cropped,1),size(tmp_cropped,2),size(imgSeq,3));
        vesselMask = false(size(tmp_cropped,1),size(tmp_cropped,2),size(imgSeq,3));
        vesselPerim = false(size(tmp_cropped,1),size(tmp_cropped,2),size(imgSeq,3));
    end

    tmp_cropped(:,:,i) = tmp_im;

    frameRate = rate;
    t0 = 0;
    frameTime = 1/frameRate;
    Radontime = ((0.5*frameTime):frameTime:(nFrames*frameTime))'  - t0;

    filterType = 'Hamming';
    % freqScaling = 1;
    outputSize = max([nRows, nCols]);
    argsIRadon = {thetaRange, filterType, outputSize};

    imgTrans = radon(tmp_im, thetaRange);
    [nRowsT, ~] = size(imgTrans);

    % Normalise each column (angle) of the radon transform
    imgNorm(:, :, i) = (imgTrans - ...
        repmat(min(imgTrans, [], 1), [nRowsT, 1]))./ ...
        repmat(max(imgTrans, [], 1) - min(imgTrans, [], 1), ...
        [nRowsT, 1]);

    % Threshold the image, fill in any holes, and extract out
    % only the largest contiguous area
    imgThresh = imgNorm(:, :, i) >= thresholdFWHM;
    imgThresh = imfill(imgThresh, 8, 'holes');
    ccT = bwconncomp(imgThresh, 4);
    ssT = regionprops(ccT);
    [~, maxIdx] = max([ssT(:).Area]);
    imgThresh = (labelmatrix(ccT) == maxIdx);

    for jAngle = 1:nAngles
        % Find the 'FWHM' edges of the transformed image
        idxEdgesTmp = [...
            find(imgThresh(:, jAngle), 1, 'first'), ...
            find(imgThresh(:, jAngle), 1, 'last')];

        % Manually threshold the transformed image using the
        % edges defined by the 'FWHM'
        imgThreshRowTmp = false(nRadii, 1);
        idxToUse = idxEdgesTmp(1) : idxEdgesTmp(2);
        imgThreshRowTmp(idxToUse) = true;
        imgThresh(:, jAngle) = imgThreshRowTmp;
        idxEdges(:, jAngle, i) = idxEdgesTmp;
    end


    doResizeInv = false;
    colsToUseInv = 1:nCols;
    rowsToUseInv = 1:nRows;
    if nRows ~= nCols
        doResizeInv = true;
        ctrOrig = floor(([nRows, nCols]+1)/2);
        ctrInv = floor((outputSize+1)/2);
        adjPixels = ctrInv - ctrOrig;
        hasRoundingProblem = all(adjPixels == 0);
        if hasRoundingProblem
            adjPixels = abs([nRows, nCols] - outputSize);
        end
        dimToAdj = find(adjPixels);
        if dimToAdj == 1
            rowsToUseInv = (1:nRows) + adjPixels(dimToAdj);
        else
            colsToUseInv = (1:nCols) + adjPixels(dimToAdj);
        end
    end

    % Invert the thresholded radon-transformed image, adjust
    % the size if necessary, and then normalise it
    imgInvTemp = iradon(imgThresh, argsIRadon{:}); %#ok<PFBNS>
    if doResizeInv
        imgInvTemp = imgInvTemp(rowsToUseInv, colsToUseInv);
    end
    imgInv(:,:,i) = imgInvTemp./max(imgInvTemp(:));

    % Threshold the inverted image, and fill in any holes in
    % the thresholded, inverted image
    imgInvThresh = imgInv(:,:,i) > thresholdInv;
    imgInvThresh = imfill(imgInvThresh, 'holes');

    % Calculate the area of the largest contiguous region
    cc = bwconncomp(imgInvThresh, connectivity);
    ss = regionprops(cc,'Area','Circularity');
    [areaPixels(i, 1), maxIdx] = max([ss(:).Area]);
    CircVec(i) = ss(maxIdx).Circularity; %Use this to check for good fit & threshold values.
    radonRVec(i) = mean(sum(imgThresh,1)); %Pixels

    % Create a binary image showing only the largest area
    % identified above
    lm = labelmatrix(cc);
    vesselMask(:, :, i) = (lm == maxIdx);
    vesselPerim(:, :, i) = bwperim(vesselMask(:,:,i)); %Perimeter pixels for plotting
    if mod(i,50)==0
        i/size(avg_ves_tmp,3)*100
    end
end
toc

areaUm2 = areaPixels.*((1/pix_um)^2); %Square um
rEq = (areaUm2/pi).^(1/2);
dEq = (4*areaUm2/pi).^(1/2);
%% %Now plot

figure
imagesc(rescale(medfilt2(imgSeq(:,:,imnumber),[medFiltRange,medFiltRange])));
daspect([1,1,1]);
%Plot radon fit
Perim = vesselPerim(:,:,imnumber);
[Px,Py] = ind2sub(size(Perim),find(Perim))
%Translate back to full image
Px = Px + rect(2) - 1;
Py = Py + rect(1) - 1;
hold on
plot(Py,Px,'r','LineStyle','none','Marker','.','MarkerSize',3);
%Plot scale bar
plot([20 50],[275 275],'k','Marker','none');
axis off;
title({'20221129WT6 PenA11 390$\mu$m depth, median filtered and TIRS fit','Scale bar 10$\mu$m'},'Interpreter','latex')

cd('X:\DataAnalysis\VesCorrPhase\Rui_2P\Plots_1_13_23')
print(gcf, '-depsc2', '2PImage_20221129_WT6_PenA11_390umDepth');

%% Plot spectral images
%1. Time series
allstats = load("X:\DataAnalysis\VesCorrPhase\Rui_2P\20221129 WT_PA_6\20221129PA6_PA11_10um_allstats.mat");
stats1 = allstats.allstats;
allstats = load("X:\DataAnalysis\VesCorrPhase\Rui_2P\20221129 WT_PA_6\20221129PA6_PA11_390um_allstats.mat");
stats2 = allstats.allstats;
time1 = stats1.time;
time2 = stats2.time;
if isequal(time1,time2)
    time2 = time2 + (1/stats2.rate)/2;
end
f1 = figure;
plot(stats1.time,stats1.RadondEq_Outl,'b'); hold on; plot(stats2.time,stats2.RadondEq_Outl,'r');
legend([depth1,' um'],[depth2,' um'],'Interpreter','latex'); ylabel('Diam (um)','Interpreter','latex'); xlabel('Time (s)','Interpreter','latex');
title('20221129WT6 PenA11 time series equivalent diameter','Interpreter','latex','FontSize',13); xlim([0 time2(end)]);
ylim([7 15])
f1.Position = [2538 901 942 207];
cd('X:\DataAnalysis\VesCorrPhase\Rui_2P\Plots_1_13_23')
print(gcf, '-depsc2', 'TimeSeriesDiam_20221129_WT6_PenA11');

%% 2.Spectra
waves = stats1.RadondEq_Outl;
waves = waves - mean(waves);
waved = stats2.RadondEq_Outl;
waved = waved - mean(waved);
times = stats1.time;
timed = stats2.time;
if sum(times - timed) == 0
    timed = timed + (1/rate)/2;
end

% figure
% plot(times,waves); hold on; plot(timed,waved);
%Interpolate each signal to get equal time points
tqs = times(1):(1/(2*rate)):times(end);
tqd = timed(1):(1/(2*rate)):timed(end);
waves_q_tmp = interp1(times,waves,tqs);
waved_q_tmp = interp1(timed,waved,tqd);
% figure
% plot(tqs,waves_q_tmp); hold on; plot(tqd,waved_q_tmp);

%Delete first wave1 value so that both start at t=0.75.
waves_q = waves_q_tmp(2:end);
tqs = tqs(2:end);
%Delete last wave1 value to make both time series the same duration
waved_q = waved_q_tmp;
waved_q(end) = [];
tqd(end) = [];
figure
plot(tqs,waves_q); hold on; plot(tqd,waved_q);

params.Fs = rate*2; %Interpolated rate is twice actual single depth rate
params.pad = 2;
params.fpass = [0 params.Fs/2]; %Hz, default is [0 Fs/2]
params.err   = [2 .05]; 
params.trialave = 0;
T = stats1.time(end);
BW = 0.02; %600s trial -> TBW =~ 7
params.tapers = [round(T*BW),round(2*T*BW-1)]; %Time-BW product and number of tapers
addpath(genpath('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\chronux_2_12'))

[S,f] = mtspectrumc(waves_q,params); %Why doesn't this have peak at 4.2 Hz? last file's spectrum did. Last file interpolation was different? (Relative time of when signal changes)
figure
plot(f,log10(S)); xlim([0 1]);
[wind,~] = ginput(2);
findf1 = round(wind(1),3);
findf2 = round(wind(2),3);
f1 = find(round(f,3)==round(findf1,3),1,'first'); %Hz lower bound
f2 = find(round(f,3)==round(findf2,3),1,'first'); %Hz upper bound
rmpath(genpath('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\chronux_2_12'))
[maxes,f_peaks] = findpeaks(S(f1:f2),f(f1:f2));
maxloc = find(maxes==max(maxes));
f_peak = f_peaks(maxloc)
freq_findx = max(find(round(f,3)==round(f_peak,3))); 
addpath(genpath('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\chronux_2_12'))

params_err = params;
params_err.err = [2,0.05];
[C,phi,S12,S1,S2,f,confC,phistd,Cerr] = coherencyc(waved_q,waves_q,params_err);

f2 = figure;
plot(f,log10(S2),'b'); xlim([0 1]); xlabel('Frequency (Hz)','Interpreter','latex'); ylabel('log10(Power)','Interpreter','latex');
% title([depth1,' um spectrum, fv = ',strfv,' Hz ',str_spectrum],'Interpreter','latex'); 
hold on
plot(f,log10(S1),'r'); xlim([0 1]); xlabel('Frequency (Hz)','Interpreter','latex'); ylabel('log10(Power)','Interpreter','latex');
title('20221129WT6 PenA11 Spectra','Interpreter','latex'); 
xline(f_peak,'-.','$f_v$','Interpreter','latex','Alpha',0.2); legend([depth1,' um'],[depth2,' um'],'Interpreter','latex');
hold off

cd('X:\DataAnalysis\VesCorrPhase\Rui_2P\Plots_1_13_23')
print(gcf, '-depsc2', 'Spectra_20221129_WT6_PenA11');

%3. Coherence
f3 = figure;
plot(f,C,'b'); hold on; plot(f,Cerr,'Color',[0,0,1,0.1]); xlim([0 rate/2]); ylim([0 1]); yline(confC,'-.','$|C|_{0.95}$','Interpreter','latex','Alpha',0.9,'LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom')
xlabel('Frequency (Hz)','Interpreter','latex'); ylabel('Magnitude Coherence','Interpreter','latex'); title('20221129WT6 PenA11 Coherence','Interpreter','latex');
xlim([0 0.5]);
cd('X:\DataAnalysis\VesCorrPhase\Rui_2P\Plots_1_13_23')
print(gcf, '-depsc2', 'Coherence_20221129_WT6_PenA11');

%4. Phase
f4 = figure;
plot(f,phi,'b'); hold on; plot(f,phi + 2*phistd','Color',[0,0,1,0.2]); plot(f,phi - 2*phistd','Color',[0,0,1,0.2]); xlim([0 0.5]); ylim([-pi pi]);
yline(0,'Alpha',0.9); xlabel('Frequency (Hz)','Interpreter','latex'); ylabel('Relative Phase (rad)','Interpreter','latex');
xline(f_peak,'-.','$f_v$','Interpreter','latex','Alpha',0.9);
%Calculate phase gradient
dx = (str2num(depth2) - str2num(depth1))/1000; %mm
phi_fv = phi(freq_findx);
k_fv = phi_fv / dx; %rad/mm
v = (2*pi*f_peak)/k_fv;
str1 = sprintf('%.2f',phi_fv);
str2 = sprintf('%.2f',abs(k_fv));
str3 = sprintf('%.2f',v);
title({'Relative Phase and 95 Percent CI',['$\phi_{fv}$ = ',str1,', $\mid k \mid$ = ',str2,' $v$ = ',str3]},'Interpreter','latex')

ylim([-pi/2,pi/2])
set(gca,'YTick',-pi/2:pi/2:pi/2)
ax = gca;
ax.YTickLabel = {'-\pi/2','0','\pi/2'}
cd('X:\DataAnalysis\VesCorrPhase\Rui_2P\Plots_1_13_23')
print(gcf, '-depsc2', 'Phase_20221129_WT6_PenA11');












%% 3/24 get vessel distance and plot correlation
clear
close all 
clc
vesmask = imread('wt_pa_5_mask.tif');
vesmask = logical(vesmask(:,:,1));
% mask2= bwskel(vesmask); %Use skeleton for distance
mask2 = vesmask; %Use total mask
figure
imagesc(mask2); %Vessel skeleton
[C,R] = ginput(1);
testdist = bwdistgeodesic(mask2,round(C),round(R));
figure; imagesc(testdist); %Check that distance transform worked.

%Click and save locations of all penetrating vessels
%For each location, calc gdistance mask, and calc vessel dist to every
%other penetrating vessel. 
%Re Calculate and plot correlation for each animal, combine.

phi_struct = load('phi_struct_20221208WT5.mat');
phi_struct = phi_struct.phi_struct;
CKFtmp = load('20221208WT5_Corr_KF_struct.mat');
Corr_KF_struct = CKFtmp.Corr_KF_struct;
CorrMat = Corr_KF_struct.CorrMat;
            
for i = 1:length(phi_struct)
    if i==1
        figure
        imagesc(mask2);
        daspect([1,1,1]);
    end
    disp(['Click PA',phi_struct(i).PA])
    pause();
    [C(i),R(i)] = ginput(1);
%     gdisttmp = bwdistgeodesic(mask2,C(i),R(i))
end
%%
mask = imread('wt_pa_5_mask.tif');
mask = mask(:,:,1);

figure; 
imagesc(mask); daspect([1,1,1]);
pause()
[cal_x,cal_y] = ginput(2); 
% prompt = "Input PA 13 x loc";
% dist_x(1) = input(prompt) %um (use this for PA 8-11)
% prompt = "Input PA 13 y loc";
% dist_y(1) = input(prompt)
% prompt = "Input PA 25 x loc";
% dist_x(2) = input(prompt)
% prompt = "Input PA 25 y loc";
% dist_y(2) = input(prompt)
% distRecord = sqrt((dist_x(2)-dist_x(1))^2 + (dist_y(2)-dist_y(1))^2); 
% distRecord = distRecord/1000; %mm

%(use this for PA 5-7)
prompt = "Insert PA 1 to 5 distance (mm) from Corr_KF_struct ";
distRecord = input(prompt);

distAiPix = sqrt((cal_x(2)-cal_x(1))^2 + (cal_y(2)-cal_y(1))^2); %Pixels in ilustrator image (units vessel-distance is in)
Ai_pix_mm = distAiPix/distRecord; %Pixels per mm in the vessel distance mask

%% Calculate vessel-distance and populate CorrMat
counter = 1;
for i=1:length(C)-1
    for j=(i+1):length(phi_struct)
        gdisttmp = bwdistgeodesic(mask2,round(C(i)),round(R(i)));
        CorrMat(counter,13) = gdisttmp(ceil(R(j)),ceil(C(j)))/Ai_pix_mm; %Populate vessel distance from PA to all other PAs
        if isnan(CorrMat(counter,13))
            CorrMat(counter,13) = gdisttmp(ceil(R(j)),floor(C(j)))/Ai_pix_mm; %Populate vessel distance from PA to all other PAs
            if isnan(CorrMat(counter,13))
                CorrMat(counter,13) = gdisttmp(floor(R(j)),floor(C(j)))/Ai_pix_mm; %Populate vessel distance from PA to all other PAs
                if isnan(CorrMat(counter,13))
                    CorrMat(counter,13) = gdisttmp(floor(R(j)),ceil(C(j)))/Ai_pix_mm; %Populate vessel distance from PA to all other PAs
                end
            end
        end
     counter = counter + 1;   
    end
   i 
end

Corr_KF_struct.CorrMat = CorrMat; %Replace CorrMat in CorrKFStruct
save('20221208WT5_Corr_KF_struct.mat','Corr_KF_struct');


%% Combine Corr Mats and plot correlation vs distance
CM1 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20221212 WT_PA_7\20221212WT7_Corr_KF_struct.mat"); CM1 = CM1.Corr_KF_struct.CorrMat;
CM2 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20221129 WT_PA_6\20221129WT6_Corr_KF_struct.mat"); CM2 = CM2.Corr_KF_struct.CorrMat;
% % CM3 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20221120 WT_PA_5\20221120WT5_Corr_KF_struct.mat"); % CM3 = CM3.Corr_KF_struct.CorrMat; %Frame rate 4.2 Hz
CM4 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20221128 WT_PA_5\20221128WT5_Corr_KF_struct.mat"); CM4 = CM4.Corr_KF_struct.CorrMat;
CM5 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20221208 WT_PA_5\20221208WT5_Corr_KF_struct.mat"); CM5 = CM5.Corr_KF_struct.CorrMat;
CM6 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20230118 WT_PA_7\20230118WT7_Corr_KF_struct.mat"); CM6 = CM6.Corr_KF_struct.CorrMat;
CM7 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20230125 WT_PA_7\20230125WT7_Corr_KF_struct.mat"); CM7 = CM7.Corr_KF_struct.CorrMat;
CM8 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20230210 WT_PA_8\20230210WT8_Corr_KF_struct.mat"); CM8 = CM8.Corr_KF_struct.CorrMat;
CM9 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20230213 WT_PA_8\20230213WT8_Corr_KF_struct.mat"); CM9 = CM9.Corr_KF_struct.CorrMat;
CM10 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20230214 WT_PA_9\20230214WT9_Corr_KF_struct.mat"); CM10 = CM10.Corr_KF_struct.CorrMat;
CM11 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20230215 WT_PA_9\20230215WT9_Corr_KF_struct.mat"); CM11 = CM11.Corr_KF_struct.CorrMat;
CM12 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20230222 WT_PA_10\20230222WT10_Corr_KF_struct.mat"); CM12 = CM12.Corr_KF_struct.CorrMat;
CM13 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20230224 WT_PA_10\20230224WT10_Corr_KF_struct.mat"); CM13 = CM13.Corr_KF_struct.CorrMat;
CM14 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20230301 WT_PA_11\20230301WT11_Corr_KF_struct.mat"); CM14 = CM14.Corr_KF_struct.CorrMat;
CM15 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20230302 WT_PA_11\20230302WT11_Corr_KF_struct.mat"); CM15 = CM15.Corr_KF_struct.CorrMat;
CM4(:,13) = []; CM5(:,13) = [];
CombinedCorrMat = [CM1;CM2;CM4;CM5;CM6;CM7;CM8;CM9;CM10;CM11;CM12;CM13;CM14;CM15];

KM1 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20221212 WT_PA_7\20221212WT7_Corr_KF_struct.mat"); KM1 = KM1.Corr_KF_struct.KFmat;
KM2 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20221129 WT_PA_6\20221129WT6_Corr_KF_struct.mat"); KM2 = KM2.Corr_KF_struct.KFmat;
% KM3 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20221120 WT_PA_5\20221120WT5_Corr_KF_struct.mat"); KM3 = KM3.Corr_KF_struct.KFmat;
KM4 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20221128 WT_PA_5\20221128WT5_Corr_KF_struct.mat"); KM4 = KM4.Corr_KF_struct.KFmat;
KM5 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20221208 WT_PA_5\20221208WT5_Corr_KF_struct.mat"); KM5 = KM5.Corr_KF_struct.KFmat;
KM6 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20230118 WT_PA_7\20230118WT7_Corr_KF_struct.mat"); KM6 = KM6.Corr_KF_struct.KFmat;
KM7 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20230125 WT_PA_7\20230125WT7_Corr_KF_struct.mat"); KM7 = KM7.Corr_KF_struct.KFmat;
KM8 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20230210 WT_PA_8\20230210WT8_Corr_KF_struct.mat"); KM8 = KM8.Corr_KF_struct.KFmat;
KM9 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20230213 WT_PA_8\20230213WT8_Corr_KF_struct.mat"); KM9 = KM9.Corr_KF_struct.KFmat;
KM10 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20230214 WT_PA_9\20230214WT9_Corr_KF_struct.mat"); KM10 = KM10.Corr_KF_struct.KFmat;
KM11 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20230215 WT_PA_9\20230215WT9_Corr_KF_struct.mat"); KM11 = KM11.Corr_KF_struct.KFmat;
KM12 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20230222 WT_PA_10\20230222WT10_Corr_KF_struct.mat"); KM12 = KM12.Corr_KF_struct.KFmat;
KM13 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20230224 WT_PA_10\20230224WT10_Corr_KF_struct.mat"); KM13 = KM13.Corr_KF_struct.KFmat;
KM14 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20230301 WT_PA_11\20230301WT11_Corr_KF_struct.mat"); KM14 = KM14.Corr_KF_struct.KFmat;
KM15 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20230302 WT_PA_11\20230302WT11_Corr_KF_struct.mat"); KM15 = KM15.Corr_KF_struct.KFmat;
CombinedKFmat = [KM1;KM2;KM4;KM5;KM6;KM7;KM8;KM9;KM10;KM11;KM12;KM13;KM14;KM15];
KFtodel = ~logical(CombinedKFmat(:,6));
CombinedKFmat(KFtodel,:) = [];
%% Stacked bar chart
%Bar graph of ups and down, stacked by animal 
X = categorical({'Down (+)','Up (-)'});
PM3 = [KM1;KM6;KM7]; PM3(PM3(:,6)~=1,:) = []; %WT7
PM2 = KM2; PM2(PM2(:,6)~=1,:) = []; %WT6
PM1 = [KM4;KM5]; PM1(PM1(:,6)~=1,:) = []; %WT5
PM4 = [KM8;KM9]; PM4(PM4(:,6)~=1,:) = []; %WT8
PM5 = [KM10;KM11]; PM5(PM5(:,6)~=1,:) = []; %WT9
PM6 = [KM12;KM13]; PM6(PM6(:,6)~=1,:) = []; %WT10
PM7 = [KM14;KM15]; PM7(PM7(:,6)~=1,:) = []; %WT11
size([PM1;PM2;PM3;PM4;PM5;PM6;PM7],1)

%Make Y with dimensions [animal1#down,animal2#down,... ;
%animal1#up,animal2#up,...]
Y = zeros(2,7);
Y(1,1) = sum(PM1(:,2)>0)
Y(1,2) = sum(PM2(:,2)>0)
Y(1,3) = sum(PM3(:,2)>0)
Y(1,4) = sum(PM4(:,2)>0)
Y(1,5) = sum(PM5(:,2)>0)
Y(1,6) = sum(PM6(:,2)>0)
Y(1,7) = sum(PM7(:,2)>0)
Y(2,1) = sum(PM1(:,2)<0)
Y(2,2) = sum(PM2(:,2)<0)
Y(2,3) = sum(PM3(:,2)<0)
Y(2,4) = sum(PM4(:,2)<0)
Y(2,5) = sum(PM5(:,2)<0)
Y(2,6) = sum(PM6(:,2)<0)
Y(2,7) = sum(PM7(:,2)<0)

figure
bar(X,Y,'stacked'); ylabel('Count','Interpreter','latex');
title('Penetrating arteriole traveling wave direction, 7 animals 138 arterioles','Interpreter','latex');
lgd = legend('WT5','WT6','WT7','WT8','WT9','WT10','WT11');
title(lgd,'Animal')
savefig('DirectionBarStacked.fig')

x = 0:138;
y = binopdf(x,138,0.5);
figure
bar(x,y,1)
xlabel('Number of successes','Interpreter','latex')
ylabel('Probability','Interpreter','latex')
title('Binomial distribution for N = 138 coin flips, P(down) = 0.5','Interpreter','latex')
savefig('BinomialDist_138Ves_P50.fig')
%% Plot top power vs bottom power, color by wave direction
%CorrMat(10, 11) is power of shallow, deep segment at fv
%Kf mat 4, 5 shallow deep
ispos = CombinedKFmat(:,2) > 0;
posplot = CombinedKFmat(ispos,:);
negplot = CombinedKFmat(~ispos,:);

figure
scatter(negplot(:,5),negplot(:,4),'filled','b')
hold on
scatter(posplot(:,5),posplot(:,4),'filled','r')

figure
scatter(negplot(:,4).^(1/2),negplot(:,5).^(1/2),'filled','b','MarkerFaceAlpha',0.4)
hold on
scatter(posplot(:,4).^(1/2),posplot(:,5).^(1/2),'filled','r','MarkerFaceAlpha',0.4)
xlabel('sqrt(Shallow power at fv)','Interpreter','latex');
ylabel('sqrt(Deep power at fv)','Interpreter','latex');
x = 0:0.1:2.5;
y = x;
hold on
plot(x,y,'k')
legend({'Outward','Inward'})
savefig('DeepvsShallow_Scatter.fig')

figure
histogram(negplot(:,5).^(1/2),'FaceColor','b','FaceAlpha',0.4,'BinWidth',0.2); hold on;
histogram(posplot(:,5).^(1/2),'FaceColor','r','FaceAlpha',0.4,'BinWidth',0.2);
xlabel('sqrt(Deep Segment Power','Interpreter','latex');
ylabel('Count','Interpreter','latex');
legend({'Outward','Inward'})

figure
histogram(negplot(:,4).^(1/2),'FaceColor','b','FaceAlpha',0.4,'BinWidth',0.2); hold on;
histogram(posplot(:,4).^(1/2),'FaceColor','r','FaceAlpha',0.4,'BinWidth',0.2);
xlabel('sqrt(Shallow Segment Power)','Interpreter','latex');
ylabel('Count','Interpreter','latex');
legend({'Outward','Inward'})
%%
%Correlation function only for significant phase differences
%Plot scatter with percent chance of observed directions
CohTodel = ~logical(CombinedCorrMat(:,12)); %Delete combinations where either vessel have not significant coherence
CombinedCorrMat(CohTodel,:) = [];

f1 = figure;
histogram(CombinedCorrMat(:,5),'BinWidth',0.1);
xlabel('Distance (mm)','Interpreter','latex');
ylabel('Count','Interpreter','latex');
title('Pair Distance Distribution','Interpreter','latex');

numSDs = 0;
binsz = 0.2; %mm
expecvec = zeros(9,3);
% CombCorrMat_tmp = CorrMat; %or Combined Corr Mat
CombCorrMat_tmp = CombinedCorrMat;
todel1 = abs(CombCorrMat_tmp(:,6)) < numSDs * CombCorrMat_tmp(:,7);
todel2 = abs(CombCorrMat_tmp(:,8)) < numSDs * CombCorrMat_tmp(:,9);
todel = or(todel1,todel2);
CombCorrMat_tmp(todel,:) = [];
for i=1:9
    inbin1 = CombCorrMat_tmp(:,5) < i*binsz;
    inbin2 = CombCorrMat_tmp(:,5) > (i-1)*binsz;
    inbin = and(inbin1,inbin2); %indices for pairs within distance bin
    
    spinvec = sign(CombCorrMat_tmp(:,3).*CombCorrMat_tmp(:,4));
    spinvec = spinvec(inbin);

    findmax = [sum(spinvec<0),sum(spinvec>0)];
    [foundmax,maxind] = max(findmax);
    if maxind == 1
        probsuccess = 0.4731; %No filtering
    elseif maxind == 2
        probsuccess = 0.5269; %No filtering
    end

%     prob(i) = 1 - binocdf(foundmax,numel(spinvec),probsuccess); %1-probability to see more than (foundmax) successes
    prob(i) = binopdf(foundmax,numel(spinvec),probsuccess); %Probability density at observed value only
    possiblevals = 1:numel(spinvec);
    pdfvals = binopdf(possiblevals,numel(spinvec),probsuccess);
    [~,pdfmaxind] = max(pdfvals);
    bpdf(i) = possiblevals(pdfmaxind); %Value with highest binopdf
    if maxind == 1
        mostlikelyavg(i) = ((numel(spinvec)-bpdf(i))-bpdf(i))/numel(spinvec);
        biolims = binoinv([0.025 0.975],numel(spinvec),probsuccess);
%         biolims = binoinv([0.05 0.95],numel(spinvec),probsuccess);
        limsavg(i,:) = ((numel(spinvec)-biolims)-biolims)./numel(spinvec);
    elseif maxind == 2
        mostlikelyavg(i) = (-(numel(spinvec)-bpdf(i))+bpdf(i))/numel(spinvec);
        biolims = binoinv([0.025 0.975],numel(spinvec),probsuccess);
%         biolims = binoinv([0.05 0.95],numel(spinvec),probsuccess);
        limsavg(i,:) = (-(numel(spinvec)-biolims)+biolims)./numel(spinvec);
    end

    meanspinvec = mean(spinvec);
    stdspinvec = std(spinvec);

    expecvec(i,1) = meanspinvec;
    expecvec(i,2) = stdspinvec;
    expecvec(i,3) = stdspinvec/sqrt(numel(spinvec));
    
end
%Make scatter from expecvec
x_dist = 0:binsz:binsz*(size(expecvec,1)-1);
x_dist = x_dist + binsz/2
y_expec = expecvec(:,1); %<sig1 * sig2>
% err = expecvec(:,3); %SDM
f2 = figure;
subplot(2,1,1);
scatter(x_dist,y_expec,'filled','blue');
hold on
s1 = scatter(x_dist,mostlikelyavg,'_r')
s1.LineWidth = 1.5;
s2 = scatter(x_dist,limsavg,'_r')
% s2.SizeData = 0.75
xlim([0 1.8]); ylim([-1 1]); yline(0);
xlabel('Distance (binned at 0.2 mm)','Interpreter','latex','FontSize',12);
ylabel('$\biggl \langle \sigma_i \sigma_j \biggr \rangle$','Interpreter','latex','FontSize',12)
grid on

title('Correlation vs euclidean distance, 7 animals, 683 pairs','Interpreter','latex','FontSize',15);

numinbin = (expecvec(:,2)./expecvec(:,3)).^2;
f3 = figure;
subplot(2,1,1);
bar(x_dist,numinbin');
xlabel('Euclidean distance bin center (mm)','Interpreter','latex');
ylabel('Vessels in bin','Interpreter','latex');

%%
%Same for vessel distance
CombinedCorrMat_VD = CombinedCorrMat;
InfTodel = isinf(CombinedCorrMat_VD(:,13)); %Delete combinations where either vessel have not significant coherence
CombinedCorrMat_VD(InfTodel,:) = [];

figure(f1); hold on;
histogram(CombinedCorrMat_VD(:,13),'BinWidth',0.1);
xlabel('Distance (mm)','Interpreter','latex');
ylabel('Count','Interpreter','latex');
title('Pair Distance Distribution','Interpreter','latex');
legend({'Euclidean','Vessel Distance'})

numSDs = 0;
binsz = 0.4; %mm
expecvec = zeros(12,3);
CombCorrMat_tmp = CombinedCorrMat_VD; %or Combined Corr Mat
for i=1:12
    inbin1 = CombCorrMat_tmp(:,13) < i*binsz; %Vessel distance is in CorrMat(i,13)
    inbin2 = CombCorrMat_tmp(:,13) > (i-1)*binsz;
    inbin = and(inbin1,inbin2); %indices for pairs within distance bin
    
    spinvec = sign(CombCorrMat_tmp(:,3).*CombCorrMat_tmp(:,4));
    spinvec = spinvec(inbin);

    findmax = [sum(spinvec<0),sum(spinvec>0)];
    [foundmax,maxind] = max(findmax);
    if maxind == 1
        probsuccess = 0.4731; %No filtering
    elseif maxind == 2
        probsuccess = 0.5269; %No filtering
    end

%     prob(i) = 1 - binocdf(foundmax,numel(spinvec),probsuccess); %1-probability to see more than (foundmax) successes
    prob(i) = binopdf(foundmax,numel(spinvec),probsuccess); %Probability density at observed value only
    possiblevals = 1:numel(spinvec);
    pdfvals = binopdf(possiblevals,numel(spinvec),probsuccess);
    [~,pdfmaxind] = max(pdfvals);
    bpdf(i) = possiblevals(pdfmaxind); %Value with highest binopdf
    if maxind == 1
        mostlikelyavg(i) = ((numel(spinvec)-bpdf(i))-bpdf(i))/numel(spinvec);
        biolims = binoinv([0.025 0.975],numel(spinvec),probsuccess);
%         biolims = binoinv([0.05 0.95],numel(spinvec),probsuccess);
        limsavg(i,:) = ((numel(spinvec)-biolims)-biolims)./numel(spinvec);
    elseif maxind == 2
        mostlikelyavg(i) = (-(numel(spinvec)-bpdf(i))+bpdf(i))/numel(spinvec);
        biolims = binoinv([0.025 0.975],numel(spinvec),probsuccess);
%         biolims = binoinv([0.05 0.95],numel(spinvec),probsuccess);
        limsavg(i,:) = (-(numel(spinvec)-biolims)+biolims)./numel(spinvec);
    end

    meanspinvec = mean(spinvec);
    stdspinvec = std(spinvec);

    expecvec(i,1) = meanspinvec;
    expecvec(i,2) = stdspinvec;
    expecvec(i,3) = stdspinvec/sqrt(numel(spinvec));
    
end
figure(f2); subplot(2,1,2)
x_dist = 0:binsz:binsz*(size(expecvec,1)-1);
x_dist = x_dist + binsz/2
y_expec = expecvec(:,1); %<sig1 * sig2>
scatter(x_dist,y_expec,'filled','k');
hold on
s1 = scatter(x_dist,mostlikelyavg,'_r')
s1.LineWidth = 1.5;
s2 = scatter(x_dist,limsavg,'_r')
% s2.SizeData = 0.75
xlim([0 4.5]); ylim([-1 1]); yline(0);
xlabel('Vessel Distance (binned at 0.5 mm)','Interpreter','latex','FontSize',12);
ylabel('$\biggl \langle \sigma_i \sigma_j \biggr \rangle$','Interpreter','latex','FontSize',12)
grid on
title('Correlation vs vessel distance, 7 animals, 644 pairs','Interpreter','latex','FontSize',15);

numinbin = (expecvec(:,2)./expecvec(:,3)).^2;
figure(f3);
subplot(2,1,2);
bar(x_dist,numinbin');
xlabel('Vessel distance bin center (mm)','Interpreter','latex');
ylabel('Vessels in bin','Interpreter','latex');

%% Correlation for euclidean and vessel distance - equal percent in each bin
CohTodel = ~logical(CombinedCorrMat(:,12)); %Delete combinations where either vessel have not significant coherence
CombinedCorrMat(CohTodel,:) = [];

f1 = figure;
histogram(CombinedCorrMat(:,5),'BinWidth',0.1);
xlabel('Distance (mm)','Interpreter','latex');
ylabel('Count','Interpreter','latex');
title('Pair Distance Distribution','Interpreter','latex');
[f,d] = ecdf(CombinedCorrMat(:,5));
figure; plot(d,f);

numSDs = 0;
numbins = 10; %10 percent of pairs in each bin
tol = (f(2)-f(1))/2;
binvec = zeros(numbins,1);
for i=1:numbins
    binvectmp = find(abs(f-(i/numbins)) < tol);
    [~,minbinvec] = min(abs(f(binvectmp)-(i/numbins)));
    binvec(i) = d(binvectmp(minbinvec));
end
binvec = [0;binvec];
expecvec = zeros(length(binvec),3);
% CombCorrMat_tmp = CorrMat; %or Combined Corr Mat
CombCorrMat_tmp = CombinedCorrMat;
todel1 = abs(CombCorrMat_tmp(:,6)) < numSDs * CombCorrMat_tmp(:,7);
todel2 = abs(CombCorrMat_tmp(:,8)) < numSDs * CombCorrMat_tmp(:,9);
todel = or(todel1,todel2);
CombCorrMat_tmp(todel,:) = [];
for i=2:size(binvec,1)
    inbin1 = CombCorrMat_tmp(:,5) < binvec(i);
    inbin2 = CombCorrMat_tmp(:,5) > binvec(i-1);
    inbin = and(inbin1,inbin2); %indices for pairs within distance bin
    x_dist(i) = median(CombCorrMat_tmp(inbin,5)); %Distance to plot is median of binned distances
    
    spinvec = sign(CombCorrMat_tmp(:,3).*CombCorrMat_tmp(:,4));
    spinvec = spinvec(inbin);

    findmax = [sum(spinvec<0),sum(spinvec>0)];
    [foundmax,maxind] = max(findmax);
    if maxind == 1
        probsuccess = 0.4731; %No filtering
    elseif maxind == 2
        probsuccess = 0.5269; %No filtering
    end

%     prob(i) = 1 - binocdf(foundmax,numel(spinvec),probsuccess); %1-probability to see more than (foundmax) successes
    prob(i) = binopdf(foundmax,numel(spinvec),probsuccess); %Probability density at observed value only
    possiblevals = 1:numel(spinvec);
    pdfvals = binopdf(possiblevals,numel(spinvec),probsuccess);
    [~,pdfmaxind] = max(pdfvals);
    bpdf(i) = possiblevals(pdfmaxind); %Value with highest binopdf
    if maxind == 1
        mostlikelyavg(i) = ((numel(spinvec)-bpdf(i))-bpdf(i))/numel(spinvec);
        biolims = binoinv([0.025 0.975],numel(spinvec),probsuccess);
%         biolims = binoinv([0.05 0.95],numel(spinvec),probsuccess);
        limsavg(i,:) = ((numel(spinvec)-biolims)-biolims)./numel(spinvec);
    elseif maxind == 2
        mostlikelyavg(i) = (-(numel(spinvec)-bpdf(i))+bpdf(i))/numel(spinvec);
        biolims = binoinv([0.025 0.975],numel(spinvec),probsuccess);
%         biolims = binoinv([0.05 0.95],numel(spinvec),probsuccess);
        limsavg(i,:) = (-(numel(spinvec)-biolims)+biolims)./numel(spinvec);
    end

    meanspinvec = mean(spinvec);
    stdspinvec = std(spinvec);

    expecvec(i,1) = meanspinvec;
    expecvec(i,2) = stdspinvec;
    expecvec(i,3) = stdspinvec/sqrt(numel(spinvec));
    numinbin(i) = numel(spinvec);
end
expecvec(1,:) = []; %First entry has nothing
prob(1) = [];
limsavg(1,:) = [];
mostlikelyavg(1) = [];
x_dist(1) = [];
numinbin(1) = [];
%Make scatter from expecvec
% x_dist = (binvec(1:end-1) + binvec(2:end))/2; %Get mean distance of every bin for plotting
%Set significance level as constant for constant number in each bin
minnum = min(numinbin);
probsuccess = 0.4731;
biolimsconst = binoinv([0.025 0.975],minnum,probsuccess);
limsavgconst = ((minnum-biolimsconst)-biolimsconst)./minnum;

y_expec = expecvec(:,1); %<sig1 * sig2>
% err = expecvec(:,3); %SDM
f2 = figure;
% subplot(2,1,1);
scatter(x_dist,y_expec,'filled','blue');
hold on
% s1 = scatter(x_dist,mostlikelyavg,'_r')
% s1.LineWidth = 1.5;
% s2 = scatter(x_dist,limsavg,'_r')
% s2 = scatter(x_dist,limsavg,'_r')
yline(limsavgconst(1),'r'); yline(limsavgconst(2),'r');
% s2.SizeData = 0.75
xlim([0 1.8]); ylim([-1 1]); yline(0);
xlabel('Distance (10Percent of all pairs in each bin)','Interpreter','latex','FontSize',12);
ylabel('$\biggl \langle \sigma_i \sigma_j \biggr \rangle$','Interpreter','latex','FontSize',12)
grid on
title('Correlation vs euclidean distance, 7 animals, 683 pairs','Interpreter','latex','FontSize',15);
%% Same for vessel distance part 2
CombinedCorrMat_VD = CombinedCorrMat;
InfTodel = isinf(CombinedCorrMat_VD(:,13)); %Delete combinations where either vessel have not significant coherence
CombinedCorrMat_VD(InfTodel,:) = [];
figure(f1); hold on;
histogram(CombinedCorrMat_VD(:,13),'BinWidth',0.1);
xlabel('Distance (mm)','Interpreter','latex');
ylabel('Count','Interpreter','latex');
title('Pair Distance Distribution','Interpreter','latex');
[f,d] = ecdf(CombinedCorrMat_VD(:,13));
figure; plot(d,f);

numSDs = 0;
numbins = 10; %10 percent of pairs in each bin
tol = (f(2)-f(1))/2;
binvec = zeros(numbins,1);
for i=1:numbins
    binvectmp = find(abs(f-(i/numbins)) < tol);
    binvec(i) = d(binvectmp);
end
binvec = [0;binvec];
expecvec = zeros(length(binvec),3);
% CombCorrMat_tmp = CorrMat; %or Combined Corr Mat
CombCorrMat_tmp = CombinedCorrMat;
todel1 = abs(CombCorrMat_tmp(:,6)) < numSDs * CombCorrMat_tmp(:,7);
todel2 = abs(CombCorrMat_tmp(:,8)) < numSDs * CombCorrMat_tmp(:,9);
todel = or(todel1,todel2);
CombCorrMat_tmp(todel,:) = [];
for i=2:size(binvec,1)
    inbin1 = CombCorrMat_tmp(:,13) < binvec(i);
    inbin2 = CombCorrMat_tmp(:,13) > binvec(i-1);
    inbin = and(inbin1,inbin2); %indices for pairs within distance bin
    x_dist(i) = median(CombCorrMat_tmp(inbin,13)); %Distance to plot is median of binned distances
    
    spinvec = sign(CombCorrMat_tmp(:,3).*CombCorrMat_tmp(:,4));
    spinvec = spinvec(inbin);

    findmax = [sum(spinvec<0),sum(spinvec>0)];
    [foundmax,maxind] = max(findmax);
    if maxind == 1
        probsuccess = 0.4731; %No filtering
    elseif maxind == 2
        probsuccess = 0.5269; %No filtering
    end

%     prob(i) = 1 - binocdf(foundmax,numel(spinvec),probsuccess); %1-probability to see more than (foundmax) successes
    prob(i) = binopdf(foundmax,numel(spinvec),probsuccess); %Probability density at observed value only
    possiblevals = 1:numel(spinvec);
    pdfvals = binopdf(possiblevals,numel(spinvec),probsuccess);
    [~,pdfmaxind] = max(pdfvals);
    bpdf(i) = possiblevals(pdfmaxind); %Value with highest binopdf
    if maxind == 1
        mostlikelyavg(i) = ((numel(spinvec)-bpdf(i))-bpdf(i))/numel(spinvec);
        biolims = binoinv([0.025 0.975],numel(spinvec),probsuccess);
%         biolims = binoinv([0.05 0.95],numel(spinvec),probsuccess);
        limsavg(i,:) = ((numel(spinvec)-biolims)-biolims)./numel(spinvec);
    elseif maxind == 2
        mostlikelyavg(i) = (-(numel(spinvec)-bpdf(i))+bpdf(i))/numel(spinvec);
        biolims = binoinv([0.025 0.975],numel(spinvec),probsuccess);
%         biolims = binoinv([0.05 0.95],numel(spinvec),probsuccess);
        limsavg(i,:) = (-(numel(spinvec)-biolims)+biolims)./numel(spinvec);
    end

    meanspinvec = mean(spinvec);
    stdspinvec = std(spinvec);

    expecvec(i,1) = meanspinvec;
    expecvec(i,2) = stdspinvec;
    expecvec(i,3) = stdspinvec/sqrt(numel(spinvec));
    
end
expecvec(1,:) = []; %First entry has nothing
prob(1) = [];
limsavg(1,:) = [];
mostlikelyavg(1) = [];
x_dist(1) = [];
%Make scatter from expecvec
% x_dist = (binvec(1:end-1) + binvec(2:end))/2; %Get mean distance of every bin for plotting
y_expec = expecvec(:,1); %<sig1 * sig2>
% err = expecvec(:,3); %SDM
figure(f2); subplot(2,1,2)
scatter(x_dist,y_expec,'filled','k');
hold on
s1 = scatter(x_dist,mostlikelyavg,'_r')
s1.LineWidth = 1.5;
s2 = scatter(x_dist,limsavg,'_r')
% s2.SizeData = 0.75
xlim([0 5]); ylim([-1 1]); yline(0);
xlabel('Distance (10Percent of all pairs in each bin)','Interpreter','latex','FontSize',12);
ylabel('$\biggl \langle \sigma_i \sigma_j \biggr \rangle$','Interpreter','latex','FontSize',12)
grid on
title('Correlation vs vessel distance, 7 animals, 644 pairs','Interpreter','latex','FontSize',15);

%% Calculate both correlations for each animal - across 2 experiment days.
%9/11 vessels measured on different days kept the same delay direction






%% Plotting
cd('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Plots_3_17_23')
%Plot 1/(dk)^2 vs. |k|
figure
scatter(abs(CombinedKFmat(:,2)),1./(CombinedKFmat(:,3).^2),'filled','MarkerFaceAlpha',1)  
ylim([0 300]);
xlabel('$|k|$','Interpreter','latex','FontSize',15);
ylabel('1/($\delta$k$)^2$','Interpreter','latex','FontSize',15);
ax = gca;
ax.TickLabelInterpreter = 'latex';
% savefig('dk2vsk_scatter.fig')
print(gcf, '-depsc2', 'dk2vsk_scatter');

% Calculate moving median for penetrating arterioles f vs k



figure
% scatter(CombinedKFmat(:,2),CombinedKFmat(:,3)./abs(CombinedKFmat(:,2)),'filled','MarkerFaceAlpha',0.4)  
scatter(CombinedKFmat(:,2),CombinedKFmat(:,3),'filled','MarkerFaceAlpha',1)  
ylim([0 2])
xline(0)
xlabel('k (rad/mm)','Interpreter','latex','FontSize',15);
% ylabel('$\delta$k/$|k|$','Interpreter','latex','FontSize',15);
ylabel('$\delta$k (rad/mm)','Interpreter','latex','FontSize',15);
title({'138 Penetrating vessel phase gradients','(+) shallow precedes deep (-) deep precedes shallow'},'Interpreter','latex');

figure
histogram(CombinedKFmat(:,2),'BinWidth',0.2,'FaceAlpha',1)
xlabel('k (rad/mm)','Interpreter','latex','FontSize',15);
ylabel('Count','Interpreter','latex','FontSize',15);    

figure
scatter(CombinedKFmat(:,3)./abs(CombinedKFmat(:,2)),CombinedKFmat(:,2),'filled','MarkerFaceAlpha',0.4)  
xlim([0 10])
yline(0)
ylabel('k (rad/mm)','Interpreter','latex','FontSize',15);
xlabel('$\delta$k/$|k|$','Interpreter','latex','FontSize',15);
title({'138 Penetrating vessel phase gradients','(+) shallow precedes deep (-) deep precedes shallow'},'Interpreter','latex');

%F vs K
ispos = CombinedKFmat(:,2) > 0;
posplot = CombinedKFmat(ispos,:);
negplot = CombinedKFmat(~ispos,:);
figure
scatter(abs(posplot(:,2)),posplot(:,1),'filled','MarkerFaceAlpha',0.1,'MarkerFaceColor','r','MarkerEdgeColor','none'); hold on;
scatter(abs(negplot(:,2)),negplot(:,1),'filled','MarkerFaceAlpha',0.1,'MarkerFaceColor','b','MarkerEdgeColor','none');
xlim([0 4]); ylim([0 0.2]);
xlabel('Magnitude K (rad/mm)','Interpreter','latex');
ylabel('Vessel peak frequency (Hz)','Interpreter','latex');
title({'Penetrating arteriole f vs k, 7 animals, 14 trials, 138 penetrating arterioles','Blue = deep precedes shallow Red = shallow precedes deep'},'Interpreter','latex','FontSize',13)
ax = gca;
ax.TickLabelInterpreter = 'latex';

figure
histogram(abs(negplot(:,2)),'BinWidth',0.25,'FaceColor','b');
hold on;
histogram(abs(posplot(:,2)),'BinWidth',0.25,'FaceColor','r');
xlim([0 4]);
xlabel('Magnitude K (rad/mm)','Interpreter','latex');
ylabel('Count','Interpreter','latex');

%Plot Power ratio versus k
kvec = CombinedKFmat(:,2);
fvec = CombinedKFmat(:,1);
PowerMetric = (CombinedKFmat(:,5)./CombinedKFmat(:,4)).^(1/2);
figure
scatter(kvec,PowerMetric,'filled','MarkerFaceAlpha',0.5);
xline(0); grid on;
xlabel('k (rad/mm)','Interpreter','latex');
ylabel('$\sqrt{Power_d/Power_s}$','Interpreter','latex');
title('Diameter power ratio vs k, Positive k = shallow preceeds deep','Interpreter','latex');
savefig('PowerRatio_vs_k.fig')

% histogram of powrer ratio, colored by k sign
ispos = sign(kvec) > 0;
isneg = ~ispos;
figure
histogram(PowerMetric(isneg),'BinWidth',0.1); 
hold on
histogram(PowerMetric(ispos),'BinWidth',0.1);
xlabel('$\sqrt{Power_d/Power_s}$','Interpreter','latex');
ylabel('Count','Interpreter','latex');
legend({'Outward','Inward'})
 
% Calculate moving median for penetrating arterioles f vs k
% Fit f vs k to moving median
kvec = abs(CombinedKFmat(:,2));
fvec = CombinedKFmat(:,1);
[f,d] = ecdf(fvec);
numbins = 10; %10 percent of pairs in each bin
% tol = (f(2)-f(1))/2;
% tol = (1/numel(f))/2;
tol = (median(diff(f)))*4;
% tol = (median(diff(f)))*6;
binvec = zeros(numbins,1);
for i=1:numbins
    binvectmp = find(abs(f-(i/numbins)) < tol);
    [isclosest,isclosestind] = min(abs(f(binvectmp)-i/numbins));
    binvec(i) = d(binvectmp(isclosestind));
end
binvec = [0;binvec];
figure; plot(d,f);
fbins = binvec';
meank = zeros(size(fbins,2)-1,1);
medk = zeros(size(fbins,2)-1,1);
sek = zeros(size(fbins,2)-1,1);
fplot = zeros(size(fbins,2)-1,1);
num_bin = zeros(size(fbins,2)-1,1);
for i=2:size(fbins,2)
    inbin1 = fvec > fbins(i-1);
    if i==size(fbins,2)
        inbin2 = fvec <= max(fvec);
    else
        inbin2 = fvec < fbins(i);
    end
    inbin = and(inbin1,inbin2);
    fvectmp = fvec(inbin);
    kvectmp = kvec(inbin);
    meank(i-1) = mean(kvectmp);
    medk(i-1) = median(kvectmp);
    sek(i-1) = std(kvectmp)/sqrt(numel(kvectmp));
    fplot(i-1) = median(fvectmp);
    num_bin(i-1) = numel(kvectmp);
end

f1 = figure;
alpha = 0.01;
scatter(kvec,fvec,'filled','MarkerFaceColor','k','MarkerEdgeColor','none','MarkerFaceAlpha',0.2);
xlim([0 3.5]); ylim([0 0.18]);
xlabel('Phase gradient magnitude (rad/mm)','Interpreter','latex');
ylabel('Vessel peak vasomotor frequency (Hz)','Interpreter','latex');
% str = sprintf('Pial arteriole f vs k, 24 animals, %.0f vessels',length(pvcomb));
% str1 = sprintf('%.2f mm minimum vessel length, T-test',minlength);
% str2 = sprintf('%.2f',alpha);
% title({str,[str1,' $\alpha$ = ',str2],'Median k values, 10pct f bins'},'Interpreter','latex')
ax = gca;
ax.TickLabelInterpreter = 'latex';
% f1.Position = [4608 45 759 665];
hold on;
% errorbar(medk,fplot,sek,'horizontal','Color','r','LineStyle','none','LineWidth',1)
scatter(medk,fplot,'filled','r');
x = 0:0.01:3.5;
y = 0.2324*x; %Fit to weighted median k values
hold on
plot(x,y,'r')

%%
% Fit f vs k 10% f bins, weighted mean k in each bin, weighted median k
kvec = abs(CombinedKFmat(:,2));
fvec = CombinedKFmat(:,1);
wts = 1./(CombinedKFmat(:,3).^2);
[f,d] = ecdf(fvec);
numbins = 10; %10 percent of pairs in each bin
% tol = (f(2)-f(1))/2;
% tol = (1/numel(f))/2;
tol = (median(diff(f)))*4;
% tol = (median(diff(f)))*6;
binvec = zeros(numbins,1);
for i=1:numbins
    binvectmp = find(abs(f-(i/numbins)) < tol);
    [isclosest,isclosestind] = min(abs(f(binvectmp)-i/numbins));
    binvec(i) = d(binvectmp(isclosestind));
end
binvec = [0;binvec];
figure; plot(d,f);
fbins = binvec';
meank = zeros(size(fbins,2)-1,1);
wtdmedk = zeros(size(fbins,2)-1,1);
sek = zeros(size(fbins,2)-1,1);
fplot = zeros(size(fbins,2)-1,1);
num_bin = zeros(size(fbins,2)-1,1);
for i=2:size(fbins,2)
    inbin1 = fvec > fbins(i-1);
    if i==size(fbins,2)
        inbin2 = fvec <= max(fvec);
    else
        inbin2 = fvec < fbins(i);
    end
    inbin = and(inbin1,inbin2);
    fvectmp = fvec(inbin);
    kvectmp = kvec(inbin);
    wtstmp = wts(inbin);

    meannum = dot(wtstmp,kvectmp);
    meandenom = sum(wtstmp);
    meank(i-1) = meannum/meandenom; %Weighted mean k

    wtstot = sum(wtstmp);
    [kvectmp,sortinds] = sort(kvectmp);
    wtstmp = wtstmp(sortinds);
    fvectmp = fvectmp(sortinds);
    wtscumsum = cumsum(wtstmp);
    issmallhalf = find(wtscumsum <= wtstot/2);
    islargehalf = find(wtscumsum >= wtstot/2);
    if isempty(intersect(issmallhalf,islargehalf))
        wtdmedk(i-1) = mean([kvectmp(issmallhalf(end)),kvectmp(islargehalf(1))]);
    else
        [indtmp,~] = intersect(issmallhalf,islargehalf);
        wtdmedk(i-1) = kvectmp(indtmp);
    end

%     medk(i-1) = median(kvectmp);
    sek(i-1) = std(kvectmp)/sqrt(numel(kvectmp));
    fplot(i-1) = median(fvectmp);
    num_bin(i-1) = numel(kvectmp);
end


f1 = figure;
alpha = 0.01;
scatter(kvec,fvec,'filled','MarkerFaceColor','k','MarkerEdgeColor','none','MarkerFaceAlpha',0.2);
xlim([0 3.5]); ylim([0 0.18]);
xlabel('Phase gradient magnitude (rad/mm)','Interpreter','latex');
ylabel('Vessel peak vasomotor frequency (Hz)','Interpreter','latex');
% str = sprintf('Pial arteriole f vs k, 24 animals, %.0f vessels',length(pvcomb));
% str1 = sprintf('%.2f mm minimum vessel length, T-test',minlength);
% str2 = sprintf('%.2f',alpha);
% title({str,[str1,' $\alpha$ = ',str2],'Median k values, 10pct f bins'},'Interpreter','latex')
ax = gca;
ax.TickLabelInterpreter = 'latex';
% f1.Position = [4608 45 759 665];
hold on;
% errorbar(medk,fplot,sek,'horizontal','Color','r','LineStyle','none','LineWidth',1)
scatter(wtdmedk,fplot,'filled','r');
% scatter(meank,fplot,'filled','r');
x = 0:0.01:3.5;
y = 0.3154*x; %Fit to weighted median k values
hold on
plot(x,y,'r')



