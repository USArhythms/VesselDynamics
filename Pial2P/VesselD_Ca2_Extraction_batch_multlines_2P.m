%VesselD_Ca2_Extraction.m
%Code for extraction and comparison of traveling signals in TxRed/GCaMP
%Updated version 5/15/23
%Version to extract diameter on muralis

%% Process to get diameter signal in matlab (5/1/23)
% cd(['/net/birdstore/Thomas_1P/Gal/',animal])
addpath('/home/jaduckwo');
roiname = 'ROI3'
animal = 'JD221223F2'
filepath = ['/net/birdstore/Thomas_1P/JD_GC_Diam2P/6_13_23/8_13_23Analysis/',animal,'/',roiname];
cd(filepath)

tmpdiam = struct(); %For saving temporary extracted diameters
% load('inds.mat');
% load('inds_50lines.mat');
load('seg_inds.mat');
origmask = imread('Mask_origsize.tif');
mask = imread('Mask.tif');
% cd(['/net/birdstore/Thomas_1P/JD_GC_Diam2P/6_13_23/',roiname]) %Load images from 6/13 analysis, nothing has changed here.
cd(filepath);
d_data = h5read([animal,'_',roiname,'_diamdat.h5'],'/frame');
cd(filepath)
d_data = reshape(d_data,[size(origmask,1),size(origmask,2),size(d_data,2)]); %X x Y x Time
disp('Done Loading')

meanimage = mean(d_data,3);

%%
crosslines = 10;
numlinegroups = length(inds)/crosslines;

linesiter = 1:crosslines:numlinegroups*crosslines;
diam = zeros(size(d_data,3),numlinegroups);
pix_mm = 4000;

cellsizes = cellfun(@length,inds,'UniformOutput',true);
maxcellsize = max(cellsizes);
indmat = NaN(maxcellsize,length(inds));
for i = 1:length(inds)
    indmat(1:length(inds{i}),i) = inds{i};
end
indmat_todel = isnan(indmat);
indmat(indmat_todel) = 1; %Temporarily, get the first pixel of the image, then delete later.

pix_mm_lines = zeros(numlinegroups,1);
for i = 1:numlinegroups
    linetmp = indmat(:,linesiter(i));
    linetmp = linetmp(~indmat_todel(:,linesiter(i)));
    [tmprow,tmpcol] = ind2sub(size(mask(:,:,1)),linetmp);
    pix_dist_line = sqrt((tmprow(end)-tmprow(1))^2 + (tmpcol(end)-tmpcol(1))^2); %Units: pixels
    pix_mm_lines(i) = length(linetmp)/((1/pix_mm)*pix_dist_line); %pixels in line per mm
end
pix_mm_lines = pix_mm_lines';

%Analysis parameters
sf = 0.3;
MedFiltVal = 9; %For image
ProfMedFiltVal = 6; %For average profile

pools = 48;
parpool(pools);

tic
parfor j = 1:size(d_data,3) %Iterate over frames
    %     for j = 1:1
    numlines = numlinegroups;
    imtmp = imresize(d_data(:,:,j),[size(mask,1),size(mask,2)],'nearest'); %Try not median filtering first
    imtmp = medfilt2(imtmp,[MedFiltVal,MedFiltVal]);
    imtmp = rescale(imtmp);

    counter = 1;
    proftmp = nan();
    mincellsize = nan(numlinegroups,1);
    for line = linesiter
        cellsizes_tmp = cellsizes(line:line+(crosslines-1));
        mincellsize(counter) = min(cellsizes_tmp);
        proftmp(1:mincellsize(counter),counter) = mean(double(imtmp(indmat(1:mincellsize(counter),line:line+(crosslines-1)))),2); %Average of 10 line profiles
        %             proftmp(1:mincellsize(counter),counter) = median(double(imtmp(indmat(1:mincellsize(counter),line:line+(crosslines-1)))),2); %Try this for trials with noisy estimates 
        counter = counter + 1;
    end
    tonan = proftmp == 0;
    proftmp = medfilt1(proftmp,ProfMedFiltVal,'omitnan');
    proftmp(tonan) = nan;

    proftmp_flip = NaN(size(proftmp));
    for i = 1:numlines
        tmpcol = proftmp(:,i);
        proftmp_flip(1:mincellsize(i),i) = flip(tmpcol(~isnan(tmpcol)));
    end
    proftmp_flip(11:end,:) = []; %Use last 10 entries to estimate minimum
    profstart = proftmp(1:10,:);
%     sf = 0.4; %FW(scalefactor)M
    
    Minvaltest = [std(profstart,[],1);std(proftmp_flip,[],1)];
    [~,Minval_isflatprof] = min(Minvaltest,[],1); %1 if start of profile is less variable, 2 if end of profile is less variable.
    Minval_isprofstart = Minval_isflatprof == 1; %Logical, which lines to use start of profile for.
    Minval = zeros(1,size(proftmp_flip,2));
    Minval(Minval_isprofstart) = median(profstart(:,Minval_isprofstart));
    Minval(~logical(Minval_isprofstart)) = median(proftmp_flip(:,~logical(Minval_isprofstart)));

%     Minval = median(proftmp_flip,1); %Min pixel intensity of each profile. Improve this by choosing which side has more constant ending profile.
    Maxval = max(proftmp,[],1,'omitnan'); %Max pixel intensity of each profile
    HMval = (Maxval - Minval)*sf + Minval;

    intpfactor = 0.01;
    %Interpolate so that we have 100x more points than starting
    %(each profile stars with silghtly different # of points, here mm/pixel will be constant)
    intpsizes = (mincellsize - 1).*(1/intpfactor)+1; %Size of interpolated profiles
    sprofintp = nan(max(intpsizes),numlines);
    for i = 1:numlines
        sprofintp(1:intpsizes(i),i) = interp1(1:mincellsize(i),proftmp(1:mincellsize(i),i),1:intpfactor:mincellsize(i));
    end

    intppix_mm = pix_mm_lines/intpfactor; %Each interpolated pixel length is (intpfactor)*pixel length

    tol = 0.005; %This does not have to be large after profile interpolation. Need to make sure vessel midpoint dip isn't within threshold.
    %         sprofintp_todel = isnan(sprofintp);
    HMval_repmat = repmat(HMval,[size(sprofintp,1),1]);
    prof_diff = abs(sprofintp - HMval_repmat);
    v = zeros(1,numlines);
    for i = 1:numlines
        small_diffs = find(prof_diff(:,i) < tol); %indices of profile distance to half-max less than tolerance

        [~,groupsep] = max(diff(small_diffs)); %last ind of group 1

        [~,group1ID] = min(abs(sprofintp(small_diffs(1:groupsep)) - HMval(i))); %Ind of closest group 1 point, within group 1
        group1ID = group1ID + small_diffs(1) - 1; %Ind of closest group 1 point, within total profile
        [~,group2ID] = min(abs(sprofintp(small_diffs(groupsep+1:end)) - HMval(i))); %Ind of closest group 2 point, within group 2
        group2ID = group2ID + small_diffs(groupsep + 1) - 1; %Ind of closest group 2 point, within total profile

        v(i) = (group2ID - group1ID)/intppix_mm(i);
    end
    diam(j,:) = v;
end
toc

%Save diameter and parameters used. 
tmpdiam.diam = diam; %Here beforepuff is col 1 puff is 2 afterpuff is 3
tmpdiam.sf = sf; %FW(scalefactor)Max 
tmpdiam.MedFiltVal = MedFiltVal;
tmpdiam.ProfMedFiltVal = ProfMedFiltVal;

% toc
disp('FWHM Done');
save('tmpdiam_2P_sf3.mat','tmpdiam')






%%
%Do outlier correction on diam
% todel= zeros(size(diam,2),1);
% figure
% for i=1:size(diam,2)
%     diamtmp = diam(:,i);
%     plot(diamtmp)
%     pause()
%     isusable = input('Is timeseries usable? 1 yes 0 no ');
%     if isusable == 1
%         [xneg,yneg] = ginput(1);
%         [xpos,ypos] = ginput(1);
%         outlneg = find(diamtmp < yneg);
%         outlpos = find(diamtmp > ypos);
%         
%         tmpbinaryneg = zeros(max(outlneg),1);
%         tmpbinaryneg(outlneg) = 1;
%         tmpbinarypos = zeros(max(outlpos),1);
%         tmpbinarypos(outlpos) = 1;
% 
%         CCneg = bwconncomp(tmpbinaryneg); 
%         CCpos = bwconncomp(tmpbinarypos);
% 
%         for j=1:numel(CCneg.PixelIdxList)
%             if CCneg.PixelIdxList{j}(end) ~= size(diam,1) && CCneg.PixelIdxList{j}(1) ~=1
%                 diamtmp(CCneg.PixelIdxList{j}) = mean([diamtmp(CCneg.PixelIdxList{j}(1)-1),diamtmp(CCneg.PixelIdxList{j}(end)+1)]);
%             elseif CCneg.PixelIdxList{j}(end) == size(diam,1) && CCneg.PixelIdxList{j}(1) ~=1
%                 diamtmp(CCneg.PixelIdxList{j}) = mean([diamtmp(CCneg.PixelIdxList{j}(1)-1),diamtmp(CCneg.PixelIdxList{j}(end))]);
%             elseif CCneg.PixelIdxList{j}(end) == size(diam,1) && CCneg.PixelIdxList{j}(1) ==1
%                 diamtmp(CCneg.PixelIdxList{j}) = mean([diamtmp(CCneg.PixelIdxList{j}(1)),diamtmp(CCneg.PixelIdxList{j}(end))]);
%             elseif CCneg.PixelIdxList{j}(end) ~= size(diam,1) && CCneg.PixelIdxList{j}(1) ==1
%                 diamtmp(CCneg.PixelIdxList{j}) = mean([diamtmp(CCneg.PixelIdxList{j}(1)),diamtmp(CCneg.PixelIdxList{j}(end)+1)]);
%             end
%         end
%         for j=1:numel(CCpos.PixelIdxList)
%             if CCpos.PixelIdxList{j}(end) ~= size(diam,1) && CCpos.PixelIdxList{j}(1) ~=1
%                 diamtmp(CCpos.PixelIdxList{j}) = mean([diamtmp(CCpos.PixelIdxList{j}(1)-1),diamtmp(CCpos.PixelIdxList{j}(end)+1)]);
%             elseif CCpos.PixelIdxList{j}(end) == size(diam,1) && CCpos.PixelIdxList{j}(1) ~=1
%                 diamtmp(CCpos.PixelIdxList{j}) = mean([diamtmp(CCpos.PixelIdxList{j}(1)-1),diamtmp(CCpos.PixelIdxList{j}(end))]);
%             elseif CCpos.PixelIdxList{j}(end) == size(diam,1) && CCpos.PixelIdxList{j}(1) ==1
%                 diamtmp(CCpos.PixelIdxList{j}) = mean([diamtmp(CCpos.PixelIdxList{j}(1)),diamtmp(CCpos.PixelIdxList{j}(end))]);
%             elseif CCpos.PixelIdxList{j}(end) ~= size(diam,1) && CCpos.PixelIdxList{j}(1) ==1
%                 diamtmp(CCpos.PixelIdxList{j}) = mean([diamtmp(CCpos.PixelIdxList{j}(1)),diamtmp(CCpos.PixelIdxList{j}(end)+1)]);
%             end
%         end
%         diam(:,i) = diamtmp;
%     else
%         todel(i) = 1;
%     end
% end
% 
% diam(:,logical(todel)) = [];
% figure
% plot(diam*1000); %um


%% Calculate cross correlations
% time = (1:size(diam,1))/toplot.rate;
% intptime = (1:0.1:size(diam,1))/toplot.rate;
% diammed = diam;
% diammed = medfilt1(diammed,5,[],1);
% meandiam = mean(diammed,1);
% diammed_mean = diammed - repmat(meandiam,size(diammed,1),1);
% pts = sum(~todel);
% % firstseg = repmat(diammed_mean(:,1),1,numpts);
% %Interpolate diameters to get smoother estimate of correlation
% diammed_mean = interp1(time,diammed_mean,intptime);
% 
% [cmat,lagmat] = xcorr(diammed_mean,'Normalized');
% 
% % figure; plot(lagmat*(intptime(2)-intptime(1)),cmat(:,1:size(diammed_mean,2))); xlim([-10 10]);
% figure; plot(lagmat*(intptime(2)-intptime(1)),cmat(:,1:pts)); xlim([-10 10]);
% intprate = 1/(intptime(2)-intptime(1)); xline(0);
% 
% %% Calculate distances
% mask = imread('Mask.tif');
% pix_mm = 1079;
% 
% dmap = bwdistgeodesic(mask,round(mean(x(1:2))),round(mean(y(1:2))));
% geodist(1) = 0;
% counter = 2;
% for i = 3:2:size(x)
%     geodist(counter) = dmap(round(mean(y(i:i+1))),round(mean(x(i:i+1))));
%     counter = counter + 1;
% end
% geodist(logical(todel)) = [];
% 
% %Calc peak lags 
% lagpts = round(intprate*5);
% zerotmp = find(lagmat == 0);
% cmat = cmat(zerotmp-lagpts:zerotmp+lagpts,1:pts); %Just take points 10s around zero
% lagmat = lagmat(zerotmp-lagpts:zerotmp+lagpts);
% 
% iterpts = 1:size(diammed_mean,2);
% counter = 1;
% for i= iterpts
%     corrtmp = cmat(:,i);
% %     plot(lagmat,corrtmp)
% %     pause()
% %     thresh = input('Type corr threshold for this location')
% %     [~,maxind(counter)] = findpeaks(corrtmp,"NPeaks",1,'MinPeakHeight',thresh);
% %     lagmat(maxind)
%     [~,maxind(counter)] = max(corrtmp);
%     counter = counter + 1;
% end
% %Manually exclude points for poor diameter extraction
% [xData, yData] = prepareCurveData(geodist/pix_mm,lagmat(maxind)/intprate);
% % Set up fittype and options.
% ft = fittype( 'poly1' );
% opts = fitoptions( 'Method', 'LinearLeastSquares' );
% opts.Lower = [-Inf -Inf];
% opts.Upper = [Inf Inf];
% [fitresult, gof] = fit(xData, yData, ft, opts );
% figure %Phase Grad vs VasoFreq
% scatter(geodist/pix_mm,lagmat(maxind)/intprate,'filled');
% hold on 
% plot(fitresult, xData, yData );
% xlabel('Distance (mm)','Interpreter','latex');
% ylabel('Lag at Peak Xcorr (s)','Interpreter','latex');
% 
% %% diamresults = struct();
% trial = 3;
% % diamresults(trial).folder = extractAfter(files(1).folder,'FileTransfer\');
% diamresults(trial).folder = extractAfter(files(1).folder,'DataAnalysis\');
% diamresults(trial).puff = 'puff';
% diamresults(trial).timeseries = diam;
% diamresults(trial).corrmat = cmat;
% diamresults(trial).distmat = geodist;
% diamresults(trial).maxlags = lagmat(maxind);
% diamresults(trial).corrv = 1/fitresult.p1; %(mm/s)
% diamresults(trial).corrR2 = gof.rsquare;
% diamresults(trial).lines_x = x;
% diamresults(trial).lines_y = y;
% diamresults(trial).todel = todel;
% time = 1:size(diam,1);
% time = time/toplot.rate;
% diamresults(trial).time = time;
% diamresults(trial).lagmat = lagmat;
% diamresults(trial).imstart = imstart;
% diamresults(trial).imend = imend;
% save('diamresults.mat','diamresults');
% 
% title([diamresults(trial).folder],'Interpreter','none');
% legend off
% savefig(['Diam Xcorr Results ',diamresults(trial).puff,'.fig'])
% 
% figure
% plot(time,mean(diam,2)*1000);
% xlabel('Time (s)','Interpreter','latex');
% ylabel('Diameter FWHM (um)','Interpreter','latex');
% title('Mean Extracted Diameter','Interpreter','latex');
% savefig(['Mean Diam ',diamresults(trial).puff,'.fig'])
% 
% %% Do the same for calcium signal waves
% stim = 'afterpuff';
% if strcmp(stim,'beforepuff');
%     trial = 1;
% elseif strcmp(stim,'puff');
%     trial = 2;
% elseif strcmp(stim,'afterpuff')
%     trial = 3;
% end
% x = diamresults(trial).lines_x;
% y = diamresults(trial).lines_y;
% wave = h5read('04.May.2023_17.15.32_JD221115M1_G05R22_dualpuff_0.1Hz_1.5xObj8xzoomsvd_afterpuff_wave.h5','/wave');
% time = (1:size(wave,2))/toplot.rate;
% intptime = (1:0.1:size(wave,2))/toplot.rate; %10x interpolation
% intprate = 1/(intptime(2)-intptime(1));
% 
% imstart = toplot.startframe;
% if imstart == 0
%     imstart = 1;
% end
% imend = toplot.endframe;
% if imstart ~= diamresults(trial).imstart %Startframe used from toplot.startframe
%     disp('Trials dont start at the same frame');
%     pause();
% end
% imiter = imstart:2:imend; %Frames taken for diameter extraction
% % if length(imiter) ~= size(wave,2)
% %     disp('diameter and wave arent same length');
% %     pause();
% %     
% %     %Beforepuff, this means wave probably has one fewer frame.
% %     wave(:,end) = []; %Check this every trial to be sure
% % end
% % Just calc results mat here, fix lengths in the section after this.
% 
% mask = imread('Mask.tif');
% pix_mm = 1079;
% dmap = bwdistgeodesic(mask,round(mean(x(1:2))),round(mean(y(1:2))));
% counter = 1;
% for i=1:2:size(x,1)
%     Caind(counter) = sub2ind(size(dmap),round(mean(y(i:i+1))),round(mean(x(i:i+1))));
%     counter = counter + 1;
% end
% for i=1:length(Caind)
%     findmaskind(i) = find(toplot.mask_ind == Caind(i));
% end
% %Delete points and use only those used for diameter extraction
% todel = diamresults(trial).todel;
% Caind(logical(todel)) = [];
% findmaskind(logical(todel)) = [];
% 
% wavenums = toplot.skel_label(findmaskind);
% wavemat = wave(wavenums,:)';
% wavematfilt = medfilt1(wavemat,5);
% wavematfilt = interp1(time,wavematfilt,intptime);
% 
% [cmat,lagmat] = xcorr(wavematfilt,'Normalized');
% pts = (numel(todel) - sum(todel));
% figure; plot(lagmat/intprate,cmat(:,1:pts)); xlim([-10 10]);
% 
% lagpts = round(intprate*10);
% zerotmp = find(lagmat == 0);
% cmat = cmat(zerotmp-lagpts:zerotmp+lagpts,1:pts); %Just take points 10s around zero
% lagmat = lagmat(zerotmp-lagpts:zerotmp+lagpts);
% iterpts = 1:size(wavematfilt,2);
% maxind = zeros(size(wavematfilt,2),1);
% counter = 1;
% for i= iterpts
%     corrtmp = cmat(:,i);
%     [~,maxind(counter)] = max(corrtmp);
%     counter = counter + 1;
% end
% 
% counter = 2;
% for i = 3:2:size(x)
%     geodist(counter) = dmap(round(mean(y(i:i+1))),round(mean(x(i:i+1))));
%     counter = counter + 1;
% end
% geodist(logical(todel)) = [];
% 
% [xData, yData] = prepareCurveData(geodist/pix_mm,lagmat(maxind)/intprate);
% % Set up fittype and options.
% ft = fittype( 'poly1' );
% opts = fitoptions( 'Method', 'LinearLeastSquares' );
% opts.Lower = [-Inf -Inf];
% opts.Upper = [Inf Inf];
% [fitresult, gof] = fit(xData, yData, ft, opts );
% figure %Phase Grad vs VasoFreq
% scatter(geodist/pix_mm,lagmat(maxind)/intprate,'filled');
% hold on
% plot(fitresult, xData, yData );
% 
% % caresults = struct();
% % trial = 1;
% caresults(trial).folder = diamresults(trial).folder;
% caresults(trial).puff = stim;
% caresults(trial).timeseries = wavemat;
% caresults(trial).timeseriesmedfilt = wavematfilt;
% caresults(trial).corrmat = cmat;
% caresults(trial).distmat = geodist;
% caresults(trial).maxlags = lagmat(maxind);
% caresults(trial).corrv = 1/fitresult.p1; %(mm/s)
% caresults(trial).corrR2 = gof.rsquare;
% caresults(trial).lines_x = x;
% caresults(trial).lines_y = y;
% caresults(trial).todel = todel;
% time = 1:size(wave,2);
% time = time/toplot.rate;
% caresults(trial).time = time;
% caresults(trial).lagmat = lagmat;
% caresults(trial).imstart = imstart;
% caresults(trial).imend = imend;
% save('caresults.mat','caresults');
% 
% %% Combine and plot summary statistics. Should use same points for diameter and Ca2+..
% cd('C:\Users\duckw\Desktop\FileTransfer');
% cafiles = dir('**\caresults.mat');
% diamfiles = dir('**\diamresults.mat');
% 
% for i=1:length(diamfiles)
%     cd(diamfiles(i).folder);
%     diammat = load('diamresults.mat');
%     diammat = diammat.diamresults;
%     camat = load('caresults.mat');
%     camat = camat.caresults;
%     if i==1
%         for j = 1:length(diammat)
%             diamv(j) = diammat(j).corrv;
%             diamr2(j) = diammat(j).corrR2;
%             cav(j) = camat(j).corrv;
%             car2(j) = camat(j).corrR2;
%         end
%     else
%         for j = 1:length(diammat)
%             diamv = [diamv,diammat(j).corrv]; 
%             diamr2 = [diamr2,diammat(j).corrR2];
%             cav = [cav,camat(j).corrv];
%             car2 = [car2,camat(j).corrR2];
%         end
%     end
% end
% cd('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\TxRed_GCaMP\5_23_23_Figs');
% 
% %Histogram of velocity (diameter and Ca2+)
% figure
% histogram(abs(diamv),'FaceColor','b','BinWidth',0.25)
% hold on
% histogram(abs(cav),'FaceColor','g','BinWidth',0.25)
% 
% %Scatter R^2 vs velocity (colored by diameter or Ca2+)
% figure
% scatter(abs(diamv),diamr2,'filled','b');
% hold on;
% scatter(abs(cav),car2,'filled','r');
% ylim([0 1]); xlim([0 4.5]); legend({'Diameter','Ca2+'});
% xlabel('Travel speed (mm/s)','Interpreter','latex');
% ylabel('Lag vs. Distance $R^2$','Interpreter','latex');
% ax = gca; ax.TickLabelInterpreter = 'latex'; ax.FontSize = 12;
% R2WtdDiamv = dot(abs(diamv),diamr2)/sum(diamr2);
% R2WtdCav = dot(abs(cav),car2)/sum(car2);
% % xline(R2WtdCav,'r'); xline(R2WtdDiamv,'b');
% r2wtdD2ndmoment = dot(diamr2,abs(diamv).^2)/sum(diamr2);
% r2wtdDuncert = sqrt(r2wtdD2ndmoment - R2WtdDiamv^2);
% r2wtdCa2ndmoment = dot(car2,abs(cav).^2)/sum(car2);
% r2wtdCauncert = sqrt(r2wtdCa2ndmoment - R2WtdCav^2);
% str = sprintf('Diam $v$ = %.1f $\pm$ %.1f, Ca $v$ = %.1f $\pm$ %.1f',R2WtdDiamv,r2wtdDuncert/sqrt(numel(diamv)),R2WtdCav,r2wtdCauncert/sqrt(numel(cav)));
% title({'$R^2$-weighted mean $\pm$ SEM',['Diam $v$ = ',num2str(round(R2WtdDiamv,2)),' $\pm$ ',num2str(round(r2wtdDuncert/sqrt(numel(diamv)),2)),', Ca$2^+$ $v$ = ',num2str(round(R2WtdCav,2)),' $\pm$ ',num2str(round(r2wtdCauncert/sqrt(numel(cav)),2))]},'Interpreter','latex');
% savefig('R2vsSpeed_Scatter_9Trials.fig');
% print(gcf, '-depsc2', 'R2vsSpeed_Scatter_9Trials');
% 
% %Calc cross-correlation of calcium and diameter (histogram of tial-wise
% %peak corr lag)
% %% Can also do this for each vessel cross section
% clear
% cd('C:\Users\duckw\Desktop\FileTransfer');
% cafiles = dir('**\caresults.mat');
% diamfiles = dir('**\diamresults.mat');
% trial = 3;
% 
% for i=1:length(diamfiles) %Iterate over beforepuff, puff, afterpuff
%     cd(diamfiles(trial).folder);
%     diammat = load('diamresults.mat');
%     diammat = diammat.diamresults;
%     camat = load('caresults.mat');
%     camat = camat.caresults;
%     toplotfiles = dir('*puff.mat');
%     stim = diammat(i).puff;
%     isfile = zeros(length(toplotfiles),1);
%     for ii = 1:length(toplotfiles)
%         if contains(toplotfiles(ii).name,['_',stim])
%             isfile(ii) = 1;
%         end
%     end
%     isfile = find(isfile);
%     load(toplotfiles(isfile).name);
% 
%     diam = mean(diammat(i).timeseries,2)*1000; %um
%     diam_mean = diam - mean(diam);
% 
%     wave = mean(camat(i).timeseries,2);
% 
%     if length(diam_mean) ~= length(wave)
%         diam_mean(end) = [];
%     end
%     %Scale wave to diameter for plotting only
%     diamamp = max(diam_mean) - min(diam_mean);
%     waveamp = max(wave) - min(wave);
%     wavescaled = wave.*(diamamp/waveamp);
% 
%     %Calc Cross correlation, spectra, phase
%     [cmat,lagmat] = xcorr(diam_mean,wave,'Normalized');
%     time = 1:length(diam_mean); time = time/toplot.rate;
%     figure; subplot(2,1,1);
%     plot(time,medfilt1(diam_mean,5),'b'); hold on; plot(time,wavescaled,'r'); xlim([0 250]);
%     subplot(2,1,2)
%     plot(lagmat/toplot.rate,cmat); xlim([-30 30]); ylim([-1 1]); xline(0);
%     [maxcorr,maxcorrind] = min(cmat); 
%     pklag = lagmat(maxcorrind)/toplot.rate;
%     title(['Lag at max correlation = ',num2str(round(pklag,2)),' s',', $|$Corr$|$ = ',num2str(round(abs(maxcorr),2))],'Interpreter','latex');
% 
%     % Save large figure with spectrogram
%     rate = toplot.rate;
%     timewave = time + 1/(2*rate); %Ca signal captured as second frame
%     %Interpolate each signal to get equal time points
%     tqD = time(1):(1/(2*rate)):time(end);
%     tqCa = timewave(1):(1/(2*rate)):timewave(end);
%     diam_q_tmp = interp1(time,diam_mean,tqD);
%     wave_q_tmp = interp1(timewave,wave,tqCa);
% %     figure; plot(tqD,diam_q_tmp,'b'); hold on; plot(tqCa,wave_q_tmp,'r');
% 
%     %Delete first wave1 value so that both start at t=0.75.
%     diam_q = diam_q_tmp(2:end);
%     tqD = tqD(2:end);
%     %Delete last wave1 value to make both time series the same duration
%     wave_q = wave_q_tmp;
%     wave_q(end) = [];
%     tqCa(end) = [];
% 
%     params.Fs = toplot.rate*2; %Interpolated rate is twice actual single depth rate
%     params.pad = 2;
%     params.fpass = [0 params.Fs/2]; %Hz, default is [0 Fs/2]
%     params.err   = [2 .05];
%     params.trialave = 0;
% %     T = stats1.time(end);
%     T = tqD(end);
%     BW = 0.02; %600s trial -> TBW =~ 7
%     params.tapers = [round(T*BW),round(2*T*BW-1)]; %Time-BW product and number of tapers
%     addpath(genpath('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\chronux_2_12'))
%     % movingwin = [100,1]; %Length of moving window and step size
%     movingwin = [60,1]; %Length of moving window and step size
%     params2 = params;
%     Tmov = movingwin(1);
%     % BW2 = 0.12;
%     BW2 = 0.2;
%     params2.Fs = rate;
%     params2.fpass = [0 rate/2];
%     params2.tapers = [round(Tmov*BW2),round(2*Tmov*BW2-1)];
% 
%     [Sdiam,fdiam] = mtspectrumc(diam_q,params);
%     [Swave,fwave] = mtspectrumc(wave_q,params);
%     if strcmp(stim,'puff')
%         f_peak = 0.05;
%         freq_findx = max(find(round(fdiam,3)==round(f_peak,3)));
%     else
%     figure
%     plot(fdiam,log10(Sdiam)); xlim([0 1]);
%     [wind,~] = ginput(2);
%     findf1 = round(wind(1),3);
%     findf2 = round(wind(2),3);
% %     f1 = find(round(fdiam,3)==round(findf1,3),1,'first'); %Hz lower bound
% %     f2 = find(round(fdiam,3)==round(findf2,3),1,'first'); %Hz upper bound
%     f1 = find(round(fdiam,2)==round(findf1,2),1,'first'); %Hz lower bound
%     f2 = find(round(fdiam,2)==round(findf2,2),1,'first'); %Hz upper bound
%     rmpath(genpath('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\chronux_2_12'))
%     [maxes,f_peaks] = findpeaks(Sdiam(f1:f2),fdiam(f1:f2));
%     maxloc = find(maxes==max(maxes));
%     f_peak = f_peaks(maxloc)
%     freq_findx = max(find(round(fdiam,3)==round(f_peak,3)));
%     addpath(genpath('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\chronux_2_12'))
%     end
% 
%     %Calc coherence and phase between the two signals
%     params_err = params;
%     params_err.err = [2,0.05];
%     [C,phi,S12,S1,S2,f,confC,phistd,Cerr] = coherencyc(diam_q,-wave_q,params_err); %Doing -wave for phase
%     freq_findx = max(find(round(f,3)==round(f_peak,3)));
%     dof = 2*round(2*T*BW-1); %Single value dof (2*tapers)
%     df = 1./((dof/2)-1);
%     p = 0.05; %Confidence level
%     confC = sqrt(1 - p.^df);
% 
%     % CALC CROSS CORRELATION FROM RAW SIGNALS
%     x2 = 1/rate*100 + 1/rate/2;
%     x1 = (1/rate)/2;
%     n = rate*(x2 - x1)+1;
%     lags_tmp = linspace(x1,x2,n);
%     lag_inds_tmp = 1:length(lags_tmp);
%     lags_neg = -lags_tmp;
%     lags1 = [flip(lags_neg),lags_tmp];
%     lag_inds = [flip(-lag_inds_tmp),lag_inds_tmp]; %For shifting vectors
%     totLags = length(lags1);
% 
%     corr1 = NaN(totLags,1);
%     for j=1:totLags
%         lag = lags1(j);
%         lag_ind = lag_inds(j);
%         %     d1 = d1_mean(210:1100);
%         %     d2 = d2_mean(210:1100);
%         %         d1 = d1_mean;
%         %         d2 = d2_mean;
%         %         d1 = diam1_filt(200:700);
%         %         d2 = diam2_filt(200:700);
%         d1 = diam_mean;
%         d2 = wave;
%         %         d1 = diam1_filt;
%         %         d2 = diam2_filt;
% 
%         N = numel(d1);
%         denom1 = std(d1)*sqrt(N-1); %For normalizaion
% 
%         %Shift d2 and do dot product
%         if lag_ind < 0 %Shift d2 to the right
%             d2((end-(abs(lag_ind)-1)):end) = [];
%             d1(1:abs(lag_ind)) = [];
%             %         corr1(i) = dot(d1,d2);
%             meand2 = mean(d2);
%             denom2 = sqrt(sum((d2 - meand2).^2));
%             corr1(j) = dot(d1,d2)/(denom1*denom2); %Normalized version
%         elseif lag_ind > 0 %Shift d2 to the left
%             d2(1:abs(lag_ind)) = [];
%             d1((end-(abs(lag_ind)-1)):end) = [];
%             %         corr1(i) = dot(d1,d2);
%             meand2 = mean(d2);
%             denom2 = sqrt(sum((d2 - meand2).^2));
%             corr1(j) = dot(d1,d2)/(denom1*denom2); %Normalized version
%         end
%     end
%     figure
%     plot(lags1,corr1)
%     xlim([-10 10]); xline(0);
%     [maxcorr,maxlag] = min(corr1);
%     maxlag = lags1(maxlag);
% 
%     % PLOT COMBINED FIGURE
%     strfv = sprintf('%.3f',f_peak);
%     str_spectrum = sprintf('BW = %.3f Hz',BW);
%     str_sg = sprintf('Wind = %.0f s, Step = %.0f s',movingwin(1),movingwin(2));
% 
%     fig = figure('units','inches','outerposition',[0 0 8.5 11]); hold on;
%     subplot(4,2,[1,2])
%     yyaxis left; plot(tqD,medfilt1(diam_q,10) + mean(diam),'b'); ylabel('Diam FWHM (um)','Interpreter','latex'); 
%     hold on; yyaxis right; plot(tqCa,wave_q,'r'); ylabel('GCaMP dF/F','Interpreter','latex'); set(gca,'YDir','reverse');
%     xlabel('Time (s)','Interpreter','latex'); xlim([0 tqD(end)]);
%     title([extractBefore(diammat(i).folder,'\Default'),' ',diammat(i).puff],'Interpreter','none','FontSize',11); 
%     ax = gca; ax.FontSize = 11;
% 
%     subplot(4,2,3)
%     plot(fdiam,log10(Sdiam/sum(Sdiam)),'b'); xlim([0 0.5]); xlabel('Frequency (Hz)','Interpreter','latex'); ylabel('log(Power)','Interpreter','latex');
%     hold on
%     plot(fwave,log10(Swave/sum(Swave)),'r'); xlim([0 0.5]); xlabel('Frequency (Hz)','Interpreter','latex'); ylabel('log(Power)','Interpreter','latex');
%     title(['Spectra, ',str_spectrum],'Interpreter','latex');
%     hold off
% 
%     subplot(4,2,4); plot(lags1,corr1,'b'); grid minor;
%     xlim([-10 10]); xline(0); xlabel('Lag (seconds)','Interpreter','latex'); ylabel('Correlation','Interpreter','latex'); ylim([-1 1]);
%     maxloc = find(corr1 == min(corr1)); title(['lag at max $|$Corr$|$= ',num2str(round(abs(lags1(maxloc)),2)),' s'],'Interpreter','latex')
%     hold on
%     plot(lags1(maxloc),corr1(maxloc),'o','MarkerSize',4,'MarkerFaceColor','r');
%     hold off
% 
%     [S2sg,t,f2] = mtspecgramc(wave,movingwin,params2);
%     ax(2) = subplot(4,2,7);
%     imagesc(t,f2,log10(S2sg')); colormap(ax(2),'jet'); colorbar; %clim([-2 1]);
%     set(gca,'YDir','normal','YLim',[0,0.5]); title('Ca signal spectrogram, 60s moving window','Interpreter','latex');xlabel('Time (s)','Interpreter','latex');ylabel('Frequency (Hz)','Interpreter','latex');
% 
%     [SDsg,t,f2] = mtspecgramc(diam_mean,movingwin,params2);
%     ax(2) = subplot(4,2,8);
%     imagesc(t,f2,log10(SDsg')); colormap(ax(2),'jet'); colorbar; %clim([-2 1]);
%     set(gca,'YDir','normal','YLim',[0,0.5]); title('Diameter signal spectrogram, 60s moving window','Interpreter','latex');xlabel('Time (s)','Interpreter','latex');ylabel('Frequency (Hz)','Interpreter','latex');
% 
%     subplot(4,2,5)
%     plot(f,C,'b'); hold on; plot(f,Cerr,'Color',[0,0,1,0.1]); ylim([0 1]); yline(confC,'-.','$|C|_{0.95}$','Interpreter','latex','Alpha',0.9,'LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom')
%     xlabel('Frequency (Hz)','Interpreter','latex'); ylabel('Magnitude Coherence','Interpreter','latex'); title('Coherence','Interpreter','latex');
%     xline(f_peak,'-.','$f_v$','Interpreter','latex','Alpha',0.9);
%     xlim([0 0.5]);
%     subplot(4,2,6)
%     plot(f,phi,'b'); hold on; plot(f,phi + 2*phistd','Color',[0,0,1,0.2]); plot(f,phi - 2*phistd','Color',[0,0,1,0.2]); xlim([0 0.5]); ylim([-pi pi]);
%     yline(0,'Alpha',0.9); xlabel('Frequency (Hz)','Interpreter','latex'); ylabel('Relative Phase (rad)','Interpreter','latex');
%     xline(f_peak,'-.','$f_v$','Interpreter','latex','Alpha',0.9);
%     phi_fv = phi(freq_findx);
%     tau = phi_fv/(2*pi*f_peak);
%     str1 = sprintf('%.2f',phi_fv);
%     str2 = sprintf('%.2f',abs(tau));
%     str3 = sprintf('%.2f',f_peak);
%     title({'Relative Phase and 95 Percent CI',['$\phi_{fv}$ = ',str1,', $f_v$ = ',str3,' $\tau$ = ',str2,' s']},'Interpreter','latex')
% 
%     print(fig,strrep(['JD',extractAfter(extractBefore(diammat(i).folder,'\Default'),'JD'),'_',diammat(i).puff],'.',','),'-dpdf','-r0')
%     
% %     print(fig,'04_May_23_JD221115M1_dualpuff_0_1Hz_Afterpuff','-dpdf','-r0','vector')
% end
% 
% %% Calc cross correlation between Ca and diam at each point
% % figure
% % subplot(2,1,1)
% % plot(diammat(1).time,diammat(1).timeseries*1000,'Color',[1,0,0,0.1]); hold on;
% % plot(diammat(1).time,mean(diammat(1).timeseries,2)*1000,'Color',[0,0,0]);
% clearvars -except diammat camat
% 
% trial = 3;
% for i=1:size(diammat(trial).timeseries,2)
%     dtmp = diammat(trial).timeseries(:,i);
%     wavetmp = -camat(trial).timeseries(:,i);
%     if length(dtmp) > length(wavetmp)
%         dtmp(end) = [];
%     end
%     dtmp = dtmp - mean(dtmp);
%     wavetmp = wavetmp - mean(wavetmp);
% 
%     [corr(:,i),lags] = xcorr(dtmp,wavetmp,'Normalized');
%     [maxcorr(i),maxlag] = max(corr(:,i));
%     lag(i) = lags(maxlag);
% end
% 
% figure
% plot(lags/toplot.rate,corr)
% figure
% scatter(diammat(trial).distmat/1079,lag/toplot.rate,'filled')
% xlabel('Distance from first seg (mm)','Interpreter','latex');
% ylabel('Time delay from GCaMP to Diameter','Interpreter','latex');
% 
% %% Calc theoretical diameter signal lag curve from Ca lags & Ca-D transfer lags
% clear
% 
% load('diamresults.mat');
% load('caresults.mat');
% load('04.May.2023_18.05.01_JD221115M1_G05R22_dualpuff_0.05Hz_1.5xObj8xzoom_ves2svd_beforepuff.mat');
% trial = 3; %1 for beforepuff etc.
% x = diamresults(trial).distmat/1079; %Distances
% time = (1:size(diamresults(trial).timeseries,1))/toplot.rate;
% intptime = (1:0.1:size(diamresults(trial).timeseries,1))/toplot.rate; %10x interpolation
% intprate = 1/(intptime(2)-intptime(1)); %Rate that was used for diameter and Ca cross correlation calcs originally
% 
% Calags = caresults(trial).maxlags/intprate; %seconds
% Dlags = diamresults(trial).maxlags/intprate; %seconds 
% 
% %Now calc theoretical Dlags using Ca-D cross-correlation information
% for i=1:size(diamresults(trial).timeseries,2)
%     dtmp = diamresults(trial).timeseries(:,i);
%     wavetmp = -caresults(trial).timeseries(:,i);
%     if length(dtmp) > length(wavetmp)
%         dtmp(end) = [];
%     end
%     dtmp = dtmp - mean(dtmp);
%     wavetmp = wavetmp - mean(wavetmp);
% 
%     [corr(:,i),lags] = xcorr(dtmp,wavetmp,'Normalized');
%     [maxcorr(i),maxlag] = max(corr(:,i));
%     lag(i) = lags(maxlag)/toplot.rate; %Seconds (lag between Ca and D at each location)
% end
% ft = fittype( 'poly1' );
% opts = fitoptions( 'Method', 'LinearLeastSquares' );
% opts.Lower = [-Inf -Inf];
% opts.Upper = [Inf Inf];
% figure
% scatter(x,lag,'filled'); xlabel('Distance (mm)','Interpreter','latex');
% ylabel('Ca-Diameter time lag (s)','Interpreter','latex');
% [xData, yData] = prepareCurveData(x,lag);
% [fitresult, gof] = fit(xData, yData, ft, opts );
% title({strrep(['JD',extractAfter(extractBefore(diamresults(trial).folder,'\Default'),'JD'),'_',diamresults(trial).puff],'.',','),['dt/dx Ca-D = ',num2str(round(fitresult.p1,2)),' s/mm']},'Interpreter','none');
% hold on; xplot = 0:0.1:0.9; yplot = fitresult.p2 + xplot.*fitresult.p1; plot(xplot,yplot,'b');
% 
% Dlags_th = lag(1) - (-Calags + lag); %See notes
% 
% 
% 
% %Plot and save summary figure
% figure
% scatter(x,Dlags,'filled') %Experimental diameter lag vs distance
% hold on;
% scatter(x,Calags,'filled','r') %Experimental Ca lag vs distance
% scatter(x,Dlags_th,'Marker',"x",'MarkerEdgeColor','r','SizeData',75)
% xlabel('Distance (mm)','Interpreter','latex');
% ylabel('Time Lag (s)','Interpreter','latex');
% %Calc velocities
% [xData, yData] = prepareCurveData(x,Dlags);
% [fitresult, gof] = fit(xData, yData, ft, opts );
% hold on; xplot = 0:0.1:0.9; yplot = fitresult.p2 + xplot.*fitresult.p1; plot(xplot,yplot,'b');
% [xData, yData] = prepareCurveData(x,Dlags_th); [fitresult_th, gof_th] = fit(xData, yData, ft, opts );
% 1/fitresult_th.p1
% hold on; xplot = 0:0.1:0.9; yplot = fitresult_th.p2 + xplot.*fitresult_th.p1; plot(xplot,yplot,'r');
% %%
% 
% [xData, yData] = prepareCurveData(x,Calags); [fitresult_ca, gof_ca] = fit(xData, yData, ft, opts );
% 
% hold on; xplot = 0:0.1:0.9; yplot = fitresult_ca.p2 + xplot.*fitresult_ca.p1; plot(xplot,yplot,'r');
% title({strrep(['JD',extractAfter(extractBefore(diamresults(trial).folder,'\Default'),'JD'),'_',diamresults(trial).puff],'.',','),['vD_Exp = ',num2str(round(1/fitresult.p1,2)),' vD_Th = ',num2str(round(1/fitresult_th.p1,2)),' vD_ca = ',num2str(round(1/fitresult_ca.p1,2)),' mm/s']},'Interpreter','none');
% legend({'Experiment','Asuume Const Ca-D lag','Use Exp. Ca-D lag'})
% % savefig([strrep(['JD',extractAfter(extractBefore(diamresults(trial).folder,'\Default'),'JD'),'_',diamresults(trial).puff],'.',','),'_DiamLagCalc.fig'])
% 



