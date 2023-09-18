%VesselD_Ca2_Extraction.m
%Code for extraction and comparison of traveling signals in TxRed/GCaMP
%Updated version 5/15/23
%Version to extract diameter on muralis

%% Process to get diameter signal in matlab (5/1/23)
clear; close all;
% data_folder = 'Y:\Data backup\20230613 JD221211F1';
data_folder = 'C:\Users\duckw\Desktop\FileTransfer';
animal = 'JD221223F2'; %For file naming

cd(data_folder);
files = dir('*001.tif');

file = 1;
namestr = extractBefore(files(file).name,'_00001.tif')
pix_um_x = 1; %1 um x resolution
pix_um_y = 4; %0.25 um y resolution
rate = input(['INPUT TRIAL FRAME RATE ',namestr,' (Hz) '])

%Load images
cd('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P')
outvars = pyrunfile("TiffDimDetector.py 'C:\Users\duckw\Desktop\FileTransfer\ROI3_00001.tif'",'test');
ims = str2double(outvars.char)

%Make vessel masks. Tried tiffreadVolume but was very slow compared to
%reading with python
maskims = round(100*rate); %100 seconds * rate (images/s). Want multiple cycles in mask
if mod(maskims,2) == 1 %If odd add 1 image for loading
    maskims = maskims + 1;
end

outvarsD = pyrunfile("TiffImReader.py 'C:\Users\duckw\Desktop\FileTransfer\ROI3_00001.tif'",'test', r1 = int16(1), r2 = int16(maskims), r3 = int16(2)); 
Dmask = squeeze(mean(double(outvarsD),1)); %D mask is average of first 2000 images.
outvarsCa = pyrunfile("TiffImReader.py 'C:\Users\duckw\Desktop\FileTransfer\ROI3_00001.tif'",'test', r1 = int16(0), r2 = int16(maskims - 1), r3 = int16(2)); 
Camask = squeeze(mean(double(outvarsCa),1));
% figure; subplot(2,1,1); imagesc(Dmask); daspect([1,4,1]); axis off; subplot(2,1,2); imagesc(Camask); daspect([1,4,1]); axis off;
fusemask = imfuse(Dmask,Camask); figure; imshow(fusemask); daspect([1,4,1]); 
cd(['Z:\Thomas_1P\JD_GC_Diam2P\6_13_23\8_13_23Analysis\',animal,'\',namestr]);
% cd(['Z:\Thomas_1P\JD_GC_Diam2P\6_13_23\7_24Analysis','\',namestr]);
savefig('Ca_Lumen_imfust.fig');
%Save fuse image after resizing
fusemask_resize = imresize(fusemask,[size(fusemask,1),size(fusemask,2)*4],'nearest'); %Try not median filtering first
figure; imshow(fusemask_resize); daspect([1,1,1]); 
savefig('Ca_Lumen_imfust_resize.fig');
%Make combined mask 
Dmask = Dmask./max(Dmask(:)); Camask = Camask./max(Camask(:)); 
Combmask = Dmask + Camask;

% cd('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Cy5.5_GCaMP_2P');
cd(['Z:\Thomas_1P\JD_GC_Diam2P\6_13_23\8_13_23Analysis\',animal,'\',namestr]);
% if ~isfolder([animal,'_',namestr])
%     mkdir([animal,'_',namestr]);
% end
% cd([animal,'_',namestr])

%Make mask from Ca + Lumen images
if ~isfile('Mask1.tif')
    figure; imshow(Combmask); daspect([1,4,1]);
    imwrite(Combmask,'Mask1.tif');
end
%Make mask from just lumen images.
if ~isfile('lumMask1.tif')
    figure; imshow(Dmask); daspect([1,4,1]);
    imwrite(Dmask,'lumMask1.tif');
end

%Plot std to help with mask creation if needed
stdmask = squeeze(std(double(outvarsD),[],1));
stdmask = stdmask./max(stdmask(:));
imwrite(stdmask,'stdmask2.tif');

%We want a Calcium mask (GCaMP wave segments), use this for diameter calcs
%as well. Could make & use lumen mask if needed.

%Use photoshop to create binary mask from lumMask1.tif. First analysis used
%combmask(Mask1.tif), now using lumen mask.
mask = logical(rgb2gray(imread('Mask.tif')));
camask = logical(rgb2gray(imread('CaMask.tif')));

%% Make vessel graph for Ca extraction
% mask1 is 0.25 x 1 um, mask is 0.25 x 0.25 resolution.

%For GCaMPintensity extraction make sure to use Ca-mask not lumen mask.
masksize = size(camask);
orig_masksize = size(Combmask);

addpath('X:\DataAnalysis\XianCode\XianCode\')
% mask_origsize = imresize(mask,[size(fusemask,1),size(fusemask,2)]);

opt.downsample_rate = 20; %Tunable parameter
opt.imdilate_disk_r = 1;
opt.min_kept_cc_num_pixel = 35;
opt.rm_neighbor_out_of_mask_Q = true;
opt.vis_mask_Q = true; % Show the cleaned up mask
[skel_neighbor_ind, ~] = fun_get_skeleton_neighbor_ind(camask, opt);
[skel_neighbor_ind_unique, tmp_unique_idx, ~] = unique(cat(2, skel_neighbor_ind{:}), 'stable');
recon_mask = false(size(camask));
recon_mask(skel_neighbor_ind_unique) = true;
num_skel_ind_unique = numel(skel_neighbor_ind_unique);
tmp_skel_label = repelem(1 : numel(skel_neighbor_ind), cellfun(@numel, skel_neighbor_ind));
skel_label_unique = tmp_skel_label(tmp_unique_idx);
toplot.skel_label = skel_label_unique;
toplot.mask_ind = skel_neighbor_ind_unique;

%Get images, average them and save matrix for wave extraction.
%rates range from 30.03 to 58.3 Hz. Average all data within 3 frames to get
%minimum rate of 10Hz.
cd('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P')
%% JUST Load everything in, resize, and average appropriate number of images!!!
tic;
tmpimsCa = pyrunfile("TiffImReader.py 'C:\Users\duckw\Desktop\FileTransfer\ROI3_00001.tif'",'test', r1 = int32(0), r2 = int32(round(ims)), r3 = int32(2)); 
Ca_data = double(tmpimsCa); clearvars tmpimsCa; %Ca_data = shiftdim(Ca_data,1); %Do this after background subtraction
toc
%% AVERAGE EVERY N frames to increase image quality while reducing frame rate

figure; imagesc(squeeze(Ca_data(100,:,:))); daspect([1,4,1]); title('Select rectangle for noise quantification');
[noisecol,noiserow] = ginput(2);
Ca_noise = Ca_data(:,round(noiserow(1)):round(noiserow(2)),round(noisecol(1)):round(noisecol(2)));
figure; histogram(Ca_noise); title('Choose Bounsd for gaussian fit');%Fit gaussian to noise
Ca_noise = Ca_noise(:);
[gausscol, ~] = ginput(2); 
Ca_noise_tofit = Ca_noise(and(Ca_noise > gausscol(1),Ca_noise < gausscol(2)));
pd = fitdist(Ca_noise_tofit,'Normal');
noise_mean = pd.mu %Noise level to be applied as offset manually.
clearvars Ca_noise_tofit

Ca_data = Ca_data - noise_mean; %Now noise is centered around zero (no photons = 0 intensity)
clearvars tmpimsD; Ca_data = shiftdim(Ca_data,1); %Make minimum value zero in Ca_data (int16)


FrameAvg = 3;
NewRate = rate/FrameAvg;
%First clip ending frames to make frames a multiple of FrameAvg
while mod(size(Ca_data,3),FrameAvg) ~= 0
    Ca_data(:,:,end) = [];
end
%Average Calcium data
Ca_data_tmp = reshape(Ca_data,[orig_masksize(1)*orig_masksize(2),size(Ca_data,3)])'; %collapse space dimension, now time x space
Ca_data_tmp = reshape(Ca_data_tmp,[FrameAvg,size(Ca_data,3)/FrameAvg*orig_masksize(1)*orig_masksize(2)]); %Each column is (FrameAvg) frames to be averaged
Ca_data_tmp = mean(Ca_data_tmp,1); %Take mean of every (FrameAvg) frames
Ca_data_tmp = reshape(Ca_data_tmp,[size(Ca_data,3)/FrameAvg,orig_masksize(1)*orig_masksize(2)])'; %Space (pixels) x time
Ca_data_tmp = reshape(Ca_data_tmp,[size(Ca_data,1),size(Ca_data,2),size(Ca_data,3)/FrameAvg]);
Ca_data = Ca_data_tmp; clearvars Ca_data_tmp;

%Convert isotropic skel_label to non-isotropic mask (0.25 x 1)um
map = zeros(size(camask));
map(toplot.mask_ind) = toplot.skel_label;
%Get set of big map indices for each small mask pixel. Set small mask pixel
%equal to mode of big map skel labels.
map_noni = zeros(size(Combmask));
for i = 1:size(Combmask,1)
    for j = 1:size(Combmask,2)
        mapcol = 4*(j-1)+1;
        map_noni(i,j) = mode(map(i,mapcol:mapcol + 3));
    end
end
%Make new skel_neighbor_ind for wave
% skel_neighbor_ind_noni = skel_neighbor_ind;
skel_neighbor_ind_noni = cell(1,1);
for i = 1:max(map_noni(:)) %Number of segments
    skel_neighbor_ind_noni{i,1} = find(map_noni == i)';
end

%CALCUATE Calcium dF/F from aveaged images
wave = fun_get_skeleton_neighbor_stat_from_image_stack(Ca_data, skel_neighbor_ind_noni, 'mean'); %Changed to mean for high mag imaging (most pixels have negligable change in fluorescence)
wave = bsxfun(@rdivide ,bsxfun(@minus, wave, mean(wave, 2)),mean(wave, 2));%calculate dF/F
time = 1:size(wave,2); time = time/NewRate; figure; plot(time,mean(wave,'omitnan')); xlabel('Time (s)','Interpreter','latex'); ylabel('GCaMP dF/F','Interpreter','latex');
% cd('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Cy5.5_GCaMP_2P'); cd([animal,'_',namestr])
cd(['Z:\Thomas_1P\JD_GC_Diam2P\6_13_23\8_13_23Analysis\',animal,'\',namestr]);
%Save GCaMP data
h5create([animal,'_',namestr,'_wave.h5'],'/wave',[size(wave,1), size(wave,2)],'Datatype','double');
h5write([animal,'_',namestr,'_wave.h5'],'/wave', wave, [1 1],[size(wave,1), size(wave,2)]);
clearvars Ca_data

%% LOAD DIAMETER IMAGES AND AVERAGE THEM
cd('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P')
tmpimsD = pyrunfile("TiffImReader.py 'C:\Users\duckw\Desktop\FileTransfer\ROI3_00001.tif'",'test', r1 = int32(1), r2 = int32(round(ims+1)), r3 = int32(2)); 
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
clearvars D_noise_tofit

D_data = D_data - noise_mean; %Now noise is centered around zero (no photons = 0 intensity)
clearvars tmpimsD; D_data = shiftdim(D_data,1); %Make minimum value zero in D_data (int16)
while mod(size(D_data,3),FrameAvg) ~= 0
    D_data(:,:,end) = [];
end
D_data_size = size(D_data);
% Average Lumen Data
D_data_tmp = reshape(D_data,[orig_masksize(1)*orig_masksize(2),D_data_size(3)])'; %collapse space dimension, now time x space
clearvars D_data
D_data_tmp = reshape(D_data_tmp,[FrameAvg,D_data_size(3)/FrameAvg*orig_masksize(1)*orig_masksize(2)]); %Each column is (FrameAvg) frames to be averaged
D_data_tmp = mean(D_data_tmp,1); %Take mean of every (FrameAvg) frames
D_data_tmp = reshape(D_data_tmp,[D_data_size(3)/FrameAvg,orig_masksize(1)*orig_masksize(2)])'; %Space (pixels) x time
D_data_tmp = reshape(D_data_tmp,[D_data_size(1),D_data_size(2),D_data_size(3)/FrameAvg]);
D_data = D_data_tmp; clearvars D_data_tmp;

%% USE CA SEGMENTS TO GET SUM OF INTENSITY OVER TIME
%Data has some negative values..we want no photons to equal 0 in our images
%For now just use min(D_data(:))
Dwave = fun_get_skeleton_neighbor_stat_from_image_stack(D_data, skel_neighbor_ind_noni, 'sum');
%Save Dwave
cd(['Z:\Thomas_1P\JD_GC_Diam2P\6_13_23\8_13_23Analysis\',animal,'\',namestr]);
h5create([animal,'_',namestr,'_lumint.h5'],'/lumint',[size(Dwave,1), size(Dwave,2)],'Datatype','double');
h5write([animal,'_',namestr,'_lumint.h5'],'/lumint', Dwave, [1 1],[size(Dwave,1), size(Dwave,2)]);

%Calc cross correlation with gcamp signal
Dint_wave = mean(Dwave); Dint_wave = Dint_wave - mean(Dint_wave);
Caint_wave = mean(wave); Caint_wave = Caint_wave - mean(Caint_wave);
[corr,lags] = xcorr(Dint_wave,Caint_wave,'Normalized');
figure; plot(lags/NewRate,corr);

%% GET FWHM AT EACH SEGMENT, USE 10 CROSS LINES AVERAGE AT EACH SEGMENT
%Calculate cross lines in the equal pixel-size image, use
%imresize(,'nearest') to resize raw data and calculate profiles.
D_data = reshape(D_data,[size(D_data,1)*size(D_data,2),size(D_data,3)]);
cd(['Z:\Thomas_1P\JD_GC_Diam2P\6_13_23\8_13_23Analysis\',animal,'\',namestr]); %Use diameter data saved here.
h5create([animal,'_',namestr,'_diamdat.h5'],'/frame',[size(D_data,1), size(D_data,2)],'Datatype','double');
h5write([animal,'_',namestr,'_diamdat.h5'],'/frame', D_data, [1 1],[size(D_data,1), size(D_data,2)]);
%% Make cross lines for later diameter FWHM extraction.
% Use skeleton from vsl_graph, this is the common skeleton, everything is
% 1-1.
[~, ~, vsl_graph] = fun_get_skeleton_neighbor_ind(mask, opt);
[skel_neighbor_ind, ~] = fun_get_skeleton_neighbor_ind(mask, opt);
linkinds_tmp = cell2mat(vsl_graph.link.cc_ind);
maskskel = zeros(size(mask));
maskskel(linkinds_tmp) = 1;

% maskskel = bwskel(mask,'MinBranchLength', round(4000*0.1)); %Min branch length 100 um


% skelinds = find(maskskel);
% [skelrow,skelcol] = ind2sub(size(mask),skelinds);

%Now get radii
maskdist = bwdist(~mask);

[noderow,nodecol] = ind2sub(size(mask),vsl_graph.node.pos_ind);
figure; imagesc(maskskel); daspect([1,1,1]); f1 = gcf; f1.Position = [686 130 1021 957]; hold on;
plot(nodecol,noderow,'ro')
keepselecting = true;
counter = 1;
while keepselecting
    [x(counter),y(counter)] = ginput(1); 
    counter = counter + 1;
    keepselecting = input('Keep Selecting Links? 1 yes 0 no ');
end
%Find link that each point belongs to
islink = zeros(length(vsl_graph.link.cc_ind),1);
for i = 1:length(x)
    ptind = sub2ind(size(mask),round(y(i)),round(x(i)));
    for j = 1:length(vsl_graph.link.cc_ind)
        tmpinds = vsl_graph.link.cc_ind{j};
        %Dilate tmpinds to expand matching area
        tmpmap = zeros(size(mask)); tmpmap(tmpinds) = 1;
        tmpinds = find(imdilate(tmpmap,strel('disk',5)));
        if any(ismember(tmpinds,ptind))
            islink(j) = 1;
        end
    end
    i
end
linktodel = ~islink;
linkindcells = vsl_graph.link.cc_ind;
linkindcells(linktodel) = [];
linkinds = cell2mat(linkindcells); %Only inds of links to keep

keptmask = double(mask);
keptmask(linkinds) = 2;
figure; imagesc(keptmask); daspect([1,1,1]); title('Click first skeleton pixel')
[x_first,y_first] = ginput(1); 

%Downsample and calculate cross lines
%Instead of making continuous lines, put 10 lines at each segment (GCaMP
%and intensity data calculated at each segment already).

dmap = bwdistgeodesic(mask,round(x_first),round(y_first));
%Sort skeleton by increasing dmap values
linkinddistvals = dmap(linkinds);
[~,sortdist] = sort(linkinddistvals); 
sortlink = linkinds(sortdist); %link pixels put in order by distance along vessel

tmpmap = zeros(size(mask));
tmpmap(sortlink) = 1;
% sortlinkexpandmap = imdilate(tmpmap,strel('disk',3));
% sortlinkexpand = find(sortlinkexpandmap);

%Get link pixels to put crossline through. USE ONE SET OF LINES AT EACH SEG

numlines = length(skel_neighbor_ind); % #segments, Draw lines for each segment
crosslines = 10;
skelinds = linkinds; %Use the remaining links as our skeleton.
[skelrow,skelcol] = ind2sub(size(mask),skelinds);
%Get mean location & skel pixel for each segment
seg_skel = zeros(numlines,1);
for i = 1:numlines
    [current_row,current_col] = ind2sub(size(mask),skel_neighbor_ind{i});
    meanrow = mean(current_row);
    meancol = mean(current_col);

    meantoskel_dist = ((meanrow - skelrow).^2 + (meancol - skelcol).^2).^(1/2);
    [~,seg_skel(i)] = min(meantoskel_dist); %Minimum distance to each skeleton point. Gives skeleton for each segment.
end

%Delete seg_skel entries if they don't appear in sortlink
seg_skel_inds = skelinds(seg_skel);
[Lia,Locb] = ismember(seg_skel_inds,sortlink); %All seg_skel_inds should be member of sortlink.

seg_skel(~Lia) = []; %Now only inds in analyzed branch remain
seg_skel_inds(~Lia) = [];

%Put seg_skel values in order
Locb(~Lia) = [];
[~,Locb_sortinds] = sort(Locb);
seg_skel = seg_skel(Locb_sortinds); %Same order as sortlink.
seg_skel_inds = seg_skel_inds(Locb_sortinds);

%Get cross slopes for every ordered skeleton pixel. Then group by segment
xplot = 1:size(mask,2);
%Iterate over downsampled link pixels
[sortlink_crossrow,sortlink_crosscol] = ind2sub(size(mask),sortlink);
slopemean_num = 5;
midslope = nan(length(sortlink),1);
for i = 1:length(sortlink)
    currentind = sortlink(i); %Current link index to put crossline
    [current_row,current_col] = ind2sub(size(mask),currentind);
    if (i - slopemean_num) > 0 && (i + slopemean_num) < length(sortlink) %we can calculate from surrounding points
        meanrow_next = mean(sortlink_crossrow((i+1):(i+slopemean_num)));
        meanrow_prev = mean(sortlink_crossrow((i-slopemean_num):(i-1)));
        meancol_next = mean(sortlink_crosscol((i+1):(i+slopemean_num)));
        meancol_prev = mean(sortlink_crosscol((i-slopemean_num):(i-1)));
        midslope(i) = (meanrow_next-meanrow_prev)/(meancol_next-meancol_prev);
    end
end
crossslope = -1./midslope;
maskradii = maskdist(sortlink);

crosslengths_radii = 2.5

crosslengths = crosslengths_radii*maskradii;

addpath '\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis'
counter = 1;
%Now we have slope for every skeleton pixel. Get ind cell array by
%collecting groups of (crosslines) at every segment pixel.
for i = 1:length(seg_skel_inds)
    skelpixtmp = seg_skel_inds(i);
    skelpix_loc = find(sortlink == skelpixtmp); %Where in the whole skeleton pixel list this segment's pixel falls
    %Get surrounding skeleton pixels
    if skelpix_loc < length(sortlink) - (crosslines-1) && skelpix_loc > (crosslines-1) %Make sure there is enough room for additional cross lines
        skelpix_surround = sortlink((skelpix_loc - (crosslines-1)):2:(skelpix_loc + (crosslines-1))); %Get surrounding skeleton pixels (skelpix_loc is the center)
        [skelpix_surround_row,skelpix_surround_col] = ind2sub(size(mask),skelpix_surround); %(x,y) coords of surrounding skeleton pixels

        singlevecs_row = skelpix_surround_row(2:end) - skelpix_surround_row(1:end-1);
        singlevecs_col = skelpix_surround_col(2:end) - skelpix_surround_col(1:end-1);
        meanvec = [mean(singlevecs_col),mean(singlevecs_row)]; %Col, row (x,y)
        meanvecslope = meanvec(2)/meanvec(1);
        meancrossslope = -1/meanvecslope;

        crossslopes = repmat(meancrossslope,[crosslines,1]);
        crosslength = crosslengths((skelpix_loc - (crosslines-1)):2:(skelpix_loc + (crosslines-1)));
        meancrosslength = mean(crosslength);

        for j = 1:length(crosslength)
            if abs(crossslopes(j)) ~= Inf
                ycross_eq = crossslopes(j).*(xplot - skelpix_surround_col(j)) + skelpix_surround_row(j);
                [xcross_subs,ycross_subs] = bresenham(xplot(1),ycross_eq(1),xplot(end),ycross_eq(end));
            else %Infinite slope -> vertical line
                ycross_subs = 1:size(mask,1);
                xcross_subs = repmat(skelpix_surround_col(j),[1,length(ycross_subs)]);
            end
            crossdists = ((ycross_subs - skelpix_surround_row(j)).^2 + (xcross_subs - skelpix_surround_col(j)).^2).^(1/2);
            keepcross = crossdists < meancrosslength;
            %Check if any pixels are outside of image
            xcross_subs = xcross_subs(keepcross); ycross_subs = ycross_subs(keepcross);
            if any(xcross_subs > size(mask,2)) || any(ycross_subs > size(mask,1)) || any(xcross_subs < 1) || any(ycross_subs < 1)
                linetodel(counter) = 1;
                counter = counter + 1;
            else
                linetodel(counter) = 0;
                inds{counter} = sub2ind(size(mask),ycross_subs,xcross_subs);
                counter = counter + 1;
            end
        end
    end
end
inds = inds(~cellfun('isempty',inds));
howuneven = mod(length(inds),crosslines);
if howuneven ~= 0
    cellsz = cell2mat(cellfun(@length,inds,'uni',false));
    if cellsz(end-(howuneven)) ~= cellsz(end-(howuneven-1)) && howuneven ~= 0
        inds((end-(howuneven - 1)):end) = [];
    elseif cellsz(howuneven+1) ~= cellsz(howuneven) && howuneven ~= 0
        inds(1:howuneven) = [];
    end
end

%Keep only unique inds
%Make array with each row is one set of lines. Then use 
indscheck = zeros(max(cell2mat(cellfun(@length,inds,'uni',false))),length(inds)/crosslines);
group = 1;
for i = 1:length(inds)
    if mod(i,crosslines) == 1 %Just look at the first line in each group
        indscheck(1:length(inds{1,i}),group) = inds{1,i};
        group = group + 1;
    end
end
indscheck = indscheck';
[~,indsia,~] = unique(indscheck,'rows'); %Inds
indstokeep = (indsia * crosslines) - (0:(crosslines-1));
indstokeep = sort(indstokeep(:));
inds = inds(indstokeep); %Only keep unique line groups

maskplot = double(logical(rgb2gray(imread('Mask.tif'))));
for i=1:length(inds)
    maskplot(inds{i}) = 2;
end
maskplot(skelinds) = 2;
figure; imagesc(maskplot); daspect([1,1,1]);

save('seg_inds.mat','inds');

% SAVE PLOTTING INFO AND EXPERIMENT INFO

vsl_graph.link.linktodel = linktodel;

toplot.sortlink = sortlink;
toplot.seg_skel_inds = seg_skel_inds;
toplot.vsl_graph = vsl_graph;
toplot.rate = rate;
toplot.NewRate = NewRate;
% toplot.Maxsf = sf;
toplot.crosslengths_radii = crosslengths_radii;
toplot.crosslines = crosslines;
toplot.frameavg = FrameAvg;
toplot.manualoffset = noise_mean;
save('toplot.mat','toplot')


%%
% time = 1:length(tmpdiam.diam(:,1));
% time = time/(47.18/3);
% figure
% for i = 1:size(tmpdiam.diam,2)
%     plot(time,medfilt1(tmpdiam.diam(:,i),5));
%     title(num2str(i)); ylim([0 0.03]);
%     pause()
% end


%%

% for trial = 1:3
%     clearvars -except puffiter trial toplotfiles diamresults tmpdiam puff filepath inds
%     warning('off','all');
%     load(toplotfiles(puffiter(trial)).name); %Load puff or afterpuff toplot file
% %     load('tmpdiam.mat');
%     cd('Default');
%     files = dir('*00.tif');
%     cd(filepath)
% 
%     imstart = toplot.startframe;
%     imend = toplot.endframe;
%     if imstart == 0
%         imstart = 1;
%     elseif mod(imstart,2) == 0
%         imstart = imstart - 1; %Start diameter timeseries one frame before the Ca timeseries
%     end
% 
%     if imend > length(files)
%         imend = length(files);
%     end
% 
%     if mod(imend,2) == 0
%         imend = imend - 1;
%     elseif mod(imend,2) == 1;
%         imend = imend - 2;
%     end
%     wavelength = (imend - imstart)/2 + 1;
% 
%     imiter = imstart:2:imend;
%     if length(imiter) ~= wavelength
%         disp('diam length ~= wave length')
%         pause;
%     end
% 
%     %Load images
%     im_cell = cell(numel(imiter), 1);
%     counter = 1;
%     numims = imiter(end);
%     tic
%     for i = imiter
% %     for i=imiter(1:100)
%         tmp_fn = [files(i).folder,'/',files(i).name];
%         im_cell{counter} = imread(tmp_fn);
%         counter = counter + 1;
%     end
%     toc
%     im_data = cat(3, im_cell{:});
%     clearvars im_cell
%     %%
% %     numlinegroups = 50;
%     numlinegroups = length(inds) - 9; %Groups of 10, do moving average in steps of one line
%     crosslines = 10;
%     linesiter = 1:crosslines:numlinegroups*crosslines;
%     diam = zeros(size(im_data,3),numlinegroups);
%     pix_mm = 1079;
% 
%     pools = 48;
%     if trial == 1
%         parpool(pools);
%     end
% 
%     cellsizes = cellfun(@length,inds,'UniformOutput',true);
%     maxcellsize = max(cellsizes);
%     indmat = NaN(maxcellsize,length(inds));
%     for i = 1:length(inds)
%         indmat(1:length(inds{i}),i) = inds{i};
%     end
%     indmat_todel = isnan(indmat);
%     indmat(indmat_todel) = 1; %Temporarily, get the first pixel of the image, then delete later.
% 
%     tic
%     parfor j = 1:size(im_data,3) %Iterate over frames
% %     for j = 1:1
%         numlines = numlinegroups;
%         imtmp = medfilt2(im_data(:,:,j),[9,9]);
% 
%         counter = 1;
%         proftmp = nan();
%         mincellsize = nan(numlinegroups,1);
%         for line = linesiter
%             cellsizes_tmp = cellsizes(line:line+(crosslines-1));
%             mincellsize(counter) = min(cellsizes_tmp);
%             proftmp(1:mincellsize(counter),counter) = mean(double(imtmp(indmat(1:mincellsize(counter),line:line+(crosslines-1)))),2);
% %             proftmp(1:mincellsize(counter),counter) = median(double(imtmp(indmat(1:mincellsize(counter),line:line+(crosslines-1)))),2);
% %             %Try this for trials with noisy estimates
%             counter = counter + 1;
%         end
% %         proftmp = double(imtmp(indmat));
% %         proftmp(indmat_todel) = nan;
%         proftmp = medfilt1(proftmp,7,'omitnan');
%         proftmp(proftmp == 0) = nan;
% 
%         proftmp_flip = NaN(size(proftmp));
%         for i = 1:numlines
%             tmpcol = proftmp(:,i);
%             proftmp_flip(1:mincellsize(i),i) = flip(tmpcol(~isnan(tmpcol)));
%         end
%         proftmp_flip(11:end,:) = []; %Use last 10 entries to estimate minimum        
% 
%         sf = 0.5; %FW(scalefactor)M
%         Minval = median(proftmp_flip,1); %Min pixel intensity of each profile
%         Maxval = max(proftmp,[],1,'omitnan'); %Max pixel intensity of each profile
%         HMval = (Maxval - Minval)*sf + Minval;
% 
%         intpfactor = 0.01;
%         %Interpolate so that we have 100x more points than starting
%         %(each profile stars with silghtly different # of points, here mm/pixel will be constant)
%         intpsizes = (mincellsize - 1).*(1/intpfactor)+1; %Size of interpolated profiles
%         sprofintp = nan(max(intpsizes),numlines);
%         for i = 1:numlines
% %             sprofintp(1:intpsizes(i),i) =
% %             interp1(1:cellsizes(i),proftmp(1:cellsizes(i),i),1:intpfactor:cellsizes(i));
% %             %Re-run ves_1 to correct this.
%             sprofintp(1:intpsizes(i),i) = interp1(1:mincellsize(i),proftmp(1:mincellsize(i),i),1:intpfactor:mincellsize(i));
%         end
% 
%         %         intppix_pix = length(sprofintp)/length(proftmp);
%         %         intppix_pix = sum(~isnan(sprofintp(:,1)))/sum(~isnan(proftmp(:,1)));
%         intppix_mm = pix_mm/intpfactor; %Each interpolated pixel length is (intpfactor)*pixel length
%         %         intppix_mm = pix_mm * intppix_pix; %Interpolated pix per mm
% 
%         tol = 1;
%         %         sprofintp_todel = isnan(sprofintp);
%         HMval_repmat = repmat(HMval,[size(sprofintp,1),1]);
%         prof_diff = abs(sprofintp - HMval_repmat);
%         v = zeros(1,numlines);
%         for i = 1:numlines
%             small_diffs = find(prof_diff(:,i) < tol); %indices of profile distance to half-max less than tolerance
% 
%             [~,groupsep] = max(diff(small_diffs)); %last ind of group 1
% 
%             [~,group1ID] = min(abs(sprofintp(small_diffs(1:groupsep)) - HMval(i))); %Ind of closest group 1 point, within group 1
%             group1ID = group1ID + small_diffs(1) - 1; %Ind of closest group 1 point, within total profile
%             [~,group2ID] = min(abs(sprofintp(small_diffs(groupsep+1:end)) - HMval(i))); %Ind of closest group 2 point, within group 2
%             group2ID = group2ID + small_diffs(groupsep + 1) - 1; %Ind of closest group 2 point, within total profile
% 
%            v(i) = (group2ID - group1ID)/intppix_mm;
%         end
%         diam(j,:) = v;
%     end
%     toc
%     tmpdiam(trial).diam = diam; %Here beforepuff is col 1 puff is 2 afterpuff is 3
% end
% % toc
% disp('FWHM Done');
% save('tmpdiam_multlines_mprof_50lines.mat','tmpdiam')










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



