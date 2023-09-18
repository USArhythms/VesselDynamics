%PlotVesNeuSegments.m

function [Mask1Crop,Mask1tmpvesCrop,NMask1Crop,NMask1tmpvesCrop] = PlotVesNeuSegments(analyzed_folder,data_folder,plotneu,link_results_struct_neu,link_results_struct)

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

cd(data_folder);
toplot = dir('*beforepuff_toplot.mat'); load(toplot.name);
ntoplot = dir('*beforepuff_neurons.mat'); 
Mask1 = imread('Mask1.tif');

f1 = figure; imshow(Mask1); title('Click Rectangle for ROI');
%Draw rectangle for roi
[col,row] = ginput(2);
pos = [col(1),row(1),diff(col),diff(row)];
hold on; rectangle('Position',pos,'EdgeColor','r');

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

ves = foundskel_lab(2);
dtmp = vesresults(1).dist_mat(:,ves);
dtmp = dtmp(~isnan(dtmp));
pix_mm = dtmp(end)/vesresults(1).link_lengths_mm(ves);
linesz = 1*pix_mm-1; %1mm scale bar
lastcol = toplot.mask_size(2);
lastrow = toplot.mask_size(1);
line([lastcol-linesz lastcol],[lastrow lastrow],'Color','white'); %Can move scalebar in illustrator

%Plot vessel segments
Mask1tmpves = repmat(Mask1,[1,1,3]);
%Get all skeleton labels and inds.
for i = 1:length(vesresults)
    skels(i) = vesresults(i).skel_lab(ves);
    if skels(i) ~= 0
        skelinds = ismember(toplot.skel_label,skels(i));
        maskinds{i} = toplot.mask_ind(skelinds);
    end
end

%Order maskinds in distance from some endpoint.
% figure; imagesc(toplot.mask);
% [col,row] = ginput(1);
startrow = link_results_struct(1).locs(ves,1);
startcol = link_results_struct(1).locs(ves,1);
geomask = bwdistgeodesic(toplot.mask,round(startcol),round(startrow));
figure; imagesc(geomask)
%Get mean locations for each maskind group
disttmp = zeros(length(maskinds),1);
for i = 1:length(maskinds)
    [rowtmp,coltmp] = ind2sub(toplot.mask_size,maskinds{i});
    meanrow = round(mean(rowtmp));
    meancol = round(mean(coltmp));
    disttmp(i) = geomask(meanrow,meancol);
end
[~,distI] = sort(disttmp);
%Sort maskinds
maskinds = maskinds(distI);

%Draw rectangle for roi
% f1 = figure; imshow(Mask1);
% [col,row] = ginput(2);
% pos = [col(1),row(1),diff(col),diff(row)];
% hold on; rectangle('Position',pos,'EdgeColor','r');
%%
vesmask = zeros(size(Mask1tmpves,1),size(Mask1tmpves,2));
numcolors = 2;
linescmap = lines(numcolors);
% linescmap = prism(numcolors);
% linescmap = (1/255*[24,69,59;123,189,0;255,255,255]); %Series of green
linescmap = (1/255*[0,0,0;255,255,255]);
for i = 1:length(maskinds)
    colornum = mod(i,numcolors)+1; %Need to order maskinds in distance from starting point
    colortmp = round(linescmap(colornum,:)*256);
%     colortmp = linescmap(colornum,:);
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
%%
figure
imshow(Mask1Crop);
hold on;
h = imshow(Mask1tmpvesCrop);
% ves_alphaCrop(ves_alphaCrop < 0) = 0;
% set(h,'AlphaData',ves_alphaCrop) %This is for making image transparent in
% matlab. If doing in adobe then keep as not transparent in matlab.
line([167 189],[205 205],'Color','white'); %Scale bar
% cd('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\AllSegments\ShortDistAnalysis\Supp8');
% print(gcf, '-depsc2','-cmyk', 'ExampleMask1_TB013120M4_20Mar170741_ves91_Segs');

if plotneu == 1
    load(ntoplot.name);

    Mask1tmpves = repmat(Mask1,[1,1,3]);
    neuinds = neuresults(1).inds(:,ves);
    numsegs = sum(~cellfun(@isempty,neuinds),1);
    neuinds(logical(cellfun(@isempty,neuinds))) = [];

    %Order maskinds in distance from some endpoint.
    map = zeros(toplot.mask_size);
    map(round(link_results_struct_neu(1).skel_locs(ves,1)),link_results_struct_neu(1).skel_locs(ves,2)) = 1;
    distmask = bwdist(map); 
%     figure;
%     imagesc(distmask)
%     Get mean locations for each maskind group
    disttmp = zeros(numsegs,1);
    map = toplot.mask;
    for i = 1:numsegs
%         [rowtmp,coltmp] = ind2sub(toplot.mask_size,neuinds{i,1});
%         meanrow = round(mean(rowtmp));
%         meancol = round(mean(coltmp));

        %Get corresponding location (closest skeleton segment pixel)
        disttmp(i) = distmask(link_results_struct_neu(i).skel_locs(ves,1),link_results_struct_neu(i).skel_locs(ves,2));
%         disttmp(i) = distmask(meanrow,meancol);

        neulocind = sub2ind(toplot.mask_size,link_results_struct_neu(i).skel_locs(ves,1),link_results_struct_neu(i).skel_locs(ves,2));
        map(neulocind) = i;
    end
%     figure; imagesc(map);

    [newdist,distI] = sort(disttmp);
    %Sort maskinds
    neuinds = neuinds(distI);

    neumask = zeros(size(Mask1tmpves,1),size(Mask1tmpves,2));
    numcolors = numsegs;
%     linescmap = lines(numcolors);
    % linescmap = prism(numcolors);
%     linescmap = colormap("winter");
    cmap_lims = [24,69,59;123,189,0]+1;
    cmap_intp = interp1([1;numsegs],cmap_lims,1:numsegs');
    linescmap = cmap_intp;
%     linescmap = (1/255*[24,69,59;123,189,0;255,255,255]); %Series of green
    
%     linescmap = (1/255*[0,0,0;255,255,255]);
%     coloriter = round((1:numsegs)/numsegs*256);
%     figure
    for i = 1:numsegs
%         colornum = mod(i,numcolors)+1;

%         colortmp = round(linescmap(colornum,:)*256);
        colortmp = round(linescmap(i,:));
%         colortmp = linescmap(colornum,:);
        for j = 1:length(neuinds{i,1})
            [tmpvescol,tmpvesrow] = ind2sub([size(Mask1,1),size(Mask1,2)],neuinds{i,1}(j));
            Mask1tmpves(tmpvescol,tmpvesrow,:) = colortmp;
            %         alphamask(tmpvescol,tmpvesrow) = alpha;
            neumask(tmpvescol,tmpvesrow) = 1;
        end
%         imshow(Mask1tmpves)
%         pause();
    end

    %     figure
    %     imshow(Mask1tmpves)

    NMask1Crop = imcrop(Mask1,pos);
    NMask1tmpvesCrop = imcrop(Mask1tmpves,pos);

    figure
    imshow(NMask1Crop);
    hold on;
    h = imshow(NMask1tmpvesCrop);
    line([167 189],[205 205],'Color','white'); %Scale bar

end


%Plot phase progressions
lrs = link_results_struct;
lrsn = link_results_struct_neu;
dist = lrs(1).dist_mat(:,ves); dist = dist(~isnan(dist));
phase = lrs(1).phi_mat(:,ves); phase = phase(~isnan(phase));
pix_mm = dist(end)/lrs(1).link_lengths_mm(ves);
dist = dist/pix_mm;
ndist = lrsn(1).n_dist_mat(:,ves); ndist = ndist(~isnan(ndist));
nphase = lrsn(1).n_phi_mat(:,ves); nphase = nphase(~isnan(nphase));
ndist = ndist/pix_mm;
figure; scatter(dist,phase,'filled','b');
hold on; scatter(ndist,nphase,36,linescmap/256,'filled');
title(num2str(ves));



end


% map = zeros(toplot.mask_size);
% figure
% for i = 1:92
%     map(neuinds{i,1}) = i;
%     imagesc(map);
%     title(num2str(i));
% %     pause()
% end



