%Use cross-correlation to detect periods of movement in 2P image stack
% framexcorr_fft_batch.m
% for running on muralis

clear; clc; close all;


animal = 'JD230306F3'
%Iterate over rois for each animal

for ii = 1:9
    clearvars -except ii animal

    roiname = ['ROI',num2str(ii)]
    % roiname = 'ROI9'
    filepath = ['/net/birdstore/Thomas_1P/JD_GC_Diam2P/6_13_23/8_13_23Analysis/',animal,'/',roiname];
    cd(filepath)
    origmask = imread('Mask_origsize.tif');
    mask = imread('Mask.tif');
    %Load images
    tic
    d_data = h5read([animal,'_',roiname,'_diamdat.h5'],'/frame');
    d_data = reshape(d_data,[size(origmask,1),size(origmask,2),size(d_data,2)]); %X x Y x Time
    toc;
    disp(['Done Loading ',roiname])
    %%

    meanim = mean(d_data,3);
    meanim = meanim - mean(meanim(:));

    pad = 0;
    Nx = size(meanim,2);
    Ny = size(meanim,1);
    nfftx=max(2^(nextpow2(Nx)+pad),Nx);
    nffty=max(2^(nextpow2(Ny)+pad),Ny);
    paddedmeanim = zeros(nffty,nfftx);
    paddedmeanim(1:Ny,1:Nx) = meanim;
    Jmeanim = fft2(paddedmeanim);

    corr1 = fftshift(real(ifft2((Jmeanim).*conj(Jmeanim))));
    %Do 2d interpolation
    origin = nfftx/2 + 1;
    cropsz = 5;
    corr1crop = corr1((origin-cropsz):(origin+cropsz),(origin-cropsz):(origin+cropsz));
    [X,Y] = meshgrid(-cropsz:cropsz);
    [Xq,Yq] = meshgrid(-cropsz:0.1:cropsz);
    corr1intp = interp2(X,Y,corr1crop,Xq,Yq,'cubic');

    [~,max1] = max(corr1intp,[],'all');
    % [locy1,locx1] = ind2sub(size(paddedmeanim),max1);
    [locy1,locx1] = ind2sub(size(corr1intp),max1);

    pools = 48;
    parpool(pools);

    tic
    parfor i = 1:size(d_data,3)

        % tic
        testim = d_data(:,:,i);
        testim = testim - mean(testim(:));
        paddedtestim = zeros(nffty,nfftx);
        paddedtestim(1:Ny,1:Nx) = testim;
        Jtestim = fft2(paddedtestim);
        corr = fftshift(real(ifft2((Jmeanim).*conj(Jtestim))));

        %Get sub-pixel maximum by doing 2d interpolation
        corrcrop = corr((origin-cropsz):(origin+cropsz),(origin-cropsz):(origin+cropsz));
        corrintp = interp2(X,Y,corrcrop,Xq,Yq,'cubic');

        [~,maxtest] = max(corrintp,[],'all');
        [locytest,locxtest] = ind2sub(size(corr1intp),maxtest);
        frameshift(i) = sqrt((locy1-locytest)^2 + (locx1-locxtest)^2);
        zeroshift(i) = corr(locy1,locx1);

        % toc
        % figure; subplot(2,1,1); imagesc(meanim); daspect([1,4,1]);
        % subplot(2,1,2); imagesc(testim); daspect([1,4,1]);

        %pad each image to next power of 2
        % Nx = size(meanim,2);
        % Ny = size(meanim,1);
        % nfftx=max(2^(nextpow2(Nx)),Nx);
        % nffty=max(2^(nextpow2(Ny)),Ny);
        % paddedmeanim = zeros(nffty,nfftx);
        % paddedmeanim(1:Ny,1:Nx) = meanim;
        % paddedtestim(1:Ny,1:Nx) = testim(end:-1:1,end:-1:1);

    end
    toc

    xcorr2dstruct = struct();
    xcorrstruct(1).frameshift = frameshift;
    xcorrstruct(1).zeroshift = zeroshift;
    save(['xcorrstruct_',roiname,'.mat'],'xcorrstruct')

end
