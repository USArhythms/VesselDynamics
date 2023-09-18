% PenCorrelationLength.m 2_25 updates

%Define range for vasomotion frequency, get vasomotion freq range based on
%dB drop in power. 
%Use this and uncertainties to get multiple estimates of delay direction
%that are more repeatable, especially for cases of weak/inconsistent
%vasomotion (no peak just shoulder)

%Quantify the correlation length for the penetrating arteriole data set
%x-axis = pairwise distance
%y-axis = +/- 1 (same direction or opposite direction) OR phi1 * phi2
%(multiplied phases, + if same direction, - if opposite direction)

%Also plot histogram of phase (positive and negative bump) put veins in
%there too

%Also plot simple histogram of up and down


% Import data from each trial
clear; clc; close all;
% analysisfolder = '\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20230302 WT_PA_11' %Has (x,y) coords
% analysisfolder = '\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20230301 WT_PA_11' %Has (x,y) coords
% analysisfolder = '\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20230224 WT_PA_10' %Has (x,y) coords
% analysisfolder = '\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20230222 WT_PA_10' %Has (x,y) coords
% analysisfolder = '\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20230215 WT_PA_9' %Has (x,y) coords
% analysisfolder = '\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20230214 WT_PA_9' %Has (x,y) coords
% analysisfolder = '\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20230213 WT_PA_8' %Has (x,y) coords
% analysisfolder = '\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20230210 WT_PA_8' %Has (x,y) coords
% analysisfolder = '\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20230125 WT_PA_7' %to do, has 2 penA measurements
% analysisfolder = '\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20230118 WT_PA_7';
% analysisfolder = '\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20221212 WT_PA_7';
% analysisfolder = '\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20221129 WT_PA_6';
% analysisfolder = '\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20221120 WT_PA_5';
% analysisfolder = '\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20221128 WT_PA_5';
analysisfolder = '\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\20221208 WT_PA_5';
cd(analysisfolder);
wdfld = imread('WidefieldView_NoLabels.tif');
wdfld = imread('WidefieldView.tif');
figure; imagesc(wdfld); daspect([1,1,1]);

% data_folder = '\\dk-server.dk.ucsd.edu\rul087\Data backup\20230202 WT_PA_11';
% data_folder = '\\dk-server.dk.ucsd.edu\rul087\Data backup\20230301 WT_PA_11';
% data_folder = '\\dk-server.dk.ucsd.edu\rul087\Data backup\20230224 WT_PA_10';
% data_folder = '\\dk-server.dk.ucsd.edu\rul087\Data backup\20230222 WT_PA_10';
% data_folder = '\\dk-server.dk.ucsd.edu\rul087\Data backup\20230215 WT_PA_9';
% data_folder = '\\dk-server.dk.ucsd.edu\rul087\Data backup\20230214 WT_PA_9';
% data_folder = '\\dk-server.dk.ucsd.edu\rul087\Data backup\20230213 WT_PA_8';
% data_folder = '\\dk-server.dk.ucsd.edu\rul087\Data backup\20230210 WT_PA_8';
% data_folder = '\\dk-server.dk.ucsd.edu\rul087\Data backup\20230125 WT_PA_7 stack';
% data_folder = '\\dk-server.dk.ucsd.edu\rul087\Data backup\20230118 WT_PA_7';
% data_folder = '\\dk-server.dk.ucsd.edu\rul087\Data backup\20221212 WT_PA_7';
% data_folder = '\\dk-server.dk.ucsd.edu\rul087\Data backup\20221129 WT_PA_6';
% data_folder = '\\dk-server.dk.ucsd.edu\rul087\Data backup\20221120 WT_PA_5';
% data_folder = '\\dk-server.dk.ucsd.edu\rul087\Data backup\20221128 WT_PA_5';
data_folder = '\\dk-server.dk.ucsd.edu\rul087\Data backup\20221208 WT_PA_5';

save_folder = '\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Pipeline\8_4_StimAnalysis';

animal = '20221208PA5'; %For file naming

cd(data_folder);
files = dir('*001.tif');
files([files.bytes] < 9e8) = []; %Only get the large files
files([files.bytes] > 1e10) = []; %Delete large files (not image sequence)
todel = logical(zeros(length(files),1));
for i=1:length(files)
    if contains(files(i).name,'test') || contains(files(i).name,'PV') || contains(files(i).name,'Pia') %|| contains(files(i).name,'stim') %|| contains(files(i).name,'20230210')
        todel(i) = true;
    end
end
files(todel) = [];

%%
stim = 'stim';
for file=1:length(files)

    %     if contains(animal,'PA7') %For this animal
    %         if file == 12
    %             namestr = extractBefore(files(file).name,'_6');
    %         elseif file == 8
    %             namestr = extractBefore(files(file).name,'_7');
    %         else
    %             namestr = extractBefore(files(file).name,'_7.6');
    %         end
    %     end
    if contains(files(file).name,'6.0x')
        namestr = extractBefore(files(file).name,'_6.0x');
    elseif contains(files(file).name,'7.6x')
        namestr = extractBefore(files(file).name,'_7.6');
    elseif contains(files(file).name,'_6x')
        namestr = extractBefore(files(file).name,'_6x');
    elseif contains(files(file).name,'_6.0x')
        namestr = extractBefore(files(file).name,'_6.0x');
    elseif contains(files(file).name,'6.4x')
        namestr = extractBefore(files(file).name,'_6.4x');
    elseif contains(files(file).name,'3.0x')
        namestr = extractBefore(files(file).name,'_3.0x');
    elseif contains(files(file).name,'5.0x')
        namestr = extractBefore(files(file).name,'_5.0x');
    elseif contains(files(file).name,'5.5x')
        namestr = extractBefore(files(file).name,'_5.5x');
    elseif contains(files(file).name,'5.6x')
        namestr = extractBefore(files(file).name,'_5.6x');
    elseif contains(files(file).name,'5.8x')
        namestr = extractBefore(files(file).name,'_5.8x');
    elseif contains(files(file).name,'5.7x')
        namestr = extractBefore(files(file).name,'_5.7x');
    elseif contains(files(file).name,'7.0x')
        namestr = extractBefore(files(file).name,'_7.0x');
    elseif contains(files(file).name,'7x')
        namestr = extractBefore(files(file).name,'_7x');
    elseif contains(files(file).name,'6.5x')
        namestr = extractBefore(files(file).name,'_6.5x');
    elseif contains(files(file).name,'4.5x')
        namestr = extractBefore(files(file).name,'_4.5x');
    else
        namestr = extractBefore(files(file).name,'_4.8x');
    end
    %     depth_str = extractAfter(namestr,'_');
    %     depth1 = extractBefore(depth_str,'_')
    %     depth2 = extractAfter(depth_str,'_')
    %     PA = extractBefore(namestr,'_')
    %     zoom = extractBetween(files(file).name,[depth2,'_'],'x_00001.tif')
    %     zoom = str2double(zoom{1})
%     if length(namestr) == 10
%         PA = extractBefore(namestr,5)
%         PA = extractAfter(PA,2)
%     elseif strcmp(animal,'20230125PA7')
%         PA = extractBefore(namestr,6)
%         PA = extractAfter(PA,3)
%     else
%         PA = extractBefore(namestr,4)
%         PA = extractAfter(PA,2)
%     end
    PA = extractBefore(namestr,'_')
    PA = extractAfter(PA,2)
    
    depth_str = extractAfter(namestr,[PA,'_']);
    depth1 = extractBefore(depth_str,'_')
    depth2 = extractAfter(depth_str,'_')

%     zoom = extractBetween(files(file).name,[depth2,'_'],'x_00001.tif')
    zoom = extractBetween(files(file).name,[depth2,'_'],'x_stim_00001.tif')
%     zoom = extractBetween(files(file).name,[depth2,'_'],'x_stim_rest_00001.tif')
    zoom = str2double(zoom{1})

    pix_um = 300/100;
    pix_um = pix_um / 7.6 * zoom

    cd(analysisfolder)
    resultfiles = dir('*allstats.mat');
    todel = zeros(length(resultfiles),1);
    for i=1:length(resultfiles)
%         if contains(resultfiles(i).name,['_',PA,'_'])
        if contains(resultfiles(i).name,['_PA',PA,'_'])
            todel(i) = 1;
        end
    end
    resultfiles(~logical(todel)) = [];


    isrest = zeros(length(resultfiles),1);
    isstim = zeros(length(resultfiles),1);
    for i=1:length(resultfiles)
%         if contains(resultfiles(i).name,'rest')
%             isrest(i) = 1;
%         end
        if contains(resultfiles(i).name,'stim')
            isstim(i) = 1;
        end
    end
%     if strcmp(stim,'rest')
%         resultfiles(~isrest) = [];
%     elseif strcmp(stim,'stim')
%         resultfiles(isrest) = [];
%     end
    resultfilestmp = resultfiles;
    resultfiles(~logical(isstim)) = [];
    resultfilestmp(logical(isstim)) = []; %Temporary, to get vessel locations from previous analysis.

    disp(resultfiles(1).name);
    disp(resultfiles(2).name);

    clearvars -except data_folder animal file files namestr depth_str depth1 depth2 PA stim pix_um rate resultfiles resultfilestmp analysisfolder phi_mat phi_struct oldfile oldpath oldphi_struct save_folder

    stats1 = load(resultfiles(1).name);
    stats2 = load(resultfiles(2).name);
    stats1 = stats1.allstats;
    stats2 = stats2.allstats;
    %     if ~contains(resultfiles(1).name,'_0um')
    if str2num(stats1.depth) > str2num(stats2.depth)
        stats1 = load(resultfiles(2).name);
        stats2 = load(resultfiles(1).name);
        stats1 = stats1.allstats;
        stats2 = stats2.allstats;
    end
    stats1tmp = load(resultfilestmp(1).name);
    stats2tmp = load(resultfilestmp(2).name);
    stats1tmp = stats1tmp.allstats;
    stats2tmp = stats2tmp.allstats;
    if str2num(stats1tmp.depth) > str2num(stats2tmp.depth)
        stats1tmp = load(resultfilestmp(2).name);
        stats2tmp = load(resultfilestmp(1).name);
        stats1tmp = stats1tmp.allstats;
        stats2tmp = stats2tmp.allstats;
    end
    stats1.vescenterx = stats1tmp.vescenterx;
    stats1.vescentery = stats1tmp.vescentery;
    stats2.vescenterx = stats2tmp.vescenterx;
    stats2.vescentery = stats2tmp.vescentery;

    %Calculate Phase
    try
        rate = stats1.rate
    catch
        rate = 7.25
    end
    waves = stats1.RadondEq_Outl;
    waves = waves - mean(waves);
    waved = stats2.RadondEq_Outl;
    waved = waved - mean(waved);
    times = stats1.time;
    timed = stats2.time;
    if sum(times - timed) == 0
        timed = timed + (1/rate)/2;
    end

    %Interpolate each signal to get equal time points
    tqs = times(1):(1/(2*rate)):times(end);
    tqd = timed(1):(1/(2*rate)):timed(end);
    waves_q_tmp = interp1(times,waves,tqs);
    waved_q_tmp = interp1(timed,waved,tqd);

    %Delete first wave1 value so that both start at t=0.75.
    waves_q = waves_q_tmp(2:end);
    tqs = tqs(2:end);
    %Delete last wave1 value to make both time series the same duration
    waved_q = waved_q_tmp;
    waved_q(end) = [];
    tqd(end) = [];

    params.Fs = rate*2; %Interpolated rate is twice actual single depth rate
    params.pad = 2;
    params.fpass = [0 params.Fs/4]; %Hz, default is [0 Fs/2]
    params.err   = [2 .05];
    params.trialave = 0;
    T = stats1.time(end);
    BW = 0.02; %600s trial -> TBW =~ 12
    params.tapers = [round(T*BW),round(2*T*BW-1)]; %Time-BW product and number of tapers
    addpath(genpath('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\chronux_2_12'))

    [S,f] = mtspectrumc(waves_q,params);
    if strcmp(stim,'rest')
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
        %         addpath(genpath('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis'))
        %         [manualf,windf,fvexpfit,fvextremafit,flinesc,manualfvec,windfvec,expfitfvec,extremafitfvec,linescfvec] = extractfv(f,S,params,waves,rate); %Use mtspectrum to extract fv repeatably
    else
        f_peak = 0.1; %Hz
        freq_findx = max(find(round(f,3)==round(f_peak,3)));
        addpath(genpath('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\chronux_2_12'))
    end
    %Calcuate phase between deep and shallow, using interpolated data
    params_err = params;
    params_err.err = [2,0.05];
    [C,phi,S12,S1,S2,f,confC,phistd,Cerr] = coherencyc(waved_q,waves_q,params_err);
    dof = 2*round(2*T*BW-1); %Single value dof (2*tapers)
    df = 1./((dof/2)-1);
    p = 0.05; %Confidence level
    confC = sqrt(1 - p.^df);

    %     manualf_findx = max(find(round(f,3)==round(manualf,3)))
    %     windf_findx = max(find(round(f,3)==round(windf,3)))
    %     fvexpfit_findx = max(find(round(f,3)==round(fvexpfit,3)))
    %     fvextremafit_findx = max(find(round(f,3)==round(fvextremafit,3)))
    %     flinesc_findx = max(find(round(f,3)==round(flinesc,3)))

    findx_05 = max(find(round(f,4)==0.05));
    findx_1 =  max(find(round(f,4)==0.1));
    findx_15 = max(find(round(f,4)==0.15));

    figure
    plot(f,phi,'b'); hold on; plot(f,phi + 2*phistd','Color',[0,0,1,0.2]); plot(f,phi - 2*phistd','Color',[0,0,1,0.2]); xlim([0 0.5]); ylim([-pi pi]);
    yline(0,'Alpha',0.9); xlabel('Frequency (Hz)','Interpreter','latex'); ylabel('Relative Phase (rad)','Interpreter','latex');
    xline(f_peak,'-.','$f_v$','Interpreter','latex','Alpha',0.9,'Color','k');
    xline(0.05,'-.','$f_v$','Interpreter','latex','Alpha',0.5,'Color','r');
    xline(0.1,'-.','$f_v$','Interpreter','latex','Alpha',0.5,'Color','r');
    xline(0.15,'-.','$f_v$','Interpreter','latex','Alpha',0.5,'Color','r');

    %     xline(windf,'-.','$f_v$','Interpreter','latex','Alpha',0.9,'Color','r');
    %     xline(fvexpfit,'-.','$f_v$','Interpreter','latex','Alpha',0.9,'Color','b');
    %     xline(fvextremafit,'-.','$f_v$','Interpreter','latex','Alpha',0.9,'Color','m');
    %     xline(flinesc,'-.','$f_v$','Interpreter','latex','Alpha',0.9,'Color','y');

    %Calculate phase gradient
    dz = (str2num(depth2) - str2num(depth1))/1000; %mm depth difference (mm)
    try
        dx = sqrt((stats1.vescentery - stats2.vescentery)^2 + (stats1.vescenterx - stats2.vescenterx)^2); %pixels, horizontal distance between vessel segments
        dx = dx/pix_um/1000; %mm
    catch
        %For animals where this wasn't recorded originally all except WT10
        cd(data_folder)
        imfiles = dir('*0001.tif')
        todel = zeros(length(imfiles),1);
        for i=1:length(imfiles)
            if ~contains(imfiles(i).name,['PA_',PA,'_'])
                todel(i) = 1;
            end
        end
        imfiles(logical(todel)) = [];
        im1 = imread([imfiles(1).folder,'\',imfiles(1).name],3);
        im2 = imread([imfiles(1).folder,'\',imfiles(1).name],4);
        figure; imagesc(im1); daspect([1,1,1]);
        [stats1.vescenterx,stats1.vescentery] = ginput(1);
        cd(analysisfolder)
        allstats = stats1;
        %         save([animal,'_PA',PA,'_',stim,'_',depth1,'um','_allstats.mat'],'allstats');
        save([animal,'_',PA,'_rest_',depth1,'um','_allstats.mat'],'allstats');
        figure; imagesc(im2); daspect([1,1,1]);
        [stats2.vescenterx,stats2.vescentery] = ginput(1);
        allstats = stats2;
        save([animal,'_',PA,'_rest_',depth2,'um','_allstats.mat'],'allstats');
        dx = sqrt((stats1.vescentery - stats2.vescentery)^2 + (stats1.vescenterx - stats2.vescenterx)^2); %pixels, horizontal distance between vessel segments
        dx = dx/pix_um/1000; %mm
    end
    ds = sqrt(dz^2 + dx^2); %mm, Lower estimate of vessel distance between two planes (not completely straight so still underestimating)

    phi_manualf = phi(freq_findx);
    k_manualf = phi_manualf / ds; %rad/mm
    v_manualf = (2*pi*f_peak)/k_manualf; %velocity mm/s, assuming no dispersion
    phi_05 = phi(findx_05);
    k_05 = phi_05 / ds; %rad/mm
    v_05 = (2*pi*0.05)/k_05; %velocity mm/s, assuming no dispersion
    phi_1 = phi(findx_1);
    k_1 = phi_1 / ds; %rad/mm
    v_1 = (2*pi*0.1)/k_1; %velocity mm/s, assuming no dispersion
    phi_15 = phi(findx_15);
    k_15 = phi_15 / ds; %rad/mm
    v_15 = (2*pi*0.15)/k_15; %velocity mm/s, assuming no dispersion

    str1 = sprintf('%.2f',phi_manualf);
    str2 = sprintf('%.2f',abs(k_manualf));
    str3 = sprintf('%.2f',v_manualf);
    title({'Relative Phase and 95 Percent CI',['$\phi_{fv}$ = ',str1,', $\mid k \mid$ = ',str2,' $v$ = ',str3]},'Interpreter','latex')

    if file == 1
        phi_struct = struct();
    end
    phi_struct(file).PA = PA;
    phi_vec_tmp = [phi_manualf;phi_05;phi_1;phi_15];
    k_vec_tmp = [k_manualf;k_05;k_1;k_15];
    phi_uncert_tmp = [phistd(freq_findx);phistd(findx_05);phistd(findx_1);phistd(findx_15)];
    fv_tmp = [f_peak;0.05;0.1;0.15];

    phi_mat = [phi_vec_tmp,k_vec_tmp,phi_uncert_tmp,fv_tmp];

    %         phi_mat(file,1) = phi_fv;
    %         phi_mat(file,2) = k_fv; %rad/mm
    %         phi_mat(file,3) = phistd(freq_findx);
    %         phi_mat(file,4) = phistd(freq_findx)/abs(phi_fv)*abs(k_fv); %k uncertainty
    %         phi_mat(file,5) = f_peak; %peak vasofreq (Hz)
    %         phi_mat(file,6) = (2*pi)/k_fv*f_peak; %Velocity (mm/s) assuming non dispersive propagation (v = w/k) vs. v = dw/dk
    %         phi_mat(file,7) = str2num(extractAfter(PA,'PA'))
    %         phi_mat(file,7) = str2num(PA)

    %         wdfld = imread('WidefieldLabels.tif');

    if contains(animal,'PA8') || contains(animal,'PA9') || contains(animal,'PA10') || contains(animal,'PA11')
        disp(['PA ',PA])
        prompt = "Enter recorded x location ";
        %             phi_mat(file,8) = input(prompt); %um relative to PA1
        phi_struct(file).xloc = input(prompt);
        prompt = "Enter recorded y location ";
        %             phi_mat(file,9) = input(prompt); %um relative to PA1
        phi_struct(file).yloc = input(prompt);
    else %Distances have been calculated previously, load these results and assign them appropriately here.
        if file == 1
            disp('Choose old phi_struct')
            disp(animal)
            [oldfile,oldpath] = uigetfile
        end
%     elseif contains(animal,'PA7') && contains(animal,'20230118')
%       
        oldphi_mat = load([oldpath,'\',oldfile]) %Here I already found distnces using the widefield image
        oldphi_struct = oldphi_mat.phi_struct;
        todel = zeros(length(oldphi_struct),1);
        for j = 1:length(oldphi_struct)
            if ~strcmp(PA,oldphi_struct(j).PA)
                todel(j) = 1;
            end
        end
        oldphi_struct(logical(todel)) = [];
        phi_struct(file).xloc = oldphi_struct(1).xloc;
        phi_struct(file).yloc = oldphi_struct(1).yloc;
%     else
%         figure; imagesc(wdfld); daspect([1,1,1])
%         [x,y] = ginput(1);
%         phi_struct(file).xloc = x;
%         phi_struct(file).yloc = y;
    end
    phi_struct(file).phi_mat = phi_mat;
    phi_struct(file).f = f;
    phi_struct(file).C = C;
    phi_struct(file).Cerr = Cerr;
    phi_struct(file).confC = confC;

    phi_struct(file).Ss_Sd = S2(freq_findx)/S1(freq_findx);
    phi_struct(file).Ss = S2(freq_findx);
    phi_struct(file).Sd = S1(freq_findx);
    phi_struct(file).confC = confC;

end

cd('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Pipeline\8_4_StimAnalysis');
save(['phi_struct_stim_',animal,'.mat'],'phi_struct');

%% Make CorrMat
cd(save_folder);
if strcmp(animal,'20230224PA10')
    load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Pipeline\7_8_23_Analysis\phi_struct_20230224PA10.mat");
    phi_struct(4) = []; %Keep PA18 measurement with better S/N ratio. They both had waves moving inward this day so for the corrmat it doesn't matter which to delete.
end

numpairs = nchoosek(length(phi_struct),2)
CorrMat = zeros(numpairs,12);
counter = 1;
for i=1:length(phi_struct)-1
    v1 = 2*pi/(phi_struct(i).phi_mat(1,2))*phi_struct(i).phi_mat(1,4); %Using v = 2*pi*f/k, negative bottom to top, positive top to bottom
    Pnum1 = str2num(phi_struct(i).PA);
    Pos1x = phi_struct(i).xloc;
    Pos1y = phi_struct(i).yloc;
    phi1 = phi_struct(i).phi_mat(1,1); %phi
    dphi1 = phi_struct(i).phi_mat(1,3); %std (phi)

    freq_findx = max(find(round(phi_struct(i).f,3)==round(phi_struct(i).phi_mat(1,4),3)));
%     Cerrfv1 = phi_struct(i).Cerr(1,freq_findx);
    Cerrfv1 = phi_struct(i).C(freq_findx); %Get coherence value at vasomotion frequency
    if isfield(phi_struct,'confC')
        confC1 = phi_struct(i).confC;
    else
        confC1 = 0.3568; %Same time length, tapers, confidence level used always 
    end

    for j=(i+1):length(phi_struct)
        v2 = 2*pi/(phi_struct(j).phi_mat(1,2))*phi_struct(j).phi_mat(1,4);
        Pnum2 = str2num(phi_struct(j).PA);
        Pos2x = phi_struct(j).xloc;
        Pos2y = phi_struct(j).yloc;

        CorrMat(counter,1) = Pnum1;
        CorrMat(counter,2) = Pnum2;
        CorrMat(counter,3) = v1;
        CorrMat(counter,4) = v2;
        %         CorrMat(counter,3) = sign(v1*v2);
        %         CorrMat(counter,4) = v1*v2;


        if contains(animal,'PA8') || contains(animal,'PA9') || contains(animal,'PA10') || contains(animal,'PA11')
            CorrMat(counter,5) = sqrt((Pos2x-Pos1x)^2 + (Pos2y-Pos1y)^2) / 1000; %Pair Wise distance (mm) Use this for WT8-12 where locations are recorded
        else
            cd(analysisfolder)
            load('wdf_pix_mm_struct.mat');
            Wdf_pix_mm = wdf_pix_mm_struct.value;
            CorrMat(counter,5) = sqrt((Pos2x-Pos1x)^2 + (Pos2y-Pos1y)^2) / Wdf_pix_mm; %Pair Wise distance (mm)
        end

        CorrMat(counter,6) = phi1; %phase 1
        CorrMat(counter,7) = dphi1; %phase 1 standard dev
        CorrMat(counter,8) = phi_struct(j).phi_mat(1,1); %phase 2
        CorrMat(counter,9) = phi_struct(j).phi_mat(1,3); %phase 2 standard dev

        CorrMat(counter,10) = phi_struct(j).Ss; %Shallow power at vasomotor frequency
        CorrMat(counter,11) = phi_struct(j).Sd; %Deep power at vasomotor frequency

        freq_findx = max(find(round(phi_struct(j).f,3)==round(phi_struct(j).phi_mat(1,4),3)));
%         Cerrfv2 = phi_struct(j).Cerr(1,freq_findx);
        Cerrfv2 = phi_struct(j).C(freq_findx);

        if isfield(phi_struct,'confC')
            confC2 = phi_struct(i).confC;
        else
            confC2 = 0.3568; %Same time length, tapers, confidence level used always
        end

        if Cerrfv1 > confC1 && Cerrfv2 > confC2
            CorrMat(counter,12) = 1; %We have significant coherence
        else
            CorrMat(counter,12) = 0; %We don't have significant coherence
        end

        counter = counter + 1;
    end

    KFmat(i,1) = phi_struct(i).phi_mat(1,4); %fv for this vessel (from top cross section)
    KFmat(i,2) = phi_struct(i).phi_mat(1,2); %k_fv for this vessel
    KFmat(i,3) = phi_struct(i).phi_mat(1,3)/abs(phi_struct(i).phi_mat(1,1))*abs(phi_struct(i).phi_mat(1,2)); %std k_fv (from std phi, assume no uncertainty in distance)
    KFmat(i,4) = phi_struct(i).Ss;
    KFmat(i,5) = phi_struct(i).Sd;

    %Get coherence at vasofreq (in phi_mat(1,4))
    freq_findx = max(find(round(phi_struct(i).f,3)==round(phi_struct(i).phi_mat(1,4),3)));
%     Cerrfv = phi_struct(i).Cerr(1,freq_findx); %Coherence lower 95% CI
    Cerrfv = phi_struct(i).C(freq_findx);% Coherence value
    if isfield(phi_struct,'confC')
        confC = phi_struct(i).confC;
    else
        confC = 0.3568; %Same time length, tapers, confidence level used always 
    end
    if Cerrfv > confC
        KFmat(i,6) = 1; %We have significant coherence
    else
        KFmat(i,6) = 0; %We don't have significant coherence
    end

    if i==length(phi_struct)-1
        KFmat(i+1,1) = phi_struct(i+1).phi_mat(1,4); %fv for this vessel (from top cross section)
        KFmat(i+1,2) = phi_struct(i+1).phi_mat(1,2); %k_fv for this vessel
        KFmat(i+1,3) = phi_struct(i+1).phi_mat(1,3)/abs(phi_struct(i+1).phi_mat(1,1))*abs(phi_struct(i+1).phi_mat(1,2)); %std k_fv (from std phi, assume no uncertainty in distance)
        KFmat(i+1,4) = phi_struct(i+1).Ss;
        KFmat(i+1,5) = phi_struct(i+1).Sd;
        freq_findx = max(find(round(phi_struct(i+1).f,3)==round(phi_struct(i+1).phi_mat(1,4),3)));
        Cerrfv = phi_struct(i+1).C(freq_findx); %Test for mean coherence > Confidence level.

        if isfield(phi_struct,'confC')
            confC = phi_struct(i+1).confC;
        else
            confC = 0.3568; %Same time length, tapers, confidence level used always
        end
        if Cerrfv > confC
            KFmat(i+1,6) = 1; %We have significant coherence
        else
            KFmat(i+1,6) = 0; %We don't have significant coherence
        end
    end
end

cd(save_folder);
if strcmp(animal,'20230224PA10') %For this trial we measured #18 twice, want to keep that information for the phi_struct and KFmat, but not correlate its direction with itself in the CorrMat.
     Corr_KF_struct = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Pipeline\7_8_23_Analysis\20230224PA10_Corr_KF_struct.mat");
     Corr_KF_struct = Corr_KF_struct.Corr_KF_struct;
     Corr_KF_struct(1).CorrMat = CorrMat; %Just save the CorrMat, KFmat and phi_struct should keep all measurements.
     save([animal,'_Corr_KF_struct.mat'],'Corr_KF_struct')
else
    Corr_KF_struct(1).phi_struct = phi_struct;
    Corr_KF_struct(1).CorrMat = CorrMat;
    Corr_KF_struct(1).KFmat = KFmat;
    save([animal,'_Corr_KF_struct_stim.mat'],'Corr_KF_struct')
end


%% Combine all
clear; close all; clc;
KM1 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Pipeline\7_8_23_Analysis\20221212PA7_Corr_KF_struct.mat"); KM1 = KM1.Corr_KF_struct.KFmat;
KM2 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Pipeline\7_8_23_Analysis\20221129PA6_Corr_KF_struct.mat"); KM2 = KM2.Corr_KF_struct.KFmat;
KM4 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Pipeline\7_8_23_Analysis\20221128PA5_Corr_KF_struct.mat"); KM4 = KM4.Corr_KF_struct.KFmat;
KM5 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Pipeline\7_8_23_Analysis\20221208PA5_Corr_KF_struct.mat"); KM5 = KM5.Corr_KF_struct.KFmat;
KM6 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Pipeline\7_8_23_Analysis\20230118PA7_Corr_KF_struct.mat"); KM6 = KM6.Corr_KF_struct.KFmat;
KM7 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Pipeline\7_8_23_Analysis\20230125PA7_Corr_KF_struct.mat"); KM7 = KM7.Corr_KF_struct.KFmat;
KM8 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Pipeline\7_8_23_Analysis\20230210PA8_Corr_KF_struct.mat"); KM8 = KM8.Corr_KF_struct.KFmat;
KM9 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Pipeline\7_8_23_Analysis\20230213PA8_Corr_KF_struct.mat"); KM9 = KM9.Corr_KF_struct.KFmat;
KM10 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Pipeline\7_8_23_Analysis\20230214PA9_Corr_KF_struct.mat"); KM10 = KM10.Corr_KF_struct.KFmat;
KM11 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Pipeline\7_8_23_Analysis\20230215PA9_Corr_KF_struct.mat"); KM11 = KM11.Corr_KF_struct.KFmat;
KM12 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Pipeline\7_8_23_Analysis\20230222PA10_Corr_KF_struct.mat"); KM12 = KM12.Corr_KF_struct.KFmat;
KM13 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Pipeline\7_8_23_Analysis\20230224PA10_Corr_KF_struct.mat"); KM13 = KM13.Corr_KF_struct.KFmat;
KM14 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Pipeline\7_8_23_Analysis\20230301PA11_Corr_KF_struct.mat"); KM14 = KM14.Corr_KF_struct.KFmat;
KM15 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Pipeline\7_8_23_Analysis\20230302PA11_Corr_KF_struct.mat"); KM15 = KM15.Corr_KF_struct.KFmat;
CombinedKFmat = [KM1;KM2;KM4;KM5;KM6;KM7;KM8;KM9;KM10;KM11;KM12;KM13;KM14;KM15];


CM1 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Pipeline\7_8_23_Analysis\20221212PA7_Corr_KF_struct.mat"); CM1 = CM1.Corr_KF_struct.CorrMat;
CM2 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Pipeline\7_8_23_Analysis\20221129PA6_Corr_KF_struct.mat"); CM2 = CM2.Corr_KF_struct.CorrMat;
CM4 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Pipeline\7_8_23_Analysis\20221128PA5_Corr_KF_struct.mat"); CM4 = CM4.Corr_KF_struct.CorrMat;
CM5 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Pipeline\7_8_23_Analysis\20221208PA5_Corr_KF_struct.mat"); CM5 = CM5.Corr_KF_struct.CorrMat;
CM6 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Pipeline\7_8_23_Analysis\20230118PA7_Corr_KF_struct.mat"); CM6 = CM6.Corr_KF_struct.CorrMat;
CM7 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Pipeline\7_8_23_Analysis\20230125PA7_Corr_KF_struct.mat"); CM7 = CM7.Corr_KF_struct.CorrMat;
CM8 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Pipeline\7_8_23_Analysis\20230210PA8_Corr_KF_struct.mat"); CM8 = CM8.Corr_KF_struct.CorrMat;
CM9 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Pipeline\7_8_23_Analysis\20230213PA8_Corr_KF_struct.mat"); CM9 = CM9.Corr_KF_struct.CorrMat;
CM10 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Pipeline\7_8_23_Analysis\20230214PA9_Corr_KF_struct.mat"); CM10 = CM10.Corr_KF_struct.CorrMat;
CM11 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Pipeline\7_8_23_Analysis\20230215PA9_Corr_KF_struct.mat"); CM11 = CM11.Corr_KF_struct.CorrMat;
CM12 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Pipeline\7_8_23_Analysis\20230222PA10_Corr_KF_struct.mat"); CM12 = CM12.Corr_KF_struct.CorrMat;
CM13 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Pipeline\7_8_23_Analysis\20230224PA10_Corr_KF_struct.mat"); CM13 = CM13.Corr_KF_struct.CorrMat; %This day we measured PA18 twice, both positive phase, so am choosing one tral for the corrmat
CM14 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Pipeline\7_8_23_Analysis\20230301PA11_Corr_KF_struct.mat"); CM14 = CM14.Corr_KF_struct.CorrMat;
CM15 = load("\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Pipeline\7_8_23_Analysis\20230302PA11_Corr_KF_struct.mat"); CM15 = CM15.Corr_KF_struct.CorrMat;
CombinedCorrMat = [CM1;CM2;CM4;CM5;CM6;CM7;CM8;CM9;CM10;CM11;CM12;CM13;CM14;CM15];

cd('\\dk-server.dk.ucsd.edu\jaduckwo\DataAnalysis\VesCorrPhase\Rui_2P\Pipeline\7_8_23_Analysis'); 
save('CombinedKFmat.mat','CombinedKFmat');
save('CombinedCorrmat_WT10_OnePA18Measurement.mat','CombinedCorrMat');




