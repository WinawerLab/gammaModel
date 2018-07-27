addpath('~/Documents/git/ecogBasicCode/')
addpath(['/Volumes/DoraBigDrive/data/visual_soc/m-files']);
addpath(genpath('~/Documents/m-files/knkutils'));
%% load all data together
%%

clear all
dataRootPath = '/Volumes/DoraBigDrive/data/visual_soc/soc_bids';
subjects = [9];
s = 1;

% analysisType = 'spectra';
% analysisType = 'spectraRERP';
% analysisType = 'spectra500';
% analysisType = 'spectraRERP500';
analysisType = 'spectra100';
% analysisType = 'spectraRERP200';
% analysisType = 'spectraShort200';
% analysisType = 'spectraLate500';
% analysisType = 'spectra1000';

if subjects(s)<10
    subj = ['0' int2str(subjects(s))];
else
    subj = int2str(subjects(s));
end

dataName = [dataRootPath '/sub-' subj '/ses-01/derivatives/ieeg/sub-' subj '_task-objects_allruns_' analysisType '.mat'];
load(dataName,'f','spectra','spectra_off','stims','runs','exclude_channels','include_channels')   

%% test fitting parameters in one electrode, one stimulus:
% channels to plot (sub-09: 68, 69):
elec = 68;

data_fft = squeeze(spectra(elec,:,:));
data_fft_off=squeeze(spectra_off(elec,:,:));

f_use4fit=[30:57 65:115 126:175 186:200];
f_sel=ismember(f,f_use4fit);

figure('Position',[0 0 1000 600])
% fit for first stim with 1 gaussian
stims_plot = unique(stims(:,2));
for k = 1:50%length(stims_plot)
    % get baseline and stimulus data
    data_base = mean(data_fft_off,1); % baseline
    data_fit = mean(data_fft(stims(:,2)==stims_plot(k),:),1); % stimuli
    
    [out_exp,bb_amp,gamma_amp,gamma_freq,gamma_width,fit_f2] = ...
        ecog_fitgamma(f,f_use4fit,data_base,data_fit);
    resamp_parms(k,1) = out_exp(1); % this is the baseline slope
    resamp_parms(k,2) = bb_amp;
    resamp_parms(k,3) = gamma_amp;
    resamp_parms(k,4) = gamma_freq;
    resamp_parms(k,5) = gamma_width;
    resamp_parms(k,6) = out_exp(2); % this is the baseline intercept

    % make figure:
    if length(stims_plot) > 5
        subplot(5,10,k),hold on
    else
        subplot(length(stims_plot),1,k),hold on
    end
    plot(log10(f),log10(data_fit),'r')
    plot(log10(f),fit_f2,'b')
    plot(log10(f),log10(data_base),'k')

    xlim([log10(20) log10(200)])
    ylim([-2.5 3.5])
    set(gca,'XTick',[log10(25) log10(50) log10(100) log10(200)])
    set(gca,'XTickLabel',{'25','50','100','200'})
    set(gca,'YTick',[])
    ylabel([int2str(k)])

end

% set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',['./figures/fit/fit_sub-' int2str(subj) '_el'  int2str(elec)])
% print('-depsc','-r300',['./figures/fit/fit_sub-' int2str(subj) '_el'  int2str(elec)])

%%
%% Fit the spectra
%%

% Choose some electodes, later do across all electrodes, overnight calculation
% electrodes = [66 67 68 69]; % sub-09

electrodes = [66 67 68 69]; 

for elec = electrodes

    stimnr_list = unique(stims(:,2));

    stims_repeat = [];
    stims_nonrepeat = [];
    cntr_rep = 0;
    cntr_nonrep = 0;
    for k = 1:length(stimnr_list)
        stim_trials = length(find(stims(:,2)==stimnr_list(k)));
        if stim_trials == 1
            cntr_nonrep = cntr_nonrep+1;
            stims_nonrepeat(cntr_nonrep) = stimnr_list(k);
        elseif stim_trials > 1
            cntr_rep = cntr_rep+1;
            stims_repeat(cntr_rep) = stimnr_list(k);
        end
    end
    
    % Repeated items - define output:
    resamp_parms_allrepeat = NaN(length(stims_repeat),7);
    resamp_parms_evenrepeat = NaN(length(stims_repeat),7);
    resamp_parms_oddrepeat = NaN(length(stims_repeat),7);
    
    % Nonrepeated items - define output:
    resamp_parms_nonrepeat = NaN(length(stims_nonrepeat),7);
    
    f_use4fit = [30:57 65:115 126:175 186:200];
    f_sel = ismember(f,f_use4fit);
    f_alpha = find(f>=8 & f<=13);
       
    % For repeated items:
    for k = 1:length(stims_repeat)
        disp(['fitting repeated stimulus ' int2str(k) ' of ' int2str(length(stims_repeat))])

        thisStims = find(stims(:,2)==stims_repeat(k));
        
        data_fft = squeeze(spectra(elec,thisStims,:));
        data_base = squeeze(mean(spectra_off(elec,:,:),2))'; 
        
        %%%% FIT ALL TRIALS
        data_fit = mean(data_fft,1); % all trials
        [out_exp,bb_amp,gamma_amp,gamma_freq,gamma_width,fit_f2] = ...
            ecog_fitgamma(f,f_use4fit,data_base,data_fit);
        resamp_parms_allrepeat(k,1)=out_exp(1); % this is the slope used in all cases
        resamp_parms_allrepeat(k,2)=bb_amp;
        resamp_parms_allrepeat(k,3)=gamma_amp;
        resamp_parms_allrepeat(k,4)=gamma_freq;
        resamp_parms_allrepeat(k,5)=gamma_width;
        resamp_parms_allrepeat(k,6)=out_exp(2); % this is the baseline intercept  
        % calculate alpha change:
        resamp_parms_allrepeat(k,7) = mean(log10(data_fit(f_alpha)) - log10(data_base(f_alpha)));
        
        %%%% FIT EVEN TRIALS
        data_fit = mean(data_fft(2:2:end,:),1); % even trials
        [out_exp,bb_amp,gamma_amp,gamma_freq,gamma_width,fit_f2] = ...
            ecog_fitgamma(f,f_use4fit,data_base,data_fit);
        resamp_parms_evenrepeat(k,1)=out_exp(1); % this is the slope used in all cases
        resamp_parms_evenrepeat(k,2)=bb_amp;
        resamp_parms_evenrepeat(k,3)=gamma_amp;
        resamp_parms_evenrepeat(k,4)=gamma_freq;
        resamp_parms_evenrepeat(k,5)=gamma_width;
        resamp_parms_evenrepeat(k,6)=out_exp(2); % this is the baseline intercept  
        % calculate alpha change:
        resamp_parms_evenrepeat(k,7) = mean(log10(data_fit(f_alpha)) - log10(data_base(f_alpha)));
        
        %%%% FIT ODD TRIALS
        data_fit = mean(data_fft(1:2:end,:),1); % off trials
        [out_exp,bb_amp,gamma_amp,gamma_freq,gamma_width,fit_f2] = ...
            ecog_fitgamma(f,f_use4fit,data_base,data_fit);
        resamp_parms_oddrepeat(k,1)=out_exp(1); % this is the slope used in all cases
        resamp_parms_oddrepeat(k,2)=bb_amp;
        resamp_parms_oddrepeat(k,3)=gamma_amp;
        resamp_parms_oddrepeat(k,4)=gamma_freq;
        resamp_parms_oddrepeat(k,5)=gamma_width;
        resamp_parms_oddrepeat(k,6)=out_exp(2); % this is the baseline intercept  
        % calculate alpha change:
        resamp_parms_oddrepeat(k,7) = mean(log10(data_fit(f_alpha)) - log10(data_base(f_alpha)));       
    end
    
    % For nonrepeated items:
    for k = 1:length(stims_nonrepeat)
        disp(['fitting nonrepeated stimulus ' int2str(k) ' of ' int2str(length(stims_nonrepeat))])

        thisStim = find(stims(:,2)==stims_nonrepeat(k));
        
        data_fit = squeeze(spectra(elec,thisStim,:))';
        data_base = squeeze(mean(spectra_off(elec,:,:),2))'; 
        
        %%%% FIT ALL TRIALS
        [out_exp,bb_amp,gamma_amp,gamma_freq,gamma_width,fit_f2] = ...
            ecog_fitgamma(f,f_use4fit,data_base,data_fit);
        resamp_parms_nonrepeat(k,1)=out_exp(1); % this is the slope used in all cases
        resamp_parms_nonrepeat(k,2)=bb_amp;
        resamp_parms_nonrepeat(k,3)=gamma_amp;
        resamp_parms_nonrepeat(k,4)=gamma_freq;
        resamp_parms_nonrepeat(k,5)=gamma_width;
        resamp_parms_nonrepeat(k,6)=out_exp(2); % this is the baseline intercept  
        % calculate alpha change:
        resamp_parms_nonrepeat(k,7) = mean(log10(data_fit(f_alpha)) - log10(data_base(f_alpha)));
    end

    [a,b] = fileparts(dataName);
    dataFitName = [a '/' b '_fitEl' int2str(elec)];
    clear a b

    save(dataFitName,'resamp_parms_allrepeat','resamp_parms_evenrepeat','resamp_parms_oddrepeat','resamp_parms_nonrepeat')
end

%% plot all stimulus conditions in a row
clear all

dataRootPath = '/Volumes/DoraBigDrive/data/visual_soc/soc_bids';
subj = '09';

%%%%% Pick an electrode:
elec = 66;

% analysisType = 'spectra200';
% analysisType = 'spectra'; % 300 ms window
analysisType = 'spectra500';

% Load data:
dataFitName = [dataRootPath '/sub-' subj '/derivatives/ieeg/sub-' subj '_task-objects_allruns_' analysisType '_fitEl' int2str(elec) '.mat'];
load(dataFitName)
%%

figure('Position',[0 0 800 300])

% repeats:
bb_base = resamp_parms_allrepeat(1,6); % from the baseline is the same for resamp_parms(:,6)
ecog_bb = squeeze(resamp_parms_allrepeat(:,2))-bb_base;
ecog_g = squeeze(resamp_parms_allrepeat(:,3));
ecog_a = squeeze(resamp_parms_allrepeat(:,7));
ylims = [...
    [min(ecog_bb) max(ecog_bb)];...
    [0 max(ecog_g)];...
    [min(ecog_a) max(ecog_a)] ];

subplot(3,5,1),hold on
bar(ecog_bb,1,'b','EdgeColor',[1 1 1]);
xlim([0 size(resamp_parms_allrepeat,1)+1]),ylim(ylims(1,:))
title(['elec ' int2str(elec) ' repeated stimuli'])

subplot(3,5,6),hold on
bar(ecog_g,1,'b','EdgeColor',[1 1 1]);
xlim([0 size(resamp_parms_allrepeat,1)+1]),ylim(ylims(2,:))

subplot(3,5,11),hold on
bar(ecog_a,1,'b','EdgeColor',[1 1 1]);
xlim([0 size(resamp_parms_allrepeat,1)+1]),ylim(ylims(3,:))

% non repeats:
bb_base = resamp_parms_nonrepeat(1,6); % from the baseline is the same for resamp_parms(:,6)
ecog_bb = squeeze(resamp_parms_nonrepeat(:,2))-bb_base;
ecog_g = squeeze(resamp_parms_nonrepeat(:,3));
ecog_a = squeeze(resamp_parms_nonrepeat(:,7));
ylims = [...
    [min(ecog_bb) max(ecog_bb)];...
    [0 max(ecog_g)];...
    [min(ecog_a) max(ecog_a)] ];

subplot(3,5,2:4),hold on
bar(ecog_bb,1,'b','EdgeColor',[1 1 1]);
xlim([0 size(resamp_parms_nonrepeat,1)+1]),ylim(ylims(1,:))
title(['elec ' int2str(elec) ' non repeated stimuli'])

subplot(3,5,7:9),hold on
bar(ecog_g,1,'b','EdgeColor',[1 1 1]);
xlim([0 size(resamp_parms_nonrepeat,1)+1]),ylim(ylims(2,:))

subplot(3,5,12:14),hold on
bar(ecog_a,1,'b','EdgeColor',[1 1 1]);
xlim([0 size(resamp_parms_nonrepeat,1)+1]),ylim(ylims(3,:))

% plot correlation repeated even and odd items:
bb_base = resamp_parms_evenrepeat(1,6); % from the baseline is the same for resamp_parms(:,6)
ecog_bb_even = squeeze(resamp_parms_evenrepeat(:,2))-bb_base;
ecog_g_even = squeeze(resamp_parms_evenrepeat(:,3));
ecog_a_even = squeeze(resamp_parms_evenrepeat(:,7));
bb_base = resamp_parms_oddrepeat(1,6); % from the baseline is the same for resamp_parms(:,6)
ecog_bb_odd = squeeze(resamp_parms_oddrepeat(:,2))-bb_base;
ecog_g_odd = squeeze(resamp_parms_oddrepeat(:,3));
ecog_a_odd = squeeze(resamp_parms_oddrepeat(:,7));

subplot(3,5,5),hold on
plot(ecog_bb_even,ecog_bb_odd,'k.')
axis tight
title(['r^2 = ' num2str(corr(ecog_bb_even,ecog_bb_odd).^2),3])
subplot(3,5,10),hold on
plot(ecog_g_even,ecog_g_odd,'k.')
axis tight
title(['r^2 = ' num2str(corr(ecog_g_even,ecog_g_odd).^2),3])
subplot(3,5,15),hold on
plot(ecog_a_even,ecog_a_odd,'k.')
axis tight
title(['r^2 = ' num2str(corr(ecog_a_even,ecog_a_odd).^2),3])

% set(gcf,'PaperPositionMode','auto')
% print('-depsc','-r300',['./figures/sub-' int2str(subj) '_' analysisType '_el' int2str(elec) '_bb_g_a'])
% print('-dpng','-r300',['./figures/sub-' int2str(subj) '_' analysisType '_el' int2str(elec) '_bb_g_a'])

%% Reliability

clear all

dataRootPath = '/Volumes/DoraBigDrive/data/visual_soc/soc_bids';
subj = '09';

electrodes = [66 67 68 69];
    
analysisTypes = {'spectraShort200','spectra100','spectra200','spectra','spectra500','spectraLate500','spectra1000'}; % window length for epoch 0-1000 ms: 200, 300 and 500

cod_bb = zeros(length(electrodes),length(analysisTypes));
cod_g = zeros(length(electrodes),length(analysisTypes));
cod_a = zeros(length(electrodes),length(analysisTypes));

for elInd = 1:length(electrodes)
    elec = electrodes(elInd);
    
    for kk = 1:length(analysisTypes)

        %%%% Choose an analysis type for Broadband:
        analysisType = analysisTypes{kk};
        %%%% Load the fitting results:
        dataFitName = [dataRootPath '/sub-' subj '/ses-01/derivatives/ieeg/sub-' subj '_task-objects_allruns_' analysisType '_fitEl' int2str(elec) '.mat'];
        load(dataFitName)
        resamp_parms = resamp_parms_allrepeat;
        %%%% Broadband estimate, one value per image:
        bb_base = resamp_parms(1,6); % from the baseline is the same for resamp_parms(:,6)
    %     ecog_bb = squeeze(resamp_parms(:,2))-bb_base;
    %     ecog_g = squeeze(resamp_parms(:,3));
    %     ecog_a = squeeze(resamp_parms(:,7));

        cod_bb(elInd,kk) = calccod(resamp_parms_evenrepeat(:,2)-bb_base,resamp_parms_oddrepeat(:,2)-bb_base);
        cod_g(elInd,kk) = calccod(resamp_parms_evenrepeat(:,3),resamp_parms_oddrepeat(:,3));
        cod_a(elInd,kk) = calccod(resamp_parms_evenrepeat(:,7),resamp_parms_oddrepeat(:,7));
    end
end

figure('Position',[0 0 400 400])

subplot(2,1,1),hold on
title('Reliability for different window length')
bar(cod_bb./100)
xlabel('electrode')
ylabel('BB COD (R)')
ylim([-0.1 1]), set(gca,'XTick',[1 2 3 4],'XTickLabel',electrodes), grid on

subplot(2,1,2),hold on
bar(cod_g./100)
xlabel('electrode')
ylabel('Gamma COD (R)')
ylim([-0.1 1]), set(gca,'XTick',[1 2 3 4],'XTickLabel',electrodes), grid on
legend('200Short','100','200','300','500','500Late','1000')

% set(gcf,'PaperPositionMode','auto')
% print('-depsc','-r300',['./figures/spectralReliability/sub-' subj '_CODdifferentAnalyses'])
% print('-dpng','-r300',['./figures/spectralReliability/sub-' subj '_CODdifferentAnalyses'])

%% contrast to noise TODO