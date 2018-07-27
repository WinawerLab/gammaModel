addpath('~/Documents/git/ecogBasicCode/')
addpath(['/Volumes/DoraBigDrive/data/visual_soc/m-files']);
addpath(genpath('~/Documents/m-files/knkutils'));

%%
%% load all data together
%%

clear all
dataRootPath = '/Volumes/DoraBigDrive/data/visual_soc/soc_bids';
subjects = [19,23,24];
s = 3;

% analysisType = 'spectra';
% analysisType = 'spectraRERP';
% analysisType = 'spectra500';
% analysisType = 'spectraRERP500';
% analysisType = 'spectra100';
analysisType = 'spectra200';
% analysisType = 'spectraRERP200';

subj = subjects(s);
dataName = [dataRootPath '/sub-' int2str(subj) '/ses-01/derivatives/ieeg/sub-' int2str(subj) '_task-soc_allruns_' analysisType '.mat'];
load(dataName,'f','spectra','spectra_off','stims','runs','exclude_channels','include_channels')   

%% test fitting parameters in one electrode, one stimulus:
% channels to plot (s1: 108 | s2: 53 54 | s3: 45 46):
elec = 109;

data_fft = squeeze(spectra(elec,:,:));
data_fft_off=squeeze(spectra_off(elec,:,:));

f_use4fit=[30:57 65:115 126:175 186:200];
f_sel=ismember(f,f_use4fit);

figure('Position',[0 0 1000 600])
% fit for first stim with 1 gaussian
stims_plot=[1:47];
for k=1:length(stims_plot)
    % get baseline and stimulus data
    data_base=mean(data_fft_off,1); % baseline
    data_fit=mean(data_fft(stims==stims_plot(k),:),1); % stimuli
    
    [out_exp,bb_amp,gamma_amp,gamma_freq,gamma_width,fit_f2] = ...
        ecog_fitgamma(f,f_use4fit,data_base,data_fit);
    resamp_parms(k,1)=out_exp(1); % this is the baseline slope
    resamp_parms(k,2)=bb_amp;
    resamp_parms(k,3)=gamma_amp;
    resamp_parms(k,4)=gamma_freq;
    resamp_parms(k,5)=gamma_width;
    resamp_parms(k,6)=out_exp(2); % this is the baseline intercept

    % make figure:
    if length(stims_plot)>5
        subplot(5,ceil(length(stims_plot)/5),k),hold on
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
%% Fit gamma/bb for electrodes that we want to analyze 
%%

% Choose some electodes, later do across all electrodes, overnight calculation
% electrodes = [107 108 109 115 120 121]; % S1
% electrodes = [53 54]; % S2
% electrodes = [45 46]; % S3

electrodes = [107 108 115]; % S1

% nr of bootstrap for resampling per stimulus
nr_boots = 100;

for jj = 1:length(electrodes)

    elec = electrodes(jj);
    % define output:
    resamp_parms = NaN(max(stims),nr_boots,7);

    data_fft = squeeze(spectra(elec,:,:));
    data_fft_off=squeeze(spectra_off(elec,:,:));

    f_use4fit = [30:57 65:115 126:175 186:200];
    f_sel = ismember(f,f_use4fit);
    f_alpha = find(f>=8 & f<=13);

    for k = 1:max(stims)
        disp(['fitting stimulus ' int2str(k) ' of ' int2str(max(stims))])
        % for the baseline, do not resamplel just average across 1290 trials
        data_base = mean(data_fft_off,1); % baseline

        for ii = 1:nr_boots
            % get stimulus data       
            data_fit = data_fft(stims==k,:); 

            % from baseline data, random sample with replacement
            trial_set = randsample(size(data_fit,1),size(data_fit,1),true);

            % average across resampled trials
            data_fit = mean(data_fit(trial_set,:),1);

            % do the fitting
            [out_exp,bb_amp,gamma_amp,gamma_freq,gamma_width,fit_f2] = ...
                ecog_fitgamma(f,f_use4fit,data_base,data_fit);
            resamp_parms(k,ii,1) = out_exp(1); % this is the slope used in all cases
            resamp_parms(k,ii,2) = bb_amp;
            resamp_parms(k,ii,3) = gamma_amp;
            resamp_parms(k,ii,4) = gamma_freq;
            resamp_parms(k,ii,5) = gamma_width;
            resamp_parms(k,ii,6) = out_exp(2); % this is the baseline intercept
            clear trial_set
            
            % calculate alpha change
            resamp_parms(k,ii,7) = mean(log10(data_fit(f_alpha)) - log10(data_base(f_alpha)));
        end
    end

    [a,b] = fileparts(dataName);
    dataFitName = [a '/' b '_fitEl' int2str(elec)];
    clear a b

    save(dataFitName,'resamp_parms')
end

%% plot all stimulus conditions in a row
clear all
dataRootPath = '/Volumes/DoraBigDrive/data/visual_soc/soc_bids';

%%%%% Pick a subject:
subjects = [19,23,24];
s = 3; subj = subjects(s);

%%%%% Pick an electrode:
% electrodes = [107 108 109 115 120 121]; % S1
% electrodes = [53 54]; % S2
% electrodes = [45 46]; % S3
elec = 46;

% Choose an analysis type:
% analysisType = 'spectra';
% analysisType = 'spectraRERP';
% analysisType = 'spectra300';
% analysisType = 'spectraRERP500';
analysisType = 'spectra200';
% analysisType = 'spectra100';
% analysisType = 'spectraRERP200';

dataFitName = [dataRootPath '/sub-' int2str(subj) '/ses-01/derivatives/ieeg/sub-' int2str(subj) '_task-soc_allruns_' analysisType '_fitEl' int2str(elec) '.mat'];

load(dataFitName,'resamp_parms')

bb_base = resamp_parms(1,1,6); % from the baseline is the same for resamp_parms(:,:,6)
ecog_bb = squeeze(median(resamp_parms(:,:,2),2))-bb_base;
ecog_bb_err=[squeeze(quantile(resamp_parms(:,:,2),.025,2)) ...
    squeeze(quantile(resamp_parms(:,:,2),.975,2))]'-bb_base;

ecog_g = squeeze(median(resamp_parms(:,:,3),2));
ecog_g_err=[squeeze(quantile(resamp_parms(:,:,3),.025,2)) ...
    squeeze(quantile(resamp_parms(:,:,3),.975,2))]';

ecog_a = squeeze(median(resamp_parms(:,:,7),2));
ecog_a_err=[squeeze(quantile(resamp_parms(:,:,7),.025,2)) ...
    squeeze(quantile(resamp_parms(:,:,7),.975,2))]';

figure('Position',[0 0 1000 100])

ylims = [...
    [min(ecog_bb_err(:)) max(ecog_bb_err(:))];...
    [0 max(ecog_g_err(:))];...
    [min(ecog_a_err(:)) max(ecog_a_err(:))] ];

subplot(1,2,1),hold on
bar(ecog_bb,1,'b','EdgeColor',[0 0 0]);
plot([1:86; 1:86],ecog_bb_err,'k');
ylim(ylims(1,:))
title(['elec ' int2str(elec)])

subplot(1,2,2),hold on
bar(ecog_g,1,'b','EdgeColor',[0 0 0]);
plot([1:86; 1:86],ecog_g_err,'k');
ylim(ylims(2,:))

% subplot(3,1,3),hold on
% bar(ecog_a,1,'b','EdgeColor',[0 0 0]);
% plot([1:86; 1:86],ecog_a_err,'k');
% ylim(ylims(3,:))

for s = 1:2
    subplot(1,2,s),hold on
    % plot stimulus cutoffs
    stim_change=[38.5 46.5 50.5 54.5 58.5 68.5 73.5 78.5 82.5];
    for k=1:length(stim_change)
        plot([stim_change(k) stim_change(k)],ylims(s,:),'Color',[.5 .5 .5],'LineWidth',2)
    end
    text([19 40 46 52 55 61 70 74 79 84],zeros(10,1)-.25,{'space','orie','grat','pl','circ','zcon','sp','zsp','coh','nm'})
    xlim([0 87])
end

set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300',['./figures/sub-' int2str(subj) '_' analysisType '_el' int2str(elec) '_bb_g_a'])
print('-dpng','-r300',['./figures/sub-' int2str(subj) '_' analysisType '_el' int2str(elec) '_bb_g_a'])

%%
%%  Analysis reliability
%%

clear all
dataRootPath = '/Volumes/DoraBigDrive/data/visual_soc/soc_bids';
subjects = [19,23,24];
s = 2;

% analysisType = 'spectra';
% analysisType = 'spectraRERP';
% analysisType = 'spectra500';
% analysisType = 'spectraRERP500';
% analysisType = 'spectra200';
% analysisType = 'spectraRERP200';
analysisType = 'spectra100';
% analysisType = 'spectra300';

subj = subjects(s);
dataName = [dataRootPath '/sub-' int2str(subj) '/ses-01/derivatives/ieeg/sub-' int2str(subj) '_task-soc_allruns_' analysisType '.mat'];
load(dataName,'f','spectra','spectra_off','stims','runs','exclude_channels','include_channels')   

if s==1
    electrodes = [107 108 109 115 120 121]; % S1
elseif s==2
    electrodes = [53 54]; % S2
elseif s==3
    electrodes = [45 46]; % S3
end

for eNr = 1:length(electrodes)
    elec = electrodes(eNr);

    data_fft = squeeze(spectra(elec,:,:));
    data_fft_off = squeeze(spectra_off(elec,:,:));

    f_use4fit = [30:57 65:115 126:175 186:200];
    f_sel = ismember(f,f_use4fit);

    resamp_parms_odd = zeros(max(stims),6);
    resamp_parms_even = zeros(max(stims),6);

    for kk = 1:max(stims)
        % get baseline and stimulus data
        data_base = mean(data_fft_off,1); % baseline
        data_fit = data_fft(stims==kk,:); % stimuli

        data_fit_odd = mean(data_fit(1:2:size(data_fit,1),:),1);
        data_fit_even = mean(data_fit(2:2:size(data_fit,1),:),1);

        [out_exp,bb_amp,gamma_amp,gamma_freq,gamma_width,fit_f2] = ...
            ecog_fitgamma(f,f_use4fit,data_base,data_fit_odd);
        resamp_parms_odd(kk,1) = out_exp(1); % this is the baseline slope
        resamp_parms_odd(kk,2) = bb_amp;
        resamp_parms_odd(kk,3) = gamma_amp;
        resamp_parms_odd(kk,4) = gamma_freq;
        resamp_parms_odd(kk,5) = gamma_width;
        resamp_parms_odd(kk,6) = out_exp(2); % this is the baseline intercept

        [out_exp,bb_amp,gamma_amp,gamma_freq,gamma_width,fit_f2] = ...
            ecog_fitgamma(f,f_use4fit,data_base,data_fit_even);
        resamp_parms_even(kk,1) = out_exp(1); % this is the baseline slope
        resamp_parms_even(kk,2) = bb_amp;
        resamp_parms_even(kk,3) = gamma_amp;
        resamp_parms_even(kk,4) = gamma_freq;
        resamp_parms_even(kk,5) = gamma_width;
        resamp_parms_even(kk,6) = out_exp(2); % th
    end
    
    [a,b] = fileparts(dataName);
    dataFitName = [a '/' b '_fitEl' int2str(elec) '_evenodd'];
    clear a b
    save(dataFitName,'resamp_parms_even','resamp_parms_odd')
end

%% plot reliability

clear all
dataRootPath = '/Volumes/DoraBigDrive/data/visual_soc/soc_bids';
subjects = [19,23,24];
s = 1;

analysisTypes ={'spectra','spectra500','spectra300','spectra200','spectra100',...
    'spectraRERP','spectraRERP500','spectraRERP200'};

subj = subjects(s);
% load(dataName,'f','spectra','spectra_off','stims','runs','exclude_channels','include_channels')   

if s==1
    electrodes = [107 108 109 115 120 121]; % S1
elseif s==2
    electrodes = [53 54]; % S2
elseif s==3
    electrodes = [45 46]; % S3
end

cod_bb = zeros(length(analysisTypes),length(electrodes));
cod_g = zeros(length(analysisTypes),length(electrodes));

for kk = 1:length(analysisTypes)
    for eNr = 1:length(electrodes)
        elec = electrodes(eNr);
        
        dataFitName = [dataRootPath '/sub-' int2str(subj) '/ses-01/derivatives/ieeg/sub-' int2str(subj) '_task-soc_allruns_' analysisTypes{kk} '_fitEl' int2str(elec) '_evenodd.mat'];
        load(dataFitName,'resamp_parms_even','resamp_parms_odd')
        
        bb_base = resamp_parms_even(1,6); % from the baseline is the same for resamp_parms_odd(:,6)
        bb_even = squeeze(resamp_parms_even(:,2))-bb_base;
        bb_odd = squeeze(resamp_parms_odd(:,2))-bb_base;
        g_even = squeeze(resamp_parms_even(:,3));
        g_odd = squeeze(resamp_parms_odd(:,3));

        cod_bb(kk,eNr) = calccod(bb_even,bb_odd);
        cod_g(kk,eNr) = calccod(g_even,g_odd);   
    end
end


figure
subplot(2,1,1),hold on
bar(cod_bb')
subplot(2,1,2),hold on
bar(cod_g')

set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300',['./figures/sub-' int2str(subj) '_reliability_bb_g_a'])
print('-dpng','-r300',['./figures/sub-' int2str(subj) '_reliability_bb_g_a'])
