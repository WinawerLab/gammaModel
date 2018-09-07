clear all

addpath('~/Documents/git/ecogBasicCode/')
addpath(['/Volumes/DoraBigDrive/data/visual_soc/m-files']);
addpath(genpath('/Users/dorahermes-miller/Documents/m-files/knkutils'));

% script to calculate average power change in bb and gamma in frequency
% ranges without fitting


%% load all data together
%%
clear all
dataRootPath = '/Volumes/DoraBigDrive/data/visual_soc/soc_bids';
subjects = [19,23,24];
s = 1;

% analysisType = 'spectra';
% analysisType = 'spectra500';
% analysisType = 'spectra300';
analysisType = 'spectra200';
% analysisType = 'spectra100';
% analysisType = 'spectraRERP';
% analysisType = 'spectraRERP500';
% analysisType = 'spectraRERP200';

subj = subjects(s);
dataName = [dataRootPath '/sub-' int2str(subj) '/ses-01/derivatives/ieeg/sub-' int2str(subj) '_task-soc_allruns_' analysisType '.mat'];
load(dataName,'f','spectra','spectra_off','stims','runs','exclude_channels','include_channels')   

%%
% channels to plot (s1: 115 | s2: 53 54 | s3: 45 46):
electrodes = [107 108 109 115 120 121]; % S1
% electrodes = [53 54]; % S2
% electrodes = [45 46]; % S3

electrodes = [109]; % S3

for k = 1:length(electrodes)
    elec = electrodes(k);
    
    data_fft = squeeze(spectra(elec,:,:));
    data_base = mean(squeeze(spectra_off(elec,:,:)),1);
    
    f_use4fit=[30:57 65:115 126:175 186:200];
    f_sel=ismember(f,f_use4fit);
    
    % calculate 1./f from baseline
    [out_exp,bb_amp,gamma_amp,gamma_freq,gamma_width,fit_f2] = ...
        ecog_fitgamma(f,f_use4fit,data_base,data_base);
    
    % subtract 1./f
    data_fft_adj = bsxfun(@minus,log10(data_fft),fit_f2');
    
    data_bb = mean(data_fft_adj(:,f>100 & f<200),2);
    data_g = mean(data_fft_adj(:,f>25 & f<80),2)-mean(data_fft_adj(:,f>100 & f<200),2);
    
    data_bb_stim = zeros(max(stims),2);
    data_g_stim = zeros(max(stims),2);
    for ss = 1:max(stims)
        data_bb_stim(ss,1) = mean(data_bb(stims==ss));
        data_bb_stim(ss,2) = std(data_bb(stims==ss))./sqrt(length(find(stims==ss)));
        data_g_stim(ss,1) = mean(data_g(stims==ss));
        data_g_stim(ss,2) = std(data_g(stims==ss))./sqrt(length(find(stims==ss)));
    end   
    
    figure('Position',[0 0 600 400])
    
    ylims = [...
        [min(data_bb_stim(:)) max(data_bb_stim(:))];...
        [min(data_g_stim(:)) max(data_g_stim(:))]];

    subplot(2,1,1),hold on
    bar(data_bb_stim(:,1),1,'b','EdgeColor',[1 1 1]);
    plot([1:size(data_bb_stim,1); [1:size(data_bb_stim,1)]],[data_bb_stim(:,1)-data_bb_stim(:,2) data_bb_stim(:,1)+data_bb_stim(:,2)]','k');
    title(['elec ' int2str(elec)])
    
    subplot(2,1,2),hold on
    bar(data_g_stim(:,1),1,'b','EdgeColor',[1 1 1]);
    plot([1:size(data_g_stim,1); [1:size(data_g_stim,1)]],[data_g_stim(:,1)-data_g_stim(:,2) data_g_stim(:,1)+data_g_stim(:,2)]','k');
    
    for s = 1:2
        subplot(2,1,s),hold on
        % plot stimulus cutoffs
        stim_change=[38.5 46.5 50.5 54.5 58.5 68.5 73.5 78.5 82.5];
        for k=1:length(stim_change)
            plot([stim_change(k) stim_change(k)],ylims(s,:),'r','MarkerSize',20)
        end
    end
    subplot(2,1,1),hold on
    text([19 40 46 52 55 61 70 74 79 84],zeros(10,1)-.25,{'space','orie','grat','pl','circ','zcon','sp','zsp','coh','nm'})
    
%     set(gcf,'PaperPositionMode','auto')
%     print('-depsc','-r300',['./figures/sub-' int2str(subj) '_' analysisType '_el' int2str(elec) '_MeanSpectbb_g_a'])
%     print('-dpng','-r300',['./figures/sub-' int2str(subj) '_' analysisType '_el' int2str(elec) '_MeanSpectbb_g_a'])

end



%% reliability
% clear all
dataRootPath = '/Volumes/DoraBigDrive/data/visual_soc/soc_bids';
subjects = [19,23,24];
s = 3;

analysisTypes ={'spectra','spectra500','spectra300','spectra200','spectra100',...
    'spectraRERP','spectraRERP500','spectraRERP200'};

if s==1
    electrodes = [107 108 109 115 120 121]; % S1
elseif s==2
    electrodes = [53 54]; % S2
elseif s==3
    electrodes = [45 46]; % S3
end

subj = subjects(s);

cod_bb = zeros(length(electrodes),length(analysisTypes));
cod_g = zeros(length(electrodes),length(analysisTypes));


for a_t = 1:length(analysisTypes)
    analysisType = analysisTypes{a_t};

    dataName = [dataRootPath '/sub-' int2str(subj) '/ses-01/derivatives/ieeg/sub-' int2str(subj) '_task-soc_allruns_' analysisType '.mat'];
    load(dataName,'f','spectra','spectra_off','stims','runs','exclude_channels','include_channels')   

    for k = 1:length(electrodes)
        elec = electrodes(k);

        data_fft = squeeze(spectra(elec,:,:));
        data_base = mean(squeeze(spectra_off(elec,:,:)),1);

        f_use4fit=[30:57 65:115 126:175 186:200];
        f_sel=ismember(f,f_use4fit);

        % calculate 1./f from baseline
        [out_exp,bb_amp,gamma_amp,gamma_freq,gamma_width,fit_f2] = ...
            ecog_fitgamma(f,f_use4fit,data_base,data_base);

        % subtract 1./f
        data_fft_adj = bsxfun(@minus,log10(data_fft),fit_f2');

        data_bb = mean(data_fft_adj(:,f>100 & f<200),2);
        data_g = mean(data_fft_adj(:,f>25 & f<80),2)-mean(data_fft_adj(:,f>100 & f<200),2);

        data_bb_stim = zeros(max(stims),2);
        data_g_stim = zeros(max(stims),2);

        for ss = 1:max(stims)
            bb_stim = data_bb(stims==ss);
            g_stim = data_g(stims==ss);
            data_bb_stim(ss,1) = mean(bb_stim(1:2:end));
            data_bb_stim(ss,2) = mean(bb_stim(2:2:end));
            data_g_stim(ss,1) = mean(g_stim(1:2:end));
            data_g_stim(ss,2) = mean(g_stim(2:2:end));
            clear bb_stim g_stim
        end   

        cod_bb(k,a_t) = calccod(data_bb_stim(:,1),data_bb_stim(:,2));
        cod_g(k,a_t) = calccod(data_g_stim(:,1),data_g_stim(:,2));
    end 
end

%%
figure
subplot(2,1,1),hold on
bar(cod_bb)
subplot(2,1,2),hold on
bar(cod_g)

set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300',['./figures/sub-' int2str(subj) '_reliability_bb_g_a_meanSpectra'])
print('-dpng','-r300',['./figures/sub-' int2str(subj) '_reliability_bb_g_a_meanSpectra'])

