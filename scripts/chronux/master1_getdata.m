clear all

addpath(['/Volumes/DoraBigDrive/data/visual_soc/m-files']);
addpath(genpath('/Users/dora/Documents/m-files/Chronux/'))

%% epoch for Chronux timefreq

% NOTE: later add exclusion bad epochs
clear all
dataRootPath = '/Volumes/DoraBigDrive/data/visual_soc/soc_bids';

subjects=[19 23 24];
s = 1;
subj = subjects(s);
dataName = dir([dataRootPath '/sub-' int2str(subj) '/ses-01/ieeg/sub-' int2str(subj) '_ses-01_task-soc_run-*_ieeg_preproc.mat']);
stimName = dir([dataRootPath '/sub-' int2str(subj) '/ses-01/ieeg/sub-' int2str(subj) '_ses-01_task-soc_run-*_events.tsv']);
nr_runs = length(dataName);

all_data_epoch = [];
all_stims = [];

for data_nr = 1:nr_runs % how many runs
    %%%% load preprocessed data
    load([dataRootPath '/sub-' int2str(subj) '/ses-01/ieeg/' dataName(data_nr).name]);
    
    %%%% load stimulus information
    stim = readtable([dataRootPath '/sub-' int2str(subj) '/ses-01/ieeg/' stimName(data_nr).name],...
        'FileType','text','Delimiter','\t');
    stims = stim.trial_type;
    
    %%%% NOTCH FILTER 
    % notch filter data at 60, 120 and 180 Hz
    data = ecog_notch(data,srate);
    
    %%%% MAKE EPOCHS
    onset_trial = round(stim.onset*srate);%from seconds to samples
    epoch_l = 1.2; % epoch length: -0.2:1 sec
    data_epoch = zeros(size(data,2),length(onset_trial),epoch_l*srate);
    for elec = 1:size(data,2)
        for l = 1:length(onset_trial)
            data_epoch(elec,l,:) = ...
                data(onset_trial(l)-.2*srate+1:onset_trial(l)+(epoch_l-.2)*srate,elec)';
        end
    end
    % define t - time vector for each epoch
    t = [1:epoch_l*srate]/srate - 0.2;    
   
    clear data
    all_data_epoch=cat(2,all_data_epoch,data_epoch);
    all_stims=cat(1,all_stims,stims);
    
    disp(['done run ' int2str(data_nr)])
end

% output we have all_data_epoch all_stims

%% check ERP
electr=find(chan_lbls==107);
s=[39:45];
% s=[10 29];
figure,hold on
% plot(t,squeeze(all_data_epoch(electr,find(all_stims==s,1),:)))
plot(t,squeeze(all_data_epoch(electr,ismember(all_stims,s),:)))
% plot(t,squeeze(mean(all_data_epoch(electr,ismember(all_stims,s),:),2)))


%% quick figure

cm1=[repmat([0 0 0],100,1)];
cm1(1:40,1)=[0.7]';
cm1(1:40,2)=[0.7:-0.6/39:0.1]';
cm1(1:40,3)=[0.7:-0.7/39:0]';
cm1(40:100,1)=[0.7:(1-0.7)/60:1]';
cm1(40:100,2)=[0.1:.9/60:1]';

cm2=[repmat([0 0 0],100,1)];
cm2(1:30,3)=[0.7]';
cm2(1:30,1)=[0.7:-0.7/29:0]';
cm2(1:30,2)=[0.7:-0.7/29:0]';
cm2(30:100,3)=[0.7:(1-0.7)/70:1]';
cm2(30:100,2)=[0:1/70:1]';
cm=[cm2(end:-1:1,:); cm1];

movingwin=[.200 .05];
params.pad=-1;
params.tapers=[3 5];
params.fpass=[0 200];
params.Fs=srate;
params.trialave=0;

electrodes=[120];

% conds_plot=[5 10 15 39 40 50 54 58 83 87];
conds_plot={5,[10 29],15,[39:46],40,50,54,58,83,87};
% conds_plot={39 40 41 42 43 44 45 46};

for l = 1:length(electrodes)
    elec=electrodes(l);
    figure('Color',[1 1 1],'Position',[0 0 800 500])
    subplot(1,1,1),hold on

    % regress erp
    % data_temp = ecog_regresserp(all_data_epoch(elec,:,:),t>-0.1 & t<0,all_stims);
    data_temp = all_data_epoch(elec,:,:);

    % calculate baseline
    data2use=squeeze(data_temp(1,all_stims==87,t>.25 & t<.5))';
    [S1b,t_tf_b,f_b]=mtspecgramc(data2use,movingwin,params);
    S1b=mean(mean(S1b,3),1);

    S1_all=zeros(21,size(f_b,2),length(conds_plot));

    %%%%% all responses:
    for k=1:length(conds_plot)
        data2use=squeeze(data_temp(1,ismember(all_stims,conds_plot{k}),:))';

        % calculate spectgram
        [S1,t_tf,f]=mtspecgramc(data2use,movingwin,params);
        t_tf=t_tf+t(1);
        % normalize wrt baseline
        S1=mean(S1,3);
        S1=S1./repmat(S1b,[size(S1,1) 1]);
        S1_all(:,:,k)=S1;
        subplot(2,5,k)
    %     imagesc(t,f,log(S1)',[-3 3])
        imagesc(t_tf,f,log10(S1)',[-1.3 1.3])
%         imagesc(t_tf,f,log10(S1)',[-1 1])
        axis xy
        colormap(cm)
    %     colormap(jet)
        set(gca,'XTick',[0 .5])
        title(['stim ' int2str(conds_plot{k})])
    end

colorbar
% set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',['./figures/ersp/gr_' subj '_el'  int2str(electr) '_cm3_regressERP'])
% print('-depsc','-r300',['./figures/ersp/gr_' subj '_el'  int2str(electr) '_cm3_regressERP'])
% close all
end


%%
% plot average gamma and bb traces
figure('Color',[1 1 1],'Position',[0 0 800 200])
subplot(1,1,1),hold on
for k=1:nr_cond
    subplot(2,10,k),hold on
    plot(t_tf,mean(S1_all(:,f>30 & f<70,k),2),'g');
    xlim([t_tf(1) t_tf(end)]),ylim([0 30])
    set(gca,'XTick',[0 .5])
    
    subplot(2,10,10+k),hold on
    plot(t_tf,mean(S1_all(:,f>100 & f<175,k),2),'b');
    xlim([t_tf(1) t_tf(end)]),ylim([0 12])
    set(gca,'XTick',[0 .5])
end

set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',['./figures/ersp/gr_' subj '_el'  int2str(electr) '_av3070_80200'])
% print('-depsc','-r300',['./figures/ersp/gr_' subj '_el'  int2str(electr) '_av3070_80200'])

%%
% plot average spectra
figure('Color',[1 1 1],'Position',[0 0 200 200])
cond_label={'white n','pink n','brown n','gr 4','gr 8','gr 16','gr 32','blank','b white n','plaid'};
cond_color={[.2 .5 1],[0 .2 1],[0 0 .8],[.5 0 0],[1 0 0],[1 .5 .5],[1 .8 .8],[0 0 0],[.5 .8 1],'g'};
cond_pat={'-','-','-','-','-','-','-','-','-','-'};

subplot(1,1,1),hold on
for k=[1:7]%[1:7 9 10]
    subplot(1,1,1),hold on
    plot(f,mean(S1_all(t_tf>0.15 & t_tf<.5,:,k),1),'Color',cond_color{k});
    ylim([0 40])
end
set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',['./figures/ersp/gr_' subj '_el'  int2str(electr) '_avspectra'])
print('-depsc','-r300',['./figures/ersp/gr_' subj '_el'  int2str(electr) '_avspectra'])


%% cohgram

movingwin=[.400 .05];
params.pad=-1;
params.tapers=[3 5];
% params.tapers=[4 7];
% params.fpass=[0 200];
params.fpass=[35 200];
params.Fs=1000;
params.trialave=1;
params.err=[2 0.05/(17*67)];
% params.err=[2 0.001];

% electrodes=[79 104];% example 97 99 / 114 99 / 115 99 but see 98 90 / 65 103??
electrodes=[112 94];
% electrodes=[112 5];

% 112 V1/V2-f
% 113 V1-f
% 94 hV4
% 95 hV4 - bad
% 99 hV4/VO1
% 100 hV4/VO1
% pretty good: 112 - 94; 112 - 100; 113 - 100

for l=1:size(electrodes,1)
electr1=electrodes(l,1);
elec1=find(chan_lbls==electr1);
electr2=electrodes(l,2);
elec2=find(chan_lbls==electr2);

figure('Color',[1 1 1],'Position',[0 0 800 500])

%%%%% all responses:
for k=1:7%1:8
    data2use1=squeeze(all_data_epoch(elec1,all_stims==k,:))';
    data2use2=squeeze(all_data_epoch(elec2,all_stims==k,:))';
    
%     % subtract ERP
%     elec_erp=zeros(size(data2use1));
%     for tr=1:size(data2use1,2)
%         elec_erp(:,tr)=data2use1(:,tr)-mean(data2use1(t>-0.1 & t<0,tr));
%     end
%     elec_erp=mean(elec_erp,2);
%     data2use1=data2use1-repmat(elec_erp,[1 size(data2use1,2)]);
%     % subtract ERP
%     elec_erp=zeros(size(data2use2));
%     for tr=1:size(data2use2,2)
%         elec_erp(:,tr)=data2use2(:,tr)-mean(data2use2(t>-0.1 & t<0,tr));
%     end
%     elec_erp=mean(elec_erp,2);
%     data2use2=data2use2-repmat(elec_erp,[1 size(data2use2,2)]);
    
    % calculate coherency without jackknife
%     [C,phi,S12,S1,S2,t_tf,f]=...
%         cohgramc(data2use1,data2use2,movingwin,params);
    % calculate coherency with jackknife
    [C,phi,S12,S1,S2,t_tf,f,confC,phistd,Cerr]=...
        cohgramc(data2use1,data2use2,movingwin,params);

    t_tf=t_tf+t(1);
    subplot(2,5,k)
    imagesc(t_tf,f,C',[0 .5]),hold on
    
    % lower confidence interval
    C_contour=squeeze(Cerr(1,:,:));
%     % only include areas where some adjacent thing is also significant
%     C_contour_test=bwlabel(C_contour);
%     for m=1:max(C_contour_test(:))
%         if length(find(C_contour_test==m))<2
%             C_contour(C_contour_test==m)=0;
%         end
%     end

    imcontour(t_tf,f,C_contour',[0.0001],'k','LineWidth',1)
    
    axis xy
    colormap(jet)

    disp([' done stim ' int2str(k)])
end

subplot(2,5,1),ylabel(['els ' int2str(electr1) ' ' int2str(electr2)])
cond_label={'white n','pink n','brown n','gr 4','gr 8','gr 16','gr 32','blank'};
for k=1:7
    subplot(2,5,k),hold on
    title(cond_label{k})
%     ylim([0 200])
    ylim([35 200])
%     plot([t_tf(1) t_tf(end)],[27 27],'k-')
%     plot([t_tf(1) t_tf(end)],[54 54],'k-')
%     plot([t_tf(1) t_tf(end)],[82 82],'k-')
end
subplot(2,5,9)
imagesc(t_tf,f,C',[0 .5]),hold on
colorbar
title('just for colorbar')
set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',['./figures/coh/' subj '_el'  int2str(electr1) '_' int2str(electr2) '_subtrERP'])
% print('-dpng','-r300',['./figures/coh/' subj '_el'  int2str(electr1) '_' int2str(electr2) '_stats1'])
% print('-depsc','-r300',['./figures/coh/' subj '_el'  int2str(electr1) '_' int2str(electr2) '_stats1'])
% print('-depsc','-r300',['./figures/coh/' subj '_el'  int2str(electr1) '_' int2str(electr2) '_stats1_colorbar'])

print('-dpng','-r300',['./figures/coh/' subj '_el'  int2str(electr1) '_' int2str(electr2) '_stats5'])
print('-depsc','-r300',['./figures/coh/' subj '_el'  int2str(electr1) '_' int2str(electr2) '_stats5'])

end

%%
%% quick figure for average and specific trial
%%

% movingwin=[.100 .01];
movingwin=[.200 .05];
params.pad=-1;
params.tapers=[3 5];
params.fpass=[0 200];
params.Fs=srate;
params.trialave=0;
load loc_colormap

electrodes=[108];%

electr=electrodes(1);
elec=find(chan_lbls==electr);

figure('Color',[1 1 1],'Position',[0 0 600 600])
subplot(1,1,1),hold on

% baseline
data2use=squeeze(all_data_epoch(elec,all_stims==8,t>.25 & t<.5))';
[S1b,t_tf_b,f_b]=mtspecgramc(data2use,movingwin,params);
S1b=mean(mean(S1b,3),1);

%%%%% all responses:
for k=1:8
    data2use=squeeze(all_data_epoch(elec,all_stims==k,:))';
    
%     % subtract ERP
%     elec_erp=zeros(size(data2use));
%     for tr=1:size(data2use,2)
%         elec_erp(:,tr)=data2use(:,tr)-mean(data2use(t>-0.1 & t<0,tr));
%     end
%     elec_erp=mean(elec_erp,2);
%     data2use=data2use-repmat(elec_erp,[1 size(data2use,2)]);
    
    % calculate spectgram
    [S1,t_tf,f]=mtspecgramc(data2use,movingwin,params);
    t_tf=t_tf-.2;
    % normalize wrt baseline
    S1=mean(S1,3);
    S1=S1./repmat(S1b,[size(S1,1) 1]);
    subplot(2,4,k)
    imagesc(t,f,log(S1)',[-3 3])
    axis xy
    colormap(cm)
end


subplot(2,4,1),ylabel(['elec ' int2str(electr)])
cond_label={'white n','pink n','brown n','gr 4','gr 8','gr 16','gr 32','blank'};
for k=1:8 
    subplot(2,4,k)
    title(cond_label{k})
end

figure('Color',[1 1 1],'Position',[0 0 700 700])
epochs_plot=find(all_stims==7);
data2use=squeeze(all_data_epoch(elec,:,:))';
% calculate spectgram
[S1,t_tf,f]=mtspecgramc(data2use,movingwin,params);
t_tf=t_tf-.2;
for k=1:length(epochs_plot)
% normalize wrt baseline
S1_sel=S1(:,:,epochs_plot(k));% pick epoch
S1_sel=S1_sel./repmat(S1b,[size(S1_sel,1) 1]);
subplot(5,6,k)
imagesc(t,f,log(S1_sel)',[-3 3])
axis xy
colormap(cm)
title(['gr32 ep ' int2str(epochs_plot(k))])
end

set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',['./figures/cases/' subj '_el'  int2str(electr)])
