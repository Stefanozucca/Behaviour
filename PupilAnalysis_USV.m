%% This Script is to analyse pupil signal from CSV files obtained from Bonsai Workflow and aligned to USV with ePhys signals


close all
clear all


%%First we need to load the table with all the data to be analyzed

Dataset_path='C:\Users\stefa\Desktop\Stez - Work\Torino2021_Postdoc\Projects\MSCA_2023\Behaviour\PupilTracking_USV';
Dataset_filename='Dataset_PupilTracking.xlsx';

PlotOption=0; %Put 1 if you want to plot single traces from all recordings of Pupil + USV
filterdata=1;%Put 1 if you want to filter the DATATABLE

DATATABLE=readtable(fullfile(Dataset_path,Dataset_filename));

%Filter table

%Filter based on days

%Look at first 3 days only (before the exposure to male)
if filterdata==1
    DATATABLE=DATATABLE(DATATABLE.DAY<4,:);
%     DATATABLE=DATATABLE(strcmp(DATATABLE.CYCLE,'DIESTRO')==1,:);
end

for thisfile=1:size(DATATABLE,1)
    FilePath=DATATABLE.FOLDER{thisfile};
    AnimalID=DATATABLE.ANIMALID{thisfile};
    PupilFile=DATATABLE.PUPILTRACKING{thisfile};
    SyncFile=DATATABLE.SYNCFILE{thisfile};
    CamCh=DATATABLE.CAMERACH{thisfile};
    SpeakerCh=DATATABLE.SPEAKERCH{thisfile};

    [DATATABLE.Pupil{thisfile},DATATABLE.USVSignal{thisfile}]=syncdata(FilePath,AnimalID,PupilFile,SyncFile,CamCh,SpeakerCh,PlotOption);
end

%% Analyse the Pupil tracking 

%First we analyse the z-scored pupil aligned at USV onset

PreWindow=60; %Number of seconds before stim onset
PostWindow=80; %Number of seconds after stim onset
InterStimInterval=20; %Seconds between two consecutive stimuli
PupilParameter=3; %Put 2 for pupil area and 3 for pupil diameter
repnum=20; %Number of repetition for shuffled data
HabMaxTime=300; %Duration of habituation in seconds

for thisfile=1:size(DATATABLE,1)
    pupil=DATATABLE.Pupil{thisfile};
    usv=DATATABLE.USVSignal{thisfile};
    
    %calculate USV acquisition freq
    usv_fr=round(1/(mean(diff(usv(:,1)))));
    pupil_fr=round(1/mean(diff(pupil(:,1))));

    %zscore the pupil
    pupil(:,PupilParameter)=zscore(pupil(:,PupilParameter));

    %Find USV onsets
    ttx=find(usv(:,2)>5);
    ttx2=diff(ttx)/usv_fr;
    USV_Onset=find(ttx2>InterStimInterval)+1;
    USV_Onset=[ttx(1);ttx(USV_Onset)];
    USV_Onset_time=usv(USV_Onset,1);

%     [~,USV_Onset_time]=findpeaks(usv(:,2),usv(:,1),'MinPeakDistance',0.1,'Threshold',1);

%     if PlotOption==1
%         figure
%         plot(usv(:,1),usv(:,2),'k');
%         hold on
%         plot(usv(USV_Onset,1),usv(USV_Onset,2),'r*')
%     end
    
    AlignedTime=-PreWindow:1/pupil_fr:PostWindow;
    Aligned_Pupil=zeros(length(USV_Onset_time),length(AlignedTime));
    %Aligned for each stimulus the pupil trace
    for thisstim=1:length(USV_Onset_time)
        ttx=min(find(pupil(:,1)>USV_Onset_time(thisstim)));
        
        p_start=ttx-PreWindow*pupil_fr;
        p_end=ttx+PostWindow*pupil_fr;

        Aligned_Pupil(thisstim,:)=pupil(p_start:p_end,PupilParameter);
        
    end
    
    %Store the data
    DATATABLE.Aligned_Pupil{thisfile}=Aligned_Pupil;

    %create shuffled data for comparison
    
    %We want to do 20 repetitions each
    Aligned_Pupil_shuffled=zeros(repnum,length(AlignedTime));
    for thisrep=1:repnum

        %Select time range
        st=PreWindow*pupil_fr;
        nd=(HabMaxTime-PostWindow)*pupil_fr;
        ttx=pupil(st:nd,1);
        onset=randi(numel(ttx))+PreWindow*pupil_fr;

        sf_start=onset-PreWindow*pupil_fr;
        sf_end=onset+PostWindow*pupil_fr;

        Aligned_Pupil_shuffled(thisrep,:)=pupil(sf_start:sf_end,PupilParameter);

    end

    %Store the data
    DATATABLE.Aligned_Pupil_shuffled{thisfile}=Aligned_Pupil_shuffled;

    %Calculate the area under the curve before and after the stimulus for
    %each trial



end


%% Analyse data for each animal

%First find the number of animals

AnimalList=unique(DATATABLE.ANIMALID);

for thisanimal=1:length(AnimalList)

    ttx=find(strcmp(DATATABLE.ANIMALID,AnimalList{thisanimal})==1);

    ttx2=find(strcmp(DATATABLE.ANIMALID,AnimalList{thisanimal})==1);


    All_trials_pupil=cat(1,DATATABLE.Aligned_Pupil{ttx});
    All_trials_pupil_shuffled=cat(1,DATATABLE.Aligned_Pupil_shuffled{ttx2});

    AVG=mean(All_trials_pupil,1);
    
    figure
    subplot(2,1,1)
    imagesc(AlignedTime,1:1:size(All_trials_pupil,1),All_trials_pupil);
    hold on
    plot([0 0],[1 size(All_trials_pupil,1)],'r-','LineWidth',2)

    subplot(2,1,2)
    boundedline(AlignedTime,mean(All_trials_pupil,1),sem(All_trials_pupil),'m')
    hold on
%     boundedline(AlignedTime,mean(All_trials_pupil_shuffled,1),sem(All_trials_pupil_shuffled),'k')
    ylims=ylim;
    plot([0 0],[ylims(1) ylims(2)],'r-','LineWidth',2);
    xlim([-PreWindow PostWindow]);

    
    AVG_pre=mean(All_trials_pupil(:,1:PreWindow*pupil_fr),2);
    AVG_post=mean(All_trials_pupil(:,PreWindow*pupil_fr:(PreWindow+PostWindow)*pupil_fr),2);
    
    if thisanimal==1
        All_AVG_pre=AVG_pre;
        All_AVG_post=AVG_post;

    else
        All_AVG_pre=[All_AVG_pre;AVG_pre];
        All_AVG_post=[All_AVG_post;AVG_post];
    end

end

[p,h]=signrank(All_AVG_pre,All_AVG_post);
boxplot_scatter_2(All_AVG_pre,All_AVG_post);
title(['p = ' num2str(p)])



%%

function [Pupil,USVSignal]=syncdata(FilePath,AnimalID,PupilFile,SyncFile,CamCh,SpeakerCh,PlotOption)
    
    %Load the PupilFile
    %First Column: Pupil Area
    %Second Column: Pupil Diameter
    %Third Column: Timestamps
    
    Pupil=readmatrix(fullfile(FilePath,AnimalID,'Pupil',PupilFile));

    %First we need to clean the pupil tracking file removing repetitions
    ttx=diff(Pupil(:,1)); %Look at the difference in timestamps
    %Find zeros and remove them
    Pupil(ttx==0,:)=[];    
    
    %Set at zero the starting time and have it as ms
    Pupil(:,1)=(Pupil(:,1)-Pupil(1,1))/10^9;
    %Set the fs for pupil;
    fs_pupil=round(1/mean(diff(Pupil(:,1))));
%     %Smooth the pupil signal
%     Pupil(:,2)=smooth(Pupil(:,2),'moving',120);
%     Pupil(:,3)=smooth(Pupil(:,3),'moving',120);

    Pupil(:,2)=medfilt1(Pupil(:,2),60);
    Pupil(:,3)=medfilt1(Pupil(:,3),60);

    %Now load the abf File to synch eye tracking with USV

    [data,si,h]=abf2load(fullfile(FilePath,AnimalID,'Ephys',SyncFile));

    %Find the channel for the speaker and for the USV
    find(contains(h.recChNames,CamCh));

    Camera=data(:,find(contains(h.recChNames,CamCh)));
    Speaker=data(:,find(contains(h.recChNames,SpeakerCh)));
    fs=1/(si*10^-6);
    Time=1/fs:1/fs:length(Camera)/fs;

    %Now we need to find the onset of the camera
    CamOnset=min(find(Camera>=median(Camera)));

    %Align the pupil signal with the USV signal and downsample the USV
    %signal
    USVSignal(:,1)=Time(CamOnset:end)-Time(CamOnset);
    USVSignal(:,2)=Speaker(CamOnset:end);

    %Now downsample the USV Signal
    ds_factor=floor(fs/fs_pupil);
    ds_USVSignal = downsample(USVSignal,ds_factor);

    USVSignal(USVSignal(:,1)>max(Pupil(:,1)),:)=[];

    if PlotOption==1
        figure
        plot(Pupil(:,1),Pupil(:,3),'k');hold on;yyaxis right; plot(ds_USVSignal(:,1),ds_USVSignal(:,2),'r-');
        %     plot(Pupil(:,1),Pupil(:,3),'k');hold on;yyaxis right; plot(USVSignal(:,1),USVSignal(:,2),'r-');
    end

   
end



