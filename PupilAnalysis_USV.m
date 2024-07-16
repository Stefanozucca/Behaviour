%% This Script is to analyse pupil signal from CSV files obtained from Bonsai Workflow and aligned to USV with ePhys signals

%%First we need to load the table with all the data to be analyzed

Dataset_path='C:\Users\stefa\Desktop\Stez - Work\Torino2021_Postdoc\Projects\MSCA_2023\Behaviour\PupilTracking_USV';
Dataset_filename='Dataset_PupilTracking.xlsx';

DATATABLE=readtable(fullfile(Dataset_path,Dataset_filename));


for thisfile=1:size(DATATABLE,1)
    FilePath=DATATABLE.FOLDER{thisfile};
    AnimalID=DATATABLE.ANIMALID{thisfile};
    PupilFile=DATATABLE.PUPILTRACKING{thisfile};
    SyncFile=DATATABLE.SYNCFILE{thisfile};
    CamCh=DATATABLE.CAMERACH{thisfile};
    SpeakerCh=DATATABLE.SPEAKERCH{thisfile};

    syncdata(FilePath,AnimalID,PupilFile,SyncFile,CamCh,SpeakerCh)
end



function syncdata(FilePath,AnimalID,PupilFile,SyncFile,CamCh,SpeakerCh)
    
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
    %Smooth the pupil signal
    Pupil(:,3)=smooth(Pupil(:,3),'moving',60);

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
    
    figure
%     plot(Pupil(:,1),Pupil(:,3),'k');hold on;yyaxis right; plot(ds_USVSignal(:,1),ds_USVSignal(:,2),'r-');
    plot(Pupil(:,1),Pupil(:,3),'k');hold on;yyaxis right; plot(USVSignal(:,1),USVSignal(:,2),'r-');


   
end

