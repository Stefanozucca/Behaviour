
%% Organise and load all the raw data

%% Load MiniscopeData

%Load Processed data
FolderPath='C:\Users\stefa\Desktop\Stez - Work\Torino2021_Postdoc\Projects\MSCA_2023\Miniscope\AN621\20240827\Miniscope';
FileName='Uncompressed\AN621_Miniscope2024-08-27T16_02_29-1b_data_processed.mat';
MiniscopeData=load(fullfile(FolderPath,FileName));


%Load csv file
FileName_csv='AN621_Miniscope_2024-08-27T16_02_29.csv';
MiniscopeTable=readtable(fullfile(FolderPath,FileName_csv));
%Remove columns not useful for the analysis
MiniscopeTable(:,3:end)=[];
%Remove the first frame as it is black
MiniscopeTable(1,:)=[];

%Now add the dFF data from the processed file
MiniscopeTable.dFF=MiniscopeData.dff';
% Fsub=F-0.95*Fneu;
% dff=(Fsub(tt,:)-mean(Fsub(tt,:),2))./mean(Fsub(tt,:),2);
% MiniscopeTable.dFF_suite2=dff';


%Rename the columns names
MiniscopeTable.Properties.VariableNames={'Time','FrameID','dFF'};

%Set the time variable
MiniscopeTable.Time=(MiniscopeTable.Time-MiniscopeTable.Time(1))/1000;


%% Load Pupil Tracking data

FolderPath='C:\Users\stefa\Desktop\Stez - Work\Torino2021_Postdoc\Projects\MSCA_2023\Miniscope\AN621\20240827\Pupil';
FileName_csv='AN621_2024-08-27T16_02_29.csv';

%load pupil infos
PupilTable=readtable(fullfile(FolderPath,FileName_csv));

%Rename variable names
PupilTable.Properties.VariableNames={'FrameID','Time','Area','Diameter','MiniscopeFrameID'};

%Fix the time
PupilTable.Time=(PupilTable.Time-PupilTable.Time(1))/1000000000;

%Now remove duplicates
PupilTable(diff(PupilTable.FrameID)==0,:)=[];

Camera_FR=1/median(unique(diff(PupilTable.Time)));

%Now clean the Pupil Signal
PupilTable.Area=medfilt1(PupilTable.Area,60);
PupilTable.Diameter=medfilt1(PupilTable.Diameter,60);

%% Load ePhys data to synchronize all files
FolderPath='C:\Users\stefa\Desktop\Stez - Work\Torino2021_Postdoc\Projects\MSCA_2023\Miniscope\AN621\20240827\Ephys';
FileName='2024_08_27_0007.abf';
CamCh='IN 5';
OdorCh='IN 7';
USVCh='IN 2';
ISI_Odor=20; %InterStimulsInterval for odors in seconds
Odor_duration=10; %Duration of odor in seconds


[data,si,h]=abf2load(fullfile(FolderPath,FileName));

Camera=squeeze(data(:,find(contains(h.recChNames,CamCh)),:)); Camera=cat(1,Camera(:));
Speaker=squeeze(data(:,find(contains(h.recChNames,USVCh)),:)); Speaker=cat(1,Speaker(:));
Odor=squeeze(data(:,find(contains(h.recChNames,OdorCh)),:)); Odor = cat(1,Odor(:));

fs=1/(si*10^-6);
Time=1/fs:1/fs:length(Camera)/fs;

%Now find camera onset to align data with Pupil and Miniscope
thr=mean(Camera)-3*std(Camera);
CamOnset=min(find(Camera<=thr));

%Find Camera Onset in time
CamOnset=CamOnset/fs;

% We have to adjust the onset based on the number of skipped frames at the
% beginning of the camera recording with bonsai
SkipTime=(PupilTable.FrameID(1)-1)/Camera_FR;
CamOnset=CamOnset+SkipTime;

%Now find Odor Onsets
thr=mean(Odor)+20*std(Odor);
[ii,OdorOnset]=findpeaks(Odor,'MinPeakHeight',thr,'MinPeakDistance',ISI_Odor*fs);
%From frames to time
OdorOnset=OdorOnset/fs;
%Adjust with Camera Onset
OdorOnset=OdorOnset-CamOnset;
OdorOffset=OdorOnset+Odor_duration;

%Adjust Miniscope Time based on 

%% Display Data

%For each unit plot the dFF and the Pupil trace
for thisunit=1:size(MiniscopeTable.dFF,2)
    dff=MiniscopeTable.dFF(:,thisunit);

    figure
    plot(PupilTable.Time,PupilTable.Area,'k');
    hold on
    ylimit=ylim;
    for ii=1:length(OdorOnset)
        plot([OdorOnset(ii) OdorOnset(ii)],[ylimit(1) ylimit(2)],'r--')
    end

    yyaxis right

    plot(MiniscopeTable.Time,dff,'m');

    xlim([0 200])
end


%% Functions

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