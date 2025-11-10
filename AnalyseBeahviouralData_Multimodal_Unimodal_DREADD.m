%Analyse behavioural data from Ethovision Tracking Output.
close all
clear all

%First get the excel file containing the list of the data
DatasetPath='C:\Users\stefa\Desktop\Stez - Work\Torino2021_Postdoc\Projects\MSCA_2023\Behaviour';
% DatasetName='Dataset_Unimodal_Multimodal.xlsx';
% DatasetName='Dataset_Unimodal_Individual.xlsx';
% DatasetName='Dataset_Unimodal_Multimodal_NeutralOdor.xlsx';
DatasetName='Dataset_Unimodal_Individual_DREADD.xlsx';
% DatasetName='Dataset_Unimodal_Multimodal_CAVcre_DREADD.xlsx';
% DatasetName='Dataset_Unimodal_Individual_CAVcre_DREADD.xlsx';
% DatasetName='Dataset_SocialPreference_CAVcre_DREADD.xlsx';
% DatasetName='Dataset_Unimodal_multimodal_DREADD.xlsx';
DATATABLE = readtable(fullfile(DatasetPath,DatasetName));

SavePath='C:\Users\stefa\Desktop\Stez - Work\Torino2021_Postdoc\Projects\MSCA_2023\Behaviour\Analysed';
SaveFileName='DREADD_Multi_vs_Uni_20251027.mat';

%Remove animals not to be included
DATATABLE(DATATABLE.TOINCLUDE==0,:)=[];

%Define the time intervals to be used for occupancy analysis
Interval=30; %Define time in seconds
MaxTime=600; %Define the maximum time
TimeIntervals=0:Interval:MaxTime;
FilterCycle=0; %Set at 1 if you want to filter dataset based on female estrous cycle
FilterCondition=0; %Set at 1 if you want to filter based on treatment condition (i.e. with dreadd vs mcherry)
Condition='MCHERRY'; %Set DREADD or MCHERRY to filter mice based on condition
ToKeep='ESTRO';
SelectNose=1; %Set at 1 if you want to use the nose point as centroid to evaluate time spent

Include_Syllables=0; %Set at 1 if you want to analyse the songs of males and align with the behaviour. This is for experiments with individual as stimuli
Audio_Duration=300; %Put here the time of the audio file used for Syllable detection (Set it equal for all animals)

%Add the option to filter based on estrous stages
if FilterCycle==1
    DATATABLE=DATATABLE(strcmp(DATATABLE.CYCLE,ToKeep)==1,:);
end

%Add the option to filter based on estrous stages
if FilterCondition==1
    DATATABLE=DATATABLE(strcmp(DATATABLE.CONDITION,Condition)==1,:);
end

%Now we cycle across animals and we load the road data

for thisanimal=1:size(DATATABLE,1)

    fprintf('\nLoading %01d of %01d ...',thisanimal,size(DATATABLE,1))

    %Get the folder path
    FolderPath=DATATABLE.FOLDER{thisanimal};

    %Get Animal Folder Name
    FolderName=DATATABLE.FILENAME{thisanimal};

    %Get the list of files for RawData and load them
    files=dir(fullfile(FolderPath,FolderName,'*.xlsx'));
    
    %Now we cycle across files to load the information from ethovision
    for thisfile=1:size(files,1)
        ttx=readtable(fullfile(files(thisfile).folder,files(thisfile).name));

        %Now we start taking all the variables
        Time=ttx.TrialTime;
        %Correct Time so all time slots start at zero
        Time=Time-Time(1);
        Tracking=ttx.XCenter;
        Tracking(:,2)=ttx.YCenter;
        Speed=ttx.Velocity;
        %In some cases the Tracking variable could be a cell instead of a
        %matrix as there might be missing values at the beginning
        if iscell(Tracking)
            Tracking=str2double(Tracking);
            Speed=str2double(Speed);
        end
        

        %Evaluate acquisition frequency
        dTime=diff(Time);
        AcquisitionFreq=round(1/mean(dTime));
        
        %Fix NaNs - There are positions where it is defined as NaNs
        %(Probably it was not tracking it correctly)

        ii=find(isnan(Tracking(:,1))==1);

        %Now Check the position of Unimodal/Multimodal speaker
        if DATATABLE.UNIMODAL(thisanimal)==1
            if SelectNose==1
                Unimodal=ttx.InZone_Speaker1_Nose_point_;
                Unimodal(:,2)=ttx.InZone_Speaker1_Zone_Nose_point_;

                Multimodal=ttx.InZone_Speaker2_Nose_point_;
                Multimodal(:,2)=ttx.InZone_Speaker2_Zone_Nose_point_;

                if ismember('InZone_Zone1_Nose_point_',ttx.Properties.VariableNames)
%                     Zone1=ttx.InZone_Zone2_Nose_point_;
                    Zone1=ttx.InZone_Zone1_Nose_point_;
                else
                    Zone1=zeros(size(Multimodal,1),1);
                end

            else

                Unimodal=ttx.InZone_Speaker1_Center_point_;
                Unimodal(:,2)=ttx.InZone_Speaker1_Zone_Center_point_;

                Multimodal=ttx.InZone_Speaker2_Center_point_;
                Multimodal(:,2)=ttx.InZone_Speaker2_Zone_Center_point_;

                if ismember('InZone_Zone1_Center_point_',ttx.Properties.VariableNames)
                    Zone1=ttx.InZone_Zone1_Center_point_;
%                     Zone1=ttx.InZone_Zone2_Center_point_;
                else
                    Zone1=zeros(size(Multimodal,1),1);
                end
            end

        else
            if SelectNose==1
                Multimodal=ttx.InZone_Speaker1_Nose_point_;
                Multimodal(:,2)=ttx.InZone_Speaker1_Zone_Nose_point_;

                Unimodal=ttx.InZone_Speaker2_Nose_point_;
                Unimodal(:,2)=ttx.InZone_Speaker2_Zone_Nose_point_;

                if ismember('InZone_Zone1_Nose_point_',ttx.Properties.VariableNames)
                    Zone1=ttx.InZone_Zone1_Nose_point_;
%                     Zone1=ttx.InZone_Zone2_Nose_point_;
                else
                    Zone1=zeros(size(Multimodal,1),1);
                end

            else
                Multimodal=ttx.InZone_Speaker1_Center_point_;
                Multimodal(:,2)=ttx.InZone_Speaker1_Zone_Center_point_;

                Unimodal=ttx.InZone_Speaker2_Center_point_;
                Unimodal(:,2)=ttx.InZone_Speaker2_Zone_Center_point_;

                if ismember('InZone_Zone1_Center_point_',ttx.Properties.VariableNames)
                    Zone1=ttx.InZone_Zone1_Center_point_;
%                     Zone1=ttx.InZone_Zone2_Center_point_;
                else
                    Zone1=zeros(size(Multimodal,1),1);
                end
            end

        end

        if iscell(Multimodal)
            Unimodal=str2double(Unimodal);
            Multimodal=str2double(Multimodal);
            Zone1=str2double(Zone1);
        end
        
        Multimodal(ii,:)=0; Multimodal(isnan(Multimodal))=0;
        Unimodal(ii,:)=0; Unimodal(isnan(Unimodal))=0;
        Zone1(ii,:)=0;

        %Now calculate for each time interval the time spent in each zone
        TimeUnimodal=zeros(length(TimeIntervals),2);
        TimeMultimodal=zeros(length(TimeIntervals),2);
        TimeZone1=zeros(length(TimeIntervals),1);
        for thistime=1:length(TimeIntervals)-1
            %Find Indeces For that time interval
            ttx=find(Time>=TimeIntervals(thistime)&Time<TimeIntervals(thistime+1));
            
            %Now calculate the time spent in each zone
            TimeUnimodal(thistime,:)=sum(Unimodal(ttx,:),1)/AcquisitionFreq;
            TimeMultimodal(thistime,:)=sum(Multimodal(ttx,:),1)/AcquisitionFreq;
%             TimeZone1(thistime,:)=sum(Zone1(ttx,:),1)/AcquisitionFreq;

        end

        %Now we add them to the DataTable depending of whether they are
        %from habituation trials or from test trials and splitting them
        %between Unimodal and Bimodal based on speaker position

        %Get the file name and check if it's habituation or test

        p=contains(files(thisfile).name,'Hab');

        if p==1 %Means it's habituation trial
            DATATABLE.ACQUISITIONFREQ{thisanimal}=AcquisitionFreq;
            DATATABLE.HAB_Time{thisanimal}=Time;
            DATATABLE.HAB_Tracking{thisanimal}=Tracking;
            DATATABLE.HAB_Speed{thisanimal}=Speed;

            DATATABLE.HAB_Unimodal{thisanimal}=Unimodal;
            DATATABLE.HAB_Multimodal{thisanimal}=Multimodal;
            DATATABLE.HAB_TimeUnimodal{thisanimal}=TimeUnimodal;
            DATATABLE.HAB_TimeMultimodal{thisanimal}=TimeMultimodal;
            DATATABLE.HAB_TimeZone{thisanimal}=TimeZone1;
        else
            DATATABLE.TEST_Time{thisanimal}=Time;
            DATATABLE.TEST_Tracking{thisanimal}=Tracking;
            DATATABLE.TEST_Speed{thisanimal}=Speed;

            DATATABLE.TEST_Unimodal{thisanimal}=Unimodal;
            DATATABLE.TEST_Multimodal{thisanimal}=Multimodal;
            DATATABLE.TEST_TimeUnimodal{thisanimal}=TimeUnimodal;
            DATATABLE.TEST_TimeMultimodal{thisanimal}=TimeMultimodal;
            DATATABLE.TEST_TimeZone{thisanimal}=TimeZone1;

        end


    end

    if Include_Syllables==1

        %%Load the file with the syllables

        files=dir(fullfile(FolderPath,FolderName,'USVs'));

        % Initialize an empty cell array to store tables
        Syllables_Table=table;
        % Loop through each file in the directory
        for i = 1:length(files)
            % Get the file name
            fileName = files(i).name;
            % Check if the file has a .csv extension
            if contains(fileName, '.csv')
                % Get the full path of the .csv file
                filePath = fullfile(FolderPath,FolderName,'USVs', fileName);

                % Load the .csv file into a table and store it in the cell array
                csvTables = readtable(filePath);


                %Check if it is the first file otherwise concatenate and correct for the
                %time
                if contains(fileName,'_1_dat')
                    Syllables_Table=csvTables;
                else
                    %Fix the time by adding the overall duration of the
                    %first audio file
                    csvTables.start=csvTables.start+Audio_Duration;
                    csvTables.xEnd=csvTables.xEnd+Audio_Duration;
                    Syllables_Table=[Syllables_Table;csvTables];
                end
            end
        end

        %Now evaluate the number of cells per Time Interval
        Syllables=zeros(length(TimeIntervals),1);
        for thisinterval=1:length(TimeIntervals)-1
            ii=find(Syllables_Table.start>=TimeIntervals(thisinterval)&Syllables_Table.start<TimeIntervals(thisinterval+1));
            if ~isempty(ii)
                Syllables(thisinterval)=length(ii);
            end
        end
        DATATABLE.Syllables{thisanimal}=Syllables;
        DATATABLE.SyllablesTime{thisanimal}=Syllables_Table.start;
        DATATABLE.SyllablesDuration{thisanimal}=Syllables_Table.duration;
        DATATABLE.SyllablesMaxFreq{thisanimal}=Syllables_Table.maxfreq;
        DATATABLE.SyllablesMaxAmp{thisanimal}=Syllables_Table.maxamp;
    end

end

%% Plot the results for each time interval for the different zones


%Evaluate the time spent in each time interval and the statistical
%difference between Unimodal and Multimodal

%First reorganize data across animals

AllUnimodal=[];
AllMultimodal=[];
AllUnimodal_Zone=[];
AllMultimodal_Zone=[];
AllZone=[];
for thisanimal=1:size(DATATABLE,1)
    iix=DATATABLE.TEST_TimeUnimodal{thisanimal};
    iix2=DATATABLE.TEST_TimeMultimodal{thisanimal};
    iix3=DATATABLE.TEST_TimeZone{thisanimal};
    AllUnimodal(thisanimal,:)=iix(:,1);
    AllUnimodal_Zone(thisanimal,:)=iix(:,2);
    AllMultimodal(thisanimal,:)=iix2(:,1);
    AllMultimodal_Zone(thisanimal,:)=iix2(:,2);

    AllZone(thisanimal,:)=iix3;
end

AllZone(isnan(AllZone))=0;

%Evaluate the cumulative time spent
AllMultimodal_cumulative=cumsum(AllMultimodal,2);
AllMultimodal_Zone_cumulative=cumsum(AllMultimodal_Zone,2);
AllUnimodal_cumulative=cumsum(AllUnimodal,2);
AllUnimodal_Zone_cumulative=cumsum(AllUnimodal_Zone,2);
AllZone_cumulative=cumsum(AllZone,2);


%Now split between DREADD and MCHERRY

DREADD_idx=find(strcmp(DATATABLE.CONDITION,'DREADD')==1);
MCHERRY_idx=find(strcmp(DATATABLE.CONDITION,'MCHERRY')==1);

AllMultimodal_cumulative_DREADD=AllMultimodal_cumulative(DREADD_idx,:);
AllMultimodal_Zone_cumulative_DREADD=AllMultimodal_Zone_cumulative(DREADD_idx,:);
AllUnimodal_cumulative_DREADD=AllUnimodal_cumulative(DREADD_idx,:);
AllUnimodal_Zone_cumulative_DREADD=AllUnimodal_Zone_cumulative(DREADD_idx,:);

AllMultimodal_cumulative_MCHERRY=AllMultimodal_cumulative(MCHERRY_idx,:);
AllMultimodal_Zone_cumulative_MCHERRY=AllMultimodal_Zone_cumulative(MCHERRY_idx,:);
AllUnimodal_cumulative_MCHERRY=AllUnimodal_cumulative(MCHERRY_idx,:);
AllUnimodal_Zone_cumulative_MCHERRY=AllUnimodal_Zone_cumulative(MCHERRY_idx,:);


%% Check for statistical Differences - 

%Start with DREADD
pVal_Time_DREADD=[];
for thistime=1:length(TimeIntervals)
    %Big Area
    Uni=AllUnimodal(DREADD_idx,thistime);
    Mul=AllMultimodal(DREADD_idx,thistime);
    pVal_Time_DREADD(thistime,1)=signrank(Uni,Mul);
    
    %Zone
    Uni=AllUnimodal_Zone(DREADD_idx,thistime);
    Mul=AllMultimodal_Zone(DREADD_idx,thistime);
    pVal_Time_DREADD(thistime,2)=signrank(Uni,Mul);

    %Cumulative Big Area
    Uni=AllUnimodal_cumulative(DREADD_idx,thistime);
    Mul=AllMultimodal_cumulative(DREADD_idx,thistime);
    pVal_Time_DREADD(thistime,3)=signrank(Uni,Mul);

    %Cumulative Zone
    Uni=AllUnimodal_Zone_cumulative(DREADD_idx,thistime);
    Mul=AllMultimodal_Zone_cumulative(DREADD_idx,thistime);
    pVal_Time_DREADD(thistime,4)=signrank(Uni,Mul);

    
end

%Repeat for MCHERRY
pVal_Time_MCHERRY=[];
for thistime=1:length(TimeIntervals)
    %Big Area
    Uni=AllUnimodal(MCHERRY_idx,thistime);
    Mul=AllMultimodal(MCHERRY_idx,thistime);
    pVal_Time_MCHERRY(thistime,1)=signrank(Uni,Mul);
    
    %Zone
    Uni=AllUnimodal_Zone(MCHERRY_idx,thistime);
    Mul=AllMultimodal_Zone(MCHERRY_idx,thistime);
    pVal_Time_MCHERRY(thistime,2)=signrank(Uni,Mul);

    %Cumulative Big Area
    Uni=AllUnimodal_cumulative(MCHERRY_idx,thistime);
    Mul=AllMultimodal_cumulative(MCHERRY_idx,thistime);
    pVal_Time_MCHERRY(thistime,3)=signrank(Uni,Mul);

    %Cumulative Zone
    Uni=AllUnimodal_Zone_cumulative(MCHERRY_idx,thistime);
    Mul=AllMultimodal_Zone_cumulative(MCHERRY_idx,thistime);
    pVal_Time_MCHERRY(thistime,4)=signrank(Uni,Mul);

    
end

%% Evaluate the preference Index (Mult-Uni)/(Mult+Uni)

for thisanimal=1:size(AllMultimodal,1)
    Mul=AllMultimodal(thisanimal,:);
    Uni=AllUnimodal(thisanimal,:);

    PrefIndex(thisanimal,:)=(Mul-Uni)./(Mul+Uni);

    Mul=AllMultimodal_Zone(thisanimal,:);
    Uni=AllUnimodal_Zone(thisanimal,:);

    PrefIndex_Zone(thisanimal,:)=(Mul-Uni)./(Mul+Uni);

    Mul=AllMultimodal_cumulative(thisanimal,:);
    Uni=AllUnimodal_cumulative(thisanimal,:);
    PrefIndex_Cumulative(thisanimal,:)=(Mul-Uni)./(Mul+Uni);

    Mul=AllMultimodal_Zone_cumulative(thisanimal,:);
    Uni=AllUnimodal_Zone_cumulative(thisanimal,:);

    PrefIndex_Zone_Cumulative(thisanimal,:)=(Mul-Uni)./(Mul+Uni);

end

PrefIndex(isnan(PrefIndex))=0;
PrefIndex_Zone(isnan(PrefIndex_Zone))=0;
PrefIndex_Cumulative(isnan(PrefIndex_Cumulative))=0;
PrefIndex_Zone_Cumulative(isnan(PrefIndex_Zone_Cumulative))=0;
%Evaluate statistical differences for preference index

for thistime=1:length(TimeIntervals)-1
    %Big Area
    pVal_PrefIndex_DREADD(thistime,1)=signrank(PrefIndex(DREADD_idx,thistime));
    pVal_PrefIndex_MCherry(thistime,1)=signrank(PrefIndex(MCHERRY_idx,thistime));

    pVal_PrefIndex_MCherry_vs_DREADD(thistime,1)=ranksum(PrefIndex(DREADD_idx,thistime),PrefIndex(MCHERRY_idx,thistime));

    %Zone
    pVal_PrefIndex_DREADD(thistime,2)=signrank(PrefIndex_Zone(DREADD_idx,thistime));
    pVal_PrefIndex_MCherry(thistime,2)=signrank(PrefIndex_Zone(MCHERRY_idx,thistime));
    pVal_PrefIndex_MCherry_vs_DREADD(thistime,2)=ranksum(PrefIndex_Zone(DREADD_idx,thistime),PrefIndex_Zone(MCHERRY_idx,thistime));

    %Cumulative Big Area
    pVal_PrefIndex_DREADD(thistime,3)=signrank(PrefIndex_Cumulative(DREADD_idx,thistime));
    pVal_PrefIndex_MCherry(thistime,3)=signrank(PrefIndex_Cumulative(MCHERRY_idx,thistime));
    pVal_PrefIndex_MCherry_vs_DREADD(thistime,3)=ranksum(PrefIndex_Cumulative(DREADD_idx,thistime),PrefIndex_Cumulative(MCHERRY_idx,thistime));


    %Zone
    pVal_PrefIndex_DREADD(thistime,4)=signrank(PrefIndex_Zone_Cumulative(DREADD_idx,thistime));
    pVal_PrefIndex_MCherry(thistime,4)=signrank(PrefIndex_Zone_Cumulative(MCHERRY_idx,thistime));
    pVal_PrefIndex_MCherry_vs_DREADD(thistime,4)=ranksum(PrefIndex_Zone_Cumulative(DREADD_idx,thistime),PrefIndex_Zone_Cumulative(MCHERRY_idx,thistime));


end

%% Plot results

%Plot data 
figure('Position',[100 100 1200 500])
subplot(2,6,1)
plot(TimeIntervals,mean(AllMultimodal(DREADD_idx,:),1),'r-');
hold on
plot(TimeIntervals,mean(AllUnimodal(DREADD_idx,:),1),'k-');
xlabel('Total Time'); ylabel('Time Spent (s)')
title('Big Zone')

subplot(2,6,2)
plot(TimeIntervals,mean(AllMultimodal_Zone(DREADD_idx,:),1),'m-');
hold on
plot(TimeIntervals,mean(AllUnimodal_Zone(DREADD_idx,:),1),'k-');
xlabel('Total Time'); ylabel('Time Spent (s)')
title('Small Zone')


subplot(2,6,3)
boundedline(TimeIntervals,mean(AllMultimodal_cumulative(DREADD_idx,:),1),nansem(AllMultimodal_cumulative(DREADD_idx,:),1),'r')
hold on
boundedline(TimeIntervals,mean(AllUnimodal_cumulative(DREADD_idx,:),1),nansem(AllUnimodal_cumulative(DREADD_idx,:),1),'k')
xlabel('Total Time'); ylabel('Cumulative Time Spent (s)')


subplot(2,6,4)
boundedline(TimeIntervals,mean(AllMultimodal_Zone_cumulative(DREADD_idx,:),1),nansem(AllMultimodal_Zone_cumulative(DREADD_idx,:),1),'m')
hold on
boundedline(TimeIntervals,mean(AllUnimodal_Zone_cumulative(DREADD_idx,:),1),nansem(AllUnimodal_Zone_cumulative(DREADD_idx,:),1),'k')
xlabel('Total Time'); ylabel('Cumulative Time Spent (s)')

subplot(2,6,4)
boundedline(TimeIntervals,mean(AllMultimodal_Zone_cumulative(DREADD_idx,:),1),nansem(AllMultimodal_Zone_cumulative(DREADD_idx,:),1),'m')
hold on
boundedline(TimeIntervals,mean(AllUnimodal_Zone_cumulative(DREADD_idx,:),1),nansem(AllUnimodal_Zone_cumulative(DREADD_idx,:),1),'k')
xlabel('Total Time'); ylabel('Cumulative Time Spent (s)')

subplot(2,6,5);
boundedline(TimeIntervals,nanmean(PrefIndex_Cumulative(DREADD_idx,:),1),nansem(PrefIndex_Cumulative(DREADD_idx,:),1),'r')
hold on
plot([TimeIntervals(1) TimeIntervals(end)],[0 0],'k--')
ylim([-1 1])
ylabel('Preference Index')
xlabel('Time Interval')
title('Big Zone')

subplot(2,6,6)
boundedline(TimeIntervals,nanmean(PrefIndex_Zone_Cumulative(DREADD_idx,:),1),nansem(PrefIndex_Zone_Cumulative(DREADD_idx,:),1),'m')
hold on
plot([TimeIntervals(1) TimeIntervals(end)],[0 0],'k--')
ylim([-1 1])
ylabel('Preference Index')
xlabel('Time Interval')
title('Small Zone')


%MCherry
subplot(2,6,7)
plot(TimeIntervals,mean(AllMultimodal(MCHERRY_idx,:),1),'r-');
hold on
plot(TimeIntervals,mean(AllUnimodal(MCHERRY_idx,:),1),'k-');
xlabel('Total Time'); ylabel('Time Spent (s)')
title('Big Zone')

subplot(2,6,8)
plot(TimeIntervals,mean(AllMultimodal_Zone(MCHERRY_idx,:),1),'m-');
hold on
plot(TimeIntervals,mean(AllUnimodal_Zone(MCHERRY_idx,:),1),'k-');
xlabel('Total Time'); ylabel('Time Spent (s)')
title('Small Zone')


subplot(2,6,9)
boundedline(TimeIntervals,mean(AllMultimodal_cumulative(MCHERRY_idx,:),1),nansem(AllMultimodal_cumulative(MCHERRY_idx,:),1),'r')
hold on
boundedline(TimeIntervals,mean(AllUnimodal_cumulative(MCHERRY_idx,:),1),nansem(AllUnimodal_cumulative(MCHERRY_idx,:),1),'k')
xlabel('Total Time'); ylabel('Cumulative Time Spent (s)')


subplot(2,6,10)
boundedline(TimeIntervals,mean(AllMultimodal_Zone_cumulative(MCHERRY_idx,:),1),nansem(AllMultimodal_Zone_cumulative(MCHERRY_idx,:),1),'m')
hold on
boundedline(TimeIntervals,mean(AllUnimodal_Zone_cumulative(MCHERRY_idx,:),1),nansem(AllUnimodal_Zone_cumulative(MCHERRY_idx,:),1),'k')
xlabel('Total Time'); ylabel('Cumulative Time Spent (s)')

subplot(2,6,11);
boundedline(TimeIntervals,nanmean(PrefIndex_Cumulative(MCHERRY_idx,:),1),nansem(PrefIndex_Cumulative(MCHERRY_idx,:),1),'r')
hold on
plot([TimeIntervals(1) TimeIntervals(end)],[0 0],'k--')
ylim([-1 1])
ylabel('Preference Index')
xlabel('Time Interval')
title('Big Zone')

subplot(2,6,12)
boundedline(TimeIntervals,nanmean(PrefIndex_Zone_Cumulative(MCHERRY_idx,:),1),nansem(PrefIndex_Zone_Cumulative(MCHERRY_idx,:),1),'m')
hold on
plot([TimeIntervals(1) TimeIntervals(end)],[0 0],'k--')
ylim([-1 1])
ylabel('Preference Index')
xlabel('Time Interval')
title('Small Zone')

% %Comparison "In Zone"
% subplot(2,6,6)
% boundedline(TimeIntervals,mean(AllZone_cumulative(DREADD_idx,:),1),nansem(AllZone_cumulative(DREADD_idx,:),1),'b','alpha')
% hold on
% boundedline(TimeIntervals,mean(AllZone_cumulative(MCHERRY_idx,:),1),nansem(AllZone_cumulative(MCHERRY_idx,:),1),'r','alpha')



%% Split and plot data based on animal cycle
%We do this only for one time interval
TimeSlot=600;
TimeInterval_idx=max(find(TimeIntervals<TimeSlot));

All_Animals_Multi_Zone_Dreadd=AllMultimodal_Zone_cumulative(DREADD_idx,TimeInterval_idx);
All_Animals_Uni_Zone_Dreadd=AllUnimodal_Zone_cumulative(DREADD_idx,TimeInterval_idx);

All_Animals_Multi_Dreadd=AllMultimodal_cumulative(DREADD_idx,TimeInterval_idx);
All_Animals_Uni_Dreadd=AllUnimodal_cumulative(DREADD_idx,TimeInterval_idx);

All_Animals_Multi_Zone_mCherry=AllMultimodal_Zone_cumulative(MCHERRY_idx,TimeInterval_idx);
All_Animals_Uni_Zone_mCherry=AllUnimodal_Zone_cumulative(MCHERRY_idx,TimeInterval_idx);

All_Animals_Multi_mCherry=AllMultimodal_cumulative(MCHERRY_idx,TimeInterval_idx);
All_Animals_Uni_mCherry=AllUnimodal_cumulative(MCHERRY_idx,TimeInterval_idx);


plotMultimodal_Unimodal(All_Animals_Uni_Zone_Dreadd,All_Animals_Multi_Zone_Dreadd);
title('DREADD Zone');ylim([0 350])

plotMultimodal_Unimodal(All_Animals_Uni_Dreadd,All_Animals_Multi_Dreadd);
title('DREADD Arm');ylim([0 350])

plotMultimodal_Unimodal(All_Animals_Uni_Zone_mCherry,All_Animals_Multi_Zone_mCherry);
title('mCherry Zone');ylim([0 350])

plotMultimodal_Unimodal(All_Animals_Uni_mCherry,All_Animals_Multi_mCherry);
title('mCherry Arm');ylim([0 350])

figure
boxplot_scatter_2(All_Animals_Uni_mCherry,All_Animals_Multi_mCherry,All_Animals_Uni_Dreadd,All_Animals_Multi_Dreadd)

figure
boxplot_scatter_2(All_Animals_Uni_Zone_mCherry,All_Animals_Multi_Zone_mCherry,All_Animals_Uni_Zone_Dreadd,All_Animals_Multi_Zone_Dreadd)


% 
% 
% 
% for thistime=1:length(TimeIntervals)
%     Var=[AllUnimodal_Zone_cumulative_MCHERRY(:,thistime);AllMultimodal_Zone_cumulative_MCHERRY(:,thistime);AllUnimodal_Zone_cumulative_DREADD(:,thistime);AllMultimodal_Zone_cumulative_DREADD(:,thistime)];
%     groups=[ones(length(MCHERRY_idx),1);2*ones(length(MCHERRY_idx),1);3*ones(length(DREADD_idx),1);4*ones(length(DREADD_idx),1) ];
%     [p(thistime),tbl,stats] = kruskalwallis(Var,groups);
%     
%     c{thistime}=multcompare(stats);

% end
% 
% for thistime=1:length(TimeIntervals)
%     Var=[AllUnimodal_cumulative_MCHERRY(:,thistime);AllMultimodal_cumulative_MCHERRY(:,thistime);AllUnimodal_cumulative_DREADD(:,thistime);AllMultimodal_cumulative_DREADD(:,thistime)];
%     groups=[ones(length(MCHERRY_idx),1);2*ones(length(MCHERRY_idx),1);3*ones(length(DREADD_idx),1);4*ones(length(DREADD_idx),1) ];
%     [p(thistime),tbl,stats] = kruskalwallis(Var,groups);
%     
%     c{thistime}=multcompare(stats);
% 
% end

%% Assess the running behaviour o the animal in test e habituation condition

RunningThr=2; %Set here the threshold for considering the animal running

Cumulative_distance=zeros(size(DATATABLE,1),2);
for thisanimal=1:size(DATATABLE,1)
    %First get the habituation speed
    Hab_Speed=DATATABLE.HAB_Speed{thisanimal};
    %Get the test speed
    Test_Speed=DATATABLE.TEST_Speed{thisanimal};
    
    %Evaluate the overall average amount of time spent running
    ttx=find(Hab_Speed>RunningThr);
    Fraction_Run_Hab(thisanimal)=length(ttx)/length(Hab_Speed);
    ttx=find(Test_Speed>RunningThr);
    Fraction_Run_Test(thisanimal)=length(ttx)/length(Test_Speed);

    %Set all stationary timepoints as NaNs
    ttx=find(Hab_Speed<=RunningThr);
    Hab_Speed(ttx)=NaN;
    ttx=find(Test_Speed<=RunningThr);
    Test_Speed(ttx)=NaN;

    %Calcualte AVG running speed per time bin
    TestTime=DATATABLE.TEST_Time{thisanimal};
    HabTime=DATATABLE.HAB_Time{thisanimal};

    for thisbin=1:length(TimeIntervals)
        if thisbin==length(TimeIntervals)
            ttx=find(HabTime>=TimeIntervals(thisbin));
            ttx2=find(TestTime>=TimeIntervals(thisbin));
        else
            ttx=find(HabTime>=TimeIntervals(thisbin)&HabTime<TimeIntervals(thisbin+1));
            ttx2=find(TestTime>=TimeIntervals(thisbin)&TestTime<TimeIntervals(thisbin+1));
        end

        AVG_HAB_Speed(thisbin,thisanimal)=mean(Hab_Speed(ttx),'omitnan');
        AVG_TEST_Speed(thisbin,thisanimal)=mean(Test_Speed(ttx2),'omitnan');
    end


    %Evaluate cumulative travelled distance during TEST
    XY=DATATABLE.TEST_Tracking{thisanimal};
    for ii=2:size(XY,1)
        distance(ii)=sqrt((XY(ii,1)-XY(ii-1,1))^2 + (XY(ii,2)-XY(ii-1,2))^2);
    end

    Cumulative_distance(thisanimal,1)=max(cumsum(distance,2));

    %Evaluate cumulative travelled distance during HAB
    XY=DATATABLE.HAB_Tracking{thisanimal};
    for ii=2:size(XY,1)
        distance(ii)=sqrt((XY(ii,1)-XY(ii-1,1))^2 + (XY(ii,2)-XY(ii-1,2))^2);
    end
    Cumulative_distance(thisanimal,2)=max(cumsum(distance,2));
end


%% Store and save analyzed data

%Generate a structure where to collect all the data

DataOutput.AllMultimodal_cumulative=AllMultimodal_cumulative;
DataOutput.AllUnimodal_cumulative=AllUnimodal_cumulative;
DataOutput.AllMultimodal_Zone_cumulative=AllMultimodal_Zone_cumulative;
DataOutput.AllMultimodal_cumulative=AllUnimodal_Zone_cumulative;

save(fullfile(SavePath,SaveFileName),'DataOutput');




%% Check Frequency

DREADDTable=DATATABLE(DREADD_idx,:);
[pValFreq_DREADD,AllFreqUni_Zone_cum_D,AllFreqMul_Zone_cum_D,All_TimeSpent_Mean_D,All_TimeSpent_D]=EvaluateFrequency(DREADDTable,TimeIntervals);

MCherryTable=DATATABLE(MCHERRY_idx,:);
[pValFreq_mCherry,AllFreqUni_Zone_cum_M,AllFreqMul_Zone_cum_M,All_TimeSpent_Mean_M,All_TimeSpent_M]=EvaluateFrequency(MCherryTable,TimeIntervals);


%% Plot Frequency Data

%Plot cumulative distribution
figure("Position",[100 100 200 250])
boundedline(TimeIntervals(1:end-1),mean(AllFreqUni_Zone_cum_M,1),nansem(AllFreqUni_Zone_cum_M,1),'k')
hold on
boundedline(TimeIntervals(1:end-1),mean(AllFreqMul_Zone_cum_M,1),nansem(AllFreqMul_Zone_cum_M,1),'m')
% xlim([0 540])
ylim([0 30])
xlabel('Total Time'); ylabel('# of visits')
set(gca,'FontName','Arial','FontSize',8,'LineWidth',1,'TickDir','out','Box','off')

%Plot cumulative distribution
figure("Position",[100 100 200 250])
boundedline(TimeIntervals(1:end-1),mean(AllFreqUni_Zone_cum_D,1),nansem(AllFreqUni_Zone_cum_D,1),'k')
hold on
boundedline(TimeIntervals(1:end-1),mean(AllFreqMul_Zone_cum_D,1),nansem(AllFreqMul_Zone_cum_D,1),'m')
% xlim([0 540])
ylim([0 30])
xlabel('Total Time'); ylabel('# of visits')
set(gca,'FontName','Arial','FontSize',8,'LineWidth',1,'TickDir','out','Box','off')

TimeBin=1:1:100;

figure
subplot(2,1,1)
histogram(cat(2,All_TimeSpent_M{:,1}),TimeBin,'FaceColor','k');
hold on
histogram(cat(2,All_TimeSpent_M{:,2}),TimeBin,'FaceColor','m');
title('mCherry')

subplot(2,1,2)
histogram(cat(2,All_TimeSpent_D{:,1}),TimeBin,'FaceColor','k');
hold on
histogram(cat(2,All_TimeSpent_D{:,2}),TimeBin,'FaceColor','m');
title('DREADD')

figure
subplot(211)
boxplot_scatter_2(All_TimeSpent_Mean_M(:,1),All_TimeSpent_Mean_M(:,2))
title('mCherry')

subplot(212)
boxplot_scatter_2(All_TimeSpent_Mean_D(:,1),All_TimeSpent_Mean_D(:,2))
return

%% Now look at the frequency (entering the zone) between Unimodal vs Multimodal


for thisanimal=1:size(DATATABLE,1)
    tUni=DATATABLE.TEST_Unimodal{thisanimal};
    tMul=DATATABLE.TEST_Multimodal{thisanimal};

    %Look at the first derivative to identify enters and exits time points

    Enter_Uni=find(diff(tUni(:,1))>0)/DATATABLE.ACQUISITIONFREQ{thisanimal};
    Enter_Uni_Zone=find(diff(tUni(:,2))>0)/DATATABLE.ACQUISITIONFREQ{thisanimal};
    Exit_Uni=find(diff(tUni(:,1))<0)/DATATABLE.ACQUISITIONFREQ{thisanimal};
    Exit_Uni_Zone=find(diff(tUni(:,2))<0)/DATATABLE.ACQUISITIONFREQ{thisanimal};

    Enter_Mul=find(diff(tMul(:,1))>0)/DATATABLE.ACQUISITIONFREQ{thisanimal};
    Enter_Mul_Zone=find(diff(tMul(:,2))>0)/DATATABLE.ACQUISITIONFREQ{thisanimal};
    Exit_Mul=find(diff(tMul(:,1))<0)/DATATABLE.ACQUISITIONFREQ{thisanimal};
    Exit_Mul_Zone=find(diff(tMul(:,2))<0)/DATATABLE.ACQUISITIONFREQ{thisanimal};


    %Store the variables in the datatable
    DATATABLE.TEST_Enter_Uni{thisanimal}=Enter_Uni;
    DATATABLE.TEST_Enter_Uni_Zone{thisanimal}=Enter_Uni_Zone;
    DATATABLE.TEST_Exit_Uni{thisanimal}=Exit_Uni;
    DATATABLE.TEST_Exit_Uni_Zone{thisanimal}=Exit_Uni_Zone;

    %Store the variables in the datatable
    DATATABLE.TEST_Enter_Mul{thisanimal}=Enter_Mul;
    DATATABLE.TEST_Enter_Mul_Zone{thisanimal}=Enter_Mul_Zone;
    DATATABLE.TEST_Exit_Mul{thisanimal}=Exit_Mul;
    DATATABLE.TEST_Exit_Mul_Zone{thisanimal}=Exit_Mul_Zone;


end

%Now evaluate and store the following measurments:
% Frequency at different time intervals
% AVG Visit Duration
% First Visit Duration
AllFreqUni=[];AllFreqUni_cum=[];AllFreqUni_Zone=[];AllFreqUni_Zone_cum=[];
AllFreqMul=[];AllFreqMul_cum=[];AllFreqMul_Zone=[];AllFreqMul_Zone_cum=[];
for thisanimal=1:size(DATATABLE,1)

    %Evaluate the entering frequency

    ii=DATATABLE.TEST_Enter_Uni{thisanimal};
    ii2=DATATABLE.TEST_Enter_Mul{thisanimal};

    iix=DATATABLE.TEST_Enter_Uni_Zone{thisanimal};
    iix2=DATATABLE.TEST_Enter_Mul_Zone{thisanimal};

    FreqUni=[];FreqMul=[];FreqUniZone=[];FreqMulZone=[];
    for thistime=1:length(TimeIntervals)-1
        ttx=find(ii>=TimeIntervals(thistime)&ii<TimeIntervals(thistime+1));
        FreqUni(thistime,1)=length(ttx);

        ttx=find(ii2>=TimeIntervals(thistime)&ii2<TimeIntervals(thistime+1));
        FreqMul(thistime,1)=length(ttx);
        
        ttx=find(iix>=TimeIntervals(thistime)&iix<TimeIntervals(thistime+1));
        FreqUniZone(thistime,1)=length(ttx);

        ttx=find(iix2>=TimeIntervals(thistime)&iix2<TimeIntervals(thistime+1));
        FreqMulZone(thistime,1)=length(ttx);

    end

        FreqUni(:,2)=cumsum(FreqUni);
        FreqMul(:,2)=cumsum(FreqMul);
        FreqUniZone(:,2)=cumsum(FreqUniZone);
        FreqMulZone(:,2)=cumsum(FreqMulZone);

    DATATABLE.TEST_FreqUni{thisanimal}=FreqUni;
    DATATABLE.TEST_FreqMul{thisanimal}=FreqMul;
    DATATABLE.TEST_FreqUniZone{thisanimal}=FreqUniZone;
    DATATABLE.TEST_FreqMulZone{thisanimal}=FreqMulZone;
    
    %Evaluate the duration of the first visit
    DATATABLE.TEST_FirstD_Uni(thisanimal)=DATATABLE.TEST_Exit_Uni{thisanimal}(1)-DATATABLE.TEST_Enter_Uni{thisanimal}(1);
    DATATABLE.TEST_FirstD_Mul(thisanimal)=DATATABLE.TEST_Exit_Mul{thisanimal}(1)-DATATABLE.TEST_Enter_Mul{thisanimal}(1);
    DATATABLE.TEST_FirstD_Uni_Zone(thisanimal)=DATATABLE.TEST_Exit_Uni_Zone{thisanimal}(1)-DATATABLE.TEST_Enter_Uni_Zone{thisanimal}(1);
    DATATABLE.TEST_FirstD_Mul_Zone(thisanimal)=DATATABLE.TEST_Exit_Mul_Zone{thisanimal}(1)-DATATABLE.TEST_Enter_Mul_Zone{thisanimal}(1);



    AllFreqUni(thisanimal,:)=FreqUni(:,1);
    AllFreqUni_cum(thisanimal,:)=FreqUni(:,2);
    AllFreqUni_Zone(thisanimal,:)=FreqUniZone(:,1);
    AllFreqUni_Zone_cum(thisanimal,:)=FreqUniZone(:,2);
    AllFreqMul(thisanimal,:)=FreqMul(:,1);
    AllFreqMul_cum(thisanimal,:)=FreqMul(:,2);
    AllFreqMul_Zone(thisanimal,:)=FreqMulZone(:,1);
    AllFreqMul_Zone_cum(thisanimal,:)=FreqMulZone(:,2);
end


for ii=1:size(AllFreqMul,2)

    pValFreq(ii,1)=signrank(AllFreqMul(:,ii),AllFreqUni(:,ii));
    pValFreq(ii,2)=signrank(AllFreqMul_Zone(:,ii),AllFreqUni_Zone(:,ii));
    pValFreq(ii,3)=signrank(AllFreqMul_cum(:,ii),AllFreqUni_cum(:,ii));
    pValFreq(ii,4)=signrank(AllFreqMul_Zone_cum(:,ii),AllFreqUni_Zone_cum(:,ii));

end


%%


for thisanimal=1:size(DATATABLE,1)
    time=DATATABLE.TEST_Time{thisanimal}; tracking=DATATABLE.TEST_Tracking{thisanimal};Syll=DATATABLE.SyllablesTime{thisanimal};
    
    pos=[];
    for thissyll=1:length(Syll)
        ttx=min(find(time>=Syll(thissyll)));

        pos(thissyll,:)=tracking(ttx,:);
    end

    figure
%     hist3([pos(:,1),pos(:,2)],'Nbins',[50 50],"EdgeColor","none",'CDataMode','auto')
    h=histogram2(pos(:,1),pos(:,2),[50 50],'DisplayStyle','tile','ShowEmptyBins','on','Normalization','probability')
    view(2)
    axis tight
    colormap jet
    caxis([0 0.1])
    colorbar

end