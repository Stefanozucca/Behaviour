%Analyse behavioural data from Ethovision Tracking Output.
close all
clear all

%First get the excel file containing the list of the data
DatasetPath='C:\Users\stefa\Desktop\Stez - Work\Torino2021_Postdoc\Projects\MSCA_2023\Behaviour';
DatasetName='Dataset_Unimodal_Multimodal.xlsx';
DATATABLE = readtable(fullfile(DatasetPath,DatasetName));

%Remove animals not to be included
DATATABLE(DATATABLE.TOINCLUDE==0,:)=[];



%Define the time intervals to be used for occupancy analysis
Interval=60; %Define time in seconds
MaxTime=600; %Define the maximum time
TimeIntervals=0:Interval:MaxTime;
FilterCycle=0; %Set at 1 if you want to filter dataset based on female estrous cycle

%Add the option to filter based on estrous stages
if FilterCycle==1
    ToKeep='DIESTRO';
    DATATABLE=DATATABLE(strcmp(DATATABLE.CYCLE,ToKeep)==1,:);
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

        %Evaluate acquisition frequency
        dTime=diff(Time);
        AcquisitionFreq=round(1/mean(dTime));
        
        %Fix NaNs - There are positions where it is defined as NaNs
        %(Probably it was not tracking it correctly)

        ii=find(isnan(Tracking(:,1))==1);

        %Now Check the position of Unimodal/Multimodal speaker
        
        if DATATABLE.UNIMODAL(thisanimal)==1
            Unimodal=ttx.InZone_Speaker1_Center_point_;
            Unimodal(:,2)=ttx.InZone_Speaker1_Zone_Center_point_;

            Multimodal=ttx.InZone_Speaker2_Center_point_;
            Multimodal(:,2)=ttx.InZone_Speaker2_Zone_Center_point_;

        else
            Multimodal=ttx.InZone_Speaker1_Center_point_;
            Multimodal(:,2)=ttx.InZone_Speaker1_Zone_Center_point_;

            Unimodal=ttx.InZone_Speaker2_Center_point_;
            Unimodal(:,2)=ttx.InZone_Speaker2_Zone_Center_point_;

        end
        
        Multimodal(ii,:)=0;
        Unimodal(ii,:)=0;

        %Now calculate for each time interval the time spent in each zone
        TimeUnimodal=zeros(length(TimeIntervals),2);
        TimeMultimodal=zeros(length(TimeIntervals),2);
        for thistime=1:length(TimeIntervals)-1
            %Find Indeces For that time interval
            ttx=find(Time>=TimeIntervals(thistime)&Time<TimeIntervals(thistime+1));
            
            %Now calculate the time spent in each zone
            TimeUnimodal(thistime,:)=sum(Unimodal(ttx,:),1)/AcquisitionFreq;
            TimeMultimodal(thistime,:)=sum(Multimodal(ttx,:),1)/AcquisitionFreq;


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
        else
            DATATABLE.TEST_Time{thisanimal}=Time;
            DATATABLE.TEST_Tracking{thisanimal}=Tracking;
            DATATABLE.TEST_Speed{thisanimal}=Speed;

            DATATABLE.TEST_Unimodal{thisanimal}=Unimodal;
            DATATABLE.TEST_Multimodal{thisanimal}=Multimodal;
            DATATABLE.TEST_TimeUnimodal{thisanimal}=TimeUnimodal;
            DATATABLE.TEST_TimeMultimodal{thisanimal}=TimeMultimodal;

        end


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
for thisanimal=1:size(DATATABLE,1)
    iix=DATATABLE.TEST_TimeUnimodal{thisanimal};
    iix2=DATATABLE.TEST_TimeMultimodal{thisanimal};
    AllUnimodal(thisanimal,:)=iix(:,1);
    AllUnimodal_Zone(thisanimal,:)=iix(:,2);
    AllMultimodal(thisanimal,:)=iix2(:,1);
    AllMultimodal_Zone(thisanimal,:)=iix2(:,2);
end


%Evaluate the cumulative time spent
AllMultimodal_cumulative=cumsum(AllMultimodal,2);
AllMultimodal_Zone_cumulative=cumsum(AllMultimodal_Zone,2);
AllUnimodal_cumulative=cumsum(AllUnimodal,2);
AllUnimodal_Zone_cumulative=cumsum(AllUnimodal_Zone,2);

figure
subplot(2,2,1)
plot(TimeIntervals,mean(AllMultimodal,1),'r-');
hold on
plot(TimeIntervals,mean(AllUnimodal,1),'k-');
xlabel('Total Time'); ylabel('Time Spent (s)')
title('Big Zone')

subplot(2,2,2)
plot(TimeIntervals,mean(AllMultimodal_Zone,1),'m-');
hold on
plot(TimeIntervals,mean(AllUnimodal_Zone,1),'k-');
xlabel('Total Time'); ylabel('Time Spent (s)')
title('Small Zone')


subplot(2,2,3)
boundedline(TimeIntervals,mean(AllMultimodal_cumulative,1),std(AllMultimodal_cumulative,1),'r','alpha')
hold on
boundedline(TimeIntervals,mean(AllUnimodal_cumulative,1),std(AllUnimodal_cumulative,1),'k','alpha')
xlabel('Total Time'); ylabel('Cumulative Time Spent (s)')


subplot(2,2,4)
boundedline(TimeIntervals,mean(AllMultimodal_Zone_cumulative,1),std(AllMultimodal_Zone_cumulative,1),'m','alpha')
hold on
boundedline(TimeIntervals,mean(AllUnimodal_Zone_cumulative,1),std(AllUnimodal_Zone_cumulative,1),'k','alpha')
xlabel('Total Time'); ylabel('Cumulative Time Spent (s)')



%% Check for statistical Differences

for thistime=1:length(TimeIntervals)
    %Big Area
    Uni=AllUnimodal(:,thistime);
    Mul=AllMultimodal(:,thistime);
    pVal_Time(thistime,1)=signrank(Uni,Mul);
    
    %Zone
    Uni=AllUnimodal_Zone(:,thistime);
    Mul=AllMultimodal_Zone(:,thistime);
    pVal_Time(thistime,2)=signrank(Uni,Mul);

    %Cumulative Big Area
    Uni=AllUnimodal_cumulative(:,thistime);
    Mul=AllMultimodal_cumulative(:,thistime);
    pVal_Time(thistime,3)=signrank(Uni,Mul);

    %Cumulative Zone
    Uni=AllUnimodal_Zone_cumulative(:,thistime);
    Mul=AllMultimodal_Zone_cumulative(:,thistime);
    pVal_Time(thistime,4)=signrank(Uni,Mul);

    
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

figure
subplot(2,2,1);
boundedline(TimeIntervals,mean(PrefIndex,1),nansem(PrefIndex,1),'r','alpha')
hold on
plot([TimeIntervals(1) TimeIntervals(end)],[0 0],'k--')
ylim([-0.5 0.5])
ylabel('Preference Index')
xlabel('Time Interval')
title('Big Zone')

subplot(2,2,2)
boundedline(TimeIntervals,mean(PrefIndex_Zone,1),nansem(PrefIndex_Zone,1),'m','alpha')
hold on
plot([TimeIntervals(1) TimeIntervals(end)],[0 0],'k--')
ylim([-0.5 0.5])
ylabel('Preference Index')
xlabel('Time Interval')
title('Small Zone')


subplot(2,2,3)
boundedline(TimeIntervals,mean(PrefIndex_Cumulative,1),nansem(PrefIndex_Cumulative,1),'r','alpha')
hold on
plot([TimeIntervals(1) TimeIntervals(end)],[0 0],'k--')
ylim([-0.5 0.5])
ylabel('Cumulative Preference Index')
xlabel('Time Interval')


subplot(2,2,4)
boundedline(TimeIntervals,mean(PrefIndex_Zone_Cumulative,1),nansem(PrefIndex_Zone_Cumulative,1),'m','alpha')
hold on
plot([TimeIntervals(1) TimeIntervals(end)],[0 0],'k--')
ylim([-0.5 0.5])
ylabel('Cumulative Preference Index')
xlabel('Time Interval')

%Evaluate statistical differences for preference index

for thistime=1:length(TimeIntervals)-1
    %Big Area
    pVal_PrefIndex(thistime,1)=signrank(PrefIndex(:,thistime));

    %Zone
    pVal_PrefIndex(thistime,2)=signrank(PrefIndex_Zone(:,thistime));

    %Cumulative Big Area
    pVal_PrefIndex(thistime,3)=signrank(PrefIndex_Cumulative(:,thistime));

    %Zone
    pVal_PrefIndex(thistime,4)=signrank(PrefIndex_Zone_Cumulative(:,thistime));

end


%% Now evaluate the preference considering only the time after having sampled both stimuli

%We just look at the 5 minutes after the sampling of both stimuli
TimeWindow=300; %Add the value in seconds of the time window you want to analyse

for thisanimal=1:size(DATATABLE,1)

    Time=DATATABLE.TEST_Time{thisanimal};
    Uni=DATATABLE.TEST_Unimodal{thisanimal};
    Mul=DATATABLE.TEST_Multimodal{thisanimal};
    
    Start=DATATABLE.TIME_SAMPLED(thisanimal);
    
    iix=find(Time>=Start&Time<=TimeWindow);

    Time_Mul(thisanimal,:)=sum(Mul(iix,:),1)/AcquisitionFreq;
    Time_Uni(thisanimal,:)=sum(Uni(iix,:),1)/AcquisitionFreq;


end


boxplot_scatter_2(Time_Uni(:,1),Time_Mul(:,1), Time_Uni(:,2),Time_Mul(:,2))
pVal(1)=signrank(Time_Uni(:,1),Time_Mul(:,1));
pVal(2)=signrank(Time_Uni(:,2),Time_Mul(:,2));

