%Analyse behavioural data from Ethovision Tracking Output.
close all
clear all

%First get the excel file containing the list of the data
DatasetPath='C:\Users\stefa\Desktop\Stez - Work\Torino2021_Postdoc\Projects\MSCA_2023\Behaviour';
% DatasetName='Dataset_Unimodal_Multimodal.xlsx';
DatasetName='Dataset_Unimodal_Individual.xlsx';
DATATABLE = readtable(fullfile(DatasetPath,DatasetName));

%Remove animals not to be included
DATATABLE(DATATABLE.TOINCLUDE==0,:)=[];



%Define the time intervals to be used for occupancy analysis
Interval=60; %Define time in seconds
MaxTime=600; %Define the maximum time
TimeIntervals=0:Interval:MaxTime;
FilterCycle=0; %Set at 1 if you want to filter dataset based on female estrous cycle

Include_Syllables=1; %Set at 1 if you want to analyse the songs of males and align with the behaviour. This is for experiments with individual as stimuli
Audio_Duration=300; %Put here the time of the audio file used for Syllable detection (Set it equal for all animals)

%Add the option to filter based on estrous stages
if FilterCycle==1
    ToKeep='ESTRO';
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
boundedline(TimeIntervals,mean(AllMultimodal_cumulative,1),nansem(AllMultimodal_cumulative,1),'r','alpha')
hold on
boundedline(TimeIntervals,mean(AllUnimodal_cumulative,1),nansem(AllUnimodal_cumulative,1),'k','alpha')
xlabel('Total Time'); ylabel('Cumulative Time Spent (s)')


subplot(2,2,4)
boundedline(TimeIntervals,mean(AllMultimodal_Zone_cumulative,1),nansem(AllMultimodal_Zone_cumulative,1),'m','alpha')
hold on
boundedline(TimeIntervals,mean(AllUnimodal_Zone_cumulative,1),nansem(AllUnimodal_Zone_cumulative,1),'k','alpha')
xlabel('Total Time'); ylabel('Cumulative Time Spent (s)')

%Evaluate Mean and Median time for Multimodal and Unimodal
AVG_Time_BigZone_Cumulative(:,1)=mean(AllMultimodal_cumulative,1);
AVG_Time_BigZone_Cumulative(:,2)=mean(AllUnimodal_cumulative,1);

AVG_Time_Zone_Cumulative(:,1)=mean(AllMultimodal_Zone_cumulative,1);
AVG_Time_Zone_Cumulative(:,2)=mean(AllUnimodal_Zone_cumulative,1);

MED_Time_BigZone_Cumulative(:,1)=median(AllMultimodal_cumulative,1);
MED_Time_BigZone_Cumulative(:,2)=median(AllUnimodal_cumulative,1);

MED_Time_Zone_Cumulative(:,1)=median(AllMultimodal_Zone_cumulative,1);
MED_Time_Zone_Cumulative(:,2)=median(AllUnimodal_Zone_cumulative,1);

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
boundedline(TimeIntervals,nanmean(PrefIndex,1),nansem(PrefIndex,1),'r','alpha')
hold on
plot([TimeIntervals(1) TimeIntervals(end)],[0 0],'k--')
ylim([-1 1])
ylabel('Preference Index')
xlabel('Time Interval')
title('Big Zone')

subplot(2,2,2)
boundedline(TimeIntervals,nanmean(PrefIndex_Zone,1),nansem(PrefIndex_Zone,1),'m','alpha')
hold on
plot([TimeIntervals(1) TimeIntervals(end)],[0 0],'k--')
ylim([-1 1])
ylabel('Preference Index')
xlabel('Time Interval')
title('Small Zone')


subplot(2,2,3)
boundedline(TimeIntervals,nanmean(PrefIndex_Cumulative,1),nansem(PrefIndex_Cumulative,1),'r','alpha')
hold on
plot([TimeIntervals(1) TimeIntervals(end)],[0 0],'k--')
ylim([-1 1])
ylabel('Cumulative Preference Index')
xlabel('Time Interval')


subplot(2,2,4)
boundedline(TimeIntervals,nanmean(PrefIndex_Zone_Cumulative,1),nansem(PrefIndex_Zone_Cumulative,1),'m','alpha')
hold on
plot([TimeIntervals(1) TimeIntervals(end)],[0 0],'k--')
ylim([-1 1])
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

%% Split and plot data based on animal cycle
%We do this only for one time interval
TimeSlot=300;
TimeInterval_idx=max(find(TimeIntervals<TimeSlot));

Estro_Multi=AllMultimodal_Zone_cumulative(strcmp(DATATABLE.CYCLE,'ESTRO'),TimeInterval_idx);
Estro_Uni=AllUnimodal_Zone_cumulative(strcmp(DATATABLE.CYCLE,'ESTRO'),TimeInterval_idx);

Diestro_Multi=AllMultimodal_Zone_cumulative(strcmp(DATATABLE.CYCLE,'DIESTRO'),TimeInterval_idx);
Diestro_Uni=AllUnimodal_Zone_cumulative(strcmp(DATATABLE.CYCLE,'DIESTRO'),TimeInterval_idx);

%Check if there is a significant difference
pVal_Estro=signrank(Estro_Uni,Estro_Multi);
pVal_Diestro=signrank(Diestro_Uni,Diestro_Multi);

plotMultimodal_Unimodal(Estro_Uni, Estro_Multi);

plotMultimodal_Unimodal(Diestro_Uni,Diestro_Multi);

%% Now, if exposed to individuals, evaluate the contribution of ultrasound vocalizations for preference

%First set condition to run this part only if USVs are present
if Include_Syllables==1

    %For each animal evaluate the # of syllables along the entire
    %behavioural test

    ii=DATATABLE.Syllables(:);
    Syllables_intime=cat(2,ii{:});
    Syllables_Fraction=Syllables_intime./sum(Syllables_intime,1);
    figure
    subplot(2,1,1)
    boundedline(TimeIntervals,nanmean(Syllables_intime,2),nansem(Syllables_intime,2));
    subplot(2,1,2)
    boundedline(TimeIntervals,nanmean(Syllables_Fraction,2),nansem(Syllables_Fraction,2));

    
    %Evaluate the number cumulative number of syllables
    Syllables_time_cumulative=cumsum(Syllables_intime,1);
    Syllables_Fraction_cumulative=cumsum(Syllables_Fraction,1);
    
    %Calculate the correlation coefficient between the time spent and the
    %syllable number for each time bin in the cumulative distribution
    USVs_Time_Cumulative_CC=zeros(length(TimeIntervals),2);

    figure
    %     plotsize=ceil(sqrt(length(TimeIntervals)));
    for ii=1:length(TimeIntervals)
        
%         var1=Syllables_time_cumulative(ii,:);
%         var2=AllMultimodal_Zone_cumulative(:,ii);
        var1=Syllables_time_cumulative(ii,:);
        var2=AllMultimodal_Zone_cumulative(:,ii);

        [a,p]=corrcoef(var1,var2);

        USVs_Time_Cumulative_CC(ii,1)=a(2,1);
        USVs_Time_Cumulative_CC(ii,2)=p(2,1);

        mdl = fitlm(var1,var2);

        subplot(6,11,ii)
        %         scatter(Syllables_time_cumulative(ii,:),AllMultimodal_Zone_cumulative(:,ii));
        plot(mdl)
        lgd = findobj('type', 'legend');
        set(lgd, 'visible', 'off')
        if p(2,1)<=0.05
        title(['CC = ' num2str(USVs_Time_Cumulative_CC(ii,1))]);
        end
    end


end


return
%% Now evaluate the preference considering only the time after having sampled both stimuli

%We just look at the 5 minutes after the sampling of both stimuli
TimeWindow=30; %Add the value in seconds of the time window you want to analyse

for thisanimal=1:size(DATATABLE,1)

    Time=DATATABLE.TEST_Time{thisanimal};
    Uni=DATATABLE.TEST_Unimodal{thisanimal};
    Mul=DATATABLE.TEST_Multimodal{thisanimal};
    
Start=DATATABLE.TIME_SAMPLED_OFF(thisanimal);
    
    iix=find(Time>=Start&Time<=(Start+TimeWindow));

    Time_Mul(thisanimal,:)=sum(Mul(iix,:),1)/AcquisitionFreq;
    Time_Uni(thisanimal,:)=sum(Uni(iix,:),1)/AcquisitionFreq;


end



pVal(1)=signrank(Time_Uni(:,1),Time_Mul(:,1));
pVal(2)=signrank(Time_Uni(:,2),Time_Mul(:,2));
boxplot_scatter_2(Time_Uni(:,1),Time_Mul(:,1), Time_Uni(:,2),Time_Mul(:,2))
title([num2str(pVal(1)) ' - ' num2str(pVal(2))])

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