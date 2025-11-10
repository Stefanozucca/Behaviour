function [pValFreq,AllFreqUni_Zone_cum,AllFreqMul_Zone_cum,All_TimeSpent_Mean,All_TimeSpent]=EvaluateFrequency(DATATABLE,TimeIntervals)

MinDuration=1;

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

    %Remove short trips due to tracking instability
    tt=find(diff(Enter_Uni)<=MinDuration);Enter_Uni(tt+1)=[];
    tt=find(diff(Enter_Uni_Zone)<=MinDuration);Enter_Uni_Zone(tt+1)=[];
    tt=find(diff(Enter_Mul)<=MinDuration);Enter_Mul(tt+1)=[];
    tt=find(diff(Enter_Mul_Zone)<=MinDuration);Enter_Mul_Zone(tt+1)=[];
    
    %Evaluate mean duration for each event.
    TimeSpent=[];
    for thisenter=1:length(Enter_Uni_Zone)-1
        ttx=min(find(Exit_Uni_Zone>=Enter_Uni_Zone(thisenter+1)));
        if ~isempty(ttx)
            TimeSpent(thisenter)=Exit_Uni_Zone(ttx)-Enter_Uni_Zone(thisenter);
        else
            TimeSpent(thisenter)=NaN;
        end
    end
    All_TimeSpent_Mean(thisanimal,1)=median(TimeSpent,'omitnan');
    All_TimeSpent{thisanimal,1}=TimeSpent;

    TimeSpent=[];
    for thisenter=1:length(Enter_Mul_Zone)-1
        ttx=min(find(Exit_Mul_Zone>=Enter_Mul_Zone(thisenter+1)));
        if ~isempty(ttx)
            TimeSpent(thisenter)=Exit_Mul_Zone(ttx)-Enter_Mul_Zone(thisenter);
        else
            TimeSpent(thisenter)=NaN;
        end
    end
    All_TimeSpent_Mean(thisanimal,2)=median(TimeSpent,'omitnan');
    All_TimeSpent{thisanimal,2}=TimeSpent;

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

%     %Evaluate the duration of the first visit
%     DATATABLE.TEST_FirstD_Uni(thisanimal)=DATATABLE.TEST_Exit_Uni{thisanimal}(1)-DATATABLE.TEST_Enter_Uni{thisanimal}(1);
%     DATATABLE.TEST_FirstD_Mul(thisanimal)=DATATABLE.TEST_Exit_Mul{thisanimal}(1)-DATATABLE.TEST_Enter_Mul{thisanimal}(1);
%     DATATABLE.TEST_FirstD_Uni_Zone(thisanimal)=DATATABLE.TEST_Exit_Uni_Zone{thisanimal}(1)-DATATABLE.TEST_Enter_Uni_Zone{thisanimal}(1);
%     DATATABLE.TEST_FirstD_Mul_Zone(thisanimal)=DATATABLE.TEST_Exit_Mul_Zone{thisanimal}(1)-DATATABLE.TEST_Enter_Mul_Zone{thisanimal}(1);
% 


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
end