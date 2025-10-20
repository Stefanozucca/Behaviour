%% This script is to analyse data from Elevated Plus Maze acquired using Bonsai

clear all
close all

FolderPath='E:\Behaviour\EPM'; %Define the folder where to find the CSV files
MetadataPath='C:\Users\stefa\Desktop\Stez - Work\Torino2021_Postdoc\Projects\MSCA_2023\Behaviour\DREADD'; %Define the folder where the metadata table is stored
MetadataFilename='EPM_MetaData.xlsx'; %Define the name of matadata file
TimeWindow1=600; %Set the time window in seconds
TimeWindow2=300; %Set the time window in seconds
Analyse_raw_data=0; %Set at 1 if you want to re-analyse the whole dataset (re-defining all ROIs) or 0 if you want to load analyzed data


if Analyse_raw_data==1

    SavePath='C:\Users\stefa\Desktop\Stez - Work\Torino2021_Postdoc\Projects\MSCA_2023\Behaviour\DREADD';
    FileName='EPM_Analyzed_20250728.mat';

    Datatable=readtable(fullfile(MetadataPath,MetadataFilename));

    %Cycle across animals

    for thisanimal=1:size(Datatable,1)
        %Get the bonsai .csv file for animal tracking

        tt=readtable(fullfile(FolderPath,Datatable.CSV_FILENAME{thisanimal}));

        %We need to clean the table removing duplicatez

        [~,b]=unique(tt.Item3);

        tt=tt(b,:);
        tt.Item4=tt.Item4/1000000000;

        AcquisitionFreq=round(1/mean(diff(tt.Item4)));

        %Get Centroid Infos
        XY=[];
        XY(:,1)=tt.Item1;
        XY(:,2)=tt.Item2;

        %Now get the starting frame from the datatable and set that as the
        %starting point of the animal tracking

        XY=XY(Datatable.START_FR(thisanimal):end,:);


        %Store XY tracking in the table
        Datatable.Tracking{thisanimal}=XY;


        vid=VideoReader(fullfile(FolderPath,Datatable.AVI_FILENAME{thisanimal}));

        this_frame=read(vid,60);

        %Define the ROI to evaluate animal occupancy
        a=figure;
        imagesc(this_frame);

        % Allow user to select two ROIs manually
        textHandle = text(size(this_frame,2)*0.05, size(this_frame,1)*0.05, 'Draw OpenArm_Left', 'Color', 'w', 'FontSize', 12, 'FontWeight', 'Bold');
        roi1 = drawpolygon('LineWidth',2,'Color','k'); % OpenArm_Left
        delete(textHandle);
        textHandle =text(size(this_frame,2)*0.05, size(this_frame,1)*0.05, 'Draw OpenArma_Right', 'Color', 'w', 'FontSize', 12, 'FontWeight', 'Bold');
        roi2 = drawpolygon('LineWidth',2,'Color','w'); %OpenArma_Right
        delete(textHandle);
        textHandle =text(size(this_frame,2)*0.05, size(this_frame,1)*0.05, 'Draw CloseArm_Up', 'Color', 'r', 'FontSize', 12, 'FontWeight', 'Bold');
        roi3 = drawpolygon('LineWidth',2,'Color','r'); %CloseArm_Up
        delete(textHandle);
        textHandle =text(size(this_frame,2)*0.05, size(this_frame,1)*0.05, 'Draw CloseArm_Down', 'Color', 'm', 'FontSize', 12, 'FontWeight', 'Bold');
        roi4 = drawpolygon('LineWidth',2,'Color','m'); %CloseArm_Down
        delete(textHandle);

        %Store the coordinates
        ROI_Coordinates={roi1.Position;roi2.Position;roi3.Position;roi4.Position};
        close (a)

        % Get ROI coordinates
        roi1Mask = inpolygon(XY(:,1), XY(:,2), ROI_Coordinates{1}(:,1), ROI_Coordinates{1}(:,2));
        roi2Mask = inpolygon(XY(:,1), XY(:,2), ROI_Coordinates{2}(:,1), ROI_Coordinates{2}(:,2));
        roi3Mask = inpolygon(XY(:,1), XY(:,2), ROI_Coordinates{3}(:,1), ROI_Coordinates{3}(:,2));
        roi4Mask = inpolygon(XY(:,1), XY(:,2), ROI_Coordinates{4}(:,1), ROI_Coordinates{4}(:,2));

        %Create a vector which has N rows (N = length
        %XY_downsampled coordinates) and 6 columns
        Occupancy=zeros(size(XY,1),4);
        Occupancy(roi1Mask==1,1)=1; %Set at 1 when the animal is inside the OpenArm_Left
        Occupancy(roi2Mask==1,2)=1; %Set at 1 when the animal is inside the OpenArm_Right
        Occupancy(roi3Mask==1,3)=1; %Set at 1 when the animal is inside the ClosedArm_Up
        Occupancy(roi4Mask==1,4)=1; %Set at 1 when the animal is inside the ClosedArm_Down

        %Store the variables
        Datatable.ROI_Coordinates{thisanimal}=ROI_Coordinates;
        Datatable.Occupancy{thisanimal}=Occupancy;

        %Evaluate time spent in each Arm
        %Define the number of frames to count
        FramNumb=TimeWindow1*AcquisitionFreq;
        TimeSpent_Arm=sum(Occupancy(1:FramNumb,:));

        Datatable.TimeSpent_10m{thisanimal}=TimeSpent_Arm;

        %Repeat for second time window
        FramNumb=TimeWindow2*AcquisitionFreq;
        TimeSpent_Arm=sum(Occupancy(1:FramNumb,:));

        Datatable.TimeSpent_5m{thisanimal}=TimeSpent_Arm;


    end

    save(fullfile(SavePath,FileName),'Datatable');

end

%% Evaluate time spent in Open vs Closed arms

%Load data if skip analysis from raw data
%Define load path
LoadPath='C:\Users\stefa\Desktop\Stez - Work\Torino2021_Postdoc\Projects\MSCA_2023\Behaviour\DREADD';
FileName='EPM_Analyzed_20250728.mat';

load(fullfile(LoadPath,FileName));


TimeSpent_10m=cell2mat(Datatable.TimeSpent_10m(:))/AcquisitionFreq;
TimeSpent_5m=cell2mat(Datatable.TimeSpent_5m(:))/AcquisitionFreq;

%Sum Time between the two open/closed arms
TimeSpent_10m(:,5)=TimeSpent_10m(:,1)+TimeSpent_10m(:,2); %Total Open Arms
TimeSpent_10m(:,6)=TimeSpent_10m(:,3)+TimeSpent_10m(:,4); %Total Closed Arms

TimeSpent_5m(:,5)=TimeSpent_5m(:,1)+TimeSpent_5m(:,2); %Total Open Arms
TimeSpent_5m(:,6)=TimeSpent_5m(:,3)+TimeSpent_5m(:,4); %Total Closed Arms


%Now split data based on conditions

dreadd=find(strcmp(Datatable.GROUP,'DREADD')==1);
mcherry=find(strcmp(Datatable.GROUP,'MCHERRY')==1);

figure
subplot(1,2,1)
boxplot_scatter_2(TimeSpent_10m(mcherry,5),TimeSpent_10m(mcherry,6),TimeSpent_10m(dreadd,5),TimeSpent_10m(dreadd,6))
ylabel('Time Spent (s)')

subplot(1,2,2)
boxplot_scatter_2(TimeSpent_5m(mcherry,5),TimeSpent_5m(mcherry,6),TimeSpent_5m(dreadd,5),TimeSpent_5m(dreadd,6))
ylabel('Time Spent (s)')

% Perform statistical comparison
%Concatenate data and define group

mcherry_o=TimeSpent_10m(mcherry,5); x1=ones(length(mcherry),1);
mcherry_c=TimeSpent_10m(mcherry,6); x2=2*ones(length(mcherry),1);
dreadd_o=TimeSpent_10m(dreadd,5); x3=3*ones(length(dreadd),1);
dreadd_c=TimeSpent_10m(dreadd,6); x4=4*ones(length(dreadd),1);

To_Comp_10m=[mcherry_o;mcherry_c;dreadd_o;dreadd_c];
Group=[x1;x2;x3;x4];

[p,tbl,stats] = kruskalwallis(To_Comp_10m,Group);

[c,m,h] = multcompare(stats);

