%% k-means-PART9 How many principal states are necessary for InFoRM Rivalry using Gratings?
% Aim: Superimpose Gaussian 6 state model with k means individual clusters
% JSK Feb 2024
clear all;
load('Analysis_ReplayMe_joystickCBR.mat', 'data'); % load data generated in PART 6

data_REPLAY =data;
clear data
load('Analysis_kmeans_PART5_joystickCBR.mat', 'data'); % load data generated in PART 6



load('Analysis_kmeans_PART5_joystickCBR.mat', 'data'); % load data generated in PART 6

close all;
n=1;
numtrials=size(data.LowContrast.Followme(1).jstkFBratio);numtrials=numtrials(1);

%% load Gaussian models and superimpose it with k means cluster during Replay Me

% % high contrast
% kmeans_in_gaussianmap_counter_highcont_Xsubj=[];group_val=[];
% for subj = 1:length(data_REPLAY.HighContrast.Replayme)%run through each subject
%     %1 step) find the optimal number of clusters
%     optimalNumcluster=round(data.HighContrast.Rivalme(1).Xsubjects.kmeans_allminX(subj),0); % these are the optimal numbers of clusters for this participant / condition
%     %2 step) run k means for that number of clusters
%     kmeans_in_gaussianmap_counter_Xtrials=[];
%     for trial=1:numtrials %run through each trial
%              kmeans_in_gaussianmap_counter(1:30)=nan; %create even length between trials/subjects as k means varies from 2-25
%         %3 step) load the Rival Me datapoint classification for each trial
%         dataPointClassification_Replayme= data_REPLAY.HighContrast.Replayme(subj).classification_Joystick(trial,:)';
%         %4 step) run k means for the optimal cluster number
%         data.HighContrast.Replayme(subj).XY(trial).matrixXY=[data.HighContrast.Replayme(subj).jstkLRratio(trial,:);data.HighContrast.Replayme(subj).jstkFBratio(trial,:) ];%create XY array
%         X=data.HighContrast.Replayme(subj).XY(trial).matrixXY;
%         X=X/2;% convert values tomrange from 0-1
%         X=X+0.5;
%         %scatter(X(1,:),X(2,:))
%         [idx,C,sumdist] = kmeans(X',optimalNumcluster,'Distance','cityblock','Display','final');
%         %5 step) run through each centroid and assign Replayme classified data to the k means data
%         for i=1:length(C)
%             centroid_temp =  X'==C(i,:); % find centroid:  logicals 1 and 1
%             temp =   dataPointClassification_Replayme(centroid_temp(:,1)==1 & centroid_temp(:,2)==1);%shows all the coordinates on which the centroid lands
%             if isempty(temp)==1 % if there is no direct corresponding values then find nearest neighbor
%                 tempC=C(i,:);
%                 for dataPoint=1:length(dataPointClassification_Replayme) % work though each data point
%                 [xLocs, yLocs]=meshgrid(linspace(0,1),linspace(0,1)); % grid for classified space
%                  distFromEachPixel=sqrt((X(1,dataPoint)-tempC(1)).^2+(X(2,dataPoint)-tempC(2)).^2); % calculate distance from each classified pixel
%                 %distFromEachPixel=sqrt((X(1, tempC(1))-xLocs).^2+(X(2, tempC(1))-yLocs).^2); % calculate distance from each classified pixel
%                 [~,closestPixel]=min(distFromEachPixel(:)); % find index of closest classified pixel
%                 temp_kmeans_in_gaussianmap_counter(dataPoint)=data.HighContrast.Followme(subj).gaussClass(closestPixel); % assign value of closest classified pixel to new data
%                 kmeans_in_gaussianmap_counter(i)=temp_kmeans_in_gaussianmap_counter(1);%pulls out just single data as all other represent the same gaussian map state
%                 end
%             else
%                 kmeans_in_gaussianmap_counter(i)=temp(1); %pulls out just single data as all other represent the same gaussian map state
%             end
%         end
% 
%         %6 step) combine data across trials
%         kmeans_in_gaussianmap_counter_Xtrials=[kmeans_in_gaussianmap_counter_Xtrials kmeans_in_gaussianmap_counter];
%     end % end of trial loop
%     % 7 step) combine across subject
%     kmeans_in_gaussianmap_counter_highcont_Xsubj=[kmeans_in_gaussianmap_counter_highcont_Xsubj ;kmeans_in_gaussianmap_counter_Xtrials];
%     group_val=[group_val; subj*ones(length(kmeans_in_gaussianmap_counter_Xtrials),1)'];
% end % end of subject loop
% figure
% % In which of the 6 gaussian states did the optimal k means cluster occur?
% swarmchart(kmeans_in_gaussianmap_counter_highcont_Xsubj',group_val','b','filled','MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.75)
% title('Replay: # of K-means clusters in each perceptual state map','FontSize',12);
% subtitle('High contrast condition','FontSize',10);
% set(gca,'XTickLabel',{'','Left','Right','Piecemeal','Equal Super','Left Super','Right Super',''});
% ylabel('Participants','FontSize',12 ,'FontWeight','bold' );
% set(gca,'TickDir','out');
% box off 
% % how many k means within each state per participant?
% for subj=1: height(kmeans_in_gaussianmap_counter_highcont_Xsubj)
% temp_EVL(subj,1)=sum(kmeans_in_gaussianmap_counter_highcont_Xsubj(subj,:)==1);
% temp_EVR(subj,1)=sum(kmeans_in_gaussianmap_counter_highcont_Xsubj(subj,:)==2);
% temp_PM(subj,1) =sum(kmeans_in_gaussianmap_counter_highcont_Xsubj(subj,:)==3);
% temp_ESI(subj,1)=sum(kmeans_in_gaussianmap_counter_highcont_Xsubj(subj,:)==4);
% temp_SIL(subj,1)=sum(kmeans_in_gaussianmap_counter_highcont_Xsubj(subj,:)==5);
% temp_SIR(subj,1)=sum(kmeans_in_gaussianmap_counter_highcont_Xsubj(subj,:)==6);
% end
% highcont_countkmeans=[temp_EVL, temp_EVR, temp_PM, temp_ESI, temp_SIL ,temp_SIR];
% 
% % plot circles with colorbar as variable
% %x 1-6 states
% %y 1 circle per participant/states
% % color map number of centroids per participants and state
% x_temp=[1:size(highcont_countkmeans,2)];
% X_temp=repmat(x_temp,height(highcont_countkmeans),1)
% y_temp=[1:28];
% Y_temp=repmat(y_temp,length(x_temp),1);
% Y_temp=Y_temp';
% mean_Xsubj_highcont_countkmeans=mean(highcont_countkmeans,2);
% mean_Xstates_highcont_countkmeans=mean(highcont_countkmeans,1);
% 
% highcont_countkmeans_withmean=[highcont_countkmeans;nan nan nan nan nan nan; round(mean_Xstates_highcont_countkmeans)]
% figure
% s=scatter(X_temp(:),Y_temp(:),168,highcont_countkmeans(:),"filled");
% xlim([0 7])
% ylim([0 30])
% % Add a colorbar
% c=colorbar
% c.Label.String = '# of K-means clusters';
% c.Label.FontWeight="bold"
% c.Label.FontSize=12
% 
% colormap turbo
% title('Replay: # of K-means clusters in each perceptual state map','FontSize',12);
% subtitle('High contrast condition','FontSize',10);
% set(gca,'XTickLabel',{'','Left','Right','Piecemeal','Equal Super','Left Super','Right Super',''});
% ylabel('Participants','FontSize',12 ,'FontWeight','bold' );
% xlabel('Perceptual state','FontSize',12 ,'FontWeight','bold' );
% 
% set(gca,'TickDir','out');
% ax = gca; %change tick fontsize
% ax.FontSize = 12;
% box off 
% 
% %low vs high
% kmeans_in_gaussianmap_counter_lowvshighcont_Xsubj=[];group_val=[];
% for subj = 1:length(data.LowvsHighContrast.Replayme)%run through each subject
%     %1 step) find the optimal number of clusters
%     optimalNumcluster=round(data.LowvsHighContrast.Rivalme(1).Xsubjects.kmeans_allminX(subj),0); % these are the optimal numbers of clusters for this participant / condition
%     %2 step) run k means for that number of clusters
%     kmeans_in_gaussianmap_counter_Xtrials=[];
%     for trial=1:numtrials %run through each trial
%              kmeans_in_gaussianmap_counter(1:30)=nan; %create even length between trials/subjects as k means varies from 2-25
%            %3 step) load the Rival Me datapoint classification for each trial
%         dataPointClassification_Replayme= data_REPLAY.LowvsHighContrast.Replayme(subj).classification_Joystick(trial,:)';
%         %4 step) run k means for the optimal cluster number
%         data.LowvsHighContrast.Replayme(subj).XY(trial).matrixXY=[data.LowvsHighContrast.Replayme(subj).jstkLRratio(trial,:);data.LowvsHighContrast.Replayme(subj).jstkFBratio(trial,:) ];%create XY array
%         X=data.LowvsHighContrast.Replayme(subj).XY(trial).matrixXY;
%         X=X/2;% convert values tomrange from 0-1
%         X=X+0.5;
%         %scatter(X(1,:),X(2,:))
%         [idx,C,sumdist] = kmeans(X',optimalNumcluster,'Distance','cityblock','Display','final');
%         %5 step) run through each centroid and assign Replayme classified data to the k means data
%         for i=1:length(C)
%             centroid_temp =  X'==C(i,:); % find centroid:  logicals 1 and 1
%             temp =   dataPointClassification_Replayme(centroid_temp(:,1)==1 & centroid_temp(:,2)==1);%shows all the coordinates on which the centroid lands
%             if isempty(temp)==1 % if there is no direct corresponding values then find nearest neighbor
%                 tempC=C(i,:);
%                 for dataPoint=1:length(dataPointClassification_Replayme) % work though each data point
%                 [xLocs, yLocs]=meshgrid(linspace(0,1),linspace(0,1)); % grid for classified space
%                  distFromEachPixel=sqrt((X(1,dataPoint)-tempC(1)).^2+(X(2,dataPoint)-tempC(2)).^2); % calculate distance from each classified pixel
%                 %distFromEachPixel=sqrt((X(1, tempC(1))-xLocs).^2+(X(2, tempC(1))-yLocs).^2); % calculate distance from each classified pixel
%                 [~,closestPixel]=min(distFromEachPixel(:)); % find index of closest classified pixel
%                 temp_kmeans_in_gaussianmap_counter(dataPoint)=data.LowvsHighContrast.Followme(subj).gaussClass(closestPixel); % assign value of closest classified pixel to new data
%                 kmeans_in_gaussianmap_counter(i)=temp_kmeans_in_gaussianmap_counter(1);%pulls out just single data as all other represent the same gaussian map state
%                 end
%             else
%                 kmeans_in_gaussianmap_counter(i)=temp(1); %pulls out just single data as all other represent the same gaussian map state
%             end
%         end
% 
%         %6 step) combine data across trials
%         kmeans_in_gaussianmap_counter_Xtrials=[kmeans_in_gaussianmap_counter_Xtrials kmeans_in_gaussianmap_counter];
%     end % end of trial loop
%     % 7 step) combine across subject
%     kmeans_in_gaussianmap_counter_lowvshighcont_Xsubj=[kmeans_in_gaussianmap_counter_lowvshighcont_Xsubj ;kmeans_in_gaussianmap_counter_Xtrials];
%     group_val=[group_val; subj*ones(length(kmeans_in_gaussianmap_counter_Xtrials),1)'];
% end % end of subject loop
% figure
% % In which of the 6 gaussian states did the optimal k means cluster occur?
% swarmchart(kmeans_in_gaussianmap_counter_lowvshighcont_Xsubj',group_val','b','filled','MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.75)
% title('Replay: # of K-means clusters in each perceptual state map','FontSize',12);
% subtitle('Low vs high contrast condition','FontSize',10);
% set(gca,'XTickLabel',{'','Left','Right','Piecemeal','Equal Super','Left Super','Right Super',''});
% ylabel('Participants','FontSize',12 ,'FontWeight','bold' );
% set(gca,'TickDir','out');
% box off 
% xlabel('Perceptual state','FontSize',12 ,'FontWeight','bold' );
% 
% % how many k means within each state per participant?
% for subj=1: height(kmeans_in_gaussianmap_counter_lowvshighcont_Xsubj)
% temp_EVL(subj,1)=sum(kmeans_in_gaussianmap_counter_lowvshighcont_Xsubj(subj,:)==1);
% temp_EVR(subj,1)=sum(kmeans_in_gaussianmap_counter_lowvshighcont_Xsubj(subj,:)==2);
% temp_PM(subj,1) =sum(kmeans_in_gaussianmap_counter_lowvshighcont_Xsubj(subj,:)==3);
% temp_ESI(subj,1)=sum(kmeans_in_gaussianmap_counter_lowvshighcont_Xsubj(subj,:)==4);
% temp_SIL(subj,1)=sum(kmeans_in_gaussianmap_counter_lowvshighcont_Xsubj(subj,:)==5);
% temp_SIR(subj,1)=sum(kmeans_in_gaussianmap_counter_lowvshighcont_Xsubj(subj,:)==6);
% end
% lowvshighcont_countkmeans=[temp_EVL temp_EVR temp_PM temp_ESI temp_SIL temp_SIR];
% 
% % low contrast
% kmeans_in_gaussianmap_counter_lowcont_Xsubj=[];group_val=[];
% for subj = 1:length(data.LowContrast.Replayme)%run through each subject
%     %1 step) find the optimal number of clusters
%     optimalNumcluster=round(data.LowContrast.Rivalme(1).Xsubjects.kmeans_allminX(subj),0); % these are the optimal numbers of clusters for this participant / condition
%     %2 step) run k means for that number of clusters
%     kmeans_in_gaussianmap_counter_Xtrials=[];
%     for trial=1:numtrials %run through each trial
%              kmeans_in_gaussianmap_counter(1:30)=nan; %create even length between trials/subjects as k means varies from 2-25
%        %3 step) load the Rival Me datapoint classification for each trial
%         dataPointClassification_Replayme= data_REPLAY.LowContrast.Replayme(subj).classification_Joystick(trial,:)';
%         %4 step) run k means for the optimal cluster number
%         data.LowContrast.Replayme(subj).XY(trial).matrixXY=[data.LowContrast.Replayme(subj).jstkLRratio(trial,:);data.LowContrast.Replayme(subj).jstkFBratio(trial,:) ];%create XY array
%         X=data.LowContrast.Replayme(subj).XY(trial).matrixXY;
%         X=X/2;% convert values tomrange from 0-1
%         X=X+0.5;
%         %scatter(X(1,:),X(2,:))
%         [idx,C,sumdist] = kmeans(X',optimalNumcluster,'Distance','cityblock','Display','final');
%         %5 step) run through each centroid and assign Replayme classified data to the k means data
%         for i=1:length(C)
%             centroid_temp =  X'==C(i,:); % find centroid:  logicals 1 and 1
%             temp =   dataPointClassification_Replayme(centroid_temp(:,1)==1 & centroid_temp(:,2)==1);%shows all the coordinates on which the centroid lands
%             if isempty(temp)==1 % if there is no direct corresponding values then find nearest neighbor
%                 tempC=C(i,:);
%                 for dataPoint=1:length(dataPointClassification_Replayme) % work though each data point
%                 [xLocs, yLocs]=meshgrid(linspace(0,1),linspace(0,1)); % grid for classified space
%                  distFromEachPixel=sqrt((X(1,dataPoint)-tempC(1)).^2+(X(2,dataPoint)-tempC(2)).^2); % calculate distance from each classified pixel
%                 %distFromEachPixel=sqrt((X(1, tempC(1))-xLocs).^2+(X(2, tempC(1))-yLocs).^2); % calculate distance from each classified pixel
%                 [~,closestPixel]=min(distFromEachPixel(:)); % find index of closest classified pixel
%                 temp_kmeans_in_gaussianmap_counter(dataPoint)=data.LowContrast.Followme(subj).gaussClass(closestPixel); % assign value of closest classified pixel to new data
%                 kmeans_in_gaussianmap_counter(i)=temp_kmeans_in_gaussianmap_counter(1);%pulls out just single data as all other represent the same gaussian map state
%                 end
%             else
%                 kmeans_in_gaussianmap_counter(i)=temp(1); %pulls out just single data as all other represent the same gaussian map state
%             end
%         end
% 
%         %6 step) combine data across trials
%         kmeans_in_gaussianmap_counter_Xtrials=[kmeans_in_gaussianmap_counter_Xtrials kmeans_in_gaussianmap_counter];
%     end % end of trial loop
%     % 7 step) combine across subject
%     kmeans_in_gaussianmap_counter_lowcont_Xsubj=[kmeans_in_gaussianmap_counter_lowcont_Xsubj ;kmeans_in_gaussianmap_counter_Xtrials];
%     group_val=[group_val; subj*ones(length(kmeans_in_gaussianmap_counter_Xtrials),1)'];
% end % end of subject loop
% figure
% % In which of the 6 gaussian states did the optimal k means cluster occur?
% swarmchart(kmeans_in_gaussianmap_counter_lowcont_Xsubj',group_val','b','filled','MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.75)
% title('Replay: # of K-means clusters in each perceptual state map','FontSize',12);
% subtitle('Low contrast condition','FontSize',10);
% set(gca,'XTickLabel',{'','Left','Right','Piecemeal','Equal Super','Left Super','Right Super',''});
% ylabel('Participants','FontSize',12 ,'FontWeight','bold' );
% set(gca,'TickDir','out');
% box off 
% xlabel('Perceptual state','FontSize',12 ,'FontWeight','bold' );
% 
% % how many k means within each state per participant?
% for subj=1: height(kmeans_in_gaussianmap_counter_lowcont_Xsubj)
% temp_EVL(subj,1)=sum(kmeans_in_gaussianmap_counter_lowcont_Xsubj(subj,:)==1);
% temp_EVR(subj,1)=sum(kmeans_in_gaussianmap_counter_lowcont_Xsubj(subj,:)==2);
% temp_PM(subj,1) =sum(kmeans_in_gaussianmap_counter_lowcont_Xsubj(subj,:)==3);
% temp_ESI(subj,1)=sum(kmeans_in_gaussianmap_counter_lowcont_Xsubj(subj,:)==4);
% temp_SIL(subj,1)=sum(kmeans_in_gaussianmap_counter_lowcont_Xsubj(subj,:)==5);
% temp_SIR(subj,1)=sum(kmeans_in_gaussianmap_counter_lowcont_Xsubj(subj,:)==6);
% end
% lowcont_countkmeans=[temp_EVL temp_EVR temp_PM temp_ESI temp_SIL temp_SIR];
% 
% 
% 



%% load Gaussian models and superimpose it with k means cluster during Rival Me
% high contrast
kmeans_in_gaussianmap_counter_highcont_Xsubj=[];group_val=[];
for subj = 1:length(data.HighContrast.Rivalme)%run through each subject
    %1 step) find the optimal number of clusters
    optimalNumcluster=round(data.HighContrast.Rivalme(1).Xsubjects.kmeans_allminX(subj),0); % these are the optimal numbers of clusters for this participant / condition
    %2 step) run k means for that number of clusters
    kmeans_in_gaussianmap_counter_Xtrials=[];
    for trial=1:numtrials %run through each trial
             kmeans_in_gaussianmap_counter(1:30)=nan; %create even length between trials/subjects as k means varies from 2-25
        %3 step) load the Rival Me datapoint classification for each trial
        dataPointClassification_RivalMe= data.HighContrast.Rivalme(subj).dataPointClassification(trial,:)';
        %4 step) run k means for the optimal cluster number
        X=data.HighContrast.Rivalme(subj).XY(trial).matrixXY;
        %scatter(X(1,:),X(2,:))
        [idx,C,sumdist] = kmeans(X',optimalNumcluster,'Distance','cityblock','Display','final');
        %5 step) run through each centroid and assign RivalMe classified data to the k means data
        for i=1:length(C)
            centroid_temp =  X'==C(i,:); % find centroid:  logicals 1 and 1
            temp =   dataPointClassification_RivalMe(centroid_temp(:,1)==1 & centroid_temp(:,2)==1);%shows all the coordinates on which the centroid lands
            if isempty(temp)==1 % if there is no direct corresponding values then find nearest neighbor
                tempC=C(i,:);
                for dataPoint=1:length(dataPointClassification_RivalMe) % work though each data point
                [xLocs, yLocs]=meshgrid(linspace(0,1),linspace(0,1)); % grid for classified space
                 distFromEachPixel=sqrt((X(1,dataPoint)-tempC(1)).^2+(X(2,dataPoint)-tempC(2)).^2); % calculate distance from each classified pixel
                %distFromEachPixel=sqrt((X(1, tempC(1))-xLocs).^2+(X(2, tempC(1))-yLocs).^2); % calculate distance from each classified pixel
                [~,closestPixel]=min(distFromEachPixel(:)); % find index of closest classified pixel
                temp_kmeans_in_gaussianmap_counter(dataPoint)=data.HighContrast.Followme(subj).gaussClass(closestPixel); % assign value of closest classified pixel to new data
                kmeans_in_gaussianmap_counter(i)=temp_kmeans_in_gaussianmap_counter(1);%pulls out just single data as all other represent the same gaussian map state
                end
            else
                kmeans_in_gaussianmap_counter(i)=temp(1); %pulls out just single data as all other represent the same gaussian map state
            end
        end

        %6 step) combine data across trials
        kmeans_in_gaussianmap_counter_Xtrials=[kmeans_in_gaussianmap_counter_Xtrials kmeans_in_gaussianmap_counter];
    end % end of trial loop
    % 7 step) combine across subject
    kmeans_in_gaussianmap_counter_highcont_Xsubj=[kmeans_in_gaussianmap_counter_highcont_Xsubj ;kmeans_in_gaussianmap_counter_Xtrials];
    group_val=[group_val; subj*ones(length(kmeans_in_gaussianmap_counter_Xtrials),1)'];
end % end of subject loop
figure
% In which of the 6 gaussian states did the optimal k means cluster occur?
swarmchart(kmeans_in_gaussianmap_counter_highcont_Xsubj',group_val','b','filled','MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.75)
% title('Rivalry: # of K-means clusters in each perceptual state map','FontSize',12);
% subtitle('High contrast condition','FontSize',10);
set(gca,'XTickLabel',{'','Left','Right','Piecemeal','Equal Super','Left Super','Right Super',''});
ax.FontSize = 18;
ylabel('Participants','FontSize',18 ,'FontWeight','bold' );
xlabel('Perceptual state','FontSize',18 ,'FontWeight','bold' );

set(gca,'TickDir','out');
box off 
% how many k means within each state per participant?
for subj=1: height(kmeans_in_gaussianmap_counter_highcont_Xsubj)
temp_EVL(subj,1)=sum(kmeans_in_gaussianmap_counter_highcont_Xsubj(subj,:)==1);
temp_EVR(subj,1)=sum(kmeans_in_gaussianmap_counter_highcont_Xsubj(subj,:)==2);
temp_PM(subj,1) =sum(kmeans_in_gaussianmap_counter_highcont_Xsubj(subj,:)==3);
temp_ESI(subj,1)=sum(kmeans_in_gaussianmap_counter_highcont_Xsubj(subj,:)==4);
temp_SIL(subj,1)=sum(kmeans_in_gaussianmap_counter_highcont_Xsubj(subj,:)==5);
temp_SIR(subj,1)=sum(kmeans_in_gaussianmap_counter_highcont_Xsubj(subj,:)==6);
end
highcont_countkmeans=[temp_EVL, temp_EVR, temp_PM, temp_ESI, temp_SIL ,temp_SIR];

% plot circles with colorbar as variable
%x 1-6 states
%y 1 circle per participant/states
% color map number of centroids per participants and state
x_temp=[1:size(highcont_countkmeans,2)];
X_temp=repmat(x_temp,height(highcont_countkmeans),1)
y_temp=[1:28];
Y_temp=repmat(y_temp,length(x_temp),1);
Y_temp=Y_temp';
mean_Xsubj_highcont_countkmeans=mean(highcont_countkmeans,2);
mean_Xstates_highcont_countkmeans=mean(highcont_countkmeans,1);

highcont_countkmeans_withmean=[highcont_countkmeans;nan nan nan nan nan nan; round(mean_Xstates_highcont_countkmeans)]
figure
hold on
s=scatter(X_temp(:),Y_temp(:),168,highcont_countkmeans(:),"filled");
% s=scatter(X_temp(:),Y_temp(:),168,highcont_countkmeans(:),"filled",'MarkerFaceAlpha',0.9,'MarkerEdgeAlpha',0.1)

xlim([0 7])
ylim([0 30])
% Add a colorbar
c=colorbar
c.Label.String = '# of K-means clusters';
c.Label.FontWeight="bold"
c.Label.FontSize=12

colormap turbo
% title('Rivalry: # of K-means clusters in each perceptual state map','FontSize',12);
% subtitle('High contrast condition','FontSize',10);
set(gca,'XTickLabel',{'','Left','Right','Piecemeal','Equal Super','Left Super','Right Super',''});
ax.FontSize = 18;
ylabel('Participants','FontSize',18 ,'FontWeight','bold' );
xlabel('Perceptual state','FontSize',18 ,'FontWeight','bold' );

set(gca,'TickDir','out');
box off 



%low vs high
kmeans_in_gaussianmap_counter_lowvshighcont_Xsubj=[];group_val=[];
for subj = 1:length(data.LowvsHighContrast.Rivalme)%run through each subject
    %1 step) find the optimal number of clusters
    optimalNumcluster=round(data.LowvsHighContrast.Rivalme(1).Xsubjects.kmeans_allminX(subj),0); % these are the optimal numbers of clusters for this participant / condition
    %2 step) run k means for that number of clusters
    kmeans_in_gaussianmap_counter_Xtrials=[];
    for trial=1:numtrials %run through each trial
             kmeans_in_gaussianmap_counter(1:30)=nan; %create even length between trials/subjects as k means varies from 2-25
        %3 step) load the Rival Me datapoint classification for each trial
        dataPointClassification_RivalMe= data.LowvsHighContrast.Rivalme(subj).dataPointClassification(trial,:)';
        %4 step) run k means for the optimal cluster number
        X=data.LowvsHighContrast.Rivalme(subj).XY(trial).matrixXY;
        %scatter(X(1,:),X(2,:))
        [idx,C,sumdist] = kmeans(X',optimalNumcluster,'Distance','cityblock','Display','final');
        %5 step) run through each centroid and assign RivalMe classified data to the k means data
        for i=1:length(C)
            centroid_temp =  X'==C(i,:); % find centroid:  logicals 1 and 1
            temp =   dataPointClassification_RivalMe(centroid_temp(:,1)==1 & centroid_temp(:,2)==1);%shows all the coordinates on which the centroid lands
            if isempty(temp)==1 % if there is no direct corresponding values then find nearest neighbor
                tempC=C(i,:);
                for dataPoint=1:length(dataPointClassification_RivalMe) % work though each data point
                [xLocs, yLocs]=meshgrid(linspace(0,1),linspace(0,1)); % grid for classified space
                 distFromEachPixel=sqrt((X(1,dataPoint)-tempC(1)).^2+(X(2,dataPoint)-tempC(2)).^2); % calculate distance from each classified pixel
                %distFromEachPixel=sqrt((X(1, tempC(1))-xLocs).^2+(X(2, tempC(1))-yLocs).^2); % calculate distance from each classified pixel
                [~,closestPixel]=min(distFromEachPixel(:)); % find index of closest classified pixel
                temp_kmeans_in_gaussianmap_counter(dataPoint)=data.LowvsHighContrast.Followme(subj).gaussClass(closestPixel); % assign value of closest classified pixel to new data
                kmeans_in_gaussianmap_counter(i)=temp_kmeans_in_gaussianmap_counter(1);%pulls out just single data as all other represent the same gaussian map state
                end
            else
                kmeans_in_gaussianmap_counter(i)=temp(1); %pulls out just single data as all other represent the same gaussian map state
            end
        end

        %6 step) combine data across trials
        kmeans_in_gaussianmap_counter_Xtrials=[kmeans_in_gaussianmap_counter_Xtrials kmeans_in_gaussianmap_counter];
    end % end of trial loop
    % 7 step) combine across subject
    kmeans_in_gaussianmap_counter_lowvshighcont_Xsubj=[kmeans_in_gaussianmap_counter_lowvshighcont_Xsubj ;kmeans_in_gaussianmap_counter_Xtrials];
    group_val=[group_val; subj*ones(length(kmeans_in_gaussianmap_counter_Xtrials),1)'];
end % end of subject loop
figure
% In which of the 6 gaussian states did the optimal k means cluster occur?
swarmchart(kmeans_in_gaussianmap_counter_lowvshighcont_Xsubj',group_val','b','filled','MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.75)
% title('Rivalry: # of K-means clusters in each perceptual state map','FontSize',12);
% subtitle('Low vs high contrast condition','FontSize',10);
set(gca,'XTickLabel',{'','Left','Right','Piecemeal','Equal Super','Left Super','Right Super',''});
ylabel('Participants','FontSize',18 ,'FontWeight','bold' );
xlabel('Perceptual state','FontSize',18 ,'FontWeight','bold' );

set(gca,'TickDir','out');
box off 

% how many k means within each state per participant?
for subj=1: height(kmeans_in_gaussianmap_counter_lowvshighcont_Xsubj)
temp_EVL(subj,1)=sum(kmeans_in_gaussianmap_counter_lowvshighcont_Xsubj(subj,:)==1);
temp_EVR(subj,1)=sum(kmeans_in_gaussianmap_counter_lowvshighcont_Xsubj(subj,:)==2);
temp_PM(subj,1) =sum(kmeans_in_gaussianmap_counter_lowvshighcont_Xsubj(subj,:)==3);
temp_ESI(subj,1)=sum(kmeans_in_gaussianmap_counter_lowvshighcont_Xsubj(subj,:)==4);
temp_SIL(subj,1)=sum(kmeans_in_gaussianmap_counter_lowvshighcont_Xsubj(subj,:)==5);
temp_SIR(subj,1)=sum(kmeans_in_gaussianmap_counter_lowvshighcont_Xsubj(subj,:)==6);
end
lowvshighcont_countkmeans=[temp_EVL temp_EVR temp_PM temp_ESI temp_SIL temp_SIR];
mean_Xsubj_lowvshighcont_countkmeans=mean(lowvshighcont_countkmeans,2);
mean_Xstates_lowvshighcont_countkmeans=mean(lowvshighcont_countkmeans,1);

% low contrast
kmeans_in_gaussianmap_counter_lowcont_Xsubj=[];group_val=[];
for subj = 1:length(data.LowContrast.Rivalme)%run through each subject
    %1 step) find the optimal number of clusters
    optimalNumcluster=round(data.LowContrast.Rivalme(1).Xsubjects.kmeans_allminX(subj),0); % these are the optimal numbers of clusters for this participant / condition
    %2 step) run k means for that number of clusters
    kmeans_in_gaussianmap_counter_Xtrials=[];
    for trial=1:numtrials %run through each trial
             kmeans_in_gaussianmap_counter(1:30)=nan; %create even length between trials/subjects as k means varies from 2-25
        %3 step) load the Rival Me datapoint classification for each trial
        dataPointClassification_RivalMe= data.LowContrast.Rivalme(subj).dataPointClassification(trial,:)';
        %4 step) run k means for the optimal cluster number
        X=data.LowContrast.Rivalme(subj).XY(trial).matrixXY;
        %scatter(X(1,:),X(2,:))
        [idx,C,sumdist] = kmeans(X',optimalNumcluster,'Distance','cityblock','Display','final');
        %5 step) run through each centroid and assign RivalMe classified data to the k means data
        for i=1:length(C)
            centroid_temp =  X'==C(i,:); % find centroid:  logicals 1 and 1
            temp =   dataPointClassification_RivalMe(centroid_temp(:,1)==1 & centroid_temp(:,2)==1);%shows all the coordinates on which the centroid lands
            if isempty(temp)==1 % if there is no direct corresponding values then find nearest neighbor
                tempC=C(i,:);
                for dataPoint=1:length(dataPointClassification_RivalMe) % work though each data point
                [xLocs, yLocs]=meshgrid(linspace(0,1),linspace(0,1)); % grid for classified space
                 distFromEachPixel=sqrt((X(1,dataPoint)-tempC(1)).^2+(X(2,dataPoint)-tempC(2)).^2); % calculate distance from each classified pixel
                %distFromEachPixel=sqrt((X(1, tempC(1))-xLocs).^2+(X(2, tempC(1))-yLocs).^2); % calculate distance from each classified pixel
                [~,closestPixel]=min(distFromEachPixel(:)); % find index of closest classified pixel
                temp_kmeans_in_gaussianmap_counter(dataPoint)=data.LowContrast.Followme(subj).gaussClass(closestPixel); % assign value of closest classified pixel to new data
                kmeans_in_gaussianmap_counter(i)=temp_kmeans_in_gaussianmap_counter(1);%pulls out just single data as all other represent the same gaussian map state
                end
            else
                kmeans_in_gaussianmap_counter(i)=temp(1); %pulls out just single data as all other represent the same gaussian map state
            end
        end

        %6 step) combine data across trials
        kmeans_in_gaussianmap_counter_Xtrials=[kmeans_in_gaussianmap_counter_Xtrials kmeans_in_gaussianmap_counter];
    end % end of trial loop
    % 7 step) combine across subject
    kmeans_in_gaussianmap_counter_lowcont_Xsubj=[kmeans_in_gaussianmap_counter_lowcont_Xsubj ;kmeans_in_gaussianmap_counter_Xtrials];
    group_val=[group_val; subj*ones(length(kmeans_in_gaussianmap_counter_Xtrials),1)'];
end % end of subject loop
figure
% In which of the 6 gaussian states did the optimal k means cluster occur?
swarmchart(kmeans_in_gaussianmap_counter_lowcont_Xsubj',group_val','b','filled','MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.75)
%title('Rivalry: # of K-means clusters in each perceptual state map','FontSize',12);
%subtitle('Low contrast condition','FontSize',10);
ax = gca; %change tick fontsize
ax.FontSize = 18;

set(gca,'XTickLabel',{'','EV-L','EV-R','PM','E-SI','L-SI','R-SI',''},'FontSize',16);

ylabel('Participants','FontSize',18 ,'FontWeight','bold' );
xlabel('Perceptual state','FontSize',18 ,'FontWeight','bold' );

set(gca,'TickDir','out');
box off 

% how many k means within each state per participant?
for subj=1: height(kmeans_in_gaussianmap_counter_lowcont_Xsubj)
temp_EVL(subj,1)=sum(kmeans_in_gaussianmap_counter_lowcont_Xsubj(subj,:)==1);
temp_EVR(subj,1)=sum(kmeans_in_gaussianmap_counter_lowcont_Xsubj(subj,:)==2);
temp_PM(subj,1) =sum(kmeans_in_gaussianmap_counter_lowcont_Xsubj(subj,:)==3);
temp_ESI(subj,1)=sum(kmeans_in_gaussianmap_counter_lowcont_Xsubj(subj,:)==4);
temp_SIL(subj,1)=sum(kmeans_in_gaussianmap_counter_lowcont_Xsubj(subj,:)==5);
temp_SIR(subj,1)=sum(kmeans_in_gaussianmap_counter_lowcont_Xsubj(subj,:)==6);
end
lowcont_countkmeans=[temp_EVL temp_EVR temp_PM temp_ESI temp_SIL temp_SIR];


mean_Xsubj_lowcont_countkmeans=mean(lowcont_countkmeans,2);
mean_Xstates_lowcont_countkmeans=mean(lowcont_countkmeans,1);

%% ANOVA 1 way- number of centroids

Y=[mean_Xsubj_lowcont_countkmeans mean_Xsubj_lowvshighcont_countkmeans mean_Xsubj_highcont_countkmeans];
[p,tbl,stats]=anova1(Y);
%% partial eta squared - Effect size for ANOVAs
% .01: Small effect size
% .06: Medium effect size
% .14 or higher: Large effect size
%Partial eta squared = SSeffect / (SSeffect + SSerror)
partial_eta_squared_phase=tbl{2,2}/(tbl{2,2} + tbl{4,2}); % partial eta^2 group
Y_mean_Xsubj_lowcont_countkmeans_meanXsubj=mean(Y(:,1));
Y_mean_Xsubj_lowcont_countkmeans_stdXsubj=std(Y(:,1));

Y_mean_Xsubj_lowvshighcont_countkmeans_meanXsubj=mean(Y(:,2));
Y_mean_Xsubj_lowvshighcont_countkmeans_stdXsubj=std(Y(:,2));

Y_mean_Xsubj_highcont_countkmeans_meanXsubj=mean(Y(:,3));
Y_mean_Xsubj_highcont_countkmeans_stdXsubj=std(Y(:,3));

%% ANOVA 2 way- number of centroids

Y=[lowcont_countkmeans lowvshighcont_countkmeans highcont_countkmeans];
g1=[ones(height(lowcont_countkmeans),1) 2*ones(height(lowcont_countkmeans),1) 3*ones(height(lowcont_countkmeans),1)...
   4* ones(height(lowcont_countkmeans),1) 5*ones(height(lowcont_countkmeans),1) 6*ones(height(lowcont_countkmeans),1)...
  ones(height(lowcont_countkmeans),1) 2*ones(height(lowcont_countkmeans),1) 3*ones(height(lowcont_countkmeans),1)...
   4* ones(height(lowcont_countkmeans),1) 5*ones(height(lowcont_countkmeans),1) 6*ones(height(lowcont_countkmeans),1)...
ones(height(lowcont_countkmeans),1) 2*ones(height(lowcont_countkmeans),1) 3*ones(height(lowcont_countkmeans),1)...
   4* ones(height(lowcont_countkmeans),1) 5*ones(height(lowcont_countkmeans),1) 6*ones(height(lowcont_countkmeans),1)];
g2=[ones(length(lowcont_countkmeans(:)),1) 2*ones(length(lowvshighcont_countkmeans(:)),1) 3*ones(length(highcont_countkmeans(:)),1)];

[p,tbl,stats]= anovan(Y(:),{g1(:), g2(:)},"Model","interaction", ...
    "Varnames",["states","contrasts"])
%% partial eta squared - Effect size for ANOVAs
% .01: Small effect size
% .06: Medium effect size
% .14 or higher: Large effect size
%Partial eta squared = SSeffect / (SSeffect + SSerror)
partial_eta_squared_state=tbl{2,2}/(tbl{2,2} + tbl{5,2}); % partial eta^2 state
partial_eta_squared_contrast=tbl{3,2}/(tbl{3,2} + tbl{5,2}); % partial eta^2 contrast

Y=[lowcont_countkmeans ;lowvshighcont_countkmeans ;highcont_countkmeans];

Y_mean_Xstates=mean(Y,1);
Y_std_Xstates=std(Y,1);
%% rel proprtion
Y_mean_Xstates_total=round(sum(Y_mean_Xstates))
Y_mean_Xstates_EL_relProp=Y_mean_Xstates(1)/Y_mean_Xstates_total*100
Y_mean_Xstates_ER_relProp=Y_mean_Xstates(2)/Y_mean_Xstates_total*100
Y_mean_Xstates_PM_relProp=Y_mean_Xstates(3)/Y_mean_Xstates_total*100
Y_mean_Xstates_ES_relProp=Y_mean_Xstates(4)/Y_mean_Xstates_total*100
Y_mean_Xstates_SIL_relProp=Y_mean_Xstates(5)/Y_mean_Xstates_total*100
Y_mean_Xstates_SIR_relProp=Y_mean_Xstates(6)/Y_mean_Xstates_total*100

%% save data
testname=['Analysis_kmeans_PART9_joystickCBR'];

while exist([testname, '.mat'],'file')~=0
    n=n+1;
    testname=[testname,num2str(n)];
end
dataFile=sprintf('%s.mat',testname);
save(dataFile, 'data')