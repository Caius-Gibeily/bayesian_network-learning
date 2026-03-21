function [L3,L4,rois] = fitTargets2(coordsTD,coordsASD,gridsize,radius,alpha,window,plotInt)

arguments
    coordsTD cell
    coordsASD cell
    gridsize {mustBeFloat}
    radius {mustBeFloat}
    alpha {mustBeFloat}
    window int16 = 20
    plotInt logical = false
end


% Collate fixation coordinates per frame
frames = cell(length(coordsTD{1}),1);
for i = 1:length(coordsTD{1})
    frames{i} = cell2mat(cellfun(@(x) x(i,:),coordsTD,'UniformOutput',false));
    frames{i} = [frames{i}; cell2mat(cellfun(@(x) x(i,:),coordsASD,'UniformOutput',false))];
end
%%
gridsize; % 20
binsy = 1:gridsize:500; % 10 pixel padding
binsx = 1:gridsize:660;
[x,y] = meshgrid(binsx, binsy);

L3 = cell(length(coordsTD{1}),1);
L3{1} = zeros(size(x));
rois = cell(length(coordsTD{1}),1);
rois{1} = zeros(size(x));
maxes = [];
maxes(1) = 1;
%alpha = 0.6

for k = 2:length(coordsTD{1})-2

    xy = [frames{k-1};frames{k};frames{k+1}];
    
    xy = xy(xy(:,1)>0,:);
 
    nns = sum(squareform(pdist(xy,"euclidean"))<=radius);
    nns = nns/max(nns);

    vq = griddata(xy(:,1),xy(:,2),nns,x,y,"cubic");
    vq(isnan(vq)) = 0;
    
    vqsmooth = imgaussfilt(vq,1.5);

    vqsmooth = zscore(vqsmooth,0,"all");
    rois{k} =vqsmooth;
    
    if k > 3
        vqsmooth = alpha * rois{k} + (1-alpha) * rois{k-1} + (1-alpha)^2 * rois{k-2} + (1-alpha)^3 * rois{k-3};
    end
    
    L3{k} = findgauspeaks(x,y,vqsmooth,20,max(maxes),unique(L3{k-1}));

    L3{k} = localAlign(L3{k-1},L3{k},0.4,max(maxes));
    maxes(k) = max(L3{k-1},[],"all");
    
    if true(plotInt)
        imagesc(L3{k}); axis equal; axis tight; colormap('parula');
        hold on
        clim([0 10])
        scatter(xy(:,1)/640*(640/gridsize),xy(:,2)/480*(480/gridsize),"filled","red")
        % 
        shg
        clf
    end
    display(k)
end

%% Post-processing
% Short-range alignment

wsize = 1:window;
L3a = L3;

for k = 1:length(coordsTD{1})-window-1
    maskWindow = L3a(wsize);
    L3a(wsize) = reclassifyTargets(maskWindow);
    wsize = wsize+1;
end
L3a(end-2) = [];
all_vals = unique(cell2mat(L3a));
new_labels = 0:length(all_vals)-1;
map = containers.Map(all_vals, new_labels);
L3a = cellfun(@(mat) arrayfun(@(x) map(x), mat), L3a, 'UniformOutput', false);


% Extract singleton targets
uCells = cell2mat(cellfun(@(x) unique(x),L3a,'UniformOutput',false));
celltab = tabulate(uCells);
singletons = celltab(celltab(:,2)<=5,1);

% Smooth masks
x_new = 1:1:640;
y_new = 1:1:480;
[X_new, Y_new] = meshgrid(x_new, y_new);
L4 = cellfun(@(m) smoothMasks(m,x,y,X_new,Y_new),L3a,'UniformOutput',false);

% Merge adjacent budded singleton and mother cell
warning('off')
for i = 1:length(singletons)
    L4 = cellfun(@(x) merge_nearest_cluster(x,singletons(i)),L4,'UniformOutput',false);
    display(i)
end
% Remove isolated singletons 
uCells2 = cell2mat(cellfun(@(x) unique(x),L4,'UniformOutput',false));
celltab2 = tabulate(uCells2);
singletonsDetached = celltab2(celltab2(:,2)<=5,1);
for i = 1:length(singletonsDetached)
    for j = 1:length(L4)
        L4{j}(L4{j}==singletonsDetached(i)) = 0;
    end
    display(i)
end

%% Temporal bounds
% 
% % 
% video = "C:\Users\cgibeil@emory.edu\OneDrive - Emory\Documents\PhD\P1-CondProbs\Core_scripts\0335PEER_movie.mov"
% vclip = VideoReader(video)
% % vframes = read(vclip);
% % 
% % cmap = jet(150); % You can replace 'jet' with 'parula', 'hsv', etc.
% % 
% % % Shuffle the colors randomly
% 
% coordsTD = getCoords(Participants2,SessionNumbers2,cliptabs(3))
% coordsASD = getCoords(ASDparticipants2,ASDsessionNumbers2,cliptabs(3))
% 
% cmap = cmap(randperm(200), :)
% figure
% for i = 2:1682
% 
%     %imshow(vframes(:,:,:,i))
%     hold on
%     coordsetTD = cell2mat(cellfun(@(x) x(i,:),coordsTD,'UniformOutput',false));
%     coordsetASD = cell2mat(cellfun(@(x) x(i,:),coordsASD,'UniformOutput',false));
%     coordsetTD(coordsetTD(:,1)==0,:) = [];
%     coordsetASD(coordsetASD(:,1)==0,:) = [];
% 
%     scatter(coordsetTD(:,1),coordsetTD(:,2),"filled","blue")
%     scatter(coordsetASD(:,1),coordsetASD(:,2),"filled","green")
% 
%     h = imagesc(targetCell{3}{i});
% 
%     alphaMap = targetCell{3}{i}~=-1;
%     set(h, 'AlphaData', alphaMap*0.8);
%     %clim([0 120])
%     clim([1 150])
%     colormap(cmap)
%     % subplot(1,2,2)
%     % imagesc(rois{i})
%     % colormap("parula")
%     % 
%     shg
%     clf
% end

% writerObj = VideoWriter('Targets_filtered-medAlign20f2.avi');
% writerObj.FrameRate = 30;
%   % set the seconds per image
% % open the video writer
% open(writerObj);
% % write the frames to the video
% for i=5:length(F)-3
%     % convert the image to a frame
%     frame = F(i) ;    
%     writeVideo(writerObj, frame);
% end
% % close the writer object
% close(writerObj);
% 
% uCells2 = cell2mat(cellfun(@(x) unique(x),L4,'UniformOutput',false));
% celltab2 = tabulate(uCells2);
% 
% lineage = zeros(length(1:max(celltab2)),1682);
% for i = 1:length(L4)
%     uVals = unique(L4{i});
%     uVals(uVals==0)=[];
%     for j = 1:length(uVals)
%         lineage(uVals(j),i) = 1;
%     end
% end
% 
% 
% 
% %% Quantify timings of saccades 
% % Proportions of TD/ASD viewers across ROIs
