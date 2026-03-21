frames = cell(1682,1);
for i = 1:1682
    frameTD = cell2mat(cellfun(@(x) x(i,:),coordsTD,'UniformOutput',false));
    frames{i} = [frameTD; cell2mat(cellfun(@(x) x(i,:),coordsASD,'UniformOutput',false))];
end
% 
% registrationTD = zeros(length(coords),1682-2);
% for i = 1:1682-3
%     L4{i}(L4{i}==0) = -1;
%     for j = 1:length(coords)
%         if coords{j}(i,1)>0
%             registrationTD(j,i) = L4{i}(coords{j}(i,2),coords{j}(i,1));
%         end
%     end
% end
% 
for i = 1:1682-2
    [countsTD{i},uValsTD{i}] = groupcounts(cell2mat(cellfun(@(x) x(i),registrationTD2,'UniformOutput',false)));
end
% 
% registrationASD = zeros(length(coordsASD),1682-2);
% for i = 1:1682-3
%     L4{i}(L4{i}==0) = -1;
%     for j = 1:length(coordsASD)
%         if coordsASD{j}(i,1)>0
%             registrationASD(j,i) = L4{i}(coordsASD{j}(i,2),coordsASD{j}(i,1));
%         end
%     end
% end

for i = 1:1682-2
    [countsASD{i},uValsASD{i}] = groupcounts(cell2mat(cellfun(@(x) x(i),registrationASD2,'UniformOutput',false)));
end
% 
totalOutTD = cell2mat(cellfun(@(x) x(1)/sum(x),countsTD,'UniformOutput',false))
totalOutASD = cell2mat(cellfun(@(x) x(1)/sum(x),countsASD,'UniformOutput',false))
boxplot([totalOutTD;totalOutASD]')

% % Total missing
totalMissTD = cell2mat(cellfun(@(x,y) x(y==0)/sum(x),countsTD,uValsTD,'UniformOutput',false))
totalMissASD = cell2mat(cellfun(@(x,y) x(y==0)/sum(x),countsASD,uValsASD,'UniformOutput',false))
boxplot([totalMissTD;totalMissASD]')

uCells2 = cell2mat(cellfun(@(x) unique(x),targets,'UniformOutput',false));
uCells2 = uCells2(~isnan(uCells2));
alltargets = unique(uCells2);
% % collapse cell counts and values
countsTDnormed = cellfun(@(x) 100*x./sum(x), countsTD,'UniformOutput',false);
cs = cell2mat(countsTDnormed(:));
vals = cell2mat(uValsTD(:));
tabTD = [vals,round(cs)];
% Average occupancy across groups
countsASDnormed = cellfun(@(x) 100*x./sum(x), countsASD,'UniformOutput',false);
cs = cell2mat(countsASDnormed(:));
vals = cell2mat(uValsASD(:));
tabASD = [vals,round(cs)];

percentagestargets = zeros(length(alltargets),2);
for i = 1:length(alltargets)
    percTD = mean(tabTD(tabTD(:,1)==alltargets(i),2));
    percASD = mean(tabASD(tabASD(:,1)==alltargets(i),2));
    percTD = 100*percTD/sum(percTD+percASD);
    percASD = 100-percTD;
    percentagestargets(i,:) = [percTD,percASD];
end
imagesc(percentagestargets)
colormap(winter)
xlabel("Group")
xticks([1,2])
xticklabels({"TD","ASD"})
L4 = targets;
C4 = L4;
for i = 1:length(L4)
    t = unique(L4{i});
    t = t(~isnan(t));
    t(t==0) = [];
    t(t==-1) = [];
    for j = 1:length(t)
        [idx,~] = find(uValsTD{i}==t(j));
        [idxASD,~] = find(uValsASD{i}==t(j));

        if ~isempty(idx)
            countTD = countsTD{i}(idx)./sum(countsTD{i});
        else 
            countTD = 0;
        end
        if ~isempty(idxASD)
            countASD = countsASD{i}(idxASD)./sum(countsASD{i});
        else
            countASD = 0;
        end
        cval = round(100*countASD/(countASD+countTD));
        C4{i}(C4{i}==t(j)) = cval;
    end
    display(i)
end
%%
video = "C:\Users\cgibeil@emory.edu\OneDrive - Emory\Documents\PhD\P1-CondProbs\Core_scripts\0335PEER_movie.mov"
vclip = VideoReader(video)
vframes = read(vclip);
prange = winter(100);

demoframes = sort(randperm(1500,10))
demoframes = [158,487,612,632,795,913,1077,1162,1230,1451]
for i = 1:length(demoframes)

    imshow(vframes(:,:,:,demoframes(i)))
    hold on

    xy = [frames{demoframes(i)}];
    xy = xy(xy(:,1)>0,:);
    coordsetTD = cell2mat(cellfun(@(x) x(demoframes(i),:),coordsTD,'UniformOutput',false));
    coordsetASD = cell2mat(cellfun(@(x) x(demoframes(i),:),coordsASD,'UniformOutput',false));
    coordsetTD(coordsetTD(:,1)==0,:) = [];
    coordsetASD(coordsetASD(:,1)==0,:) = [];
    

    scatter(coordsetTD(:,1),coordsetTD(:,2),"filled","blue")
    scatter(coordsetASD(:,1),coordsetASD(:,2),"filled","green")
    h = imagesc(targets{demoframes(i)});
    colormap(cmap)
    alphaMap = targets{demoframes(i)}~=0;
    set(h, 'AlphaData', alphaMap*0.7);
    clim([1 200])
    
    colormap(cmap)
    text(20,50,strcat("Frame: ",string(demoframes(i))),"FontWeight","bold","Color","white")
    
    F(i) = getframe(gcf);
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    saveas(gcf,strcat("cellularPipeline_demoFrame_",string(demoframes(i)),".svg"))
    clf
end

for i = 1:length(demoframes)
    imagesc(densityMaps{demoframes(i)})
    colormap("parula")
    hold on
    coordsetTD = cell2mat(cellfun(@(x) x(demoframes(i),:),coordsTD,'UniformOutput',false));
    coordsetASD = cell2mat(cellfun(@(x) x(demoframes(i),:),coordsASD,'UniformOutput',false));
    coordsetTD(coordsetTD(:,1)==0,:) = [];
    coordsetASD(coordsetASD(:,1)==0,:) = [];
    
    scatter(coordsetTD(:,1)*36/480,coordsetTD(:,2)*36/480,10,"filled","blue")
    scatter(coordsetASD(:,1)*48/640,coordsetASD(:,2)*48/640,10,"filled","green")
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    saveas(gcf,strcat("cellularPipeline_demoDensityMap_",string(demoframes(i)),".svg"))

    clf
end
cmap = hsv(200);
cmap = cmap(randperm(200), :)
cmap(1,:) = [0 0 0]
for i = 1:length(demoframes)
    imagesc(L3{demoframes(i)})
    colormap(cmap)
    clim([1 200])
    hold on
    coordsetTD = cell2mat(cellfun(@(x) x(demoframes(i),:),coordsTD,'UniformOutput',false));
    coordsetASD = cell2mat(cellfun(@(x) x(demoframes(i),:),coordsASD,'UniformOutput',false));
    coordsetTD(coordsetTD(:,1)==0,:) = [];
    coordsetASD(coordsetASD(:,1)==0,:) = [];
    
    scatter(coordsetTD(:,1)*36/480,coordsetTD(:,2)*36/480,"filled","blue")
    scatter(coordsetASD(:,1)*48/640,coordsetASD(:,2)*48/640,"filled","green")
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    saveas(gcf,strcat("cellularPipeline_demoEarlyMasks_",string(demoframes(i)),".svg"))
    clf
end

cmap = hsv(200);
cmap = cmap(randperm(200), :)
for i = 1:length(demoframes)
    imagesc(targets{demoframes(i)})
    colormap(cmap)
    clim([1 200])
    hold on
    coordsetTD = cell2mat(cellfun(@(x) x(demoframes(i),:),coordsTD,'UniformOutput',false));
    coordsetASD = cell2mat(cellfun(@(x) x(demoframes(i),:),coordsASD,'UniformOutput',false));
    coordsetTD(coordsetTD(:,1)==0,:) = [];
    coordsetASD(coordsetASD(:,1)==0,:) = [];
    
    scatter(coordsetTD(:,1),coordsetTD(:,2),"filled","blue")
    scatter(coordsetASD(:,1),coordsetASD(:,2),"filled","green")
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    saveas(gcf,strcat("cellularPipeline_demoFinalMasks_",string(demoframes(i)),".svg"))
    clf
end

writerObj = VideoWriter('Targets-propGroup-TD_blueNew.avi');
writerObj.FrameRate = 30;
  % set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for i=5:length(F)-3
    % convert the image to a frame
    frame = F(i) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);