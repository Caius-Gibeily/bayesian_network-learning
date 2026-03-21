function plotMovies(movieDir, coordsTD, coordsASD, plotOpts)

arguments
    movieDir
    coordsTD
    coordsASD
    plotOpts struct
end

if ~isfield(plotOpts,"epochs"),         plotOpts.epochs = []; end
if ~isfield(plotOpts,"targetCells"),    plotOpts.targetCells = []; end
if ~isfield(plotOpts,"targetIdentity"), plotOpts.targetIdentity = []; end
if ~isfield(plotOpts,"bydx"),            plotOpts.bydx = 1; end
if ~isfield(plotOpts,"pColorList"),      plotOpts.pColorList = "greenblue"; end
if ~isfield(plotOpts,"pSize"),           plotOpts.pSize = 1; end
if ~isfield(plotOpts,"tColorList"),      plotOpts.tColorList = "parula"; end
if ~isfield(plotOpts,"saveDir")
    plotOpts.saveDir = "sVideo_Clip" + ...
        string(plotOpts.clip) + ...
        "-annotated.avi";
else
     plotOpts.saveDir = plotOpts.saveDir + "sVideo_Clip" + ...
        string(plotOpts.clip) + ...
        "-annotated.avi";
end
display(plotOpts.saveDir)

if plotOpts.bydx == 0
    plotOpts.pColorList = "states";
else
    plotOpts.pColorList = "greenblue";
end
      


%video = %strcat("C:\Users\cgibeil@emory.edu\OneDrive - Emory\Documents\PhD\P1-CondProbs\Core_scripts\Cellular_targets_pipeline\clips\","0",string(clips(f)),"PEER_movie.mov");
vclip = VideoReader(movieDir);
vframes = read(vclip);
if plotOpts.pColorList == "greenblue"
    cmap =  greengreyblue_colormap(100);
elseif plotOpts.pColorList == "states"
    cmap = [
        1 0 0;
        0 0 1;
        0 0.5 0;
        1 0.647 0];
    %cmap = hsv(5);
end

%coordsTD = getCoords(Participants2,SessionNumbers2,clips(f));
%coordsASD = getCoords(ASDparticipants2,ASDsessionNumbers2,clips(f));

entDur = 60;
F = struct;
F.cdata = []; 
F.colormap = [];

imshow(imresize(vframes(:,:,:,1),[480 640]))
tmp = getframe(gcf);
frameSize = size(tmp.cdata);

for j = 1:length(plotOpts.epochs)-1
    if ~isscalar(plotOpts.pSize)
        infEpoch = plotOpts.pSize(plotOpts.pSize(:,2)==j,:);
        if isempty(infEpoch)
            F = insertEpoch(F,j,plotOpts.epochs(j+1)-plotOpts.epochs(j),entDur,frameSize,1);
        else
            F = insertEpoch(F,j,plotOpts.epochs(j+1)-plotOpts.epochs(j),entDur,frameSize,0);
        end
    else
         F = insertEpoch(F,j,plotOpts.epochs(j+1)-plotOpts.epochs(j),entDur,frameSize,0);
    end

    for i = plotOpts.epochs(j):plotOpts.epochs(j+1)
        
        imshow(imresize(vframes(:,:,:,i),[480 640]))
        hold on
        
        coordsetTD = cell2mat(cellfun(@(x) x(i,:),coordsTD,'UniformOutput',false));
        coordsetASD = cell2mat(cellfun(@(x) x(i,:),coordsASD,'UniformOutput',false));
        if ~isempty(plotOpts.clusterArrayTD)
            [~, ~, statesTD] = unique(plotOpts.clusterArrayTD(:,j));
            [~, ~, statesASD] = unique(plotOpts.clusterArrayASD(:,j));
        end

        if ~isscalar(plotOpts.pSize)
              
            if isempty(infEpoch)
                pSizeTD = 0.2;
                pSizeASD = 0.2;
            else
                 
                 [~, loc] = ismember(statesTD, infEpoch(:,1)); 
                 pSizeTD = infEpoch(loc,3);

                 [~, loc] = ismember(statesASD, infEpoch(:,1)); 
                 pSizeASD = infEpoch(loc,3);
            end
            
            %pSizeTD = plotOpts.pSize(plotOpts.pSize(:,2)==j & plotOpts.pSize(:,1) == plotOpts.clusterArrayTD(:,j),3);
            %pSizeASD = plotOpts.pSize(plotOpts.pSize(:,2)==j & plotOpts.pSize(:,1) == plotOpts.clusterArrayTD(:,j),4);
        else
            pSizeTD = 0.2;
            pSizeASD = 0.2;
        end
        if plotOpts.bydx == 1
            scatter(coordsetTD(:,1),coordsetTD(:,2),pSizeTD*100,"x","filled","MarkerFaceColor",cmap(end,:))
            scatter(coordsetASD(:,1),coordsetASD(:,2),pSizeASD*100,".","filled","MarkerFaceColor",cmap(1,:))
            colormap(cmap)
            axis off
        elseif plotOpts.bydx == 0 & ~isempty(plotOpts.clusterArrayTD)
            scatter(coordsetTD(:,1),coordsetTD(:,2),"x","SizeData",pSizeTD*100,"CData",cmap(statesTD,:))
            scatter(coordsetASD(:,1),coordsetASD(:,2),"+","SizeData",pSizeASD*100,"CData",cmap(statesASD,:))
            colormap(cmap)
            axis off
        end
        
        if ~isempty(plotOpts.targetCells)
            if ~isempty(plotOpts.targetIdentity)
                targetPlots = changem(plotOpts.targetCells{i}, ...
                    plotOpts.targetIdentity(:,2), ...
                    plotOpts.targetIdentity(:,1));

                h = imagesc(targetPlots);
                alphaMap = plotOpts.targetCells{i}~=-1;
                
                set(h, 'AlphaData', alphaMap*0.3);
                colormap(parula(5))
            else
                [~,~,uTargets] = unique(unique(plotOpts.targetCells{i}));
                targetPlots = changem(plotOpts.targetCells{i},uTargets, ...
                    unique(plotOpts.targetCells{i}));
                targetPlots(plotOpts.targetCells{i}==-1) = -1;

                h = imagesc(targetPlots);
                alphaMap = plotOpts.targetCells{i}~=-1;
                
                set(h, 'AlphaData', alphaMap*0.3);
                colormap(gray(5))
            end

            clim([1,5])
            axis off
        end
        xlim([0 640])
        ylim([0 480])
        F(end+1) = getframe(gcf);
        clf
        display(i)
    end
  
end
  
writerObj = VideoWriter(plotOpts.saveDir);
writerObj.FrameRate = 30;

open(writerObj);
F = F(2:end-1)
for n=1:length(F)
    frame = F(n) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);
clear F

function F = insertEpoch(F, eNum, eDur, entDur,frameSize,jointEdge)

    clf
    ax = axes('Units','normalized','Position',[0 0 1 1]);
    
    for c = 1:entDur
        imagesc(ax, ones(frameSize(1),frameSize(2)))   % MATCH movie size
        colormap(ax,'gray')
        
        if jointEdge == 0
            text(ax, frameSize(2)/2, frameSize(1)/2, ...
                {["Epoch: " + eNum], ...
                 ["Duration: " + round(eDur/30,2) + " s"]}, ...
                'Color','w','FontSize',14,'FontWeight','bold', ...
                'HorizontalAlignment','center')
        else
            text(ax, frameSize(2)/2, frameSize(1)/2, ...
                {["Epoch: " + eNum], ...
                 ["Duration: " + round(eDur/30,2) + " s"], ...
                 ["No joint edges"]}, ...
                'Color','w','FontSize',14,'FontWeight','bold', ...
                'HorizontalAlignment','center')
        end
    
        axis(ax,'off')
        F(end+1) = getframe(gcf);
    end
end


end


    