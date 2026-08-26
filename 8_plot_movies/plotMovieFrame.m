function plotMovieFrame(movieDir, coordsTD, coordsASD, plotOpts, frame)

arguments
    movieDir
    coordsTD
    coordsASD
    plotOpts struct
    frame
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
        "-annotated.png";
else
     plotOpts.saveDir = plotOpts.saveDir + "sVideo_Clip_" + ...
        string(frame) + "_" + string(plotOpts.state) + ...
        string(plotOpts.clip) + ...
        "-annotated.png";
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


imshow(imresize(vframes(:,:,:,frame),[480 640]))
hold on


coordsetTD = cell2mat(cellfun(@(x) x(frame,:),coordsTD,'UniformOutput',false));
coordsetASD = cell2mat(cellfun(@(x) x(frame,:),coordsASD,'UniformOutput',false));
if ~isempty(plotOpts.clusterArrayTD)
    id = find(frame >= plotOpts.epochs, 1, 'last');
    display(id)
    [~, ~, statesTD] = unique(plotOpts.clusterArrayTD(:,id));
    [~, ~, statesASD] = unique(plotOpts.clusterArrayASD(:,id));
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
    scatter(coordsetTD(statesTD==plotOpts.state,1),coordsetTD(statesTD==plotOpts.state,2),"x","SizeData",pSizeTD*1000,"CData",cmap(statesTD(statesTD==plotOpts.state),:))
    scatter(coordsetASD(statesASD==plotOpts.state,1),coordsetASD(statesASD==plotOpts.state,2),"+","SizeData",pSizeASD*1000,"CData",cmap(statesASD(statesASD==plotOpts.state),:))
    colormap(cmap)
    axis off
end

if ~isempty(plotOpts.targetCells)
    if ~isempty(plotOpts.targetIdentity)
        targetPlots = changem(plotOpts.targetCells{frame}, ...
            plotOpts.targetIdentity(:,2), ...
            plotOpts.targetIdentity(:,1));

        h = imagesc(targetPlots);
        alphaMap = plotOpts.targetCells{frame}~=-1;
        
        set(h, 'AlphaData', alphaMap*0.3);
        colormap(parula(5))
    else
        [~,~,uTargets] = unique(unique(plotOpts.targetCells{frame}));
        targetPlots = changem(plotOpts.targetCells{frame},uTargets, ...
            unique(plotOpts.targetCells{frame}));
        targetPlots(plotOpts.targetCells{frame}==-1) = -1;

        h = imagesc(targetPlots);
        alphaMap = plotOpts.targetCells{frame}~=-1;
        
        set(h, 'AlphaData', alphaMap*0.3);
        colormap(gray(5))
    end

    clim([1,5])
    axis off
end
xlim([0 640])
ylim([0 480])

exportgraphics(gcf, plotOpts.saveDir, 'Resolution', 300)

end