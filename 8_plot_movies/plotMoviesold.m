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
if ~isfield(plotOpts,"outDir")
    plotOpts.outDir = "sVideo_Clip" + ...
        regexp(movieDir,"\d{4}","match","once") + ...
        "-targets_test.avi";
end


if plotOpts.bydx == 0
    plotOpts.pColorList = "hsv";
else
    plotOpts.pColorList = "greenblue";
end
      


%video = %strcat("C:\Users\cgibeil@emory.edu\OneDrive - Emory\Documents\PhD\P1-CondProbs\Core_scripts\Cellular_targets_pipeline\clips\","0",string(clips(f)),"PEER_movie.mov");
vclip = VideoReader(movieDir);
vframes = read(vclip);
if plotOpts.pColorList == "greenblue"
    cmap =  greengreyblue_colormap(100);
elseif plotOpts.pColorList == "hsv"
    cmap = hsv(5);
end

%coordsTD = getCoords(Participants2,SessionNumbers2,clips(f));
%coordsASD = getCoords(ASDparticipants2,ASDsessionNumbers2,clips(f));

entDur = 60;
F = struct;
F.cdata = []; 
F.colormap = [];

for j = 1:length(plotOpts.epochs)-1   
    F = insertEpoch(F,j,plotOpts.epochs(j+1)-plotOpts.epochs(j),[480 640],entDur);
    for i = plotOpts.epochs(j):plotOpts.epochs(j+1)
        
        imshow(vframes(:,:,:,i))
        hold on
        
        coordsetTD = cell2mat(cellfun(@(x) x(i,:),coordsTD,'UniformOutput',false));
        coordsetASD = cell2mat(cellfun(@(x) x(i,:),coordsASD,'UniformOutput',false));
        if plotOpts.bydx == 1
            scatter(coordsetTD(:,1),coordsetTD(:,2),plotOpts.pSize*100,"x","filled","MarkerFaceColor",cmap(end,:))
            scatter(coordsetASD(:,1),coordsetASD(:,2),plotOpts.pSize*100,"x","filled","MarkerFaceColor",cmap(1,:))
            colormap(cmap)
            axis off
        elseif ~isempty(plotOpts.clusterArrayTD)
            scatter(coordsetTD(:,1),coordsetTD(:,2),"x","SizeData",plotOpts.pSize*100,"CData",cmap(plotOpts.clusterArrayTD(:,j),:))
            scatter(coordsetASD(:,1),coordsetASD(:,2),"x","SizeData",plotOpts.pSize*100,"CData",cmap(plotOpts.clusterArrayASD(:,j),:))
            colormap(cmap)
            axis off
        end
        
        if ~isempty(plotOpts.targetCells)
            if ~isempty(plotOpts.targetIdentity)
                targetPlots = changem(plotOpts.targetCells{i}, ...
                    plotOpts.targetIdentity(:,2), ...
                    plotOpts.targetIdentity(:,1));
            else
                targetPlots = plotOpts.targetCells{i};
            end

            h = imagesc(targetPlots);
            alphaMap = plotOpts.targetCells{i}~=-1;
            
            set(h, 'AlphaData', alphaMap*0.3);
            colormap("turbo")
            clim([1 5])
            colorbar
            axis off
        end

        F(end+1) = getframe(gcf);
        clf
    end
  
end
  
writerObj = VideoWriter(plotOpts.outDir);
writerObj.FrameRate = 30;

open(writerObj);
F = F(2:end)
for n=1:length(F)
    frame = F(n) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);
clear F


function F = insertEpoch(F,eNum,eDur,dim,entDur)
    
    for c = 1:entDur
        epoframe = ones(dim);
        hold on
        imagesc(epoframe)
        colormap("gray")
        
        text(dim(1)/2,dim(2)/2, {[strcat("Epoch: ",string(eNum))], ...
            [strcat("Duration: ",string(round(eDur/30,2)),"s")]}, ...
            'Color', 'white', 'FontSize', 14, 'FontWeight', 'bold');
        xlim([0 dim(2)])
        ylim([0 dim(1)])
        axis off
        F(end+1) = getframe(gcf);

    end
end

end


    