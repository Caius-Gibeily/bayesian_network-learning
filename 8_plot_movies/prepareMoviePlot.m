function prepareMoviePlot(data,clip,movieDir,plotStates,bydx,targetCells,targetIdentity,pSize,saveDir)

arguments
    data
    clip
    movieDir
    plotStates = 0
    bydx = 1
    targetCells = []
    targetIdentity = []
    pSize = []
    saveDir = []
    
end


clips = [335 340 348 351 352];
clipNum = find(clips==clip);


coordsTD = data{clipNum}.coordsTD;
coordsASD = data{clipNum}.coordsASD;
if ~isempty(targetCells)
    plotOpts.targetCells = targetCells;
end
if ~isempty(targetIdentity)
    plotOpts.targetIdentity = targetIdentity;
end
if ~isempty(pSize)
    plotOpts.pSize = pSize;
end
if bydx ~= 0 
    plotOpts.bydx = 1;
else
    plotOpts.bydx = 0;
end
if plotStates ~= 0 
    plotOpts.clusterArrayTD = data{clipNum}.clusterArrayTD;
    plotOpts.clusterArrayASD = data{clipNum}.clusterArrayASD;
end
if ~isempty(saveDir)
    plotOpts.saveDir = saveDir;
end
plotOpts.epochs = data{clipNum}.epochs;
plotOpts.clip = clip;
plotMovies(movieDir,coordsTD,coordsASD,plotOpts)


end