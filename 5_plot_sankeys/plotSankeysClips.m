function plotSankeysClips(data,savedir,avgProp,avgProbe,clips)

for c = 1:length(clips)
    
    [ssTD, propsTD,probesTD,countsTD,~,durationsTD] = formatSankey2(data{c}.clusterArrayTD, ...
        avgProp{c}.TD,avgProbe{c}.TD,data{c}.epochs); 
    [ssASD, propsASD,probesASD,countsASD,~,durationsASD] =  formatSankey2(data{c}.clusterArrayASD, ...
        avgProp{c}.ASD,avgProbe{c}.ASD,data{c}.epochs);

    [diffTD,diffASD] = subtractMatchingPairs(ssTD,ssASD);

    % TD
    figure("Units","normalized","Position",[0.05,0.2,0.35,0.4])
    plotSankey(ssTD,propsTD,probesTD,durationsTD,cell2mat(diffTD(:,3)),cell2mat(countsTD), ...
        "pColorList","parula","diffCol","greenblue")
    saveas(gcf,strcat(savedir, "viewPathsClustered-TD_", string(clips(c)), ".svg"))
    
    % ASD
    figure("Units","normalized","Position",[0.05,0.2,0.35,0.4])
    plotSankey(ssASD,propsASD,probesASD,durationsASD, cell2mat(diffASD(:,3)),cell2mat(countsASD), ...
        "pColorList","parula","diffCol","greenblue")
    saveas(gcf,strcat(savedir, "viewPathsClustered-ASD_", string(clips(c)), ".svg"))

end

end