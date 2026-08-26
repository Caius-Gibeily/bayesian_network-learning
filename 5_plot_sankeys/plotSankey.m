function plotSankey(ss, props, probes, durations, diffs, counts, nodeOpt, colOpt)

arguments
    ss                                  % [source, target, weight]
    props                            % property labels
    probes                            % probe labels
    durations
    diffs                               % edge differences
    counts                             % node heights

    nodeOpt.nodeList = []                     % optional
    nodeOpt.nodeLayers = []                   % optional

    colOpt.pColorList = "parula"              % optional
    colOpt.diffCol = "greenblue"                   % optional
end

rng(2356)

edgeCols = [1 1 1]

% -------------------------------------------------------------------------
% Determine if optional layer/node info exists
% -------------------------------------------------------------------------
hasLayers = ~isempty(nodeOpt.nodeLayers);
hasList   = ~isempty(nodeOpt.nodeList);

% ========================================================================
% ========== CASE 1: No nodeList / nodeLayers provided ==================
% ========================================================================
if ~hasLayers 

    %figure('Name','sankey demo2','Units','normalized','Position',[.05,.2,.5,.56])

    SK = SSankey( ...
        ss(:,1), ss(:,2), ss(:,3), ...
        'NodeList', nodeOpt.nodeList, ...
        'EdgeDiff', diffs, ...
        'NodeHeight', counts);

    SK.LayerOrder       = 'reverse';
    SK.ColorList        = repmat(edgeCols,15,1);
    SK.RenderingMethod  = 'right';
    SK.Align            = 'center';
    SK.LabelLocation    = 'right';


    if colOpt.pColorList == "parula"
        SK.pColorList = parula(5);
        SK.values = 1:5;
        SK.EdgeDiffCmap = greengreyblue_colormap(61);

    elseif colOpt.pColorList == "rand"
        SK.pColorList = rand(208, 3);
        SK.values = -3:200;
    end

    if colOpt.diffCol == "greenblue"
        SK.EdgeDiffCmap = greengreyblue_colormap(61);

    elseif colOpt.diffCol == "redgreen"
        SK.EdgeDiffCmap = redgreygreen_colormap(61);
    end

    % if colOpt.edgeCol == "convdiv"
    %     SK.ColorList = redgreygreen_colormap();
    % elseif colOpt.edgeCol == "dx"
    %     SK.ColorList = repmat(edgeCols,15,1);
    % end

    FontCell = {'FontSize',15,'FontName','Times New Roman', ...
                'HorizontalAlignment','center','VerticalAlignment','bottom'};

    %SK.BlockScale = 0.3;
    SK.BlockScale = durations*10;
    epochScales = unique(SK.BlockScale);

    SK.draw(props, probes)

    % Set borders
    for i = 1:SK.BN
        SK.setBlock(i,'EdgeColor',[0,0,0],'LineWidth',1)
    end

    % Layer labels
    chars = 'A':'Z';
    for i = 1:SK.LN
        text(i + epochScales(i)/2, min(min(SK.LayerPos(:,3:4))), ...
            [chars(i)], FontCell{:})
    end
  
else

    SK = SSankey( ...
        ss(:,1), ss(:,2), ss(:,3), ...
        'EdgeDiff', diffs, ...
        'NodeHeight', counts, ...
        'NodeList', nodeOpt.nodeList, ...
        'Layer', nodeOpt.nodeLayers);

    SK.LayerOrder       = 'reverse';
    SK.ColorList        = repmat(edgeCols,15,1);
    SK.ColorList = [1 1 1];
    SK.RenderingMethod  = 'right';
    SK.Align            = 'center';
    SK.LabelLocation    = 'right';
    

    if colOpt.pColorList == "parula"
        SK.pColorList = parula(5);
        SK.values = 1:5;
        SK.EdgeDiffCmap = greengreyblue_colormap(61);

    elseif colOpt.pColorList == "rand"
        SK.pColorList = rand(208, 3);
        SK.values = -3:200;
    end

    if colOpt.diffCol == "greenblue"
        SK.EdgeDiffCmap = greengreyblue_colormap(61);

    elseif colOpt.diffCol == "redgreen"
        SK.EdgeDiffCmap = redgreygreen_colormap(61);
    end    
    nLayers = range(nodeOpt.nodeLayers);
    SK.BlockScale = 0.3 * nLayers;
    SK.draw(props, probes)

    % Set borders
    for i = 1:SK.BN
        SK.setBlock(i,'EdgeColor',[0,0,0],'LineWidth',1.5)
    end

    % Natural layer labeling
    layers = unique(nodeOpt.nodeLayers);

    
    
    FontCell = {'FontSize',15,'FontName','Times New Roman', ...
                'HorizontalAlignment','center','VerticalAlignment','bottom'};
    
    chars = 'A':'Z';
    for i = 1:length(layers)
        text(layers(i) + SK.BlockScale/2, ...
            min(min(SK.LayerPos(:,3:4))), [chars(i), FontCell{:}])
    end



end

end
