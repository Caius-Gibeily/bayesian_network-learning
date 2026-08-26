classdef SSankey < handle
% Copyright (c) 2023-2025, Zhaoxu Liu / slandarer
% =========================================================================
% @author : slandarer
% 公众号  : slandarer随笔
% 知乎    : slandarer
% =========================================================================
% # update 2.0.0(2024-02-04)
% see natureSankeyDemo1.m
%
% + 层向右对齐(Align layers to the right)
%   try : obj.LayerOrder='reverse';
%
% + 单独调整每层间隙大小(Adjust the Sep size of each layer separately)
%   try : obj.Sep=[.2,.06,.05,.07,.07,.08,.15];
% =========================================================================
% # update 3.0.0(2024-04-15)
% see sankeyDemo9.m sankeyDemo10.m sankeyDemo11.m
% 
% + 通过邻接矩阵创建桑基图(Creating a Sankey diagram through adjacency matrix)
%   method 1 :
%     SK=SSankey([],[],[],'AdjMat',adjMat);
%   method 2 :
%     SK=SSankey([],[],[],'NodeList',nodeList,'AdjMat',adjMat)
%   method 3 :
%     SK=SSankey([],[],[]);
%     SK.AdjMat=adjMat;
% 
%   try : 
%     adjMat=zeros(10,10);
%     layerNum=[3,3,2,2];
%     layerInd=cumsum([0,layerNum]);
%     for i=1:length(layerInd)-2
%         adjMat(layerInd(i)+1:layerInd(i+1),layerInd(i+1)+1:layerInd(i+2))=randi([1,6],[layerNum([i,i+1])]);
%     end
%     disp(adjMat)
%     SK=SSankey([],[],[],'NodeList',nodeList,'AdjMat',adjMat);
%     SK.draw()
%
% + 每层情况可被设置(Each layer state can be set)
%   try : obj.Layer = [1,1,1, 2,2,2, 3,3, 4,4,...];
% 
% + 每个节点可在x方向上位移(Each node can be displaced in the x-direction)
%   try : obj.moveBlockX(n,dx)
% =========================================================================
% # update 3.1.0(2024-05-15)
% see sankeyDemo12.m sankeyDemo13.m
% + 为链接添加显示数值的文本(Display value labels for each link)
%   try : SK.ValueLabelLocation='left';
% =========================================================================
% # update 4.0.0(2024-05-17)
% see  sankeyDemo14.m sankeyDemo15.m
% + 增添节点及链接(Add node and link)
%   try : obj.addNode(name,layer)
%   try : obj.addLink(source,target,value)
% =========================================================================
% # version 5.0.0
% + 左键添加数据提示框，右键隐藏高亮 
%   Left-click to add data tooltip, right-click to hide highlight

    properties
        Source;Target;Value;
        SourceInd;TargetInd;
        Layer;LayerPos;MovePos;LayerOrder='normal';
        AdjMat;BoolMat;
        RenderingMethod='interp'   % 'left'/'right'/'interp'/'map'/'simple'
        LabelLocation='left'       % 'left'/'right'/'top'/'center'/'bottom'
        ValueLabelLocation='none'  % 'left'/'right'/'center'/'none'
        ValueLabelFormat=@(X)num2str(X);
        Align='center'             % 'up'/'down'/'center'
        BlockScale=0.05;           %  BlockScale>0 ! !
        Sep=0.05;                  %  Sep>=0 ! !
        NodeList={};
        %              Alpha                               text format
        dataTipFormat = {1, 'Source:', 'Target:', 'Value:', 'auto'}
        ColorList=[[65,140,240;252,180,65;224,64,10;5,100,146;191,191,191;26,59,105;255,227,130;18,156,221;
                    202,107,75;0,92,219;243,210,136;80,99,129;241,185,168;224,131,10;120,147,190;65,140,240;252,180,65;224,64,10;5,100,146;191,191,191;26,59,105;255,227,130;18,156,221;
                    202,107,75;0,92,219;243,210,136;80,99,129;241,185,168;224,131,10;120,147,190]./255;
                   [127,91,93;187,128,110;197,173,143;59,71,111;104,95,126;76,103,86;112,112,124;
                    72,39,24;197,119,106;160,126,88;238,208,146;127,91,93;187,128,110;197,173,143;59,71,111;104,95,126;76,103,86;112,112,124;
                    72,39,24;197,119,106;160,126,88;238,208,146]./255];
        BlockHdl;LinkHdl;LabelHdl;ValueLabelHdl;ax;Parent;
        BN;LN;VN;TotalLen;SepLen;
        arginList={'RenderingMethod','LabelLocation','ValueLabelLocation','BlockScale','Layer'...
                   'Sep','Align','ColorList','Parent','NodeList','AdjMat','EdgeDiff','NodeHeight','BlockScale','pColorList','EdgeDiffCmap','values'}
        %values = -1:6
        %pColorList = parula(8)
        %values = -3:200
        values = []
        pColorList = [] %rand(208, 3);
        % --- NEW properties to support per-node stacked-proportion rendering
        NodeProps = {};    % cell array: each cell is numeric vector of proportions per node
        NodeProbes = {};   % cell array: each cell is probe indices mapping to colors, same length as NodeProps{n}
        NodeHeight = []; 
        EdgeDiff = [];     % numeric vector, one value per link
        %EdgeDiffCmap = [linspace(0, 0, 61)', ... % Red channel (stays 0)
        %           linspace(0,1, 61)', ... % Green channel (1 to 0)
        %           linspace(1,0, 61)'];    % Blue channel (0 to 1)
        %EdgeDiffCmap = redgreygreen_colormap(61);
        %EdgeDiffCmap = greengreyblue_colormap(61);
        EdgeDiffCmap = []
    end
% 构造函数 =================================================================
    methods
        function obj=SSankey(varargin)
            % 获取基本数据 -------------------------------------------------
            if ~isempty(varargin) && isa(varargin{1},'matlab.graphics.axis.Axes')
                obj.ax=varargin{1};varargin(1)=[]; 
            else  
            end
            if ~isempty(varargin)
                obj.Source=varargin{1};
            else
                obj.Source={};
            end
            if length(varargin)>=2
                obj.Target=varargin{2};
            else
                obj.Target={};
            end
            if length(varargin)>=3
                obj.Value=varargin{3};
            else
                obj.Value={};
            end
            varargin(1:min(3,length(varargin)))=[];
            % 获取其他信息 -------------------------------------------------
            for i=1:2:(length(varargin)-1)
                tid=ismember(obj.arginList,varargin{i});
                if any(tid)
                obj.(obj.arginList{tid})=varargin{i+1};
                end
            end
            if isempty(obj.ax)&&(~isempty(obj.Parent)),obj.ax=obj.Parent;end
            if isempty(obj.ax),obj.ax=gca;end
            obj.ax.NextPlot='add';
            % 基本数据预处理 -----------------------------------------------
            if isempty(obj.NodeList)
                if isempty(obj.Source)
                    if ~isempty(obj.AdjMat)
                        obj.NodeList=compose('node%d',1:size(obj.AdjMat,1));
                    end
                else
                    obj.NodeList=[obj.Source;obj.Target];
                    obj.NodeList=unique(obj.NodeList,'stable');
                end
            end
            obj.BN=length(obj.NodeList);
            if length(obj.NodeList)>size(obj.ColorList,1)
                obj.ColorList=[obj.ColorList;rand(length(obj.NodeList),3).*.7];
            end
            obj.MovePos=zeros(obj.BN,4);
            % ensure BlockHdl is initialized as cell for multi-segment nodes
            obj.BlockHdl = cell(1,obj.BN);
            
            % obj.VN=length(obj.Value);
            % 坐标区域基础设置 ---------------------------------------------
            obj.ax.YDir='reverse';
            obj.ax.XColor='none';
            obj.ax.YColor='none';

            
        end
% 绘图函数 =================================================================
        function draw(obj, props, probes)
            % store props and probes for later updates (refresh/move)
            if nargin < 2 || isempty(props)
                obj.NodeProps = cell(1,obj.BN);
            else
                obj.NodeProps = props;
            end
            if nargin < 3 || isempty(probes)
                obj.NodeProbes = cell(1,obj.BN);
            else
                obj.NodeProbes = probes;
            end

            % 生成整体邻接矩阵 ---------------------------------------------
            obj.getAdjMat()
            % help SSankey
            obj.BoolMat=abs(obj.AdjMat)>0;
            if any(any(obj.BoolMat+obj.BoolMat.'==2))
                warning('Currently, bidirectional flow sankey diagram plotting is not supported.')
            end
            obj.VN=sum(sum(obj.BoolMat));
            % 计算每个对象位于的层、每层方块长度、每个方块位置 ----------------
            if isempty(obj.Layer)
                obj.getLayer()
            end
            obj.getLayerPos()
            % Draw connections  -----------------------------------------------------
            for i=1:obj.VN
                obj.drawLink(i)
            end
            % Draw rectangles  -----------------------------------------------------
            for i=1:obj.BN
                obj.drawNode(i, obj.NodeProps, obj.NodeProbes)
            end
            % -------------------------------------------------------------
            axis tight;
        end
% =========================================================================
        function setBlock(obj,n,varargin)
            % apply set(...) to all segment patches of node n
            segs = {};
            if iscell(obj.BlockHdl) && numel(obj.BlockHdl) >= n
                segs = obj.BlockHdl{n};
            elseif ~isempty(obj.BlockHdl) && numel(obj.BlockHdl) >= n
                segs = obj.BlockHdl(n);
            end
            for s = 1:numel(segs)
                try set(segs(s),varargin{:}); catch; end
            end
        end
        function setLink(obj,n,varargin)
            set(obj.LinkHdl(n),varargin{:})
        end
        function setLabel(obj,n,varargin)
            set(obj.LabelHdl(n),varargin{:})
        end
        function setValueLabel(obj,n,varargin)
            set(obj.ValueLabelHdl(n),varargin{:})
        end
% =========================================================================
        function addLink(obj,S,T,V)
            obj.getAdjMat()
            if isempty(obj.BlockHdl)
                obj.AdjMat(S,T)=obj.AdjMat(S,T)+abs(V);
            else
                if obj.AdjMat(S,T)==0
                    obj.AdjMat(S,T)=obj.AdjMat(S,T)+abs(V);
                    obj.getLayerPos()
                    [M,N]=find(obj.AdjMat~=0);
                    obj.drawLink(find(M==S&N==T))
                else
                    obj.AdjMat(S,T)=obj.AdjMat(S,T)+abs(V);
                    obj.getLayerPos()
                end
                % disp(obj.AdjMat)
                obj.refresh()
            end
        end
        function addNode(obj,name,layer)
            obj.getAdjMat()
            obj.AdjMat(end+1,:)=0;obj.AdjMat(:,end+1)=0;
            if nargin<2
                obj.NodeList{end+1}=compose('node%d',size(obj.AdjMat,1));
            else
                obj.NodeList{end+1}=name;
            end
            obj.BN=length(obj.NodeList);
            obj.BoolMat=abs(obj.AdjMat)>0;
            if any(any(obj.BoolMat+obj.BoolMat.'==2))
                warning('Currently, bidirectional flow sankey diagram plotting is not supported.')
            end
            obj.VN=sum(sum(obj.BoolMat));
            if isempty(obj.Layer)
                obj.getLayer()
                if nargin<3,obj.Layer(end)=max(obj.Layer);else,obj.Layer(end)=layer;end
            else
                if nargin<3,obj.Layer(end+1)=max(obj.Layer);else,obj.Layer(end+1)=layer;end
            end
            obj.dst(end+1,:)=rand(1,3).*.7;
            obj.MovePos(end+1,:)=0;
            % expand BlockHdl, NodeProps, NodeProbes
            if ~iscell(obj.BlockHdl)
                tmp = cell(1,obj.BN-1);
                for ii = 1:min(numel(obj.BlockHdl),numel(tmp))
                    tmp{ii} = obj.BlockHdl(ii);
                end
                obj.BlockHdl = tmp;
            else
                obj.BlockHdl{end+1} = [];
            end
            if numel(obj.NodeProps) < obj.BN
                obj.NodeProps{obj.BN} = [];
            end
            if numel(obj.NodeProbes) < obj.BN
                obj.NodeProbes{obj.BN} = [];
            end
            % -------------------------------------------------------------
            if isempty(obj.BlockHdl)
            else
                obj.getLayerPos()
                obj.drawNode(length(obj.NodeList), obj.NodeProps, obj.NodeProbes)
                N=find(obj.Layer==obj.Layer(end));
                for n=1:length(N)
                    obj.moveBlock(N(n))
                end
            end
        end
% =========================================================================
        function refresh(obj)
            tLayerPos = obj.MovePos + obj.LayerPos;
            display(obj.LayerPos)
            obj.BoolMat = abs(obj.AdjMat) > 0;
            if any(any(obj.BoolMat + obj.BoolMat.' == 2))
                warning('Currently, bidirectional flow sankey diagram plotting is not supported.')
            end
            obj.VN = sum(sum(obj.BoolMat));

            % update block patches (now possibly multiple segments per node)
            for n = 1:obj.BN
                % update label position
                switch obj.LabelLocation
                    case 'right', set(obj.LabelHdl(n),'Position',[tLayerPos(n,2),mean(tLayerPos(n,[3,4]))]);
                    case 'left',  set(obj.LabelHdl(n),'Position',[tLayerPos(n,1),mean(tLayerPos(n,[3,4]))]);
                    case 'top',   set(obj.LabelHdl(n),'Position',[mean(tLayerPos(n,[1,2])),tLayerPos(n,3)]);
                    case 'center',set(obj.LabelHdl(n),'Position',[mean(tLayerPos(n,[1,2])),mean(tLayerPos(n,[3,4]))]);
                    case 'bottom',set(obj.LabelHdl(n),'Position',[mean(tLayerPos(n,[1,2])),tLayerPos(n,4)]);
                end

                % find segment handles for node n
                segs = {};
                if iscell(obj.BlockHdl) && numel(obj.BlockHdl) >= n
                    segs = obj.BlockHdl{n};
                elseif ~isempty(obj.BlockHdl) && numel(obj.BlockHdl) >= n
                    segs = obj.BlockHdl(n);
                end

                if ~isempty(segs)
                    % get stored proportions (if available), fallback to equal segments
                    if ~isempty(obj.NodeProps) && numel(obj.NodeProps) >= n && ~isempty(obj.NodeProps{n})
                        p = double(obj.NodeProps{n}(:));
                        if sum(p) > eps
                            p = p / sum(p);
                        else
                            p = ones(size(p))/numel(p);
                        end
                    else
                        p = ones(1,numel(segs))/numel(segs);
                    end

                    x = tLayerPos(n,[1,2,2,1]);
                    y_current = tLayerPos(n,3);
                    len = tLayerPos(n,4) - tLayerPos(n,3);

                    for m = 1:numel(segs)
                        y1 = y_current;
                        y2 = y_current + p(m) * len;
                        yy = [y1, y1, y2, y2];
                        try
                            set(segs(m),'XData',x,'YData',yy);
                        catch
                            % ignore if handle deleted or invalid
                        end
                        y_current = y2;
                    end
                end
            end

            [obj.SourceInd,obj.TargetInd] = find(obj.AdjMat~=0);
            for n = 1:obj.VN
                tSource = obj.SourceInd(n);
                tTarget = obj.TargetInd(n);

                % compute source/target vertical offsets normalized inside node extents
                src_top = tLayerPos(tSource,3);
                src_bottom = tLayerPos(tSource,4);
                src_h = src_bottom - src_top;
                tgt_top = tLayerPos(tTarget,3);
                tgt_bottom = tLayerPos(tTarget,4);
                tgt_h = tgt_bottom - tgt_top;

                rowOut = sum(obj.AdjMat(tSource,:));
                if rowOut <= eps
                    tS1 = src_top + src_h/2;
                    tS2 = tS1;
                else
                    tS1 = src_top + (sum(obj.AdjMat(tSource,1:(tTarget-1))) / rowOut) * src_h;
                    tS2 = tS1 + (obj.AdjMat(tSource,tTarget) / rowOut) * src_h;
                end

                colIn = sum(obj.AdjMat(:,tTarget));
                if colIn <= eps
                    tT1 = tgt_top + tgt_h/2;
                    tT2 = tT1;
                else
                    tT1 = tgt_top + (sum(obj.AdjMat(1:(tSource-1),tTarget)) / colIn) * tgt_h;
                    tT2 = tT1 + (obj.AdjMat(tSource,tTarget) / colIn) * tgt_h;
                end

                tX = [tLayerPos(tSource,1), tLayerPos(tSource,2), tLayerPos(tTarget,1), tLayerPos(tTarget,2)];
                qX = linspace(tLayerPos(tSource,1), tLayerPos(tTarget,2), 200);
                qT = linspace(0,1,50);
                qY1 = interp1(tX, [tS1, tS1, tT1, tT1], qX, 'pchip');
                qY2 = interp1(tX, [tS2, tS2, tT2, tT2], qX, 'pchip');
                YY = qY1 .* (qT'.*0+1) + (qY2 - qY1) .* (qT');

                try
                    set(obj.LinkHdl(n),'YData',YY,'XData',qX);
                    set(obj.ValueLabelHdl(n),'String',[' ',obj.ValueLabelFormat(obj.AdjMat(obj.SourceInd(n),obj.TargetInd(n)))]);
                catch
                    % handle missing link handles gracefully
                end

                switch obj.ValueLabelLocation
                    case 'left'
                        set(obj.ValueLabelHdl(n),'Position',[tLayerPos(tSource,2),(tS1+tS2)/2]);
                    case 'right'
                        set(obj.ValueLabelHdl(n),'Position',[tLayerPos(tTarget,1),(tT1+tT2)/2]);
                    case 'center'
                        set(obj.ValueLabelHdl(n),'Position',[(tLayerPos(tSource,2)+tLayerPos(tTarget,1))/2, (tS1+tS2+tT1+tT2)/4]);
                    case 'none'
                        set(obj.ValueLabelHdl(n),'Position',[tLayerPos(tSource,2),(tS1+tS2)/2]);
                end
            end
        end

        function drawLink(obj,n)
            % Draw a single link (n) using normalized offsets inside source/target nodes
            [obj.SourceInd,obj.TargetInd] = find(obj.AdjMat~=0);
            tSource = obj.SourceInd(n);
            tTarget = obj.TargetInd(n);

            % node vertical extent and height
            src_top = obj.LayerPos(tSource,3);
            src_bottom = obj.LayerPos(tSource,4);
            src_h = src_bottom - src_top;
            tgt_top = obj.LayerPos(tTarget,3);
            tgt_bottom = obj.LayerPos(tTarget,4);
            tgt_h = tgt_bottom - tgt_top;

            % compute cumulative offsets normalized to node height
            rowOut = sum(obj.AdjMat(tSource,:));
            if rowOut <= eps
                % no outgoing mass: place link at middle of node
                tS1 = src_top + src_h/2;
                tS2 = tS1;
            else
                leftCum = sum(obj.AdjMat(tSource,1:(tTarget-1)));
                tS1 = src_top + (leftCum / rowOut) * src_h;
                tS2 = tS1 + (obj.AdjMat(tSource,tTarget) / rowOut) * src_h;
            end

            colIn = sum(obj.AdjMat(:,tTarget));
            if colIn <= eps
                tT1 = tgt_top + tgt_h/2;
                tT2 = tT1;
            else
                upCum = sum(obj.AdjMat(1:(tSource-1),tTarget));
                tT1 = tgt_top + (upCum / colIn) * tgt_h;
                tT2 = tT1 + (obj.AdjMat(tSource,tTarget) / colIn) * tgt_h;
            end

            % prepare X and Y for the interpolation
            tX = [obj.LayerPos(tSource,1), obj.LayerPos(tSource,2), obj.LayerPos(tTarget,1), obj.LayerPos(tTarget,2)];
            if abs(tX(1)-tX(3)) < eps
                warning('Currently, flow between the same layer is not supported.')
            end
            qX = linspace(obj.LayerPos(tSource,1), obj.LayerPos(tTarget,2), 200);
            qT = linspace(0,1,50);

            qY1 = interp1(tX, [tS1, tS1, tT1, tT1], qX, 'pchip');
            qY2 = interp1(tX, [tS2, tS2, tT2, tT2], qX, 'pchip');

            XX = repmat(qX, [50,1]);
            YY = qY1 .* (qT'.*0+1) + (qY2 - qY1) .* (qT');


            if ~isempty(obj.EdgeDiff)
                range = -30:30;
                d = obj.EdgeDiff(n);
                v = round(d*100);
                if v > 30 
                    v = 30;
                elseif v < -30
                    v = -30;
                end

                    
                idx = find(range==v);
                %dnorm = (d - min(obj.EdgeDiff)) / (max(obj.EdgeDiff) - min(obj.EdgeDiff) + eps);
            
                % pick RGB from colormap
                %idx = max(1, round(dnorm * size(obj.EdgeDiffCmap,1)));
                col = obj.EdgeDiffCmap(idx,:);
                MeshC = ones(50,200,3);
                MeshC(:,:,1) = col(1);
                MeshC(:,:,2) = col(2);
                MeshC(:,:,3) = col(3);
            else
                % fallback to original behaviour
                MeshC = ones(50,200,3);
                switch obj.RenderingMethod
                    case 'left'
                        MeshC(:,:,1) = MeshC(:,:,1) .* obj.ColorList(tSource,1);
                        MeshC(:,:,2) = MeshC(:,:,2) .* obj.ColorList(tSource,2);
                        MeshC(:,:,3) = MeshC(:,:,3) .* obj.ColorList(tSource,3);
                    case 'right'
                        MeshC(:,:,1) = MeshC(:,:,1) .* obj.ColorList(tTarget,1);
                        MeshC(:,:,2) = MeshC(:,:,2) .* obj.ColorList(tTarget,2);
                        MeshC(:,:,3) = MeshC(:,:,3) .* obj.ColorList(tTarget,3);
                    case 'interp'
                        MeshC(:,:,1) = repmat(linspace(obj.ColorList(tSource,1),obj.ColorList(tTarget,1),200),[50,1]);
                        MeshC(:,:,2) = repmat(linspace(obj.ColorList(tSource,2),obj.ColorList(tTarget,2),200),[50,1]);
                        MeshC(:,:,3) = repmat(linspace(obj.ColorList(tSource,3),obj.ColorList(tTarget,3),200),[50,1]);
                    case 'map'
                        % if you want to keep 'map' add it here
                end
            end

            
            tLinkHdl = surf(obj.ax, XX, YY, XX.*0, 'EdgeColor', 'none', 'FaceAlpha', .3, 'CData', MeshC, ...
                'UserData', n, 'ButtonDownFcn', @obj.onLinkClick);
            obj.LinkHdl = [obj.LinkHdl(1:n-1), tLinkHdl, obj.LinkHdl(n:end)];

            % Value label placement
            switch obj.ValueLabelLocation
                case 'left'
                    tValueLabelHdl = text(obj.LayerPos(tSource,2), (tS1+tS2)/2, [' ', obj.ValueLabelFormat(obj.AdjMat(tSource,tTarget))], ...
                        'FontSize',12,'FontName','Times New Roman','HorizontalAlignment','left');
                case 'right'
                    tValueLabelHdl = text(obj.LayerPos(tTarget,1), (tT1+tT2)/2, [obj.ValueLabelFormat(obj.AdjMat(tSource,tTarget)),' '], ...
                        'FontSize',12,'FontName','Times New Roman','HorizontalAlignment','right');
                case 'center'
                    tValueLabelHdl = text((obj.LayerPos(tSource,2)+obj.LayerPos(tTarget,1))/2, (tS1+tS2+tT1+tT2)/4, ...
                        obj.ValueLabelFormat(obj.AdjMat(tSource,tTarget)), 'FontSize',12,'FontName','Times New Roman','HorizontalAlignment','center');
                case 'none'
                    tValueLabelHdl = text(obj.LayerPos(tSource,2), (tS1+tS2)/2, [' ',obj.ValueLabelFormat(obj.AdjMat(tSource,tTarget))], ...
                        'FontSize',12,'FontName','Times New Roman','HorizontalAlignment','left', 'Visible','off');
            end
            obj.ValueLabelHdl = [obj.ValueLabelHdl(1:n-1), tValueLabelHdl, obj.ValueLabelHdl(n:end)];
        end

        function drawNode(obj,n, props, probes)
            % Draw node n as stacked horizontal bands whose heights
            % are given by props{n} (should sum to 1, but we normalize).
            %
            % props and probes are accepted as inputs for initial draw,
            % but the function also uses obj.NodeProps/NodeProbes (saved).

            if nargin < 3 || isempty(props)
                p = [];
            else
                % props passed as full cell array; extract entry n safely
                if iscell(props) && numel(props) >= n
                    p = props{n};
                else
                    p = [];
                end
            end
            if nargin < 4 || isempty(probes)
                pr = [];
            else
                if iscell(probes) && numel(probes) >= n
                    pr = probes{n};
                else
                    pr = [];
                end
            end

            % fallback to saved props/probes if empty
            if isempty(p)
                if ~isempty(obj.NodeProps) && numel(obj.NodeProps) >= n
                    p = obj.NodeProps{n};
                else
                    p = [];
                end
            end
            if isempty(pr)
                if ~isempty(obj.NodeProbes) && numel(obj.NodeProbes) >= n
                    pr = obj.NodeProbes{n};
                else
                    pr = [];
                end
            end

            % defensive: ensure numeric and non-negative
            if isempty(p)
                p = [];
            else
                p = double(p(:));
                p(p<0) = 0;
            end

            % if all zeros or empty, make single full segment
            if isempty(p) || all(p==0)
                p = 1;
                pr = 1;
            end

            % normalize so they sum to 1 (if they are proportions)
            total = sum(p);
            if total <= eps
                p = ones(size(p))/numel(p);
            else
                p = p / total;
            end

            % block vertical extent
            y_top = obj.LayerPos(n,3);   % top (start)
            y_bottom = obj.LayerPos(n,4);% bottom (end)
            len = y_bottom - y_top;      % note YDir='reverse' earlier; this produces correct signed length

            % We'll stack from top (y_top) downward adding p(i)*len each step.
            y_current = y_top;

            % create array to collect handles for this node
            segHandles = gobjects(1,numel(p));
            x = obj.LayerPos(n,[1,2,2,1]); % constant x for all segments

            for m = 1:numel(p)
                seg_h = p(m) * len;
                y1 = y_current;
                y2 = y_current + seg_h;

                % fill expects x and y vectors of same length
                xx = x;
                yy = [y1, y1, y2, y2];

  
                colidx = [];
                if ~isempty(pr) 
                    if pr(m) == 0
                        colidx = 0;
                    else
                        colidx = find(obj.values == pr(m));
                    end
                end
                if colidx == 0
                    col = [0 0 0];
                else
                    col = obj.pColorList(colidx,:);
                end

                segHandles(m) = fill(obj.ax, xx, yy, col, 'EdgeColor', 'none');

                % advance
                y_current = y2;
            end

            % store handles in BlockHdl as a cell entry (so one cell per node)
            if ~iscell(obj.BlockHdl)
                % convert to cell preserving any existing scalar entries
                old = obj.BlockHdl;
                tmp = cell(1,max(obj.BN,numel(old)));
                for i = 1:numel(old)
                    tmp{i} = old(i);
                end
                obj.BlockHdl = tmp;
            end
            % ensure BlockHdl long enough
            if numel(obj.BlockHdl) < n
                obj.BlockHdl{n} = segHandles;
            else
                % delete old handles first (if any)
                try
                    oldh = obj.BlockHdl{n};
                    if ~isempty(oldh)
                        delete(oldh(ishandle(oldh)));
                    end
                catch
                end
                obj.BlockHdl{n} = segHandles;
            end

            % draw label (unchanged)

            switch obj.LabelLocation
                case 'right'
                    obj.LabelHdl(n)=text(obj.ax,obj.LayerPos(n,2),mean(obj.LayerPos(n,[3,4])),...
                        ['  ',regexp(obj.NodeList{n},"^[0-9]+","match")],'FontSize',15,'FontName','Times New Roman','HorizontalAlignment','left');
                case 'left'
                    obj.LabelHdl(n)=text(obj.ax,obj.LayerPos(n,1),mean(obj.LayerPos(n,[3,4])),...
                        [obj.NodeList{n},' '],'FontSize',15,'FontName','Times New Roman','HorizontalAlignment','right');
                case 'top'
                    obj.LabelHdl(n)=text(obj.ax,mean(obj.LayerPos(n,[1,2])),obj.LayerPos(n,3),...
                        obj.NodeList{n},'FontSize',15,'FontName','Times New Roman','HorizontalAlignment','center','VerticalAlignment','bottom');
                case 'center'
                    obj.LabelHdl(n)=text(obj.ax,mean(obj.LayerPos(n,[1,2])),mean(obj.LayerPos(n,[3,4])),...
                        obj.NodeList{n},'FontSize',15,'FontName','Times New Roman','HorizontalAlignment','center');
                case 'bottom'
                    obj.LabelHdl(n)=text(obj.ax,mean(obj.LayerPos(n,[1,2])),obj.LayerPos(n,4),...
                        obj.NodeList{n},'FontSize',15,'FontName','Times New Roman','HorizontalAlignment','center','VerticalAlignment','top');
            end
        end
% =========================================================================
        function getAdjMat(obj)
            if isempty(obj.AdjMat)
                obj.AdjMat=zeros(obj.BN,obj.BN);
                for i=1:length(obj.Source)
                    obj.SourceInd(i)=find(strcmp(obj.Source{i},obj.NodeList));
                    obj.TargetInd(i)=find(strcmp(obj.Target{i},obj.NodeList));
                    obj.AdjMat(obj.SourceInd(i),obj.TargetInd(i))=obj.Value{i};
                end
            end
        end
        function getLayer(obj)
            if strcmp(obj.LayerOrder,'normal')
                obj.Layer=zeros(obj.BN,1);
                obj.Layer(sum(obj.BoolMat,1)==0)=1;
                startMat=diag(obj.Layer);
                for i=1:(obj.BN-1)
                    tLayer=(sum(startMat*obj.BoolMat^i,1)>0).*(i+1);
                    obj.Layer=max([obj.Layer,tLayer'],[],2);
                end
            else
                obj.Layer=zeros(obj.BN,1);
                obj.Layer(sum(obj.BoolMat,2)==0)=-1;
                startMat=diag(obj.Layer);
                for i=1:(obj.BN-1)
                    tLayer=(sum(startMat*(obj.BoolMat.')^i,1)<0).*(-i-1);
                    obj.Layer=min([obj.Layer,tLayer'],[],2);
                end
                obj.Layer=obj.Layer-min(obj.Layer)+1;
            end
        end
        function getLayerPos(obj)
            % Ensure Layer is column
            obj.Layer = obj.Layer(:);
            obj.LN = max(obj.Layer);
        
            % NodeHeight is BN×1 vector
            nodeH  = obj.NodeHeight(:);
        
            % sep is scalar; sepLen must match nodeH
            sep         = max(1, obj.Sep);
            sepLen      = 0.5 * nodeH;   % same size as nodeH
        
            % Output container
            obj.LayerPos = zeros(obj.BN,4); % [x1 x2 y_top y_bottom]
        
            % === LAYER-BY-LAYER VERTICAL PLACEMENT ===
            for i = 1:obj.LN
                tBlockInd = find(obj.Layer == i);
                nNodes = numel(tBlockInd);
        
                % Extract per-node heights and gaps
                nodeH_layer  = nodeH(tBlockInd);     % n×1
                sep_layer    = sepLen(tBlockInd);    % n×1
        
                % Compute starts using cumulative sum (safe, monotonic)
                % For node k:  start = sum_{j<k} (nodeH_j + sep_j)
                if nNodes == 1
                    starts = 0;
                else
                    starts = cumsum([0; nodeH_layer(1:end-1) + sep_layer(1:end-1)]);
                end
        
                ends = starts + nodeH_layer;
        
                % Store in output
                obj.LayerPos(tBlockInd,3) = starts;
                obj.LayerPos(tBlockInd,4) = ends;
            end
        
            % === HORIZONTAL POSITIONS ===
            obj.LayerPos(:,1) = obj.Layer;
            obj.LayerPos(:,2) = obj.Layer + obj.BlockScale';
        
            % === ALIGNMENT ===
            tMinY = min(obj.LayerPos(:,3));
            tMaxY = max(obj.LayerPos(:,4));
        
            for i = 1:obj.LN
                tBlockInd = find(obj.Layer == i);
                pos3 = obj.LayerPos(tBlockInd,3);
                pos4 = obj.LayerPos(tBlockInd,4);
        
                switch obj.Align
                    case 'up'
                        % no shift
        
                    case 'down'
                        shift = (tMaxY - max(pos4));
                        obj.LayerPos(tBlockInd,3) = pos3 + shift;
                        obj.LayerPos(tBlockInd,4) = pos4 + shift;
        
                    case 'center'
                        globalMid = (tMinY + tMaxY) / 2;
                        layerMid  = (min(pos3) + max(pos4)) / 2;
                        shift     = globalMid - layerMid;
                        obj.LayerPos(tBlockInd,3) = pos3 + shift;
                        obj.LayerPos(tBlockInd,4) = pos4 + shift;
                end
            end
        end


% =========================================================================
        function moveBlock(obj,n)
            tLayerPos=obj.MovePos+obj.LayerPos;

            % update block segment patches for this node n
            segs = {};
            if iscell(obj.BlockHdl) && numel(obj.BlockHdl) >= n
                segs = obj.BlockHdl{n};
            elseif ~isempty(obj.BlockHdl) && numel(obj.BlockHdl) >= n
                segs = obj.BlockHdl(n);
            end

            if ~isempty(segs)
                if ~isempty(obj.NodeProps) && numel(obj.NodeProps) >= n && ~isempty(obj.NodeProps{n})
                    p = double(obj.NodeProps{n}(:));
                    if sum(p) > eps
                        p = p / sum(p);
                    else
                        p = ones(size(p))/numel(p);
                    end
                else
                    p = ones(1,numel(segs))/numel(segs);
                end

                x = tLayerPos(n,[1,2,2,1]);
                y_current = tLayerPos(n,3);
                len = tLayerPos(n,4) - tLayerPos(n,3);

                for m = 1:numel(segs)
                    y1 = y_current;
                    y2 = y_current + p(m) * len;
                    yy = [y1, y1, y2, y2];
                    try
                        set(segs(m),'XData',x,'YData',yy);
                    catch
                    end
                    y_current = y2;
                end
            end

            switch obj.LabelLocation
                case 'right',set(obj.LabelHdl(n),'Position',[tLayerPos(n,2),mean(tLayerPos(n,[3,4]))]);
                case 'left',set(obj.LabelHdl(n),'Position',[tLayerPos(n,1),mean(tLayerPos(n,[3,4]))]);
                case 'top',set(obj.LabelHdl(n),'Position',[mean(tLayerPos(n,[1,2])),tLayerPos(n,3)]);
                case 'center',set(obj.LabelHdl(n),'Position',[mean(tLayerPos(n,[1,2])),mean(tLayerPos(n,[3,4]))]);
                case 'bottom',set(obj.LabelHdl(n),'Position',[mean(tLayerPos(n,[1,2])),tLayerPos(n,4)]);
            end
            for i=1:obj.VN
                tSource=obj.SourceInd(i);
                tTarget=obj.TargetInd(i);
                if tSource==n||tTarget==n
                    tS1=sum(obj.AdjMat(tSource,1:(tTarget-1)))+tLayerPos(tSource,3);
                    tS2=sum(obj.AdjMat(tSource,1:tTarget))+tLayerPos(tSource,3);
                    tT1=sum(obj.AdjMat(1:(tSource-1),tTarget))+tLayerPos(tTarget,3);
                    tT2=sum(obj.AdjMat(1:tSource,tTarget))+tLayerPos(tTarget,3);
                    if isempty(tS1),tS1=0;end
                    if isempty(tT1),tT1=0;end
                    tX=[tLayerPos(tSource,1),tLayerPos(tSource,2),tLayerPos(tTarget,1),tLayerPos(tTarget,2)];
                    qX=linspace(tLayerPos(tSource,1),tLayerPos(tTarget,2),200);qT=linspace(0,1,50);
                    qY1=interp1(tX,[tS1,tS1,tT1,tT1],qX,'pchip');
                    qY2=interp1(tX,[tS2,tS2,tT2,tT2],qX,'pchip');
                    YY=qY1.*(qT'.*0+1)+(qY2-qY1).*(qT');
                    set(obj.LinkHdl(i),'YData',YY,'XData',qX);
                    switch obj.ValueLabelLocation
                        case 'left'
                            set(obj.ValueLabelHdl(i),'Position',[tLayerPos(tSource,2),tS1/2+tS2/2]);
                        case 'right'
                            set(obj.ValueLabelHdl(i),'Position',[tLayerPos(tTarget,1),tT1/2+tT2/2]);
                        case 'center'
                             set(obj.ValueLabelHdl(i),'Position',[tLayerPos(tSource,2)/2+tLayerPos(tTarget,1)/2,tS1/4+tS2/4+tT1/4+tT2/4]);
                        case 'none'
                            set(obj.ValueLabelHdl(i),'Position',[tLayerPos(tSource,2),tS1/2+tS2/2]);
                    end
                end
            end
        end
        function moveBlockX(obj,n,dx)
            obj.MovePos(n,[1,2])=obj.MovePos(n,[1,2])+dx;
            obj.moveBlock(n)
        end
        function moveBlockY(obj,n,dy)
            obj.MovePos(n,[3,4])=obj.MovePos(n,[3,4])-dy;
            obj.moveBlock(n)
        end
        function onLinkClick(obj, src, event)
            if ~verLessThan('matlab', '9.7')
            if event.Button == 1
                src.FaceAlpha = obj.dataTipFormat{1};
                datatip(src, event.IntersectionPoint(1), event.IntersectionPoint(2));
                n = src.UserData;
                src.DataTipTemplate.DataTipRows(1) = ...
                dataTipTextRow(obj.dataTipFormat{2}, repmat(obj.NodeList(obj.SourceInd(n)), length(src.XData), length(src.YData)));
                src.DataTipTemplate.DataTipRows(2) = ...
                dataTipTextRow(obj.dataTipFormat{3}, repmat(obj.NodeList(obj.TargetInd(n)), length(src.XData), length(src.YData)));
                src.DataTipTemplate.DataTipRows(3) = ...
                dataTipTextRow(obj.dataTipFormat{4}, repmat(obj.AdjMat(obj.SourceInd(n),obj.TargetInd(n)), ...
                [length(src.XData), length(src.YData)]), obj.dataTipFormat{5});
            else
                src.FaceAlpha = 0.3;
            end
            end
        end
    end
% Copyright (c) 2023-2025, Zhaoxu Liu / slandarer
% =========================================================================
% @author : slandarer
% 公众号  : slandarer随笔
% 知乎    : slandarer
% -------------------------------------------------------------------------
end
