function [indivRegions2] = smoothTargetsWindow(indivRegions,cutoffVal,winsize)

frames = length(indivRegions{1});
indivRegions2 = cell(length(indivRegions),1);

for participant = 1:length(indivRegions)
    indivRegions2{participant} = zeros(frames,1);
    pad = zeros(winsize,1)+indivRegions{participant}(1,1);
    indivRegions_padded = cat(1,pad,indivRegions{participant}(:,1));
    for frame = winsize+1:frames+winsize
        
        %dat = indivRegions{participant}(frame-winsize:frame,1)+1;
        dat = indivRegions_padded(frame-winsize:frame);
        tabulated = tabulate(dat);
        tabulated = tabulated(tabulated(:,1)>0,:);
        [maxVal, ROInd] = max(tabulated(:, 2));
        prop = maxVal/sum(tabulated(:,2));

        
        if prop >= cutoffVal
            indivRegions2{participant}(frame-winsize) = tabulated(ROInd,1);
        elseif prop < cutoffVal
            if frame ~= frames+winsize
                if indivRegions_padded(frame) ~= indivRegions_padded(frame-1) && indivRegions_padded(frame-1) == indivRegions_padded(frame+1)
                    indivRegions_padded(frame) = indivRegions_padded(frame-1);
                elseif indivRegions_padded(frame) ~= indivRegions_padded(frame-1) && indivRegions_padded(frame) ~= indivRegions_padded(frame+1)
                    indivRegions_padded(frame) = indivRegions_padded(frame-1);
                elseif indivRegions_padded(frame) ~= indivRegions_padded(frame-1) && indivRegions_padded(frame) == indivRegions_padded(frame+1)
                    indivRegions_padded(frame) = indivRegions_padded(frame+1);
                end
            end
            indivRegions2{participant}(frame-winsize) = indivRegions_padded(frame);
        end
    end

end
