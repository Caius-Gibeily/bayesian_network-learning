function exportPatches(clip, outdir, viddir, targetCell, sampleInt)


video = VideoReader(strcat(viddir,"0",string(clip),"PEER_movie.mov"));
masks = targetCell; 

uVals = cellfun(@(x) unique(x(x>0)),masks,'UniformOutput',false);
targetset = unique(cell2mat(uVals));
firstframe = zeros(length(targetset),3);
firstframe(:,1) = targetset;
for i = 1:length(uVals)
    vals = uVals{i};
    for j = 1:length(vals)
        firstframe(firstframe(:,1)==vals(j),3) = i;
        if firstframe(firstframe(:,1)==vals(j),2)==0
           firstframe(firstframe(:,1)==vals(j),2) = i;
        end
    end
end

for i = 1:size(firstframe,1)

    for k = firstframe(i,2):sampleInt:firstframe(i,3)
        mask = masks{k};
        frame = read(video, k);
        
        mask(mask~=firstframe(i,1)) = 0;
        mask(mask~=0) = 1;

        [rowids, colids] = find(mask == 1);

        overlay = im2double(frame);

        patch = overlay(min(rowids):max(rowids), ...
            min(colids):max(colids),:);
        
        imwrite(patch, strcat(outdir,string(clip), "_", string(firstframe(i)), "-", string(k), ".png"))
    end
    display(strcat("Target: ",string(firstframe(i,1))))
end
end
