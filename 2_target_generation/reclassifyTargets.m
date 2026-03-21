function maskWindow = reclassifyTargets(maskWindow)
wsize = length(maskWindow);
uCells = cell2mat(cellfun(@(x) unique(x),maskWindow,'UniformOutput',false));
uCells = unique(uCells);
uCells(uCells==0) = [];

overlap = zeros(length(uCells));
n = zeros(length(uCells));
for i = 1:wsize
    for j = 1:size(overlap,1)
        if any(maskWindow{i}==uCells(j),"all")
            mask1 = maskWindow{i} == uCells(j);

            for k = 1:wsize
                for l = 1:size(overlap,2)
                    if any(maskWindow{k}==uCells(l),"all")
                        mask2 = maskWindow{k} == uCells(l);
                        n(j,l) = n(j,l) + 1;
                        tt = sum(mask1 == 1 & mask2 == 1,"all");
                        overlap(j,l) = overlap(j,l) * (n(j,l)-1)/n(j,l) + (((tt)/sum(mask2,"all"))*((tt)/sum(mask1,"all")))/n(j,l);
                    end
                end
            end
        end
    end

end

othreshold = 0.65;
for i = 1:length(uCells)
    if ~all(overlap(i,:)==0)
        [~,idx] = find(overlap(i,:)>=othreshold);
        mergeSet = uCells(idx);

        for j = 1:length(mergeSet)
            posIDs = [];
            c = 1;
            for k = 1:wsize
                maskWindow{k}(maskWindow{k}==mergeSet(j)) = uCells(i);
                if ~isempty(maskWindow{k}(maskWindow{k}==uCells(i)))
                    posIDs(c) = k;
                    c = c + 1;
                end
                overlap(idx(j),:) = repelem(0,size(overlap,2));
            end
            % for l = 2:length(posIDs)
            %     if posIDs(l)-posIDs(l-1)>1
            %         for m = 1:posIDs(l)-posIDs(l-1)
            %             maskWindow{posIDs(l-1)+m}(maskWindow{posIDs(l-1)}==uCells(i)) = uCells(i);
            %         end
            %     end
            % end


        end
    end
end
        