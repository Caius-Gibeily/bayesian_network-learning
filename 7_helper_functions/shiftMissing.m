function meanmiss = shiftMissing(indivArray,nshifts)


arguments
    indivArray
    nshifts = 1000
end

meanmiss = zeros([nshifts,size(indivArray,2)]);
for n = 1:nshifts
    shiftedArray = zeros(size(indivArray));
    for p = 1:size(indivArray,1)
        shiftedArray(p,:) = circshift(indivArray(p,:),randi([0,size(indivArray,2)]));
    end
    meanmiss(n,:) = mean(shiftedArray==0,1);
end


end