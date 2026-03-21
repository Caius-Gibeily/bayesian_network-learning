function [clipindices,cliptab] = getClips(Participants,SessionNumbers,thresh)

participantData = cell(size(Participants,1),1); % each row is one participant's matfile for this clip
% Calculate the number of substrings

clipindices = cell(length(participantData),1);
for i = 1:length(participantData)
    wholeMatfile = load(strcat('C:\Users\cgibeil@emory.edu\OneDrive - Emory\Documents\PhD\Projects\P1-CondProbs\matfiles\',Participants(i,:),'_',common_getpaddednumber(SessionNumbers(i),2),'/',Participants(i,:),'_',common_getpaddednumber(SessionNumbers(i),2),'_data.mat'));
    matfileData = wholeMatfile.(strcat(Participants(i,:),'_data'));
    prefBinData = wholeMatfile.(strcat(Participants(i,:),'_prefbin'));
    clipindices{i} = prefBinData.WhichClips;
end
clipall = cell2mat(clipindices);
[n,cliptab] = groupcounts(clipall);

figure
[n,inds] = sort(n,"descend");
cliptab = cliptab(inds);
cliptab(n<thresh) = [];
n(n<thresh) = [];
bar(string(cliptab),n)
shg

