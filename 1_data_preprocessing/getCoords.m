function [coords,ParticipantsAdj] = getCoords(Participants,SessionNumbers,ClipNum)

participantData = cell(size(Participants,1),1); % each row is one participant's matfile for this clip
% Calculate the number of substrings
presentlist = zeros(size(Participants,1),1);
for i = 1:length(participantData)
    wholeMatfile = load(strcat('C:\Users\cgibeil@emory.edu\OneDrive - Emory\Documents\PhD\Projects\P1-CondProbs\matfiles\',Participants(i,:),'_',common_getpaddednumber(SessionNumbers(i),2),'/',Participants(i,:),'_',common_getpaddednumber(SessionNumbers(i),2),'_data.mat'));
    matfileData = wholeMatfile.(strcat(Participants(i,:),'_data'));
    prefBinData = wholeMatfile.(strcat(Participants(i,:),'_prefbin'));
    clipindex = prefBinData.WhichClips==ClipNum;
    if ~all(clipindex==0)
        participantData{i} = matfileData{clipindex};
        presentlist(i) = 1;
    end
end
presentlist = logical(presentlist);
Participants = string(Participants);
ParticipantsAdj = Participants(presentlist,:);
participantData = participantData(presentlist,:);
coords = cell(length(participantData),1);


for c = 1:length(participantData)
    fprintf("\nParticipant "+c)
    coords{c} = zeros(length(participantData{1}),2);
    
    scorrected = participantData{c}(:,7:8);
    
    for frame = 1:length(scorrected)
        
        
        x = scorrected(frame,1);
        y = scorrected(frame,2);

       if participantData{c}(frame,2)==1 && participantData{c}(frame,13)>=0
                coords{c}(frame,:) = [x y];

       elseif  participantData{c}(frame,2)==1 && participantData{c}(frame,13)==-1
          coords{c}(frame,:) = [-2 -2]; % saccade
       elseif  participantData{c}(frame,2)==1 && participantData{c}(frame,13)==-3
          coords{c}(frame,:) = [-3 -3]; % blink
       elseif participantData{c}(frame,13)==-2
           coords{c}(frame,:) = [0 0]; % missing
       end % modified 15/11
           
            
        
    end
    
end

end