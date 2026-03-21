function sampleBoots = bootstrapParticipantArray(participantArray,nboots)

arguments
    participantArray
    nboots = 1000

end
sampleBoots = cell(nboots,1);
for i = 1:nboots
    sampleBoots{i} = datasample(participantArray,size(participantArray,1));
end