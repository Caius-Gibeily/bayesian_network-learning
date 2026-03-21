function [indivarrayTD,indivarrayASD] = processViewSeq(coordsTD,coordsASD,targetCell)
    % recoding
    registrationTD = getRegistrations(targetCell,coordsTD);
    registrationASD = getRegistrations(targetCell,coordsASD);
    
    % Smoothing
    registrationTD2 = smoothTargetsWindow(registrationTD,0.6,10);
    registrationASD2 = smoothTargetsWindow(registrationASD,0.6,10);
    
    % Vectorize
    indivarrayTD = zeros(length(coordsTD),length(targetCell)-2);
    indivarrayASD = zeros(length(coordsASD),length(targetCell)-2);
    
    for i = 1:size(coordsTD,1)
        indivarrayTD(i,:) = registrationTD2{i};
    end
    for i = 1:size(coordsASD,1)
        indivarrayASD(i,:) = registrationASD2{i};
    end
    

end