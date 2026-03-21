function [PaddedNumberAsString]=common_getpaddednumber(UnpaddedNumberAsNumber,HowManyPlaces)
%
% 	HELP: 			pads a number to 2,3,4,5, or 6 places, returns as string

UnpaddedNumberAsString = num2str(UnpaddedNumberAsNumber);

LengthUnpaddedNumberAsString = length(UnpaddedNumberAsString);

if HowManyPlaces == 2
    LengthUnpaddedNumberAsString = LengthUnpaddedNumberAsString + 4;
elseif HowManyPlaces == 3
    LengthUnpaddedNumberAsString = LengthUnpaddedNumberAsString + 3;
elseif HowManyPlaces == 4
    LengthUnpaddedNumberAsString = LengthUnpaddedNumberAsString + 2;
elseif HowManyPlaces == 5
    LengthUnpaddedNumberAsString = LengthUnpaddedNumberAsString + 1;
elseif HowManyPlaces == 6
    LengthUnpaddedNumberAsString = LengthUnpaddedNumberAsString;
else
    disp('*******');
    disp('**ERROR: You cannot try to pad to that many places.');
    disp('*******');
    return

end

switch LengthUnpaddedNumberAsString
    case 1
        PaddedNumberAsString=['00000' UnpaddedNumberAsString];
    case 2
        PaddedNumberAsString=['0000' UnpaddedNumberAsString];
    case 3
        PaddedNumberAsString=['000' UnpaddedNumberAsString];
    case 4
        PaddedNumberAsString=['00' UnpaddedNumberAsString];
    case 5
        PaddedNumberAsString=['0' UnpaddedNumberAsString];
    case 6
        PaddedNumberAsString=['' UnpaddedNumberAsString];
end
