function registration = getRegistrations(targets,coords)

registration = cell(length(coords),1);
for i = 1:length(coords)
    registration{i} = zeros(length(targets)-2,1);

    for j = 1:length(targets)-2

        if coords{i}(j,1)>0
            registration{i}(j) = targets{j}(coords{i}(j,2),coords{i}(j,1));
            if registration{i}(j) == -1
                registration{i}(j) = findNearestNonBackground(coords{i}(j,1),coords{i}(j,2),targets{j},25);
            end
        elseif coords{i}(j,1) == 0; registration{i}(j) = 0;
        elseif coords{i}(j,1) < 0; registration{i}(j) = coords{i}(j,1);
        else
            display("hello")
        end
    end
end

end