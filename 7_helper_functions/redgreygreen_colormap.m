function cmap = redgreygreen_colormap(n)
  

    anchors = [
        1 0 0;        
        0.5 0.5 0.5;  
        0 1 0        
    ];

    % Interpolate to n colors
    cmap = interp1([1 n/2 n], anchors, 1:n, 'linear');
end