function cmap = rgb(m)
    persistent rgb;

    if (nargin < 1)
        m = size(get(gcf,'colormap'),1);
    end

    if (isempty(rgb))
        rgb = ...
        [255   0   0;
           0 255   0;
           0   0 255];
    end;

    cmap = zeros([m 3]);

    % red 
    cmap(:,1) = interp1(0:1/2:1,rgb(:,1),0:1/(m-1):1)/255;

    % green
    cmap(:,2) = interp1(0:1/2:1,rgb(:,2),0:1/(m-1):1)/255;

    % blue
    cmap(:,3) = interp1(0:1/2:1,rgb(:,3),0:1/(m-1):1)/255;
end  
