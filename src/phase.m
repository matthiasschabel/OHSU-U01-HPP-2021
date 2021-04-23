% PHASE circular color map
%   PHASE(M,C) returns an M-by-3 matrix containing a "phase" colormap built using colormap C.
%   PHASE(M) returns an M-by-3 matrix containing a "phase" colormap.
%   PHASE, by itself, is the same length as the current colormap.
function cmap = phase(m,c)
    if (nargin < 1)
       m = size(get(gcf,'colormap'),1); 
    end

    if (nargin < 2)
        c = jet(m);
    end

%     cmap = (c(1:2:m,:)+c(2:2:m,:))/2;
    cmap = [c; flipud(c)];
end
