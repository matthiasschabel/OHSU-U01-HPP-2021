% STERN color map
%   STERN(M) returns an M-by-3 matrix containing a "stern" colormap.
%   STERN, by itself, is the same length as the current colormap.
function cmap = stern(m)
    if (nargin < 1)
       m = size(get(gcf,'colormap'),1); 
    end;

    cmap = zeros([m 3]);

    fst = 1;
    pt1 = m/16;
    pt2 = m/4;
    pt3 = m/2;
    pt4 = 3*m/4;
    lst = m;

    % red 
    cmap(fst:pt1,1) = (0:pt1-fst)./(pt1-fst);
    cmap(pt1+1:pt2,1) = cmap(pt1,1)-(0:pt2-pt1-1)./(pt2-pt1-1);
    cmap(pt2+1:lst,1) = (0:lst-pt2-1)./(lst-pt2-1);

    % green
    cmap(fst:lst,2) = (0:lst-fst)./(lst-fst);

    % blue
    cmap(fst:pt3,3) = (0:pt3-fst)./(pt3-fst);
    cmap(pt3+1:pt4,3) = cmap(pt3,3)-((0:pt4-pt3-1)./(pt4-pt3-1));
    cmap(pt4+1:lst,3) = (0:lst-pt4-1)./(lst-pt4-1);
end
