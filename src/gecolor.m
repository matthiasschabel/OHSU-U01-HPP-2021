% GECOLOR color map with "gecolor" colormap from OsiriX (http://www.osirix-viewer.com)
%   GECOLOR(M) returns an M-by-3 matrix containing a "gecolor" colormap.
%   GECOLOR, by itself, is the same length as the current colormap.
function cmap = gecolor(m)
    if (nargin < 1) 
        m = size(get(gcf,'colormap'),1); 
    end;

    cmap = zeros([m 3]);

    fst = 1;
    lo = m/4;
    mid = m/2;
    hi = 3*m/4;
    lst = m;

    % red 
    cmap(fst:lo,1) = 0.0;
    cmap(lo+1:hi,1) = (0:hi-lo-1)./(hi-lo-1);
    cmap(hi+1:lst,1) = 1.0;

    % green
    cmap(fst:lo,2) = 0.5*((0:lo-1)./(lo-1));
    cmap(lo+1:mid,2) = cmap(lo,2)-0.5*((0:mid-lo-1)./(mid-lo-1));
    cmap(mid+1:lst,2) = (0:lst-mid-1)./(lst-mid-1);

    % blue
    cmap(1:mid,3) = (0:mid-1)./(mid-1);
    cmap(mid+1:hi,3) = cmap(mid,3)-((0:hi-mid-1)./(hi-mid-1));
    cmap(hi+1:lst,3) = (0:lst-hi-1)./(lst-hi-1);
end     
