% DELTA color map with black for small values, bright red for large +, bright blue for large -.
%   DELTA(M) returns an M-by-3 matrix containing a "delta" colormap.
%   DELTA, by itself, is the same length as the current colormap.
function cmap = delta(m)
    if (nargin < 1)
        m = size(get(gcf,'colormap'),1);
    end

    mtwo = floor(m/2);
    cr = [(0:(1/(mtwo-1)):1)' zeros([mtwo 1]) zeros([mtwo 1])];
    cb = [zeros([mtwo 1]) zeros([mtwo 1]) (0:(1/(mtwo-1)):1)'];
    cmap = [flipud(cb); cr;];
end
