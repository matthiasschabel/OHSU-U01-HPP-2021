function ret = nanquantile(y,q,dim)
% FORMAT: y = nanquantile(x,q,dim)
% 
%    Quantiles, ignoring NaNs
%
%    y = nanquantile(x,q) returns the quantiles or quantiles listed in q, 
%    treating NaNs as missing values.
%
%    For vectors, nanquantile(x,q) is a vector of q quantiles of the non-NaN
%    elements of x. For matrices, nanquantile(x,q) is a set of row
%    vectors containing the q quantiles of the non-NaN elements in each
%    column of x. For N-dimensional arrays, nanquantile operates along the
%    first dimension of x.
%
%    y = nanquantile(x,q,dim) returns the q quantiles along the dimension
%    dim of x

% -------------------------------------------------------------------------
%    author:      Matthias Schabel
%    affiliation: UCAIR, University of Utah, Salt Lake City, UT 84108-1218
%    email:       matthias.schabel@hsc.utah.edu
%    
%    Version 1.0, 2008/12/08

if (nargin < 3) 
    dim = 1; 
end;

if (dim > ndims(y))
    error('quantile :: dim exceeds dimension of y');
end;

if (~isvector(q))
    error('quantile :: quantiles must be specified as a vector');
end;

if (any(~isfinite(q)))
    error('quantile :: quantiles must be in [0,1]');
end;

% if (isrowvector(y)) y = y'; end;

x = permute(y,[dim 1:(dim-1) dim+1:ndims(y)]);
sz = size(x);
x = reshape(x,[sz(1) prod(sz(2:length(sz)))]);

% ret = zeros([1 sz(2) length(q)]);
ret = zeros([sz(2) length(q)]);

for i=1:size(x,2)
%     ret(1,i,:) = nanquantile_impl(x(:,i),q);
    ret(i,:) = nanquantile_impl(x(:,i),q);
end;

ret = reshape(ret,[sz(2:end) size(ret,2)]);
    
return;


% find quantiles for each column of reduced data set
function ret = nanquantile_impl(y,q)

y = flatten(y);

ret = ones(size(q))*NaN;

% find values that are not NaN/Inf
y = y(find(isfinite(y)));

% if no values that are not NaN/Inf, return NaNs for quantiles
if (numel(y) == 0) 
    return; 
end;

% sort values
if (~issorted(y))
    y=sort(y);
end;

% iterate through desired quantiles
for i=1:length(q)
    qv = q(i)*(length(y)-1);
    ql = floor(qv);
    qh = ceil(qv);

    if (ql == qh)
        ret(i) = y(ql+1);
    else     
        y1 = y(ql+1);
        y2 = y(qh+1);

        ret(i) = interp1([ql qh],[y1 y2],qv);
    end;
end;

return;
