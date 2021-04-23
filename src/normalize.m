function out = normalize(in,dim)
    if (nargin < 2)
        if (isrowvector(in))
            dim = 2;
        else
            dim = 1;
        end
    end
    
    n = sum(in,dim);
    
    p = [dim 1:dim-1 dim+1:ndims(in)];
    
    n = permute(n,p);
    
    n = repmat(n,[size(in,dim) 1]);
    
    n = ipermute(n,p);
    
    out = in./n;
end