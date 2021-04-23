function val = isrowvector(v)

val = isvector(v) && size(v,1) == 1;

end