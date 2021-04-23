function rng = atob(a,b,N)

if (N < 1) error('atob :: N must be at least 1'); end;

if (a > b) 
    invert = true;
    a = -a;
    b = -b;
else
    invert = false;
end;

if (N == 1)
    rng = [(a+b)/2];
else
    rng = a:(b-a)/(N-1):b;
end;

if (invert)
    rng = -rng;
end;

return; 
