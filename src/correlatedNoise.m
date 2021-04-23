function out = correlatedNoise(N,covar)

if (size(covar,1) ~= size(covar,2))
    error('DCELAB::correlatedNoise : covariance matrix must be square');
end;

numberOfVariables = size(covar,1);

noise = randn([N numberOfVariables]);

cholfac = chol(covar);

out = noise*cholfac;

end