function out = lowerStruct(in)
    if (~isstruct(in))
        out = in;
        return;
    else
        f = fieldnames(in);
        
        % check for case-insensitive duplicates
        fci = cell(size(f));
        
        for i=1:length(f)
            fci{i} = lower(f{i});
        end
        
        if (length(unique(fci)) ~= length(fci))
            error('DCELAB::lowerStruct : case-insensitive fieldnames are not unique');
        end
        
        out = struct();
        
        for i=1:length(f)
            currentFieldname = f{i};
            
            out.(lower(currentFieldname)) = in.(currentFieldname);
        end
    end
end