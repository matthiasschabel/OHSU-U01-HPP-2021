% merge two structs, with fields in primary overriding fields found in secondary 
function out = mergeStructs(primary,secondary)
    if (isempty(secondary))
        secondary = struct();
    end
    
    if (~isempty(primary))
        f = fields(primary);

        for i=1:length(f)
            secondary.(f{i}) = primary.(f{i});
        end
    end
    
    out = secondary;
end
