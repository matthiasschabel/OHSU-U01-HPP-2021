function [args,nvpargs] = parseArguments(argin)
    args = {};
    nvpargs = struct();

    i = 1;
    j = 1;

    while i<=length(argin)
        if (isa(argin{i},'char'))
            if (i==length(argin))
    %             error(['DCELAB::parseArguments : argument ' argin{i} ' missing value']);
                args{j} = argin{i};
                return;
%             elseif (isa(argin{i+1},'char'))
%                 args{j} = argin{i};
%                 i = i+1;
%                 continue;
            end

            nvpargs.(argin{i}) = argin{i+1};

            i = i+2;
        else
            args{j} = argin{i};
            i = i+1;
            j = j+1;
        end
    end
return

