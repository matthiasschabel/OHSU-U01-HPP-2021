function [plotHandles,out] = compareModelFits(xx,yy,masks,model,guess,ww,labels)
    if (nargin < 7) 
        labels = [];
    end
    
    if (nargin < 6)
        ww = [];
    end
    
    markerSize = 6;
    lineWidth = 3;
    
    if (~iscell(masks))
        masks = { masks };
    end
    
    if (length(masks)<4)
        cmap = [0 0 1;1 0 0;0 1 0];
    else
        cmap = Colormap.rgb(length(masks));
    end
    
    symbols = {'o','o','o','^','v','*'};
    symbols = {'.','.','.','^','v','*'};
    
    clf;
    
    for idx=1:length(masks)
        msk = masks{idx};
        x = xx(msk)';
        y = yy(msk)';
        
        [~,sidx] = sort(x);
        
        x = x(sidx);
        y = y(sidx);
        
        if (~isempty(labels))
            lb = labels(msk);
            lb = lb(sidx);
        end
        
        if (~isempty(ww))
            w = ww(msk)';
            w = w(sidx);
        else
            w = [];
        end
        
        fit = nlfit(model,x,y,w,guess);
        
        plotPredictionIntervals = false;
        plotPredictionIntervals = true;
        
        if (idx <= length(symbols))
            sym = symbols{idx};
        else
            sym = symbols{end};
        end
        
        hold on;
        
        if (plotPredictionIntervals)
            ph = plot(x,y,sym,x,fit.fit,'-',x,fit.fitlo,'--',x,fit.fithi,'--',x,fit.predlo,'-.',x,fit.predhi,'-.');
        else
            ph = plot(x,y,sym,x,fit.fit,'-',x,fit.fitlo,'--',x,fit.fithi,'--');
        end  
        
        if (~isempty(labels))
            deltax = 0.1;
            deltay = 1;
            
            for jdx=1:length(lb)
                text(x(jdx)+deltax,y(jdx)+deltay,lb{jdx},'FontSize',8,'Color',cmap(idx,:));
            end
        end
        
        hold off;
        
        ph(1).Color = cmap(idx,:);
        ph(1).MarkerSize = markerSize;
        
        if (plotPredictionIntervals)
            maxidx = 6;
        else
            maxidx = 4;
        end
        
        for jdx=2:maxidx
            ph(jdx).Color = cmap(idx,:);
             
            if (jdx < 2)
                ph(jdx).LineWidth = lineWidth;
            elseif (jdx >=2 && jdx < 5)
                ph(jdx).LineWidth = lineWidth-1;
            else
                ph(jdx).LineWidth = lineWidth-2;
            end            
        end
        
        plotHandles{idx} = ph;
        out(idx) = fit;
    end
end
