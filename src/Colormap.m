% Colormap class
%
%   Colormap.xxx(m) returns a colormap of type 'xxx' of size m. If output
%       argument is not used, the current colormap is also set.
% 
%   Copyright 2011 Matthias Christian Schabel (matthias @ stanfordalumni . org)
%   Oregon Health & Science University
%   Advanced Imaging Research Center
%   3181 SW Sam Jackson Park Road L452
%   Portland, OR 97239
%   
classdef Colormap
    properties (Constant)
        defaultColormap = 'nih';
        defaultSize = 256;
    end
    
    methods (Static)
        function names = supportedColormaps
            names = {'delta','gecolor','nih','phase','placenta','stern','ucla','red','green','blue','autumn',...
                     'bipolar','bone','colorcube','cool','copper','flag','gray','hot','hsv','jet',...
                     'lines','pink','prism','rgb','spring','summer','thermal','white','winter'};
        end
        
        function cmap = generate(startRGB,endRGB,sz,method)
            if (nargin < 3)
                sz = 128;
            end
            
            if (nargin < 4 || strcmpi(method,'linear'))
                cmap = [interp1([1 sz],[startRGB(1) endRGB(1)],1:sz); ...
                        interp1([1 sz],[startRGB(2) endRGB(2)],1:sz); ...
                        interp1([1 sz],[startRGB(3) endRGB(3)],1:sz)]';
            elseif (ischar(method))
                cmap = [interp1([1 sz],[startRGB(1) endRGB(1)],1:sz,method); ...
                        interp1([1 sz],[startRGB(2) endRGB(2)],1:sz,method); ...
                        interp1([1 sz],[startRGB(3) endRGB(3)],1:sz,method)]';
            elseif (isa(method,'function_handle'))
                cmap = method(startRGB,endRGB,sz);
            end
            
            return;
        end
        
        function cmap = getCurrentColormap
            cmap = [];
            
            currentFigure = get(0,'CurrentFigure');
            
            if (~isempty(currentFigure))
                try
                    cmap = get(currentFigure,'Colormap');
                catch e
                    cmap = [];
                end
            end
                
            if (isempty(cmap))
                cmap = Colormap.get(Colormap.defaultColormap,Colormap.defaultSize);
            end
        end
        
        function setColormap(cmap,nargout)
            if (nargin < 2) nargout = 0; end
            
            currentFigure = get(0,'CurrentFigure');
            
            if (nargout == 0 && ~isempty(currentFigure))
                set(currentFigure,'Colormap',cmap);
            end
        end
        
        % get named colormap
        function cmap = get(name,m)   
            if (nargin < 2 || isempty(m)) 
                m = size(get(gcf,'colormap'),1); 
            end
            
            cmap = eval([name '(' num2str(m) ');']);
            
            Colormap.setColormap(cmap,nargout);
        end
        
        % DELTA color map with black for small values, bright red for large +, bright blue for large -.
        %   DELTA(M) returns an M-by-3 matrix containing a "delta" colormap.
        %   DELTA, by itself, is the same length as the current colormap.
        function cmap = delta(m)
            if (nargin < 1)
                m = size(get(gcf,'colormap'),1);
            end

            cmap = delta(m);
             
            Colormap.setColormap(cmap,nargout);
        end
        
        % GECOLOR color map with "gecolor" colormap from OsiriX (http://www.osirix-viewer.com)
        function cmap = gecolor(m)
            if (nargin < 1)
                m = size(Colormap.getCurrentColormap,1);
            end

            cmap = gecolor(m);
             
            Colormap.setColormap(cmap,nargout);
        end     
        
        % NIH color map with "nih" colormap from OsiriX (http://www.osirix-viewer.com)
        %   NIH(M) returns an M-by-3 matrix containing a "nih" colormap.
        %   NIH, by itself, is the same length as the current colormap.
        function cmap = nih(m)
            if (nargin < 1)
                m = size(Colormap.getCurrentColormap,1);
            end

            cmap = nih(m);
             
            Colormap.setColormap(cmap,nargout);
        end
        
        % PHASE circular color map
        %   PHASE(M,C) returns an M-by-3 matrix containing a "phase" colormap built using colormap C.
        %   PHASE(M) returns an M-by-3 matrix containing a "phase" colormap.
        %   PHASE, by itself, is the same length as the current colormap.
        function cmap = phase(m,c)
            if (nargin < 2)
                c = Colormap.getCurrentColormap;
            end
            
            if (nargin < 1)
                m = size(c,1);
            end

            cmap = phase(m,c);
            
            Colormap.setColormap(cmap,nargout);
        end
        
        function cmap = placenta(m)
            if (nargin < 1)
                m = size(Colormap.getCurrentColormap,1);
            end
            
            cmap = [Colormap.ucla(ceil(m/2)); Colormap.generate([1 0 0],[1 1 1],m-ceil(m/2))];
        end
        
        % STERN color map
        function cmap = stern(m)
            if (nargin < 1)
                m = size(Colormap.getCurrentColormap,1);
            end

            cmap = stern(m);
            
            Colormap.setColormap(cmap,nargout);
        end
        
        % UCLA color map with "ucla" colormap from OsiriX (http://www.osirix-viewer.com)
        function cmap = ucla(m)
            if (nargin < 1)
                m = size(Colormap.getCurrentColormap,1);
            end

            cmap = ucla(m);
            
            Colormap.setColormap(cmap,nargout);
        end  
        
        % black to red
        function cmap = red(m)
            if (nargin < 1)
                m = size(Colormap.getCurrentColormap,1);
            end
            
            cmap = zeros([m 3]);
            
            cmap(:,1) = linspace(0,1,m);
            
            Colormap.setColormap(cmap,nargout);
        end
        
        % black to green
        function cmap = green(m)
            if (nargin < 1)
                m = size(Colormap.getCurrentColormap,1);
            end
            
            cmap = zeros([m 3]);
            
            cmap(:,2) = linspace(0,1,m);
            
            Colormap.setColormap(cmap,nargout);
        end
        
        % black to blue
        function cmap = blue(m)
            if (nargin < 1)
                m = size(Colormap.getCurrentColormap,1);
            end
            
            cmap = zeros([m 3]);
            
            cmap(:,3) = linspace(0,1,m);
            
            Colormap.setColormap(cmap,nargout);
        end
        
        % MATLAB built-in colormaps
        function cmap = autumn(m)
            if (nargin < 1)
                m = size(Colormap.getCurrentColormap,1);
            end

            cmap = autumn(m);
            
            Colormap.setColormap(cmap,nargout);
        end

        function cmap = bipolar(m)
            if (nargin < 1)
                m = size(Colormap.getCurrentColormap,1);
            end

            cmap = bipolar(m);
             
            Colormap.setColormap(cmap,nargout);
        end
       
        function cmap = bone(m)
            if (nargin < 1)
                m = size(Colormap.getCurrentColormap,1);
            end

            cmap = bone(m);
            
            Colormap.setColormap(cmap,nargout);
        end
        
        function cmap = colorcube(m)
            if (nargin < 1)
                m = size(Colormap.getCurrentColormap,1);
            end

            cmap = colorcube(m);
            
            Colormap.setColormap(cmap,nargout);
        end
         
        function cmap = cool(m)
            if (nargin < 1)
                m = size(Colormap.getCurrentColormap,1);
            end

            cmap = cool(m);
            
            Colormap.setColormap(cmap,nargout);
        end
       
        function cmap = copper(m)
            if (nargin < 1)
                m = size(Colormap.getCurrentColormap,1);
            end

            cmap = copper(m);
            
            Colormap.setColormap(cmap,nargout);
        end
        
        function cmap = flag(m)
            if (nargin < 1)
                m = size(Colormap.getCurrentColormap,1);
            end

            cmap = flag(m);
            
            Colormap.setColormap(cmap,nargout);
        end
        
        function cmap = gray(m)
            if (nargin < 1)
                m = size(Colormap.getCurrentColormap,1);
            end

            cmap = gray(m);
            
            Colormap.setColormap(cmap,nargout);
        end
        
        function cmap = hot(m)
            if (nargin < 1)
                m = size(Colormap.getCurrentColormap,1);
            end

            cmap = hot(m);
            
            Colormap.setColormap(cmap,nargout);
        end
         
        function cmap = hsv(m)
            if (nargin < 1)
                m = size(Colormap.getCurrentColormap,1);
            end

            cmap = hsv(m);
            
            Colormap.setColormap(cmap,nargout);
        end
        
        function cmap = jet(m)
            if (nargin < 1)
                m = size(Colormap.getCurrentColormap,1);
            end

            cmap = jet(m);
            
            Colormap.setColormap(cmap,nargout);
        end
       
        function cmap = lines(m)
            if (nargin < 1)
                m = size(Colormap.getCurrentColormap,1);
            end

            cmap = lines(m);
            
            Colormap.setColormap(cmap,nargout);
        end
        
        function cmap = pink(m)
            if (nargin < 1)
                m = size(Colormap.getCurrentColormap,1);
            end

            cmap = pink(m);
            
            Colormap.setColormap(cmap,nargout);
        end
        
        function cmap = prism(m)
            if (nargin < 1)
                m = size(Colormap.getCurrentColormap,1);
            end

            cmap = prism(m);
             
            Colormap.setColormap(cmap,nargout);
        end

        function cmap = rgb(m)
            if (nargin < 1)
                m = size(Colormap.getCurrentColormap,1);
            end

            cmap = rgb(m);
             
            Colormap.setColormap(cmap,nargout);
        end  
        
        function cmap = spring(m)
            if (nargin < 1)
                m = size(Colormap.getCurrentColormap,1);
            end

            cmap = spring(m);
            
            Colormap.setColormap(cmap,nargout);
        end
        
        function cmap = summer(m)
            if (nargin < 1)
                m = size(Colormap.getCurrentColormap,1);
            end

            cmap = summer(m);
             
            Colormap.setColormap(cmap,nargout);
       end
        
        function cmap = thermal(m)
            if (nargin < 1)
                m = size(Colormap.getCurrentColormap,1);
            end

            cmap = thermal(m);
             
            Colormap.setColormap(cmap,nargout);
        end
        
        function cmap = white(m)
            if (nargin < 1)
                m = size(Colormap.getCurrentColormap,1);
            end

            cmap = white(m);
             
            Colormap.setColormap(cmap,nargout);
       end
        
        function cmap = winter(m)
            if (nargin < 1)
                m = size(Colormap.getCurrentColormap,1);
            end

            cmap = winter(m);
            
            Colormap.setColormap(cmap,nargout);
        end     
    end
end
