function out = nlfit(varargin)
    [args,options] = parseArguments(varargin);
    
    f = args{1};
    x = args{2};
    y = args{3};
    w = args{4};
    guess = args{5};
    
    defaultOptions = struct('fit_alpha',0.05,...
                            'prediction_alpha',0.01:0.01:0.99);
                        
    options = mergeStructs(lowerStruct(options),lowerStruct(defaultOptions));
    
    prediction_alpha = options.prediction_alpha;
    
    out.f = f;
    out.model = f;
    out.x = x;
    out.y = y;
    out.w = w;
    out.guess = guess;
    out.options = options;
    
    try
        x = x(:);
        y = y(:);
        w = w(:);
        
        hasWeights = ~isempty(w);

        if (hasWeights && ~isequal(size(w),size(y))) 
            error('DCELAB::nlfit : y and u must have same dimensions');
        end

        try
            displayPlot = strcmpi(options.DisplayPlot,'on');
        catch e
            displayPlot = false;
        end

%         if (hasWeights)
%             % weight fitting function and measured data
%             fw = @(par,x)(f(par,x)./u);
%             yw = y./u;
%         else
%             fw = f;
%             yw = y;
%         end
% 
%         [out.beta,out.res,out.J,out.cov,out.mse] = nlinfit(x,yw,fw,guess);

        if (hasWeights)
            % weight fitting function and measured data
            [out.beta,out.res,out.J,out.cov,out.mse] = nlinfit(x,y,f,guess,'Weights',w);

            % confidence interval of the parameters
            out.CI = nlparci(out.beta,out.res,'Jacobian',out.J);

            % confidence interval of the estimation
            [out.fit,delta] = nlpredci(f,x,out.beta,out.res,'Covar',out.cov,'Weights',w);

            % prediction interval of the estimation
            [~,deltaPred] = nlpredci(f,x,out.beta,out.res,'Covar',out.cov,'Weights',w,'PredOpt','observation');
            
            % generate prediction intervals for 0<alpha<1 for use in ROC analysis 
            idx=1;
            
            for idx=1:length(prediction_alpha) 
                alpha = prediction_alpha(idx);
                [~,d(:,idx)] = nlpredci(f,x,out.beta,out.res,'Covar',out.cov,'Weights',w,'Alpha',alpha,'PredOpt','observation'); 
                idx=idx+1; 
            end
        else
            [out.beta,out.res,out.J,out.cov,out.mse] = nlinfit(x,y,f,guess);

            % confidence interval of the parameters
            out.CI = nlparci(out.beta,out.res,'Jacobian',out.J);

            % confidence interval of the estimation
            [out.fit,delta] = nlpredci(f,x,out.beta,out.res,'Covar',out.cov);
            
            % prediction interval of the estimation
            [~,deltaPred] = nlpredci(f,x,out.beta,out.res,'Covar',out.cov,'PredOpt','observation');
            
            % generate prediction intervals for 0<alpha<1 for use in ROC analysis 
            idx=1;
            
            for idx=1:length(prediction_alpha) 
                alpha = prediction_alpha(idx);
                [~,d(:,idx)] = nlpredci(f,x,out.beta,out.res,'Covar',out.cov,'Alpha',alpha,'PredOpt','observation'); 
                idx=idx+1; 
            end
        end
        
        out.p = out.beta';

%         if (hasWeights)
%             out.res = out.res.*u(:);
%             out.J = out.J.*repmat(u(:),[1 length(out.beta)]);
%         end

        out.fitlo = out.fit - delta;
        out.fithi = out.fit + delta;
        
        out.predlo = out.fit - deltaPred;
        out.predhi = out.fit + deltaPred;
        
        out.predloalpha = repmat(out.fit,[1 size(d,2)]) - d;
        out.predhialpha = repmat(out.fit,[1 size(d,2)]) + d;
    catch e
        out = struct('beta',[],'res',[],'J',[],'cov',[],'mse',[],'p',[],'CI',[],'fit',[],'fitlo',[],'fithi',[]);
    end
end
