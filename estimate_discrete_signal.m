function [varargout] = estimate_discrete_signal (threshold, varargin)

nchannel = length(varargin);
varargout = cell(1,nchannel);

for i=1:nchannel
    varargout{i} = 0;
    ch = varargin{i};
    if length(ch) > 1
        [chngpts,~] = findchangepts(ch,'MinThreshold',threshold,'Statistic','mean');
        chngpts = chngpts -1;
        if isrow(chngpts)
            chngpts = [chngpts length(ch)];
        else
            chngpts = [chngpts; length(ch)];
        end
        numchngpts = length(chngpts);
        
        est=[];
        temp=1;
        for j= 1: numchngpts
            if isrow(chngpts)
                est = [est mean(ch(temp:1+chngpts(j)-1))*ones(1, 1+chngpts(j)-temp)];
            else
                est = [est; mean(ch(temp:1+chngpts(j)-1))*ones(1+chngpts(j)-temp, 1)];
            end
            temp = 1 + chngpts(j);
        end       
    end
    varargout{i} = est;
end