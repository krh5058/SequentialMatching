function data_record(fid, entry_type, record_type)

% 2-5-12
% Ken Hwang
% Scherf Lab, SLEIC, Dept. of Psych., PSU

% Function to perform data recording
% Input:
%   - fid = diary file index
%   - entry_type = Either any string to dictate header creation or data
%   cell formatted in either short or long
%   - record_type = 'summary'
%
% Ex: data_record(fid, 'h', 'summary')
% Ex: data_record(fid, data_tmp, 'summary')

% Function handle for fprintf
f = @(x,y)fprintf(fid,x,y);

% Check whether header print or data entry
if strcmp(entry_type, 'h')
    
    if strcmp(record_type,'summary')
        h_str = {'GlobalMeanRT','CleanMeanRT-lower','CleanMeanRT-upper','GoodMeanRT','Accuracy','D-prime';'%s,','%s,','%s,','%s,','%s,','%s\n'};
    end % End if
    
    % Fprintf headers
    cellfun(f,h_str(2,:),h_str(1,:));
    
elseif iscell(entry_type)
    
    if strcmp(record_type,'summary') 
        resp_desc = {'%f,','%f,','%f,','%f,','%f,','%f\n'}; % Cell for response descriptors (comma-delimited)
        % 1 = Global Mean RT
        % 2 = Clean Mean RT (Lower bound)
        % 3 = Clean Mean RT (Upper bound)
        % 4 = Good (Not Outliers and correct) Mean RT
        % 5 = Accuracy
        % 6 = D-Prime
    end % END CHECK
    
    % Fprintf data
    cellfun(f,resp_desc,entry_type);
    
end % End classification check    
    