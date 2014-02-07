function sm

% 2-10-12
% Ken Hwang
% Scherf Lab, SLEIC, Dept. of Psych., PSU

% Sequential matching task
% Presentation of three lists in random order (Faces, Greebles, Vehicles)
% Displays first stim, then mask, then second stim
% Records match/no-match response
% Data output of trial data and sumamry statistics for task subsets
%
% Usage: sm

%% 1: Start-up

monitor = monitor_set;

% Directory of this script
file_str = mfilename('fullpath');
[file_dir,~,~] = fileparts(file_str);

% Flags
end_flag = 0;

% Timepoint cell
tcell = {'baseline','post'};

% UI to select task
[sel,v] = listdlg('PromptString', 'Select a task:',...
'SelectionMode', 'single',...
'ListString',tcell);

if v
   
    tstr = tcell{sel}; % Specify task string
    
else % v == 0

    error('User Cancelled')
    
end % End if: v

% List load
load([file_dir filesep 'SM_practice.mat']);
load([file_dir filesep 'SM_' tstr '.mat']);

% Picture load
s = who([tstr '*']);
n_list = length(s);
s_prac = who('practice*'); % Assuming same number of practice conditions in SM_practice.mat as actual conditions

[x,y] = size(eval(s{1})); % Assuming all list sizes are the same
[x_prac,y_prac] = size(eval(s_prac{1})); % Assuming all practice list sizes are the same

list_cell = cell([x (y+2) n_list]); % Initializing list_cell
prac_cell = cell([x_prac (y_prac+2) n_list]); % Initializing list_cell
test_order = zeros([x, n_list]); % Initializing pres_order
prac_order = zeros([x_prac, n_list]); % Initializing prac_order
liststr = cell([1 3]); % Preallocating liststr
pracstr = cell([1 3]); % Preallocating pracstr

for i = 1:n_list % For each list
    
    list = eval(s{i}); % Evaluate list
    prac_list = eval(s_prac{i}); % Evaluate practice list
    listmatch = regexp(s{i},[tstr '_'],'split');
    pracmatch = regexp(s_prac{i},'_','split');
    liststr{i} = listmatch{2};
    pracstr{i} = pracmatch{2};
    
    % List check
    if ~strcmp(pracstr{i},liststr{i})
        disp(pracstr{i})
        disp(liststr{i})
        error('Conditions in practice list and testing list are not exact (case-sensitive).')
    end
        
    for j = 1:x % For each presentation
        
        list_cell{j,1,i} = imread([file_dir filesep 'Pictures' filesep liststr{i} filesep list{j,1}]); % Read in first picture file
        list_cell{j,2,i} = imread([file_dir filesep 'Pictures' filesep liststr{i} filesep list{j,2}]); % Read in second picture file
        list_cell{j,3,i} = list{j,1}; % Noting first stim name
        list_cell{j,4,i} = list{j,2}; % Noting second stim name
        list_cell{j,5,i} = list{j,3}; % Noting correct response
        list_cell{j,6,i} = list{j,4}; % Noting sub-description
    
    end % End for: j = 1:x
    
    for jj = 1:x_prac % For each practice presentation
        
        prac_cell{jj,1,i} = imread([file_dir filesep 'Pictures' filesep pracstr{i} filesep prac_list{jj,1}]); % Read in first picture file
        prac_cell{jj,2,i} = imread([file_dir filesep 'Pictures' filesep pracstr{i} filesep prac_list{jj,2}]); % Read in second picture file
        prac_cell{jj,3,i} = prac_list{jj,1}; % Noting first stim name
        prac_cell{jj,4,i} = prac_list{jj,2}; % Noting second stim name
        prac_cell{jj,5,i} = prac_list{jj,3}; % Noting correct response
        prac_cell{jj,6,i} = prac_list{jj,4}; % Noting sub-description
    
    end % End for: jj = 1:x_prac
    
    test_order(:,i) = Shuffle(1:x); % Order trials
    prac_order(:,i) = Shuffle(1:x_prac); % Order practice trials
    
end % End for: i = 1:n_list

% Randomized Trials
prompt1 = {'Faces:', 'Greebles:', 'Vehicles:'}; % Assuming these three conditions
dlg_title1 = 'Specify order (1,2,3)';    
num_lines1 = 1;
def1 = {'1', '2','3'};
options.Resize = 'on';
orderset = inputdlg(prompt1, dlg_title1, num_lines1, def1, options);
list_order = str2double(orderset); % Convert to str

% Loading mask
mask = imread([file_dir filesep 'Pictures' filesep 'mask.jpg']);

% Keyboard setup
KbName('UnifyKeyNames')
SameKey= KbName('q');
DifferentKey= KbName('p');
escapeKey = KbName('Escape');

%% End 1

%% 2: Trial run
% ---------- Window Setup ---------

% Hides the mouse cursor
HideCursor;
ListenChar(2)

% ***** This removes SyncTest Screen.  Double check this. *****
Screen('Preference', 'SkipSyncTests', 2); % Skips time sync tests

% Opens a graphics window on the main monitor (screen 0).  If you have
% multiple monitors connected to your computer, then you can specify
% a different monitor by supplying a different number in the second
% argument to OpenWindow, e.g. Screen('OpenWindow', 2).

[monitor.display_window, monitor.rect] = Screen('OpenWindow', monitor.whichScreen, monitor.white, [] , 32, 2); % Open Screen
Screen('FillRect',monitor.display_window,monitor.white);
Screen('Flip',monitor.display_window);

% Text formatting
Screen('TextSize',monitor.display_window,20);
Screen('TextFont',monitor.display_window,'Helvetica');
Screen('TextStyle',monitor.display_window,0);  

% User types ID
ID = Ask(monitor.display_window,'SubjectID:    ', [], monitor.white, 'GetChar', 'left', 'center');

% Making appropriate folder
mkdir([file_dir filesep 'data'], ID);

% For each list
for i2 = 1:n_list
    
    task_ind = find(list_order==i2); % Finding task index
    
    end_flag = 0; % Reset
    
    % Beginning instructions screen
    RestrictKeysForKbCheck(44); % Restrict to space

    % Starting screen
    DrawFormattedText(monitor.display_window,'Now you will go through a short block of practice.\n\n Press the Green button if the 2 pictures are the same.\n\n Press the Red button if they are different.\n\n Press the space bar when your ready.\n\n Good luck!','center', 'center', monitor.black);
    Screen('Flip',monitor.display_window);

    KbStrokeWait;
    
    for ii2 = 1:2 % For practice and test
    
        if ii2 == 1
           
            x_trial = x_prac;
            pres_order = prac_order;
            pres_cell = prac_cell;
            f_str = 'practice';
            textstr = 'You finished the practice!\n\n Remember to press the Green button if the pictures are the same and the Red button if they are different.\n\n Do you have any questions?\n\n Good luck!';
            
        else
            
            x_trial = x;
            pres_order = test_order;
            pres_cell = list_cell;
            f_str = tstr;
            textstr = 'Good job! You did great!';            
            
        end
        
        RestrictKeysForKbCheck([19 20 41]); % Restrict to q, p, and escape

        % Output prep
        filename = [file_dir filesep 'data' filesep ID filesep datestr(clock,30) '_' f_str '_' liststr{task_ind}];
        diaryfile = [filename '_diary']; % Diary file
        diary(diaryfile);
        fid = fopen([filename '.csv'],'a'); % Make csv

        fprintf(fid,'%s,','Trial');
        fprintf(fid,'%s,','FirstStim');
        fprintf(fid,'%s,','SecondStim');
        fprintf(fid,'%s,','Category');
        fprintf(fid,'%s,','Description');
        fprintf(fid,'%s,','ActualResponse');
        fprintf(fid,'%s,','SubjectResponse');
        fprintf(fid,'%s,','RT');
        fprintf(fid,'%s,','Accuracy');
        fprintf(fid,'%s,','Hit');
        fprintf(fid,'%s\n','FA');

        % For each trial
        for j2 = 1:x_trial

            % Displaying first stim
            pic1 = Screen('MakeTexture',monitor.display_window,pres_cell{pres_order(j2,task_ind),1,task_ind});
            Screen('DrawTexture',monitor.display_window,pic1);
            Screen('Flip',monitor.display_window);

            WaitSecs(.5); % Wait .5 seconds

            % Displaying mask
            masktex = Screen('MakeTexture',monitor.display_window,mask);
            Screen('DrawTexture',monitor.display_window,masktex);
            Screen('Flip',monitor.display_window);

            WaitSecs(.3); % Wait 300 milliseconds

            % Displaying second stim
            pic2 = Screen('MakeTexture',monitor.display_window,pres_cell{pres_order(j2,task_ind),2,task_ind});
            Screen('DrawTexture',monitor.display_window,pic2);
            Screen('Flip',monitor.display_window);

            % Clear keyCode
            keyCode(1:256) = 0;

            % Record start time for this trial
            RT = 0; % Resetting RT
            trial_start= GetSecs;

            keyIsDown = 0; % Resetting keyIsDown

            while keyIsDown == 0 && ((GetSecs - trial_start) <= 8) % Wait for keyboard press or 8 second wait

                [keyIsDown, secs, keyCode, ~] = KbCheck; % Check keyboard

                % Clearing if multiple keys are pressed
                if length(find(keyCode)) > 1
                    keyIsDown = 0;
                    secs = [];
                    keyCode(1:256) = 0; % Clear keyCode
                end

                RT= (secs-trial_start)*1000;

            end % End while

            % Gray Screen
            Screen('FillRect',monitor.display_window,monitor.white);
            Screen('Flip',monitor.display_window);

            WaitSecs(.5); % Wait 500 milliseconds

            % Codes what button was pressed
            if keyCode(SameKey)
                response='q';
            elseif keyCode(DifferentKey)
                response='p';
            elseif keyCode(escapeKey)
                response = '0';
                end_flag = 1;
                break;
            end % End if

            % Accuracy calculation
            if strcmp(response,pres_cell{pres_order(j2,task_ind),5,task_ind})

                accuracy=1; % Accuracy is 1

                % Hit/FA calculation
                if strcmp(response, 'q')

                    hit = 1;
                    FA = 0;

                else % Not 'q'

                    hit = 0;
                    FA = 0;

                end % End if: strcmp(response, 'q')

            else % if incorrect

                accuracy=0; % Accuracy is 0

                % Hit/FA calculation
                if strcmp(response, 'q')

                    hit = 0;
                    FA = 1;

                else % Not 'q'

                    hit = 0;
                    FA = 0;

                end % End if: strcmp(response, 'q')

            end % End if: pres_cell{pres_order(task_ind,j2),1,task_ind}

            % Data print
            fprintf(fid,'%3i,', j2);
            fprintf(fid,'%s,', pres_cell{pres_order(j2,task_ind),3,task_ind});
            fprintf(fid,'%s,', pres_cell{pres_order(j2,task_ind),4,task_ind});
            fprintf(fid,'%s,', liststr{task_ind});
            fprintf(fid,'%s,', pres_cell{pres_order(j2,task_ind),6,task_ind});
            fprintf(fid,'%s,', pres_cell{pres_order(j2,task_ind),5,task_ind});
            fprintf(fid,'%s,', response);
            fprintf(fid,'%3.2f,', RT);
            fprintf(fid,'%i,', accuracy);
            fprintf(fid,'%i,', hit);
            fprintf(fid,'%i\n', FA);

            % Texture close
            Screen('Close',pic1);
            Screen('Close',masktex);
            Screen('Close',pic2);

        end % End trial for

        % Data formatting
        fid = fopen([filename '.csv']); % File identifier for raw data
        datacell = textscan(fid,'%f%s%s%s%s%s%s%f%f%f%f','HeaderLines',1','Delimiter',','); % Read data
        fid2 = fopen([filename '_summary.csv'],'a'); % Starting summary file

        group_name = unique(datacell{5}); % Group identifiers
        group_name{end+1} = '\w*\w'; % Wildcard pattern expression for overall stats
        ind1 = find(strcmp('Identical', group_name)); % Finding Identical group (will be first index)
        d_order = [ind1 find(1:length(group_name)~=ind1)]; % Rearranging so 'Identical' is performed first

        for i3 = 1:length(d_order) % For each required stats calculation

            group_ind = find(cellfun(@(yyy)(~isempty(regexp(yyy,['(' group_name{d_order(i3)} ')'],'match'))), datacell{5})); % Respective indices

            RTvect = datacell{8}(group_ind); % RT vector (for group)
            AccInd = datacell{9}(group_ind)==1; % Accuracy Index (for group)
            ActualResponse = datacell{6}(group_ind); % Actual responses (for group)
            GroupTrials = datacell{1}(group_ind); % Trials (for group)
            HitVector = datacell{10}(group_ind); % Hits (for group)
            FAVector = datacell{11}(group_ind); % FAs (for group)

            inrange = ones([length(RTvect) 1]); % Preallocating inrange vector (0 if outlier, 1 if not)
            goodRT = zeros([length(RTvect) 1]); % Preallocating good RT vector (correct & within range)

            globalRT = mean(RTvect(AccInd)); % Mean reaction time for correct trials
            stdRT = std(RTvect(AccInd));  % Standard deviation of reaction time for correct trials
            cleanRT(1) = globalRT - 2*stdRT; % 2 Standard deviations under mean
            cleanRT(2) = globalRT + 2*stdRT; % 2 Standard deviations over mean

            % For each reaction time *** Including incorrect ***
            for j = 1:length(RTvect)

                if RTvect(j) < cleanRT(1) || RTvect(j) > cleanRT(2) % If outside of range

                    inrange(j) = 0; % Mark as outlier (0)

                else % If within range

                    if AccInd(j) == 1 % If also a correct trial

                        goodRT(j) = 1; % Mark as good trial

                    end % End if

                end % End if

            end % End for

            signal = length(find(cellfun(@(yy)(strcmp(yy,'q')),ActualResponse))); % Number of signal in task
            noise = size(GroupTrials,1) - signal; % Rest of trials are noise ('p')

            goodmeanRT = mean(RTvect(goodRT==1)); % Mean of good trials    
            acc_rate = length(find(AccInd))/length(AccInd); % Number of 1's in accuracy vector over total

            if i3 == 1 || i3 == length(d_order) % Hit rate calculation for 'Identical'(assuming first) or Overall only (assuming last)
                Hit_Rate = (length(find(HitVector)) + .5)/(signal + 1); % Number of 1's in Hit vector over total signal
            end % End if: i3 == 1 || i3 == length(d_order)

            FA_Rate = (length(find(FAVector))+.5)/(noise + 1); % Number of 1's in FA vector over total noise
            d_prime = norminv(Hit_Rate,0,1) - norminv(FA_Rate,0,1); % Calculating d_prime

            % Group Hit_Rate is taken from the 'Identical' iteration, except in the
            % Overall iteration, where Hit_Rate is calculated from the entire
            % data set

            d_cell = {globalRT, cleanRT(1), cleanRT(2), goodmeanRT, acc_rate, d_prime}; % Output cell

            % Writing onto data summary
            if i3 == length(d_order)
                fprintf(fid2,'%s\n', 'Overall'); % "Overall" header
            else % Individual groups
                fprintf(fid2,'%s\n', group_name{d_order(i3)}); % Header
            end % End if: i3 == length(d_order)

            data_record(fid2,'h','summary') % Creating headers
            data_record(fid2,d_cell,'summary') % Writing data

        end % End for: i3 = 1:(length(group_name)+1)

        % If escape was hit
        if end_flag
            break;
        end % End if: end_flag

        RestrictKeysForKbCheck(44); % Restrict to space

        % Ending screen
        Screen('FillRect',monitor.display_window,monitor.white);
        DrawFormattedText(monitor.display_window,textstr,'center', (monitor.center_H - 75), monitor.black);
        DrawFormattedText(monitor.display_window,'Press the space bar to continue.','center', (monitor.center_H + 75), monitor.black);
        Screen('Flip',monitor.display_window);

        KbStrokeWait; % Wait for keyboard stroke

    end % End for: ii2 = 1:2 
    
    % Ending break screen
    RestrictKeysForKbCheck(44); % Restrict to space
    
    Screen('FillRect',monitor.display_window,monitor.white);
    DrawFormattedText(monitor.display_window,'Please take a break.','center', 'center', monitor.black);
    Screen('Flip',monitor.display_window);

    KbStrokeWait; % Wait for keyboard stroke    
    
end % End list for

diary off;
fclose('all');
Screen('CloseAll');
ShowCursor;
ListenChar(0)
%% End 2

end % End primary function