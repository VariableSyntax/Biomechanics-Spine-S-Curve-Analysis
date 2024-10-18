function load_disp_edits()

% Written by Jalen Winfield (Research Engineer Co-op, Fall/Winter 2022)
% 03/14/2022

% Written by Josh McGuckin (Associate Research Engineer, Summer 2022)
% 07/31/2023

% Written by Izzy Lachcik (Associate Research Engineer, Winter 2024)
% 03/01/2024

% Edited by Josh McGuckin (Research Engineer, Summer 2024)
% 09/17/2024

%-----General Information-----%

% --- Version --- %
% v2.6


% --- About --- %
% The following script analyzes angular displacement from the Optotrak and
% the load applied during 6DOF testing and plots the hysteresis curve(s) of
% the last cycle. This script also calculates the magnitude of the neutral 
% zone the neutral zone stiffnesses, the magnitude of the elastic zones,
% and the elastic zone stiffnesses. Graphs and summary data tables are
% saved.

% This version of the code has been modified to allow for feature selection
% of NZ and EZ


% --- Abbreviations --- %
% AR = Axial Rotational
% EZ = Elastic Zone
% EZS = Elastic Zone Stiffness
% FE = Flexion-Extension
% LB = Lateral Bending
% NZ = Neutral Zone
% NZB = Neutral Zone Boundary
% NZS = Neutral Zone Stiffness

%% INITIALIZING ANALYSIS - ALLOWS USER TO CHOOSE WHAT FILES TO ANALYZE
% --- Collecting all 6DOF and Optotrak files --- %
clear all
clc
try
    z=1;
    while z==1
    folderdir = uigetdir(path, 'Select Folder with Load Data');
    % folder containing load data
    sixdofdir = dir(fullfile(folderdir,'*.xlsx*')); 

    if ~isempty(sixdofdir) % Extract only Load Data files (i.e.'Book')
        for p = 1:length(sixdofdir) % Extract all .xlsx files that have
            % 'Books' in the filename
            if ~isempty(strfind(sixdofdir(p).name, 'Book'))
                Books(p) = sixdofdir(p);
            end 
        end
        if exist('Books','var')
            sixdofdir = Books'; % Set new sixdordir to include only 
            % .xlsx files that have 'Book' in their filename
        end 
    end 

    while isempty(sixdofdir)
        w = questdlg(['Cannot continue because there is no Load Data files ' ...
                    'within the selected Load Data folder. ' ...
                    'Please select a folder with the appropriate .xlsx ' ...
                    'files containing the Load Data.'], ...
                    'Warning: No Load Data Found','OK','OK');
        folderdir = uigetdir(path, 'Select Folder with Load Data');
        sixdofdir = dir(fullfile(folderdir,'*.xlsx*'));

        if ~isempty(sixdofdir) % Extract only Load Data files (i.e.'Book')
            for p = 1:length(sixdofdir) % Extract all .xlsx files that have
                % 'Books' in the filename
                if ~isempty(strfind(sixdofdir(p).name, 'Book'))
                    Books(p) = sixdofdir(p);
                end 
            end
            if exist('Books','var')
                sixdofdir = Books'; % Set new sixdordir to include only 
                % .xlsx files that have 'Book' in their filename
            end 
        end 

    end    

    folderdir2 = uigetdir(path, 'Select Folder With Optotrak Data');
    optodir = dir(fullfile(folderdir2,'*_cal.xls*')); % folder containing 
    % angular displacement data (calculated relatives)

    while isempty(optodir)
        w = questdlg(['Cannot continue because there is no Optotrak Data' ...
            ' files within the selected Optotrak Data folder. Please ' ...
            'select a folder with the appropriate .xlsx files containing ' ...
            'the Optotrak Data.'], 'Warning: No Optotrak Data Found','OK', ...
            'OK');
        folderdir2 = uigetdir(path, 'Select Folder With Optotrak Data');
        optodir = dir(fullfile(folderdir2,'*_cal.xls*')); 
    end    

    loadfiles = struct2cell(sixdofdir);
    name_loadfiles = loadfiles(1,:);
    anglefiles = struct2cell(optodir);
    name_anglefiles = anglefiles(1,:);

%% Initializing variables for special cases
% In circumstances where there files used for analysis do not come in a 
% multiple of 3 or there is not a pair for every load/displacement file,
% the code will try to circumvent the issue(s). In normal cases the 
% variables skip, a, s1, and s2 are set to 0 but are then set to 1 
% depending on the case. 
    
    skip = 0;
    s1 = 0;
    s2 = 0;

% --- Creating folder for plots --- %
    newfolderdir = dir(fullfile(folderdir, ...
        'Angular Displacement VS Load Plots*'));
    newFolderName = 'Angular Displacement VS Load Plots';
    mkdir(folderdir, newFolderName);

% ---Case Scenario 1: Uneven number of data files --- %
% Sometimes data will be captured for the Optotrak during a 6DOF test
% while an excel sheet for the load data is not produced or vice versa. 
% In the case that there is missing pair for load/displacement files, 
% the user will first be prompted to reorganize the data.
% Next the user is prompted to which file is in excess and does not have
% a match/pair; selected files will be skipped in the analysis but the 
% expected order will be preserved (i.e., if a FE file is skipped, the
% next file analysis is for LB).

% --- Accounting for missmatched number of files between data sets --- %
    if length(sixdofdir)<length(optodir) % less load data files than 
        % Optotrak files
        fig = uifigure;
        selection = uiconfirm(fig,...
            {['The number of Optotrak Data files does not ' ...
            'match the number of Load Data files.'], ' '...
            'Would you like to omit certain Optotrak Data files?'},...
            'Reorganize Data','Icon','warning','Options',{'YES','NO'});
        close(fig)
        switch selection
            case 'YES'
                [skip,~] = listdlg('PromptString',...
                    {['Choose Optotrak Data files that do' ...
                    ' not have a matching Load Data file'],''},...
                    'ListString',name_anglefiles,'ListSize',[200 300]);
                %optodir(skip) = []; % Exlcude skipped Optotrak files
                if isempty(skip)
                    g = questdlg(['Cannot continue unless the Load Data ' ...
                    'files match with the corresponding Optotrak Data ' ...
                    'files. Make sure there are equal numbers of files in ' ...
                    'their respective folders and that these files ' ...
                    'correspond to the same bending test.'],['Warning:' ...
                    ' Cannot Continue'],'OK','OK');
                    close(fig)
                    break
                end
            case 'NO' % if no is chosen, script stops and cannot continue
                g = questdlg(['Cannot continue unless the Load Data ' ...
                    'files match with the corresponding Optotrak Data ' ...
                    'files. Make sure there are equal numbers of files in ' ...
                    'their respective folders and that these files ' ...
                    'correspond to the same bending test.'],['Warning:' ...
                    ' Cannot Continue'],'OK','OK');
                close(fig)
                break
        end
    elseif length(optodir)<length(sixdofdir) % if less Optotrak files than 
        % load data files
        fig2 = uifigure;
        selection2 = uiconfirm(fig2,...
            {['The number of Load Data files does not match the number of' ...
            ' Optotrak Data files.'], ' '...
            'Would you like to omit certain Load Data files?'},...
            'Reorganize Data','Icon','warning','Options',{'YES','NO'});
        close(fig2)
        switch selection2
            case 'YES'
                [skip,TF] = listdlg('PromptString',...
                    {'Choose Load Data files that do not have a matching Optotrak Data file.' ...
                ''},'ListString',name_loadfiles,'ListSize',[200 300]);
                %sixdofdir(skip) = []; % Exlcude skipped load files 
            case 'NO' % if no is chosen, script stops and cannot continue
                g = questdlg(['Cannot continue unless the Load Data ' ...
                    'files match with the corresponding Optotrak Data ' ...
                    'files. Make sure there are equal numbers of files in ' ...
                    'their respective folders and that these files ' ...
                    'correspond to the same bending test'],['Warning:' ...
                    ' Cannot Continue'],'OK','OK');
                close(fig2)
                break
        end
    end

% --- Case Scenario 2: Missing a triplet in data files --- %
% It is important to note that the code utilzes a grand index matrix 
% to index each file for analysis. Since this code is built primarily
% for use with Globus's 6DOF machine, there is always a round of testing
% that consists of FE, LB, and AR; and then it repeats for each cycle of
% testing (I'm writing this in 2022, so if anything changes, that's on 
% you). In the case that the number of files selected for analysis is 
% not a multiple of three, the script can still run since it will 'pad'
% the grand index with 0s. 

% ***Note: the code runs but it is not smart enough to know where the
% gap is on it's own and will always 'pad' the matrix at the end

% --- Accounting for missing triplet in load data set --- %
    if rem(length(sixdofdir),3)==1 % missing 2/3 load files 
        % (i.e., the number of load files is 2 less than the next highest 
        % multiple of 3)
        % Prompt user to select where the missing data is in the file list
        [missing,TF] = listdlg('PromptString',...
                    {'Choose the Load Data files that comes after ' ...
                'the missing Load Data files. If the first FE Load ' ...
                'Data file is missing, select the first file from the ' ...
                'top. If the last LB and AR Load Data file is missing, press ' ...
                'Cancel or close the window '},'ListString',name_loadfiles, ...
                'ListSize', [250 300]);
        grand_index = [1:3:length(sixdofdir);2:3:length(sixdofdir),0;...
            3:3:length(sixdofdir),0];
        % Fill in padded zeros
        grand_index(grand_index == 0) = find((grand_index == 0));
        % Nullify index where there is missing Load Data files (FE,LB,AR)
        if isempty(missing) % If the user selects the last Load Data file
            grand_index(end) = 0;   % Nullify the last LB/AR files
        else
            for d = 1:length(missing) 
                if missing(d) == 1 % If the user selects the first file
                    % Nullify the first FE/LB files
                    grand_index(missing:missing+1) = 0; 
                else % Otherwise nullify the files in front of the selected 
                    grand_index(missing-2:missing-1) = 0;
                end
            end 
        end
        a = 2;
    elseif rem(length(sixdofdir),3)==2 % missing 1/3 load file 
        % (i.e., the number of load files is 1 less than the next highest
        % multiple of 3)
        [missing,TF] = listdlg('PromptString',...
                    {'Choose the Load Data files that comes after ', ...
                    'the missing Load Data files. If the first FE Load ' ...
                'Data file is missing, select the first file from the ' ...
                'top. If the last AR Load Data file is missing, delete ' ...
                'the last Optotrak file.'},'ListString',name_loadfiles, ...
                'ListSize', [250 300]);
        grand_index = [1:3:length(sixdofdir);2:3:length(sixdofdir);...
            3:3:length(sixdofdir),0];
        % Fill in padded zeros
        grand_index(end) = grand_index(end-1) + 1;
        % Nullify index where there is missing Load Data files (FE,LB,AR)
        if missing == 1
            grand_index(missing) = 0;
        else
            grand_index(missing-1) = 0;
        end 
        a = 1;
    else % No missing files; the user made it easy on you :)
        grand_index = [1:3:length(sixdofdir);2:3:length(sixdofdir);...
            3:3:length(sixdofdir)];
        a = 0;
    end

    [numrow,numcol] = size(grand_index);
    grand_index2 = reshape(grand_index,[numcol numrow]);
    if ~isempty(grand_index2(grand_index2 == 0)) % Pad any zeros
        grand_index2(grand_index2 == 0) = find(grand_index2 == 0);
    end
%     if grand_index2(end-1)==0 % missing 2/3 load files
%         grand_index2(end-1) = grand_index2(end-2) + 1;
%     end
%     if grand_index2(end)==0 % missing 1/3 load file
%         grand_index2(end) = grand_index2(end-1) + 1;
%     end
    grand_index3 = grand_index2';

%% Defining groups
% When the grand index is created, the number of groups for analysis is 
% automatically calculated (determined by how many triplets are present 
% in the files used for analysis). The test category names are 
% auto-populated with generic names; However, you can easily change these
% names and they will appear in the summary sheet (as long as MATLAB likes
% the name you gave it).

% --- Creating group/specimen names --- %
    groups = cell(numcol,1);
    for g = 1:numcol % auto-populates generic names for number of groups 
        % present in dataset
        groups{g} = sprintf('Group %d' , g);
    end
    dlgtitle = 'Enter Test Category Names';
    dims = [1 35];
    groups = inputdlg(groups,dlgtitle,dims,groups); % allows user to create
    % their own group names

    % If user presses cancel instead, the defaul groups names will be used
    % (i.e. Group 1, Group 2, Group 3...)
    if isempty(groups)
        groups = cell(numcol,1);
        for g = 1:numcol 
            groups{g} = sprintf('Group %d' , g);
        end 
    end 

%% File Selection
% If the user feels so inclined to save some time and skip some files, 
% they can. The user will be prompted to choose from a list of excel 
% sheets from the folder containing the load data (by default you should 
% probably click 'Select All'):
% Files not selected will not have plots or summary data produced 
% (except the first file which is needed to initialize analysis).

% --- Selecting file for analysis --- %
    [indx,~] = listdlg('PromptString','Select Files for Analysis',...
        'ListString',name_loadfiles); % by default, 
    % You should probably select all
    % If the user presses 'Cancel', the code will run with all Load Data 
    % files as selected
    if isempty(indx)   
        indx = 1:length(name_loadfiles);
    end

    x = grand_index(indx+a);
    x1 = ismember(grand_index,x); % find where the selected files are in
    % the Grand Matrix
    grand_index(x1==0) = 0; % excludes unselected data
    Mgrand_index = grand_index;

    % Progress bar
    f = waitbar(0, 'Starting');

%% EXTRACTING DATA FROM FILES FOR ANALYSIS
    %for q = 1:length(sixdofdir)
    for q = 1:length(indx)    
        if grand_index(indx(q))==0 && q~=1
            continue
        end
        pause(1) % code needs a second, got a lot of pressure on its shoulders
        if (ismember(indx(q),skip) && length(sixdofdir)<length(optodir)) || s1 == 1 % less load 
            % data files than Optotrak files
            loadfilepth = fullfile(folderdir,sixdofdir(indx(q)).name);
            data = readtable(loadfilepth);
            data.Time_s_ = round(data.Time_s_,2);
            anglefilepth = fullfile(folderdir2,optodir(indx(q)+numel(skip)).name);
            [row,col] = find(Mgrand_index==indx(q)+numel(skip));
            s1 = 1;
        elseif (ismember(indx(q),skip) && length(sixdofdir)>length(optodir)) || s2 ==1 % less 
            % Optotrak files than load data files
            loadfilepth = fullfile(folderdir,sixdofdir(indx(q)+numel(skip)).name);
            data = readtable(loadfilepth);
            data.Time_s_ = round(data.Time_s_,2);
            anglefilepth = fullfile(folderdir2,optodir(indx(q)).name);
            [row,col] = find(grand_index==indx(q)+numel(skip));
            s2 = 1;
        else
            loadfilepth = fullfile(folderdir,sixdofdir(indx(q)).name);
            data = readtable(loadfilepth);
            data.Time_s_ = round(data.Time_s_,2);
            anglefilepth = fullfile(folderdir2,optodir(indx(q)).name);
            [row,col] = find(grand_index==indx(q));
        end
        
        txtfilepth = strrep(anglefilepth, '.xls', '.txt');
        movefile(anglefilepth, txtfilepth,'f');
        
        try
            txtfilestruct = importdata(txtfilepth);
    %         txtfilestruct = importdata(anglefilepth);
            warning('off','last')
            movefile(txtfilepth, anglefilepth)
            data2 = txtfilestruct.data;
            if numel(data2)<6
                txtfilestruct = importdata(anglefilepth);
                data2 = txtfilestruct.data;
            end 
        catch % Built in redundancy in case user edits the .xls file
            txtfilestruct = importdata(anglefilepth);
%             warning('off','last')
%             movefile(txtfilepth, anglefilepth)
            data2 = txtfilestruct.data;
        end

        data2(:,1) = round(data2(:,1),2); % round Optotrak data to the 2nd
        % decimal place since 6DOF reports to the 2nd decimal place
     
        colheaders = txtfilestruct.textdata{:};
        colheaders = strsplit(colheaders,'\t');

        if size(colheaders,2) < 2
            colheaders = txtfilestruct.textdata;
            colheaders = erase(colheaders, '"' );
        end 
        
        % Finds all columns (i.e., relatives) that have FE data
        FE_index = cellfun(@(x) endsWith(x,'[Ry]'),colheaders, ...
            'Uniformoutput',false);
        FE_indexMAT = cell2mat(FE_index);
        data_FE = find(FE_indexMAT);
    
        % Finds all columns (i.e., relatives) that have LB data
        LB_index = cellfun(@(x) endsWith(x,'[Rz]'),colheaders, ...
            'Uniformoutput',false);
        LB_indexMAT = cell2mat(LB_index);
        data_LB = find(LB_indexMAT);
    
        % Finds all columns (i.e., relatives) that have AR data
        AR_index = cellfun(@(x) endsWith(x,'[Rx]'),colheaders, ...
            'Uniformoutput',false);
        AR_indexMAT = cell2mat(AR_index);
        data_AR = find(AR_indexMAT);
    
        FELBAR_index = [data_FE;data_LB;data_AR];
    
        load = table2array(data(2:end,row+5));
        dof_time = data.Time_s_(2:end);
    
        if q==1
            summarydata = zeros(size(grand_index,2)*3,11, ...
                size(FELBAR_index,2)); % initializes
            % summary data
            relativeNames = cell(1,size(FELBAR_index,2));
        end
    
        % EXTRACTING DATA FOR ALL RELATIVES
        for j = 1:size(FELBAR_index,2) % iterates through data (FE>LB>AR)
            Angle = data2(:,FELBAR_index(row,j));
            relatives = colheaders{FELBAR_index(row,j)};
            relativeNames{j} = relatives(1:19);
            if endsWith(relatives, 'deg [Ry]')
                graph_name = replace(relatives, 'deg [Ry]', 'FE ');
            elseif endsWith(relatives, 'deg [Rz]')
                graph_name = replace(relatives, 'deg [Rz]', 'LB ');
            elseif endsWith(relatives, 'deg [Rx]')
                graph_name = replace(relatives, 'deg [Rx]', 'AR ');
            end
            fgraph_name = join([graph_name, groups{col}]);
        
%             if row ~= 1
%                 Angle = -Angle; % flips data for LB and AR since the test logs 
%                 % motion in the negative direction first
%             end
    
            opto_time = data2(:,1);

% --- How is the data Aligned? --- %
% The Optotrak is usually started a few seconds before the 6DOF, which
% necessitates shifting the data. To account for the time shift, the peaks
% in angular displacement are matched with the peaks in applied load. 
% The findpeaks function locates the peaks by finding 3 main peaks which
% have a minimum separation in time of about 1/3 of the length of the test
% (hence, I/3.5--where I is the length of the test). You may or may not
% know this, but the sampling frequency of the Optotrak is different from
% the 6DOF. Moreover, the 6DOF records points less frequently than the
% Optotrak and requires some Optotrak data to be cut (i.e., how can you
% capture load/displacement data if there is no load data at that time).
% The code indexes only the points in the Optotrak data that share a time
% point during loadinf cycles captured by the 6DOF (this is completed 
% after the data is shifted).
  
            % SHIFTING OPTOTRAK DATA FOR ALIGNMENT WITH LOADING/UNLOADING 
            % CYCLES
%             if j == 1
                sectioner = 3.5;
                [~,I] = max(dof_time);
% %                 [smoothpks,smoothlocs] = findpeaks(load,'MinPeakDistance',...
% %                     I/sectioner,'NPeaks',3); % 3 main peaks for load
                [~,I2] = max(opto_time);
% %                 [smoothpks2,smoothlocs2] = findpeaks(Angle,'MinPeakDistance',...
% %                     I2/sectioner,'NPeaks',3); % 3 main peaks for displacement
% %                 [smoothvals,smoothlocs3] = findpeaks(-Angle,'MinPeakDistance',...
% %                     I2/sectioner,'NPeaks',3); % 3 main valleys for displacement
    
%% Added by Izzy to replace findpeaks with localmax
                [smoothpks,smoothlocs]=islocalmax(load, 'MinSeparation', I/sectioner);

                NPeaks = 3;
                if sum(smoothpks) > NPeaks
                    [~, sortedIndices] = sort(load(smoothpks), 'descend'); % Sort peak heights in descending order
                    peakIndices = find(smoothpks); % Get indices of all local maxima
                    topNPeakIndices = peakIndices(sortedIndices(1:NPeaks)); % Get indices of the top N peaks
                    smoothpks = load(topNPeakIndices); % Extract peak values
                    smoothlocs = topNPeakIndices; % Indices of the top N peaks
                    smoothlocs= sort(smoothlocs,'ascend');
                else
                    % If there are fewer than NPeaks, consider all local maxima as peaks
                    smoothlocs = find(smoothpks);
                    smoothpks = load(smoothpks);
                end

                [smoothpks2,smoothlocs2]=islocalmax(Angle,'MinSeparation',length(Angle)/10);
%                
                NPeaks = 3;
                if sum(smoothpks2) >= NPeaks
                    [~, sortedIndices] = sort(Angle(smoothpks2), 'descend'); % Sort peak heights in descending order
                    peakIndices = find(smoothpks2); % Get indices of all local maxima
                    topNPeakIndices = peakIndices(sortedIndices(1:NPeaks)); % Get indices of the top N peaks
                    smoothpks2 = Angle(topNPeakIndices); % Extract peak values
                    smoothlocs2 = topNPeakIndices; % Indices of the top N peaks
                    smoothlocs2=sort(smoothlocs2,'ascend');
                    
                else
                    % If there are fewer than NPeaks, consider all local maxima as peaks
                    smoothlocs2 = find(smoothpks2);
                    smoothpks2 = Angle(smoothpks2);
                end

                [smoothvals,smoothlocs3]=islocalmax(-Angle,'MinSeparation',length(Angle)/10);
%                
                NPeaks = 4;
                if sum(smoothvals) >= NPeaks
                    [~, sortedIndices] = sort(-Angle(smoothvals), 'descend'); % Sort peak heights in descending order
                    peakIndices = find(smoothvals); % Get indices of all local maxima
                    topNPeakIndices = peakIndices(sortedIndices(1:NPeaks)); % Get indices of the top N peaks
                    smoothvals = -Angle(topNPeakIndices); % Extract peak values
                    smoothlocs3 = (topNPeakIndices); % Indices of the top N peaks
                    smoothlocs3=sort(smoothlocs3,'ascend');
                else
                    % If there are fewer than NPeaks, consider all local maxima as peaks
                    smoothlocs3 = find(smoothvals);
                    smoothvals = -1*-Angle(smoothvals); %,might not need to times by -1
                end

                %%

                dof_st = dof_time(smoothlocs(1)); % find the first peak load
                opto_TstLoc = find(opto_time==dof_st); % finds the location in 
                % the Optotrak time-data matrix where the time of the first 
                % peak occurs in the load cycles
                opto_PKaLoc1 = smoothlocs2(1); % location in the displacement 
                % matrix where the first peak in displacement occurs
                opto_start = opto_PKaLoc1 - opto_TstLoc;

                % If the Angle data starts with a negative value, flip data
                if row > 1 && opto_start <= 0 || row > 1 && smoothlocs3(1) < smoothlocs2(1)   
%                     while opto_start <= 0
                        try 
%                             questdlg(['The program encountered errors when ' ...
%                                 'looking for peaks in OptoTrak data. ' ...
%                                 'Attempting to flip the OptoTrak data ' ...
%                                 '(Negative values become positive) ' ...
%                                 'in order to correct.'], ...
% 	                             'Error Finding OptoTrak Peaks','OK','OK');
                            Angle = -Angle;
                            [~,smoothlocs2] = findpeaks(Angle, ...
                                'MinPeakDistance',I2/sectioner,'NPeaks',3); 
                            % 3 main peaks for displacement
                            opto_PKaLoc1 = smoothlocs2(1); 
                            % location in the displacement matrix where 
                            % the first peak in displacement occurs
                            opto_start = opto_PKaLoc1 - opto_TstLoc;
%                             sectioner = sectioner - 0.05;
                        catch
%                             sectioner = 3.5;
                            break
                        end 
%                     end
                 end
                ls_opto_time = round(opto_time-opto_time(opto_start),2);
                aindex = ismember(ls_opto_time,dof_time); % finds all time 
                % points where displacement is captured at the same frequency
                % as load
                opto_time_f = ls_opto_time(aindex);
                
%             end 
            s_angle = Angle(aindex);
            angle_f = s_angle-s_angle(1); % zeros the displacement

% --- The plots --- %
% The data is normally saved as plots. Displacement (measured in degrees)
% is plotted along the y-axis, Load (the moment/torque, measured in Nm) 
% is plotted along the x-a
% xis. The red section indicates the region with
% the least amount of resistance (i.e. NZ) for both postive and negative
% loading. The dashed black lines show the angle boundaries of the neutral 
% zone. The purple/magenta section represents the region of highest 
% resistance (i.e., EZS) for both positive and negative loading. The solid
% black line show the slope that represents NZS. The green line shows the
% slope that represents EZS for positive and negative loading. 

% There are some circumstances when the displacement data will be so close
% in magnitude that the NZ cannot be determined mathematically.
% In these such instances, the Raw Data plots will be saved. The top plot 
% shows the last cycle and is in red to indicate there will be a problem 
% calculating measurments (will show up as a row of 0s in the summary
% data sheet) The bottom plot shows data for the Optotrak and the 6DOF, 
% simultaneously, for the full test duration.


% PLOTTING DATA
            curfig = figure('visible','off');
            response2 = 0;
            loadingNZ = [];
            unloadingNZ = [];
                   
            smANGLE = smoothdata(angle_f,'movmean',15); % smooths data to allow
            % for clear vizualization for distinguishing between NZ and EZ 
            
            load_a = abs(load);

            try 
                [smoothtroughs,smoothlocs1] = findpeaks(-load_a, ...
                'MinPeakDistance',numel(load_a)/12,...
                'MinPeakHeight',-0.1); % finds location of zero (0) load
                last_cycleSTART = smoothlocs1(end-2);

                zero_load_index = -0.1;
                while length(smoothlocs1) < 6
                    zero_load_index = zero_load_index - 0.1;
                    [smoothtroughs,smoothlocs1] = findpeaks(-load_a, ...
                    'MinPeakDistance',numel(load_a)/12,...
                    'MinPeakHeight',zero_load_index); % finds location of zero (0) load
                    last_cycleSTART = smoothlocs1(end-2);
                end
            catch
                try 
                    [smoothtroughs,smoothlocs1] = findpeaks(-load_a, ...
                    'MinPeakDistance',numel(load_a)/12,...
                    'MinPeakHeight',-0.3); % finds location of zero (0) load
                    last_cycleSTART = smoothlocs1(end-2);
                catch
                    [smoothtroughs,smoothlocs1] = findpeaks(-load_a, ...
                    'MinPeakDistance',numel(load_a)/12,...
                    'MinPeakHeight',-0.5); % finds location of zero (0) load
                    last_cycleSTART = smoothlocs1(end-2);
                end
            end
            
            % olast_cycleLOAD = load(last_cycleSTART:end);
            % olast_cycleANGLE = smANGLE(last_cycleSTART:end);
            last_cycleLOAD = load(last_cycleSTART:smoothlocs1(end));
            last_cycleANGLE = smANGLE(last_cycleSTART:smoothlocs1(end));
            last_cycleTIME = dof_time(last_cycleSTART:smoothlocs1(end));
            
            % Save raw data to .xlsx
            T = table(opto_time_f, load, angle_f, 'VariableNames', ...
                {'Time (s)','6DOF Load (Nm)', ...
                'OptoTrak Displacement (deg)'});
             writetable(T, fullfile(folderdir,newFolderName, ...
                 join([fgraph_name ' - Raw Data.xlsx'])),'Sheet',1, ...
                 'Range','A1')
             T2 = table(last_cycleTIME, last_cycleLOAD, last_cycleANGLE, ...
                 'VariableNames',{'Time (s)','Last Cycle 6DOF Load (Nm)', ...
                 'Last Cycle OptoTrak Displacement (deg)'});
             writetable(T2, fullfile(folderdir,newFolderName, ...
                 join([fgraph_name ' - Raw Data.xlsx'])),'Sheet',2, ...
                 'Range','A1')

            % Plot raw data 
            rawdata_fig = figure('visible','on');
            subplot(2,1,1);
            plot(last_cycleLOAD, last_cycleANGLE,'r');
            title(fgraph_name);
            subplot(2,1,2);
            plot(opto_time_f, angle_f);
            hold on
            plot(dof_time, load,'r');
            plot(last_cycleTIME, last_cycleANGLE)
            plot(last_cycleTIME, last_cycleLOAD)
            lgd = legend('Relative Angle','Load',...
                'Last Cycle: Relative Angle', 'Last Cycle: Load');
            lgd.FontSize = 5;

            % Saves fig and .jpg to new subfolder
            saveas(rawdata_fig,fullfile(folderdir,newFolderName, ...
                       join([fgraph_name ' - Raw Data.jpg'])),'jpg'); 
            saveas(rawdata_fig,fullfile(folderdir,newFolderName, ...
                   join([fgraph_name ' - Raw Data.fig'])),'fig');

            cur_axis = axis; % get current axis to then reset 
                    % it after zooming in

            pause(0.5);
            w = questdlg(['Zoom into the bottom plot to view proper' ...
                ' alignment of the last cycle of OptoTrak data with the' ...
                ' last cycle of load data. Then press Enter to continue.'], ...
                    'Verify Last Cycle Alignment','OK','OK');

            % ALLOW USER TO ZOOM
            buttonwait = 0;
            while ~buttonwait
                buttonwait = waitforbuttonpress;
                if ~strcmp(get(gcf,'CurrentKey'),'return')
                    buttonwait = 0;
                    if ishandle(w)
                        close(w);
                    end
                    pause(0.5);
                end
            end

            % Ask if user would like to manually select features

            answer1 = questdlg(['Would you like to manually select ' ...
                 'the last cycle of the load and OptoTrak data?'], ...
	                'Manual Data Selection', ...
	                'YES','NO','NO');

            % Handle response 1
            switch answer1
                case 'YES'
                     response1 = 1;
                case 'NO'
                     response1 = 0;
            end 

            data_selection_fig = figure('visible','off');
            while response1 == 1 && ishandle(data_selection_fig)
                try
                    ls_opto_time = round(opto_time,2);
%                     aindex2 = ismember(ls_opto_time2,dof_time); % finds all time 
                    aindex =  rem(ls_opto_time, 0.1) == 0;
                    % points where displacement is captured at the same frequency
                    % as load
                    opto_time_f = ls_opto_time(aindex);
%                     opto_time_f2 = ls_opto_time2;
                    s_angle = Angle(aindex);
                    angle_f = s_angle-s_angle(1); % zeros the displacement
%                     angle_f2 = Angle - Angle(1); % zeros the displacement
                    smANGLE = smoothdata(angle_f,'movmean',15);

                    figure(data_selection_fig);
                    set(data_selection_fig, 'visible','on')
                    plot(opto_time_f, smANGLE,'o-')
                    hold on
                    grid on
                    plot(linspace(min(opto_time_f), ...
                        max(opto_time_f),10),zeros(10,1),'k--')
                    xlabel('Time (s)');
                    ylabel('Relative Angle (°)');
                    lgd = legend('Relative Angle','X-axis');
                    lgd.FontSize = 10;
                    cur_axis = axis; % get current axis to then reset 
                    % it after zooming in
                    title(['Press Enter. Then select the start of the ' ...
                        'last cycle of the OptoTrak data:']);
                    
                    % ALLOW USER TO ZOOM 
                    buttonwait = 0;
                    while ~buttonwait
                        buttonwait = waitforbuttonpress;
                        if ~strcmp(get(gcf,'CurrentKey'),'return')
                            buttonwait = 0;
                        end

                        % Flip angle data if it starts with a valley
%                         if strcmp(get(gcf,'CurrentKey'),'-')
%                             clf(data_selection_fig)
%                             plot(opto_time_f, smANGLE,'o-')
%                             hold on
%                             grid on
%                             plot(linspace(min(opto_time_f), ...
%                                 max(opto_time_f),10),zeros(10,1),'k--')
%                             xlabel('Time (s)');
%                             ylabel('Relative Angle (°)');
%                             lgd = legend('Relative Angle','X-axis');
%                             lgd.FontSize = 10;
%                             cur_axis = axis; % get current axis to then reset 
%                             % it after zooming in
%                             title(['Press Enter. Then select the start of the ' ...
%                                 'last cycle of the OptoTrak data:']);
%                         end
                    end
    
                    pause(1)
                    pass = 0;
                    while pass == 0
                        try 
                            rect = getrect(data_selection_fig);
                            x1 = rect(1);
                            y1 = rect(2);
                            x2 = rect(3) + x1;
                            y2 = rect(4) + y1;
                        catch
                            break;
                        end 
                        
                        opto_start = find((opto_time_f >= x1) & ... 
                            (opto_time_f <= x2) & (smANGLE >= y1) & ...
                            (smANGLE <= y2));
                        if isempty(opto_start)
                            questdlg(['Please be sure to draw a box ' ...
                            'with points in it.'], ...
                            'Warning: No Data Selected','OK','OK');
                            pass = 0;
                        else
                            opto_start_idx = abs(smANGLE(opto_start)) == min(abs(smANGLE(opto_start)));
                            opto_start = opto_start(opto_start_idx);
                            % Plot the selected OptoTrak starting point
                            plot(opto_time_f(opto_start), smANGLE(opto_start),'m*');
                            lgd = legend('Relative Angle', 'X-axis', ...
                                'OptoTrak Last Cycle Start');
                            pass = 1;
                        end
                    end 
    
                    title(['Press Enter. Then select the end of the '...
                        'last cycle of the OptoTrak data:']);
    
                    % ALLOW USER TO ZOOM
                    buttonwait = 0;
                    while ~buttonwait


                        buttonwait = waitforbuttonpress;
                        if ~strcmp(get(gcf,'CurrentKey'),'return')
                            buttonwait = 0;
                        end
                    end
    
                    pause(1)
                    pass = 0;
                    while pass == 0
                        try 
                            rect = getrect(data_selection_fig);
                            x1 = rect(1);
                            y1 = rect(2);
                            x2 = rect(3) + x1;
                            y2 = rect(4) + y1;
                        catch
                            break;
                        end 
        
                        opto_end = find((opto_time_f >= x1) & ...
                            (opto_time_f <= x2) & (smANGLE >= y1) & ...
                            (smANGLE <= y2));
                        
                        if isempty(opto_end)
                            questdlg(['Please be sure to draw a box ' ...
                            'with points in it.'], ...
                            'Warning: No Data Selected','OK','OK');
                            pass = 0;
                        else
                            opto_end_idx = abs(smANGLE(opto_end)) == min( ...
                                abs(smANGLE(opto_end)));
                            opto_end = opto_end(opto_end_idx);
    
                            if opto_end < opto_start || opto_end == opto_start
                                g = questdlg(['The last cycle of the Optotrak data' ...
                                    ' could not be found. Please select a new ' ...
                                    'start and end point or close out of the figure.'], ...
                                    ['Invalid: End of 3rd cycle cannot come before' ...
                                    'the beginning of the 3rd cycle.'],'OK','OK');
                                pass = 0;
        
                            else
                            % Plot the selected OptoTrak starting point
                                pass = 1;
                                last_cycleTIME = opto_time_f(opto_start:opto_end);
                                last_cycleANGLE = smANGLE(opto_start:opto_end);
                                Angle_peak_idx = find(last_cycleANGLE == max( ...
                                    last_cycleANGLE)) + opto_start;
                                
                                plot(opto_time_f(opto_end), smANGLE(opto_end),'k*');
                                plot(last_cycleTIME, last_cycleANGLE,'ro-');
                                plot(opto_time_f(Angle_peak_idx), ...
                                    smANGLE(Angle_peak_idx), 'k*-')
                                lgd = legend('Relative Angle', 'X-axis',...
                                    'OptoTrak Last Cycle Start', ...
                                    'OptoTrak Last Cycle End', ...
                                    'OptoTrak Last Cycle',...
                                    'OptoTrak Last Cycle Peak');
        
                                pause(2)
                            end
                        end 
                    end 
                    
                    clf(data_selection_fig)
                    figure(data_selection_fig);
                    plot(dof_time, load,'o-')
                    hold on
                    grid on
                    plot(linspace(min(dof_time), ...
                        max(dof_time),10),zeros(10,1),'k--')
                    xlabel('Time (s)');
                    ylabel('Load (Nm)');
                    lgd = legend('Load', 'X-axis');
                    lgd.FontSize = 10;
                    cur_axis = axis; % get current axis to then reset 
                    % it after zooming in
                    title(['Press Enter. Then select the peak of the ' ...
                        'last cycle of the load data:']);
                    
                    % ALLOW USER TO ZOOM
                    buttonwait = 0;
                    while ~buttonwait
                        
                        buttonwait = waitforbuttonpress;
                        if ~strcmp(get(gcf,'CurrentKey'),'return')
                            buttonwait = 0;
                        end
                    end
    
                    pause(1)

                    pass = 0;
                    while pass == 0
                        try 
                            rect = getrect(data_selection_fig);
                            x1 = rect(1);
                            y1 = rect(2);
                            x2 = rect(3) + x1;
                            y2 = rect(4) + y1;
                        catch
                            break;
                        end

                        load_peak = find((dof_time >= x1) & (dof_time <= x2) & ...
                            (load >= y1) & (load <= y2));

                        if isempty(load_peak)
                            questdlg(['Please be sure to draw a box ' ...
                            'with points in it.'], ...
                            'Warning: No Data Selected','OK','OK');
                            pass = 0;
                        else
                            % Find the peak of the last cycle of the load
                            % data then use the length of the last cycle
                            % angle data to back calculate the start-end
                            % points of the last cycle load data. Both
                            % vectors should be the same length and aligned
                            % based on the timing of their peaks
                            load_peak_idx = find(load == max(load(load_peak)));
                            load_start = load_peak_idx(end) - ( ...
                                Angle_peak_idx - opto_start);
                            load_end = load_peak_idx(end) + ( ...
                                opto_end - Angle_peak_idx);


%                             load_start_idx = abs(load(load_start)) == min( ...
%                                 abs(load(load_start)));
%                             load_start = load_start(load_start_idx);
                            % Plot the selected OptoTrak starting point
                            plot(dof_time(load_start), load(load_start),'m*')
                            lgd = legend('Load', 'X-axis', ...
                                'Load Last Cycle Start');
                            pass = 1;
                        end
                    end 
    
%                     pause(1)
%                     pass = 0;
%                     while pass == 0
                    % Plot the selected OptoTrak starting point
%                         load_end = load_start + numel(last_cycleANGLE)-1;
                        if load_end > numel(dof_time)
                                % If there is less load data then OptoTrak
                                % data, chop off the overhang OptoTrak data
                               load_end = numel(dof_time);
                               last_cycleANGLE = last_cycleANGLE(1:numel(load_start:load_end));
                        end
                        last_cycleTIME = dof_time(load_start:load_end);
                        last_cycleLOAD = load(load_start:load_end);
                        
                        plot(dof_time(load_end), load(load_end),'k*');
                        plot(last_cycleTIME,last_cycleLOAD,'ro-');
                        lgd = legend('Load', 'X-axis',...
                                'Load Last Cycle Start', ...
                                'Load Last Cycle End', ...
                                'Load Last Cycle');
    
                        pause(2);
                        set(data_selection_fig,'visible','off');
                        clf(data_selection_fig);
%                     end 
    
                    % Plot raw data 
                    clf(rawdata_fig)
                    figure(rawdata_fig);
                    subplot(2,1,1);
                    plot(last_cycleLOAD, last_cycleANGLE,'r');
                    title(fgraph_name);
                    subplot(2,1,2);
                    plot(opto_time_f, angle_f);
                    hold on
                    plot(dof_time, load,'r');
                    plot(last_cycleTIME, last_cycleANGLE)
                    plot(last_cycleTIME, last_cycleLOAD)
                    lgd = legend('Relative Angle','Load',...
                        'Last Cycle: Relative Angle', 'Last Cycle: Load');
                    lgd.FontSize = 5;
                   
                    subanswer = questdlg(['Are you satisfied with the ' ...
                            'selected last cycle of the OptoTrak data?'], ...
                            'Manual Data Selection', ...
                            'YES','NO','NO');
                        
                        switch subanswer
                            case 'YES'
                                response1 = 0;
                            case 'NO'
                                response1 = 1;
                        end 
    
                    % Saves fig and .jpg to new subfolder
                    saveas(rawdata_fig,fullfile(folderdir,newFolderName, ...
                               join([fgraph_name ' - Raw Data.jpg'])),'jpg'); 
                    saveas(rawdata_fig,fullfile(folderdir,newFolderName, ...
                           join([fgraph_name ' - Raw Data.fig'])),'fig');
%                     close(rawdata_fig);
                catch
                    response1 = 0;
                    break
                end 
            end 
            
            close(data_selection_fig)
            figure(rawdata_fig);
            answer2 = questdlg(['Would you like to draw the region of ' ...
                 'interest (ROI) of the neutral zone?'], ...
	                'Manual Neutral Zone Selection', ...
	                'YES','NO','NO');
            close(rawdata_fig);

            % Handle response 2
            switch answer2
                case 'YES'
                    response2 = 1;
                case 'NO'
                    response2 = 0;
            end 

            pos_NZB = smANGLE(smoothlocs1(end-1));
            neg_NZB = smANGLE(last_cycleSTART);
            NZ_index = last_cycleANGLE<=pos_NZB & last_cycleANGLE>=neg_NZB;
            NZlocs = find(NZ_index); 

            clf(curfig)
              
            if isempty(NZlocs) || abs(pos_NZB-neg_NZB) < 0.1 || response2 == 1
                % 1.) Display the problematic hysteresis curve
                while response2 == 1
                     try
%                         clf(curfig);
                        figure(curfig);
                        plot(last_cycleLOAD,last_cycleANGLE,'o-')
                        hold on
                        plot(zeros(10,1), linspace(min(last_cycleANGLE), ...
                            max(last_cycleANGLE),10),'k--')
                        xlabel('Load (Nm)');
                        ylabel('Relative Angle (°)');
                        cur_axis = axis; % get current axis to then reset 
                        % it after zooming in
                        hold on
                        grid on
                        title('Press Enter. Then select the (+) Neutral Zone Boundary:');
                        
                        % ALLOW USER TO ZOOM
                        buttonwait = 0;
                        while ~buttonwait
                            buttonwait = waitforbuttonpress;
                            if ~strcmp(get(gcf,'CurrentKey'),'return')
                                buttonwait = 0;
                            end
                        end
                     
                % 2.) Prompt the user to manually select a ROI for the (+) NZB:
                %     ([x1, y1, x2, y2] = getrect(curfig))
                        pause(1)
                        pos_pass = 0;
                        clf(curfig);
                        
                        while pos_pass == 0
                            figure(curfig);
                            plot(last_cycleLOAD,last_cycleANGLE,'o-')
                            hold on
                            plot(zeros(10,1), linspace(min(last_cycleANGLE), ...
                                max(last_cycleANGLE),10),'k--')
                            xlabel('Load (Nm)');
                            ylabel('Relative Angle (°)');
                            cur_axis = axis; % get current axis to then reset 
                            % it after zooming in
                            hold on
                            grid on
                            title('Press Enter. Then select the (+) Neutral Zone Boundary:');

                            try
                                rect_pos = getrect(curfig); % [xmin ymin width height]
                                x1_pos = rect_pos(1);
                                y1_pos = rect_pos(2);
                                x2_pos = rect_pos(3) + x1_pos;
                                y2_pos = rect_pos(4) + y1_pos;
                                %display([x1_pos,y1_pos,x2_pos,y2_pos])
                            catch
                                pos_pass = 1;
                                break;
                            end 
        
                    % 3.) In the ROI for the (+) NZB, find the angles (y's) 
                    % where the load (x's) is closest to zero (i.e. lowest 
                    % absolute value)
                            pos_NZB_range = last_cycleLOAD>=x1_pos ...
                                & last_cycleLOAD<=x2_pos ...
                                & last_cycleANGLE>=y1_pos ...
                                & last_cycleANGLE<=y2_pos;
                            if sum(pos_NZB_range) == 0 
                                g = questdlg(['The NZ cannot be determined from the' ...
                                    ' ROI selected. Please ' ...
                                    'draw new ROIs for the NZBs or close out of' ...
                                    ' the next figure.'],['Warning: Cannot ' ...
                                    'Identify Neutral Zone'],'OK','OK');
                                pos_pass = 0;
                            elseif (x1_pos > 0 && x2_pos > 0) || ( ...
                                x1_pos < 0 && x2_pos < 0)  
                                g = questdlg(['The NZ cannot be determined from the' ...
                                    ' ROI selected. Please ' ...
                                    'draw a new ROI for the (+) NZB that contains' ...
                                    ' the zero-load point (x-axis) or close out of' ...
                                    ' the next figure.'],['Warning: Cannot ' ...
                                    'Indentiy Zero-load Point'],'OK','OK');
                                pos_pass = 0;
                            else
                                pos_pass = 1;
                                posNZB_load_range = last_cycleLOAD(pos_NZB_range);
                                posNZB_angle_range = last_cycleANGLE(pos_NZB_range);                
                                pos_NZB = posNZB_angle_range(abs( ...
                                    posNZB_load_range)==min(abs(posNZB_load_range)));
                                if numel(pos_NZB) > 1
                                    pos_NZB = max(pos_NZB);
                                end

                       % 4.) Plot the (+) NZB:         
                                plot(linspace(min(last_cycleLOAD),max( ...
                                    last_cycleLOAD),20), pos_NZB*ones(1,20), 'k--')
                           end
                        end 
                    
                        title('Press Enter. Then select the (-) Neutral Zone Boundary:')
                        
                        % ALLOW USER TO ZOOM
                        buttonwait = 0;
                        while ~buttonwait
                            buttonwait = waitforbuttonpress;
                            if ~strcmp(get(gcf,'CurrentKey'),'return')
                                buttonwait = 0;
                            end
                        end
                        pause(1)
                        
                      % 5.) Prompt the user to manually select a ROI for the (-) NZB:
                        neg_pass = 0;
                        while neg_pass == 0 
                            try
                                rect_neg = getrect(curfig);
                                x1_neg = rect_neg(1);
                                y1_neg = rect_neg(2);
                                x2_neg = rect_neg(3) + x1_neg;
                                y2_neg = rect_neg(4) + y1_neg;
                                %display([x1_neg,y1_neg,x2_neg,y2_neg])
                            catch
                                neg_pass = 1;
                                break;
                            end 
        
                    % 6.) In the ROI for the (-) NZB, find the angles (y's) where 
                    % the load (x's) is closest to zero (i.e. lowest absolute
                    % value)
                            neg_NZB_range = last_cycleLOAD>=x1_neg ...
                                & last_cycleLOAD<=x2_neg ...
                                & last_cycleANGLE>=y1_neg ...
                                & last_cycleANGLE<=y2_neg;
                                
                            if sum(neg_NZB_range) == 0 
                                g = questdlg(['The NZ cannot be determined from the' ...
                                    ' ROI selected. Please ' ...
                                    'draw new ROIs for the NZBs or close out of' ...
                                    ' the next figure.'],['Warning: Cannot ' ...
                                    'Identify Neutral Zone'],'OK','OK');
                                neg_pass = 0;
                             elseif (x1_neg > 0 && x2_neg > 0) || ( ...
                                x1_neg < 0 && x2_neg < 0)  
                                g = questdlg(['The NZ cannot be determined from the' ...
                                    ' ROI selected. Please ' ...
                                    'draw a new ROI for the (+) NZB that contains' ...
                                    ' the zero-load point (x-axis) or close out of' ...
                                    ' the next figure.'],['Warning: Cannot ' ...
                                    'Indentiy Zero-load Point'],'OK','OK');
                                neg_pass = 0;
                            else
                                neg_pass = 1;
                                negNZB_load_range = last_cycleLOAD(neg_NZB_range);
                                negNZB_angle_range = last_cycleANGLE(neg_NZB_range);                
                                neg_NZB = negNZB_angle_range(abs( ...
                                    negNZB_load_range)==min(abs(negNZB_load_range)));

                                if numel(neg_NZB) > 1
                                    neg_NZB = min(neg_NZB);
                                end
                                
                                set(curfig,'visible','on')
                                plot(linspace(min(last_cycleLOAD),max( ...
                                    last_cycleLOAD),20), neg_NZB*ones(1,20), 'k--')
                                pause(1)
                                set(curfig,'visible','off')
                            end 
                        end
                        % Find the NZ between the NZ bounds 
                        NZ_index = last_cycleANGLE<pos_NZB ...
                            & last_cycleANGLE>neg_NZB;
                        NZlocs = find(NZ_index);
                     
                        % If the user selects an ROI that does not include (+) and
                        % (-) load, prompt them to reselect the ROI
                        if isempty(NZlocs) 
                            g = questdlg(['The NZ cannot be determined from the' ...
                                ' ROI selected. Please ' ...
                                'draw new ROIs for the NZBs or close out of' ...
                                ' the next figure.'],['Warning: Cannot ' ...
                                'Identify Neutral Zone'],'OK','OK');
                            response2 = 1;
                        else
                            response2 = 0;
        
                            % Separate NZ; one for positive loading and one for
                            % negative loading
                            loadingNZ = NZlocs(last_cycleLOAD(NZlocs)>0);
                            % Drop any overhang in the loading NZ
                            % This can happen when the 6DOF motion simulator moves
                            % back into position and some (+) angle is collected as
                            % seen from the Optotrack

%%      Prompt user to select the end of the loading NZ

%                             set(curfig,'visible','on')
                            figure(curfig);
%                             title(['Press Enter. Then Select the ' ...
%                                 'intersection between the ' ...
%                                 '(+) Neutral Zone Boundary and the loading ' ...
%                                 'Neutral Zone:'])
%                         
%                             % ALLOW USER TO ZOOM
%                             buttonwait = 0;
%                             while ~buttonwait
%                                 buttonwait = waitforbuttonpress;
%                                 if ~strcmp(get(gcf,'CurrentKey'),'return')
%                                     buttonwait = 0;
%                                 end
%                             end
%                             pause(1)
%                             pos_pass = 0;
%                         while pos_pass == 0 
%                             try
%                                 rect_pos = getrect(curfig);
%                                 x1_pos = rect_pos(1);
%                                 y1_pos = rect_pos(2);
%                                 x2_pos = rect_pos(3) + x1_pos;
%                                 y2_pos = rect_pos(4) + y1_pos;
%                                 %display([x1_neg,y1_neg,x2_neg,y2_neg])
%                             catch
%                                 pos_pass = 1;
%                                 break;
%                             end 
% 
%                             loadingNZ_end_range = last_cycleLOAD>=x1_pos ...
%                                  & last_cycleLOAD<=x2_pos ...
%                                 & last_cycleANGLE>=y1_pos ...
%                                  & last_cycleANGLE<=y2_pos ...
%                                 & y1_pos <= pos_NZB;
%                             
%                             if sum(loadingNZ_end_range) == 0 
%                                 g = questdlg(['The loading NZ cannot be determined from the' ...
%                                     ' ROI selected. Please ' ...
%                                     'draw new ROIs or close out of' ...
%                                     ' the next figure.'],['Warning: Cannot ' ...
%                                     'Identify the loading Neutral Zone'], ...
%                                     'OK','OK');
%                                 pos_pass = 0;
%                             else
%                                 pos_pass = 1;
%                                 loadingNZ_end_idx = find(loadingNZ_end_range~=0);
                                loadingNZ_end_idx = [find( ...
                                    last_cycleANGLE == min( ...
                                    last_cycleANGLE)):numel( ...
                                    last_cycleANGLE), 1:find( ...
                                    last_cycleANGLE == max( ...
                                    last_cycleANGLE))];
                                loadingNZ = loadingNZ_end_idx(last_cycleANGLE( ...
                                    loadingNZ_end_idx)<= pos_NZB & last_cycleANGLE(...
                                    loadingNZ_end_idx) >= neg_NZB);

                                if numel(loadingNZ) <= 20
                                    loadingNZ = [loadingNZ_end_idx(find( ...
                                        loadingNZ_end_idx == loadingNZ(1))-3),...
                                        loadingNZ_end_idx(find(...
                                    loadingNZ_end_idx == loadingNZ(1))-2), ...
                                        loadingNZ_end_idx(find(...
                                    loadingNZ_end_idx == loadingNZ(1))-1), ...
                                        loadingNZ, (loadingNZ(end)+1):(...
                                        loadingNZ(end)+3)];
                                end
                                
                                set(curfig,'visible','on')
                                plot(last_cycleLOAD(loadingNZ'), ...
                                    last_cycleANGLE(loadingNZ'), 'ro-')
%                             end
%                         end 

%%      Prompt user to select the end of the unloading NZ

%                             set(curfig,'visible','on')
                            figure(curfig);
%                             title(['Press Enter. Then Select the ' ...
%                                 'intersection between the ' ...
%                                 '(-) Neutral Zone Boundary and the unloading ' ...
%                                 'Neutral Zone:'])
%                         
%                             % ALLOW USER TO ZOOM
%                             buttonwait = 0;
%                             while ~buttonwait
%                                 buttonwait = waitforbuttonpress;
%                                 if ~strcmp(get(gcf,'CurrentKey'),'return')
%                                     buttonwait = 0;
%                                 end
%                             end
%                             pause(2)
%                             neg_pass = 0;
%                         while neg_pass == 0 
%                             try
%                                 rect_neg = getrect(curfig);
%                                 x1_neg = rect_neg(1);
%                                 y1_neg = rect_neg(2);
%                                 x2_neg = rect_neg(3) + x1_neg;
%                                 y2_neg = rect_neg(4) + y1_neg;
%                                 %display([x1_neg,y1_neg,x2_neg,y2_neg])
%                             catch
%                                 neg_pass = 1;
%                                 break;
%                             end 
% 
%                             unloadingNZ_end_range = last_cycleLOAD>=x1_neg ...
%                                  & last_cycleLOAD<=x2_neg ...
%                                 & last_cycleANGLE>=y1_neg ...
%                                  & last_cycleANGLE<=y2_neg ...
%                                  & y2_neg >= neg_NZB;
%                             
%                             if sum(unloadingNZ_end_range) == 0 
%                                 g = questdlg(['The unloading NZ cannot be determined from the' ...
%                                     ' ROI selected. Please ' ...
%                                     'draw new ROIs or close out of' ...
%                                     ' the next figure.'],['Warning: Cannot ' ...
%                                     'Identify the unloading Neutral Zone'], ...
%                                     'OK','OK');
%                                 neg_pass = 0;
%                             else
%                                 neg_pass = 1;
%                                 unloadingNZ_end_idx = find(unloadingNZ_end_range~=0);
                                unloadingNZ_end_idx = find( ...
                                last_cycleANGLE == max(...
                                last_cycleANGLE)):find(...
                                last_cycleANGLE == min(last_cycleANGLE));
                                 unloadingNZ = unloadingNZ_end_idx(last_cycleANGLE( ...
                                    unloadingNZ_end_idx) >= neg_NZB & last_cycleANGLE( ...
                                    unloadingNZ_end_idx) <= pos_NZB);

                                if numel(unloadingNZ) <= 20
                                    unloadingNZ = [unloadingNZ_end_idx(find( ...
                                        unloadingNZ_end_idx == unloadingNZ(1))-3),...
                                        unloadingNZ_end_idx(find(...
                                    unloadingNZ_end_idx == unloadingNZ(1))-2), ...
                                        unloadingNZ_end_idx(find(...
                                    unloadingNZ_end_idx == unloadingNZ(1))-1), ...
                                        unloadingNZ, (unloadingNZ(end)+1):(...
                                        unloadingNZ(end)+3)];
                                end

                                set(curfig,'visible','on');
                                plot(last_cycleLOAD(unloadingNZ), ...
                                    last_cycleANGLE(unloadingNZ), 'ro-')
                                pause(2)
%                             end
%                         end 

                        clf(curfig)
                        close(curfig)
                        end
                    
                    catch
                        response2 = 0;
                        break 
                    end
                end 
            end 

            % --- Plotting hysteresis curve (curfig) --- %
            if ~ishandle(curfig)
                curfig = figure('visible','on');
            end
            figure(curfig);
            clf(curfig);
            plot(last_cycleLOAD,last_cycleANGLE);
            title(fgraph_name);
            xlabel('Load (Nm)');
            ylabel('Relative Angle (°)');
            hold on
                
            % --- Plotting for NZ --- %
            % According to literature, the NZ is defined in many ways, but an 
            % effective method is the Zero Load method, especially for larger 
            % specimens such as humans. Essentially, the section on the hysteresis 
            % curve bounded by applied loads of 0Nm for positive loading and negative
            % loading is defined as the NZ. 
    
            plot(last_cycleLOAD(last_cycleANGLE == pos_NZB), ...
                last_cycleANGLE(last_cycleANGLE == pos_NZB),'r*');
            plot(last_cycleLOAD(last_cycleANGLE == neg_NZB), ...
                last_cycleANGLE(last_cycleANGLE == neg_NZB),'r*');
            yline(neg_NZB,'--','-NZ');
            yline(pos_NZB,'--','+NZ');

            saveas(curfig,fullfile(folderdir,newFolderName, ...
                join([fgraph_name '.jpg'])),'jpg'); 
            saveas(curfig,fullfile(folderdir,newFolderName, ...
                join([fgraph_name '.fig'])),'fig'); 
    
            if ~isempty(NZlocs) && response2 == 0
                try 
                    if isempty(loadingNZ) || isempty(unloadingNZ) 
                        % Separate NZ; one for positive loading and one for
                        % negative loading
   
                          loadingNZ = find(NZlocs < find(last_cycleANGLE == pos_NZB)); 
                          %unloadingNZ = find(NZlocs > find(last_cycleANGLE == pos_NZB)); 
                          unloadingNZ = NZlocs(loadingNZ(end)+1:end);
                          unloadingNZ = unloadingNZ(1:find(last_cycleLOAD( ...
                              unloadingNZ) == min(last_cycleLOAD(unloadingNZ))));
                          
                    end 

                catch
                    saveas(curfig,fullfile(folderdir,newFolderName, ...
                        join([fgraph_name '.jpg'])),'jpg'); 
                    saveas(curfig,fullfile(folderdir,newFolderName, ...
                        join([fgraph_name '.fig'])),'fig'); 
                    
                    clf(curfig)
                    close(curfig)
                end

                try
                    % --- Plot red lines on loading and unloading NZ
                    plot(last_cycleLOAD(loadingNZ),...
                        last_cycleANGLE(loadingNZ),'r','LineWidth',1.5);
                    plot(last_cycleLOAD(unloadingNZ),...
                        last_cycleANGLE(unloadingNZ),'r','LineWidth',1.5);

                     % --- calcualte area of hysteresis curve
                    hysteresis = polyarea(last_cycleLOAD,last_cycleANGLE);
                    %fill(last_cycleLOAD,last_cycleANGLE,'yellow');
                    fill(last_cycleLOAD,last_cycleANGLE,[1.000,0.9765,0.8118])

                    % --- Determining NZS --- %
                    posNZS_bf = polyfit(last_cycleLOAD(loadingNZ), ...
                        last_cycleANGLE(loadingNZ),1);
                    plot(last_cycleLOAD,posNZS_bf( ...
                        1)*last_cycleLOAD+posNZS_bf(2),'black');
                    negNZS_bf = polyfit(last_cycleLOAD(unloadingNZ), ...
                        last_cycleANGLE(unloadingNZ),1);
                    plot(last_cycleLOAD,negNZS_bf( ...
                        1)*last_cycleLOAD+negNZS_bf(2),'black');
                catch
                    saveas(curfig,fullfile(folderdir,newFolderName, ...
                        join([fgraph_name '.jpg'])),'jpg'); 
                    saveas(curfig,fullfile(folderdir,newFolderName, ...
                        join([fgraph_name '.fig'])),'fig'); 

                    summarydata(grand_index2(col,row),:,j) = [max( ...
                        last_cycleANGLE),min(last_cycleANGLE),...
                    max(last_cycleANGLE)-min(last_cycleANGLE),...
                    abs(pos_NZB-neg_NZB), 0, 0, 0, 0, 0, 0, 0];

                    clf(curfig)
                    close(curfig)
                end
               
    % Jalen: I do not know of any literature or some sort of consesus on how
    % to accurately or repeatedly determine the EZ. As such, this is the most
    % noteworthy limitation when it comes to using this script for data 
    % analysis. Currently, the kmeans function is utilized to section off the
    % EZ in the data (see code below). I'll be honest, this is close to being
    % an educated guess and should be improved if possible.
    
            % --- Plotting for EZ --- %
               try  
                    %[~,load_maxloc] = max(last_cycleANGLE);
                    [~,load_maxloc] = max(last_cycleLOAD);
                    
                    % Find the start of the postive EZ as where the data 
                    % in the last cycle, outside of the NZ, starts after 
                    % the max load 
                    if find(pos_NZB==last_cycleANGLE) > load_maxloc
                        % Select the range of data between the max load
                        % and the NZ as the positive EZ range
%                         pos_EZ_range = load_maxloc:find( ...
%                             last_cycleANGLE>pos_NZB, 1, 'last');
                        pos_EZ_range = find(last_cycleANGLE>pos_NZB, ...
                        1, 'first'):load_maxloc;
%                         group_select_pos = 1;
                        group_select_pos = 3;

                    else  
                        pos_EZ_range = find( ...
                            pos_NZB==last_cycleANGLE):load_maxloc;
                        group_select_pos = length(pos_EZ_range);
                    end 
                
                    % We then use K-means to separate the positive EZ into
                    % 3 distinct regions (ideally of different stiffnesses)
                    % K-mean separates the data (load & angle) by squared 
                    % euclidean distance (closest to respective averages 
                    % of the separated data) *Great spot for future edits*
              
                    % pos_EZgroups = kmeans([last_cycleLOAD( ...
                    %     pos_EZ_range)';last_cycleLOAD(pos_EZ_range)'],3);
                        % NEW CODE HERE (3/9/2023)
%                      pos_EZgroups = kmeans(last_cycleLOAD( ...
%                      pos_EZgroups = kmeans(last_cycleANGLE( ...
%                          pos_EZ_range),3,'Start',[min(last_cycleLOAD( ...
%                          pos_EZ_range)); (min(last_cycleLOAD( ...
%                          pos_EZ_range)) + max(last_cycleLOAD( ...
%                          pos_EZ_range)))/2; max(last_cycleLOAD( ...
%                          pos_EZ_range))]);

                     pos=[last_cycleLOAD(pos_EZ_range) last_cycleANGLE(pos_EZ_range)];
     
                        %replacing kmeans

                    P=idivide(int16(length(pos)),3,"round"); %dividing the positive group into 3 sections
                    Pos_group1=pos(1:P,:);
                    Pos_group2=pos(P+1:P*2,:);
                    Pos_group3=pos((P*2)+1:end,:);
                   
                    
                    pos_slopes = zeros(1,length(pos(:,1))-2);

                    for i = 1:length(pos(:,1))-2
                        pos_slopes(i) =(pos(i+2,2)-pos(i,2)) / (pos(i+2,1)-pos(i,1));
                    end
                    pos_EZS_bf= max(pos_slopes);


                    % Plot the positive EZ on top of the hysteresis curve
%                     plot(last_cycleLOAD( ...
%                         pos_EZ_range(pos_EZgroups==pos_EZgroups( ...
%                         group_select_pos))), ...
%                         last_cycleANGLE(pos_EZ_range( ...
%                         pos_EZgroups==pos_EZgroups(group_select_pos))),...
%                         'm','LineWidth',1.5);

                     plot(Pos_group3(:,1),Pos_group3(:,2),...
                        'm','LineWidth',1.5);
    
                    %[~,load_minloc] = min(last_cycleANGLE);
                    [~,load_minloc] = min(last_cycleLOAD);
                    neg_EZ_range_idx = find(last_cycleANGLE<neg_NZB);
                   
                    % As long as the NZ boundary is after the peak height,
                    % the negative EZ start where the data in the last 
                    % cycle comes after the min load
                    if find(pos_NZB==last_cycleANGLE) > load_maxloc
%                         neg_EZ_range_idx = neg_EZ_range_idx( ...
%                             neg_EZ_range_idx < load_minloc);
%                         neg_EZ_range = max(neg_EZ_range_idx):load_minloc;
%                         group_select_neg = 1;
                        neg_EZ_range = find(last_cycleANGLE<neg_NZB, ...
                        1, 'first'):load_minloc;
                        if last_cycleANGLE(1) < neg_NZB
                            neg_EZ_range = neg_EZ_range(find( ...
                                last_cycleANGLE==pos_NZB):end);
                        end 
                        group_select_neg = 3;
                    else % If, for whatever reason, the positive NZ 
                        % boundary comes before the max load, find the 
                        % start of the negative EZ as the data that comes 
                        % before the min load
%                         neg_EZ_range_idx = neg_EZ_range_idx( ...
%                             neg_EZ_range_idx < load_minloc);
%                         neg_EZ_range = min(neg_EZ_range_idx):load_minloc;
                        neg_EZ_range = find(last_cycleANGLE<neg_NZB, ...
                        1, 'first'):load_minloc;
                        if last_cycleANGLE(1) < neg_NZB
                            neg_EZ_range = neg_EZ_range(find( ...
                                last_cycleANGLE==neg_NZB):end);
                        end 
                        group_select_neg = length(neg_EZ_range);
                    end 
    
                    % neg_EZgroups = kmeans([last_cycleLOAD( ...
                    %     neg_EZ_range)';last_cycleLOAD(neg_EZ_range)'],3);
%                     neg_EZgroups = kmeans(last_cycleLOAD( ...
%                     neg_EZgroups = kmeans(last_cycleANGLE( ...
%                         neg_EZ_range),3,'Start',[min(last_cycleLOAD( ...
%                          neg_EZ_range)); (min(last_cycleLOAD( ...
%                          neg_EZ_range)) + max(last_cycleLOAD( ...
%                          neg_EZ_range)))/2; max(last_cycleLOAD( ...
%                          neg_EZ_range))]);

                    neg=[last_cycleLOAD(neg_EZ_range) last_cycleANGLE(neg_EZ_range)];

                    N=idivide(int16(length(neg)),3,"round"); %dividing the positive group into 3 sections
                    Neg_group1=neg(1:N,:);
                    Neg_group2=neg(N+1:N*2,:);
                    Neg_group3=neg((N*2)+1:end,:);

                    for i = 1:length(neg(:,1))-2
                        neg_slopes(i) =(neg(i+2,2)-neg(i,2)) / (neg(i+2,1)-neg(i,1));
                    end
                    neg_EZS_bf= max(neg_slopes);



%                     plot(last_cycleLOAD(neg_EZ_range( ...
%                         neg_EZgroups==neg_EZgroups(group_select_neg))), ...
%                         last_cycleANGLE(neg_EZ_range( ...
%                         neg_EZgroups==neg_EZgroups(group_select_neg))),...
%                         'm','LineWidth',1.5);
%                     hold on

                    plot(Neg_group3(:,1),Neg_group3(:,2),...
                        'm','LineWidth',1.5);
                    hold on
                    
               catch
                   hold off 

                   saveas(curfig,fullfile(folderdir,newFolderName, ...
                        join([fgraph_name '.jpg'])),'jpg'); 
                   saveas(curfig,fullfile(folderdir,newFolderName, ...
                        join([fgraph_name '.fig'])),'fig'); 
                    
                    summarydata(grand_index2(col,row),:,j) = [max( ...
                        last_cycleANGLE),min(last_cycleANGLE),...
                    max(last_cycleANGLE)-min(last_cycleANGLE),...
                    abs(pos_NZB-neg_NZB), 1/posNZS_bf(1),1/negNZS_bf(1),...
                    0, 0, 0, 0, hysteresis];

                    clf(curfig);
                    close(curfig);
               end 

                    % --- Determining EZS --- %
               try
%                     posEZS_bf = polyfit(last_cycleLOAD(pos_EZ_range( ...
%                         pos_EZgroups==pos_EZgroups(group_select_pos))), ...
%                     last_cycleANGLE(pos_EZ_range( ...
%                     pos_EZgroups==pos_EZgroups(group_select_pos))),1);
%                     plot(last_cycleLOAD,posEZS_bf( ...
%                         1)*last_cycleLOAD+posEZS_bf(2),'g');

                    posEZS_bf = polyfit(Pos_group3(:,1),Pos_group3(:,2),1);
                    
                    plot(last_cycleLOAD,posEZS_bf( ...
                        1)*last_cycleLOAD+posEZS_bf(2),'g');

                
%                     negEZS_bf = polyfit(last_cycleLOAD( ...
%                         neg_EZ_range(neg_EZgroups==neg_EZgroups( ...
%                         group_select_neg))), ...
%                         last_cycleANGLE(neg_EZ_range( ...
%                         neg_EZgroups==neg_EZgroups(group_select_neg))),1);
%                     plot(last_cycleLOAD,negEZS_bf( ...
%                         1)*last_cycleLOAD+negEZS_bf(2),'g');

                    negEZS_bf = polyfit(Neg_group3(:,1),Neg_group3(:,2),1);
                    
                    plot(last_cycleLOAD,negEZS_bf( ...
                        1)*last_cycleLOAD+negEZS_bf(2),'g');
    
                    ylim([min(last_cycleANGLE)-0.5*(max( ...
                        last_cycleANGLE)-min(last_cycleANGLE)) ...
                    max(last_cycleANGLE)+0.5*(max(last_cycleANGLE)-min( ...
                    last_cycleANGLE))])

                    % --- Saving Figures --- %
                    saveas(curfig,fullfile(folderdir,newFolderName, ...
                        join([fgraph_name '.jpg'])),'jpg'); 
                    saveas(curfig,fullfile(folderdir,newFolderName, ...
                        join([fgraph_name '.fig'])),'fig'); 
                    hold off
    
                    summarydata(grand_index2(col,row),:,j) = [max( ...
                        last_cycleANGLE),min(last_cycleANGLE),...
                    max(last_cycleANGLE)-min(last_cycleANGLE),...
                    abs(pos_NZB-neg_NZB),...
                    1/posNZS_bf(1),1/negNZS_bf(1),...
                    abs(last_cycleANGLE(pos_EZ_range(1))-last_cycleANGLE( ...
                    pos_EZ_range(end))), abs(last_cycleANGLE( ...
                    neg_EZ_range(1))-last_cycleANGLE(neg_EZ_range(end)))...
                    ,abs(1/posEZS_bf(1)),abs(1/negEZS_bf(1)), hysteresis];

                    clf(curfig);
                    close(curfig);
                       
               catch         
%                    set(curfig,'visible','on'); 
                   figure(curfig);
                   plot(last_cycleLOAD,last_cycleANGLE);
                    title(fgraph_name);
                    xlabel('Load (Nm)');
                    ylabel('Relative Angle (°)');
                    hold on
                    plot(last_cycleLOAD(last_cycleANGLE == pos_NZB), ...
                        last_cycleANGLE(last_cycleANGLE == pos_NZB),'r*'); 
                    plot(last_cycleLOAD(last_cycleANGLE == neg_NZB), ...
                        last_cycleANGLE(last_cycleANGLE == neg_NZB),'r*');
                    yline(neg_NZB,'--','-NZ');
                    yline(pos_NZB,'--','+NZ');
                    plot(last_cycleLOAD(unloadingNZ),...
                        last_cycleANGLE(unloadingNZ),'r','LineWidth',1.5);
                    plot(last_cycleLOAD(loadingNZ),...
                        last_cycleANGLE(loadingNZ),'r','LineWidth',1.5);
                    
                     ylim([min(last_cycleANGLE)-0.5*(max( ...
                        last_cycleANGLE)-min(last_cycleANGLE)) ...
                    max(last_cycleANGLE)+0.5*(max(last_cycleANGLE)-min( ...
                    last_cycleANGLE))])

                    saveas(curfig,fullfile(folderdir,newFolderName, ...
                        join([fgraph_name '.jpg'])),'jpg');
                    saveas(curfig,fullfile(folderdir,newFolderName, ...
                        join([fgraph_name '.fig'])),'fig');
                    hold off
    
                    summarydata(grand_index2(col,row),:,j) = [max( ...
                        last_cycleANGLE),min(last_cycleANGLE),...
                    max(last_cycleANGLE)-min(last_cycleANGLE),...
                    abs(pos_NZB-neg_NZB), 1/posNZS_bf(1),1/negNZS_bf(1),...
                    abs(last_cycleANGLE(pos_EZ_range(1))-last_cycleANGLE( ...
                    pos_EZ_range(end))), abs(last_cycleANGLE( ...
                    neg_EZ_range(1))-last_cycleANGLE(neg_EZ_range(end))),...
                    0, 0, hysteresis];

                    clf(curfig);
                    close(curfig);
                end
            end 

            if ishandle(curfig)
                close(curfig)
            elseif ishandle(rawdata_fig)
                close(rawdata_fig)
            end
        
        end
        
        % Update progress bar
        waitbar(q/length(indx), f, sprintf('Progress: %d %%', floor(q/ ...
            length(indx)*100)));
    end

    close(f)
    
%% The summary data
% The summary data is presented with each relative separated onto a
% different sheet. Any file skipped because of mismatched pairing or 
% containg 'messy' displacement data will produce 0 for each measurement.
% Any files skipped, in the absence of the above errors, will not appear
% in the summary. 

% ***Note: if the script is run again with the same data the excel sheet
% is not completely overwritten, only the rows in which changes remained.
% To produce a fresh, new excel sheet, delete the old one or change the
% name of the file.

% CREATING SUMMARY DATA SHEETS
    groups_FE = cellfun(@(x) join([x, ' (FE)']),groups,'UniformOutput',false);
    groups_LB = cellfun(@(x) join([x, ' (LB)']),groups,'UniformOutput',false);
    groups_AR = cellfun(@(x) join([x, ' (AR)']),groups,'UniformOutput',false);
    grand_groups = [groups_FE;groups_LB;groups_AR];
    var = {'MaxAngle (°)','MinAngle (°)','ROM (°)','NZmagnitude (°)', ...
        'posNZS (Nm/°)','negNZS (Nm/°)','posEZmagnitude (°)', ...
        'negEZmagnitude (°)','posEZS (Nm/°)','negEZS (Nm/°)', 'Hysteresis Area (°Nm)'};
    % Paper describing use of area under hysteresis curve: Kelbl M, 
    % Kocis J, Vesely R, Florian Z, Navrat T, Vosynek P. Biomechanical 
    % Testing of Spinal Segment Fixed by Arcofix System on the Swine Spine.
    % Asian Spine J. 2015 Aug;9(4):503-10. doi: 10.4184/asj.2015.9.4.503. 
    % Epub 2015 Jul 28. PMID: 26240706; PMCID: PMC4522437.

    summarysheets = num2cell(summarydata);
    xlsfullfile = fullfile(folderdir2,'Specimen Summary - Hystereis Curve.xls');
    for G = 1:size(summarysheets,3)
        T = cell2table(summarysheets(grand_index3(indx),:,G),...
            'RowNames',grand_groups(grand_index3(indx)),...
            'VariableNames',var);
        writetable(T,xlsfullfile,'Sheet',relativeNames{G},'WriteRowNames',true);
    end

    % ENDING CODE
    fig2 = uifigure;
    selection3 = uiconfirm(fig2,...
        {['The summary Excel file and subfolder containing the generated ' ...
        'plots have been made and placed in the folder containing the load' ...
        ' data.'],' ', 'Would you like to run the script again?'},...
        'Success!','Icon','success','Options',{'Done','Restart'});
    switch selection3
        case 'Done'
            close(fig2)
            z = 0;
            return
        case 'Restart'
            clear indx
            close(fig2)
            z = 1;
    end
    end 

catch
    warndlg(['Program closed unexpectedly. Please make sure data is ' ...
        'in the proper format and the correct folders are selected.' ...
        'Otherwise continue as usual.'], 'Warning: Program Ran Into Error');
    try
        close(all)
    catch
    end
end

%% Limiations
% This script is made for use with the Optotrak and 6DOF used by the 
% biomechanics team at Globus, files should be organized as such: 
% FE,LB, AR, repeat... The way MATLAB organizes files is not a typical
% alph-numeric ordering and can affect ordering of load data files, 
% for instance, Book1 should be changed to Book01 before running script
% since Book10 will be placed after Book1 (e.g., Book1,Book10,Book11...
% Book2,Book20,Book21... instead it should be Book01,Book02,Book03...
% Book10,Book11,Book12...). The EZ and its subsequent calculations do not
% use the most reliable methods (see Plotting for EZ in code). It is easy
% to see that if you capture the magnitude of the NZ and both the 
% magnitudes of the positive EZ and negative EZ, the sum is the range of 
% motion (ROM); however, the difference in sampling frequency between the
% Optotrak and 6DOF means that some Optotrak points will be removed. 
% As such, this should not be used to replace typical ROM analysis since 
% definite peaks can be missed. On the other hand it can be used as a good
% estimate (calculations with this script are typically 1-1.3 degrees less
% than the actual ROM as measured by the Optotrak).
