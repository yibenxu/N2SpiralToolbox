
%% isolate N2 epoch data

main_folder = pwd; % define the main_folder for the toolbox 
cd ([main_folder,'/eeglab2024.0'])
addpath(genpath(pwd))

for baseline_or_treatment = 1:2
for sub_id = 1:9
    
    clearvars -except sub_id baseline_or_treatment
    disp(['sub ',num2str(sub_id),' baseline_or_treatment ',num2str(baseline_or_treatment)])
subject_list = [{'01SC'} {'03RA'} {'06AC'} {'08RG'} {'09AS'} {'11NR'} {'12SM'} {'13PR'} {'15DC'} ];
subject = cell2mat(subject_list(sub_id));

if baseline_or_treatment == 1
    foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/',subject,'/Baseline/'];
elseif baseline_or_treatment == 2
    foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/',subject,'/Treatment/'];
end

cd(foldername)
name = dir(foldername);

% find target filename
target_filename = ['bc1pass_FFT.set']; 
for i2 = 1:size(name,1)
filename_list = name(i2).name;
    TF = contains(filename_list, target_filename);
    if TF == 1
       filename = filename_list;
       break
    end
end
if TF == 0 
    disp(['input data doesnt exist, sub ',num2str(sub_id),' baseline_or_treatment ',num2str(baseline_or_treatment)])
    continue
end

%% load full-night HdEEG data
EEG = pop_loadset(filename,foldername);

%% extract channel locations

[labels,x,y] = topoplot([],EEG.chanlocs,'style','blank','electrodes','labelpoint','chaninfo',EEG.chaninfo);
channel_locations_2d = [ x(:) y(:)];

% rotate the channel configeration by 90 degrees to make it face up
% Create rotation matrix
theta = 90; % to rotate 90 counterclockwise
R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
% Rotate your point(s)
channel_locations_2d_rotate = [];
for ichan = 1:size(channel_locations_2d,1)
point = channel_locations_2d(ichan,:)'; % arbitrarily selected
rotpoint = R*point;
channel_locations_2d_rotate(ichan,:) = rotpoint;
end

channel_locations_3d = [EEG.chanlocs(1:end).X ; EEG.chanlocs(1:end).Y ; EEG.chanlocs(1:end).Z]';
bad_channel = EEG.bad_channels;

% remove bad channels ******************
bad_channel_accu = bad_channel{1};
labels_num = zeros(size(labels,1),1); % original channel label
for i2 = 1:size(labels,1)
temp1 = labels(i2,:);
labels_num(i2) = str2double(temp1);
end
channel_locations_3d_adj = channel_locations_3d; % 3d coordinate after removing bad channels
channel_locations_2d_rotate_adj = channel_locations_2d_rotate;
labels_num_adj = labels_num; % labels after removing bad channels
[C,ia,ib] = intersect(labels_num,bad_channel_accu); % only remove bad channels that still exist
bad_channel_accu_valid = [];
if nansum(bad_channel_accu(:))~=0 && nansum(C(:))~=0
bad_channel_accu_valid = bad_channel_accu(ib);    
channel_locations_3d_adj(ia,:) = [];
channel_locations_2d_rotate_adj(ia,:) = [];
labels_num_adj(ia,:) = [];
end
% filter channels by their 3D height, anything below the ear is removed
channel_locations_3d_Zfilt = channel_locations_3d_adj;
channel_id_invalid_idx = find(channel_locations_3d_Zfilt(:,3)<-2);
channel_id_invalid = labels_num_adj(channel_id_invalid_idx);
channel_id_valid_idx = find(channel_locations_3d_Zfilt(:,3)>=-2); % valid channel id with z-axis>-2
channel_id_valid = labels_num_adj(channel_id_valid_idx);
[C2, ia2, working_channel_id_valid_idx] = intersect(channel_id_valid,labels_num); % ensure remaining channels are both working and valid

% 2D channel map from EEGLAB **************************
channel_locations_2d_rotate_valid = channel_locations_2d_rotate(working_channel_id_valid_idx,:);
upscale = 30; % 30=20x21;
channel_locations_2d_rotate_valid_rescale = channel_locations_2d_rotate_valid.*upscale;

min_x = nanmin(channel_locations_2d_rotate_valid_rescale(:,1));
min_y = nanmin(channel_locations_2d_rotate_valid_rescale(:,2));

channel_locations_2d_rotate_valid_rescale_adj = zeros(size(channel_locations_2d_rotate_valid_rescale));
channel_locations_2d_rotate_valid_rescale_adj(:,1) = channel_locations_2d_rotate_valid_rescale(:,1) - min_x + 1; %0.5
channel_locations_2d_rotate_valid_rescale_adj(:,2) = channel_locations_2d_rotate_valid_rescale(:,2) - min_y + 1; %0.5

channel_locations_2d_rotate_valid_rescale_adj_int = round(channel_locations_2d_rotate_valid_rescale_adj);
row_max = nanmax(channel_locations_2d_rotate_valid_rescale_adj_int(:,2));
col_max = nanmax(channel_locations_2d_rotate_valid_rescale_adj_int(:,1));

figure(2)
scatter(channel_locations_2d_rotate_valid_rescale_adj(:,1),channel_locations_2d_rotate_valid_rescale_adj(:,2))
xlim([1,col_max])
ylim([1,row_max])

channel_locations_2d_rotate = channel_locations_2d_rotate_valid_rescale_adj;
% 2D channel map from EEGLAB **************************

    if baseline_or_treatment == 1
        foldername = ['/import/taiji1/yixu4976/HD_EEG/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Baseline/'];
        filename = ['Channel_Location_2D_EEGLAB_20x21_',subject,'_baseline.mat'];
    elseif baseline_or_treatment == 2
        foldername = ['/import/taiji1/yixu4976/HD_EEG/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Treatment/'];
        filename = ['Channel_Location_2D_EEGLAB_20x21_',subject,'_treatment.mat'];
    end
    
    save([foldername,filename],'working_channel_id_valid_idx','channel_locations_2d_rotate','channel_locations_3d','labels','bad_channel','upscale','-v7.3');

%% isolate N2 stage data
stage_name_list{1} = {'Wake'};
stage_name_list{2} = {'N1'};
stage_name_list{3} = {'N2'};
stage_name_list{4} = {'N3'};
stage_name_list{5} = {'REM'};

EEG_stagenames = [];
EEG_stagenames = cellfun(@char, EEG.stagenames, 'UniformOutput', false); % convert to cell array with string

EEG_stageid = zeros(size(EEG.stagenames,2),1);
for stage_id = 1:5
    stage_name = stage_name_list{stage_id};
    stage_name = cell2mat(stage_name);
    index = cellfun(@(x) isequal(x,stage_name), EEG_stagenames); % index is in logic values
    EEG_stageid(index) = stage_id;
end

%% isolate N2 epoch data

EEG_stageid_N2 = EEG_stageid;
EEG_stageid_N2(EEG_stageid_N2~=3) = 0; % 3=N2;

% identify EEG epochs
EEG_stageid_N2_t_dif = EEG_stageid_N2(2:end)-EEG_stageid_N2(1:end-1);
epoch_boundary = find(EEG_stageid_N2_t_dif);
epoch_start = epoch_boundary(1:2:end);
epoch_end = epoch_boundary(2:2:end);
% 
% EEG signal data
data_N2_epochs = [];
for iepoch = 1:size(epoch_start,1)
    if epoch_start(iepoch)>size(EEG.data,2) || epoch_end(iepoch)>size(EEG.data,2) 
        continue
    end
    data_1epoch = EEG.data(:,epoch_start(iepoch):epoch_end(iepoch));
    data_N2_epochs{iepoch} = data_1epoch;
end
channel_locations = EEG.chanlocs;
bad_channel = EEG.bad_channels;
x_arousal = EEG.x_arousal;

% arousal data 
arousal_N2_epochs = [];
for iepoch = 1:size(epoch_start,1)
    if epoch_start(iepoch)>size(EEG.data,2) || epoch_end(iepoch)>size(EEG.data,2) 
        continue
    end 
    x_arousal_1epoch = x_arousal(epoch_start(iepoch):epoch_end(iepoch));
    arousal_N2_epochs{iepoch} = x_arousal_1epoch;
end
channel_csc_montage = EEG.csc_montage;

% save N2 epoch data
if baseline_or_treatment == 1
    foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Baseline/'];
    filename = ['N2_epoch_data_',subject,'_baseline.mat'];
elseif baseline_or_treatment == 2
    foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Treatment/'];
    filename = ['N2_epoch_data_',subject,'_treatment.mat'];
end
save([foldername,filename],'arousal_N2_epochs','data_N2_epochs','epoch_start','epoch_end','channel_locations','bad_channel','x_arousal','channel_csc_montage','-v7.3');
% save arousal data
if baseline_or_treatment == 1
    foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Baseline/'];
    filename = ['N2_arousal_data_',subject,'_baseline.mat'];
elseif baseline_or_treatment == 2
    foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Treatment/'];
    filename = ['N2_arousal_data_',subject,'_treatment.mat'];
end
save([foldername,filename],'arousal_N2_epochs','epoch_start','epoch_end','channel_locations','bad_channel','x_arousal','channel_csc_montage','-v7.3');

end
end

%% temporal bandpass filter

for baseline_or_treatment = 1:2
for sub_id = 1:9
    clearvars -except sub_id baseline_or_treatment
    subject_list = [{'01SC'} {'03RA'} {'06AC'} {'08RG'} {'09AS'} {'11NR'} {'12SM'} {'13PR'} {'15DC'} ];
    subject = cell2mat(subject_list(sub_id));
    disp(['sub ',num2str(sub_id),' baseline_or_treatment ',num2str(baseline_or_treatment)])

    if baseline_or_treatment == 1
        foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Baseline/'];
        filename = ['N2_epoch_data_',subject,'_baseline.mat'];
    elseif baseline_or_treatment == 2
        foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Treatment/'];
        filename = ['N2_epoch_data_',subject,'_treatment.mat'];
    end
    cd(foldername)

    if exist(filename)
        load(filename)
    else
        disp(['input data doesnt exist, sub ',num2str(sub_id),' baseline_or_treatment ',num2str(baseline_or_treatment)])
        continue
    end

    data_bandpa_epochs = [];
    for iepoch = 1:size(data_N2_epochs,2)
        dataIn = double(data_N2_epochs{iepoch});
        
        % % get rid of electric lines
        sampleFs = 0.002;
        fLowNorm = 49/(1/sampleFs)*2 ;
        fHighNorm = 51/(1/sampleFs)*2 ;
        filterOrder = 2 ;
        dataIn_50HzNotchFilt = zeros(size(dataIn)) ;
        % temporal stop filter for each channel to remove electric lines
        tic
        for ichan = 1:size(dataIn,1)
                data_1chan = dataIn(ichan,:);
                if nansum(data_1chan(:)) ~=0   
                   [b,a] = butter(filterOrder,[fLowNorm fHighNorm],'stop');
                   dataIn_50HzNotchFilt(ichan,:) = filter(b,a,data_1chan(:)) ;
                end
        end
        % temporal filter 11-15Hz
        dataIn = dataIn_50HzNotchFilt;
        flow = [11]; % first cutoff
        fhigh = [15]; % second cutoff
        N = 2;  %filter order
        Fs = 500;
        h = fdesign.bandpass('N,F3dB1,F3dB2',N,flow,fhigh,Fs); % design a filter specification object
        Hd = design(h,'butter'); % apply design method to filter specification object
        set(Hd,'arithmetic','double'); %??

        SOS = Hd.sosMatrix; % store filter coefficients in SOS as in [ sections x coefficients] (4x6) matrix
        G = Hd.ScaleValues; % store filter scale values in sections (4x1) matrix
        data_bandpa = nan(size(dataIn));  

        for ichan = 1:size(dataIn,1)
            data_1chan = dataIn(ichan,:);
            if nansum(data_1chan(:))~=0
                data_bandpa(ichan,:) = filtfilt(SOS,G,data_1chan(:)); % apply zero-phase bandpass filter (in both forword and backward directions to avoid phase distortion) 
            end
        end
        data_bandpa_epochs{iepoch} = data_bandpa;
    end
    % save N2 bandpass filtered data
    if baseline_or_treatment == 1
        foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Baseline/'];
        filename = ['N2_epoch_data_50HzNotchFilt_bandpass11To15Hz_',subject,'_baseline.mat'];
    elseif baseline_or_treatment == 2
        foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Treatment/'];
        filename = ['N2_epoch_data_50HzNotchFilt_bandpass11To15Hz_',subject,'_treatment.mat'];
    end
    save([foldername,filename],'data_bandpa_epochs','channel_locations','bad_channel','x_arousal','channel_csc_montage','-v7.3');
    disp(['data saved, sub ',num2str(sub_id),' baseline_or_treatment ',num2str(baseline_or_treatment)])

end

            
end

%% spindle epoch detection

for baseline_or_treatment = 1:2
for sub_id = 1:9
    disp(['baseline or treatment = ',num2str(baseline_or_treatment),' sub = ',num2str(sub_id)])
    clearvars -except sub_id baseline_or_treatment
    subject_list = [{'01SC'} {'03RA'} {'06AC'} {'08RG'} {'09AS'} {'11NR'} {'12SM'} {'13PR'} {'15DC'} ];
    subject = cell2mat(subject_list(sub_id));
    if baseline_or_treatment == 1
    foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Baseline/'];
    filename = ['N2_epoch_data_50HzNotchFilt_bandpass11To15Hz_',subject,'_baseline.mat'];
    elseif baseline_or_treatment == 2
    foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Treatment/'];
    filename = ['N2_epoch_data_50HzNotchFilt_bandpass11To15Hz_',subject,'_treatment.mat'];
    end
    cd(foldername)
    if exist(filename)
        load(filename)
    else
        continue
    end
        
    % first temporally smooth the SIGMA signel with root mean square (RMS) using a moving window of 200ms
    window_size = 100;
    spindle_template_epochs = [];
    for iepoch = 1:size(data_bandpa_epochs,2)
        data_bandpa_1epoch = data_bandpa_epochs{iepoch};
        epoch_duration = size(data_bandpa_1epoch,2);
        channel_count = size(data_bandpa_1epoch,1);
        data_bandpa_1epoch_RMS = zeros(size(data_bandpa_1epoch));
        for t = 1+window_size/2:epoch_duration-window_size/2
            temp1 = data_bandpa_1epoch(:,t-window_size/2:t+window_size/2); % select a window
            data_bandpa_1epoch_RMS(1:size(temp1,1),t) = sqrt(nanmean((temp1.^2),2));
        end
%         data_bandpa_1epoch_RMS = movmean(data_bandpa_1epoch_RMS,100,2);

        % find the mean std of all channels
        data_bandpa_1epoch_std = [];
        for ichan = 1:size(data_bandpa_1epoch,1)
            temp1 = data_bandpa_1epoch(ichan,:);
            data_bandpa_1epoch_std(ichan,:) = nanstd(temp1(:));
        end
        data_bandpa_1epoch_stdavg = nanmean(data_bandpa_1epoch_std(:));
    
        % use amplitude threshold to filter the RMS signal (> 1.5SD)
        std_threshold = 1.5; % 1.5SD+
        data_bandpa_1epoch_RMS_stdnorm = data_bandpa_1epoch_RMS./data_bandpa_1epoch_stdavg;
        data_bandpa_1epoch_RMS_stdnorm_filt = data_bandpa_1epoch_RMS_stdnorm;
        data_bandpa_1epoch_RMS_stdnorm_filt(data_bandpa_1epoch_RMS_stdnorm_filt<=1.5) = 0;
        
        % use duration threshold to filter the RMS signal (>0.5s && < 3s)
        % via pattern_detection_v4, we convert the 1D time series to 3D space, by assuming the time (t) dimension as the
        % x-axis in a 2D space, artificially increase the y-axis to 250 by
        % dupicating the time (or x-axis), t-axis of the 3D space is also
        % artificially increased to 2 by duplicating x and y-axis
        duration_threshold_min = 250; % 500ms+
        duration_threshold_max = 1500; % 3000ms-
        spindle_template = zeros(channel_count,epoch_duration);
%         tic
        for ichan = 1:channel_count % we do this for each channel
                temp1 = data_bandpa_1epoch_RMS_stdnorm_filt(ichan,:);
                if nansum(temp1(:))==0
                   continue 
                end
                
                start_count = 0;
                end_count = 0;
                spindle_start = [];
                spindle_end = [];
                for t = 1:epoch_duration-1
                    temp1_t = temp1(t);
                    temp1_tplus1 = temp1(t+1);
                    if temp1_t*temp1_tplus1==0 && temp1_t==0 && temp1_tplus1~=0
                       start_count = start_count+1;
                       spindle_start(start_count) = t+1;
                    elseif temp1_t*temp1_tplus1==0 && temp1_t~=0 && temp1_tplus1==0
                       end_count = end_count+1;
                       spindle_end(end_count) = t+1;                       
                    end
                end
                if nansum(spindle_start)== 0 || nansum(spindle_end(:))==0
                    continue
                end
                if size(spindle_start(:),1) == size(spindle_end(:),1) && spindle_start(1)<spindle_end(1) % % start/end/start..../start/end
                   spindle_duration = spindle_end - spindle_start;
                   spindle_duration_start_end = [spindle_duration' spindle_start' spindle_end'];                    
                elseif size(spindle_start(:),1) == size(spindle_end(:),1) && spindle_start(1)>spindle_end(1) && spindle_start(end)>spindle_end(end) % end/start/end..../start/end/start
                   spindle_duration = spindle_end(2:end) - spindle_start(1:end-1);
                   spindle_duration_start_end = [spindle_duration' spindle_start(1:end-1)' spindle_end(2:end)'];                    
                elseif size(spindle_start(:),1) > size(spindle_end(:),1) && spindle_start(1)<spindle_end(1) % start/end/start..../start/end/start
                   spindle_duration = spindle_end - spindle_start(1:end-1);
                   spindle_duration_start_end = [spindle_duration' spindle_start(1:end-1)' spindle_end'];                    
                elseif size(spindle_start(:),1) < size(spindle_end(:),1) && spindle_start(1)<spindle_end(2) % end/start..../start/end/start/end
                   spindle_duration = spindle_end(2:end) - spindle_start;
                   spindle_duration_start_end = [spindle_duration' spindle_start' spindle_end(2:end)'];
                end
                channel_id_invalid = find(spindle_duration<250);
                spindle_duration_start_end(channel_id_invalid,:) = nan;
                channel_id_invalid = find(spindle_duration>1500);
                if isempty(channel_id_invalid)~=1
                    spindle_duration_start_end(channel_id_invalid,:) = nan;
                end
                temp2 = spindle_duration_start_end(:,1);
                channel_id_invalid = find(isnan(temp2));
                spindle_duration_start_end(channel_id_invalid,:) = [];
                
                if nansum(spindle_duration_start_end(:))==0 % skip loop if no spindle found
                    continue
                end
                
                for ipatt = 1:size(spindle_duration_start_end,1)
                    spindle_start = spindle_duration_start_end(ipatt,2);
                    spindle_end = spindle_duration_start_end(ipatt,3);
                    spindle_template(ichan,spindle_start:spindle_end) = data_bandpa_1epoch_RMS_stdnorm(ichan,spindle_start:spindle_end);
                end
                
        end
%         toc
        spindle_template_epochs{iepoch} = spindle_template;
    end
    % save spindle epochs via N2 bandpass filtered data
    if baseline_or_treatment == 1
    foldername = ['/import/silo3/yixu4976/HD EEG/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Baseline/'];
    filename = ['N2_epoch_data_bandpass11To15Hz_SpindleDetectedMolle2011_',subject,'_baseline.mat'];
    elseif baseline_or_treatment == 2
    foldername = ['/import/silo3/yixu4976/HD EEG/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Treatment/'];
    filename = ['N2_epoch_data_bandpass11To15Hz_SpindleDetectedMolle2011_',subject,'_treatment.mat'];
    end
    save([foldername,filename],'spindle_template_epochs','channel_locations','bad_channel','x_arousal','channel_csc_montage','-v7.3');

    
end

end

%% find start and end of spindle epochs
for baseline_or_treatment = 1:2

for sub_id = 1:9
    clearvars -except sub_id baseline_or_treatment
    subject_list = [{'01SC'} {'03RA'} {'06AC'} {'08RG'} {'09AS'} {'11NR'} {'12SM'} {'13PR'} {'15DC'} ];
    subject = cell2mat(subject_list(sub_id));

        % load spindle data
       if baseline_or_treatment == 1
        foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Baseline/After 50Hz Notch Filter/'];
        filename = ['N2_epoch_data_bandpass11To15Hz_50HzNotchFilt_SpindleDetectedMolle2011_',subject,'_baseline.mat'];
        elseif baseline_or_treatment == 2
        foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Treatment/After 50Hz Notch Filter'];
        filename = ['N2_epoch_data_bandpass11To15Hz_50HzNotchFilt_SpindleDetectedMolle2011_',subject,'_treatment.mat'];
        end
        cd(foldername)
        if exist(filename)
            load(filename)
        else
            continue
        end
        disp(['baseline or treatment = ',num2str(baseline_or_treatment),' sub = ',num2str(sub_id)])

        spindle_epochs_start_end_duration_iepoch_1sub = []; 
        for iepoch = 1:size(spindle_template_epochs_start_end,2)
            spindle_template_epochs_start_end_1N2epoch = spindle_template_epochs_start_end{iepoch};
            spindle_template_epochs_start_end_1N2epoch_duration = spindle_template_epochs_start_end_1N2epoch(:,3)-spindle_template_epochs_start_end_1N2epoch(:,2)+1;
            spindle_template_epochs_start_end_duration_iepoch = [spindle_template_epochs_start_end_1N2epoch(:,2:3) spindle_template_epochs_start_end_1N2epoch_duration iepoch.*ones(size(spindle_template_epochs_start_end_1N2epoch,1),1) ];
            spindle_epochs_start_end_duration_iepoch_1sub = [spindle_epochs_start_end_duration_iepoch_1sub;spindle_template_epochs_start_end_duration_iepoch ];
        end
    
    if isempty(spindle_epochs_start_end_duration_iepoch_1sub)==1
        continue
    end
    if baseline_or_treatment == 1
    foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Baseline/After 50Hz Notch Filter/'];
    filename = ['N2_spindle_epochs_accu_Molle2011_',subject,'_baseline.mat'];
    elseif baseline_or_treatment == 2
    foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Treatment/After 50Hz Notch Filter/'];
    filename = ['N2_spindle_epochs_accu_Molle2011_',subject,'_treatment.mat'];
    end
    save([foldername,filename],'spindle_epochs_start_end_duration_iepoch_1sub','-v7.3');
    disp(['file saved, baseline or treatment = ',num2str(baseline_or_treatment),' sub = ',num2str(sub_id)])

end
end

%% isolate data from spindle epochs

for baseline_or_treatment = 1:2

for sub_id = 1:9
    disp(['baseline or treatment = ',num2str(baseline_or_treatment),' sub = ',num2str(sub_id)])
    clearvars -except sub_id baseline_or_treatment
    subject_list = [{'01SC'} {'03RA'} {'06AC'} {'08RG'} {'09AS'} {'11NR'} {'12SM'} {'13PR'} {'15DC'} ];
    subject = cell2mat(subject_list(sub_id));

    % load spindle epoch time
    if baseline_or_treatment == 1
    foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Baseline/After 50Hz Notch Filter/'];
    filename = ['N2_spindle_epochs_accu_Molle2011_',subject,'_baseline.mat'];
    elseif baseline_or_treatment == 2
    foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Treatment/After 50Hz Notch Filter/'];
    filename = ['N2_spindle_epochs_accu_Molle2011_',subject,'_treatment.mat'];
    end
    cd(foldername)
    if exist(filename)
        load(filename)
    else
        disp(['input doesnt exist, baseline or treatment = ',num2str(baseline_or_treatment),' sub = ',num2str(sub_id)])    
        continue
    end
    
    % load arousal signal
    if baseline_or_treatment == 1
        foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Baseline/'];
        filename = ['N2_arousal_data_',subject,'_baseline.mat'];
    elseif baseline_or_treatment == 2
        foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Treatment/'];
        filename = ['N2_arousal_data_',subject,'_treatment.mat'];
    end
    cd(foldername)
    if exist(filename)
        load(filename)
        disp(['data loaded, baseline or treatment = ',num2str(baseline_or_treatment),' sub = ',num2str(sub_id)])
    else
        continue
    end 
    
    % load band-passed signal 
    if baseline_or_treatment == 1
    foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Baseline/After 50Hz Notch Filter/'];
    filename = ['N2_epoch_data_50HzNotchFilt_bandpass11To15Hz_',subject,'_baseline.mat'];
    elseif baseline_or_treatment == 2
    foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Treatment/After 50Hz Notch Filter/'];
    filename = ['N2_epoch_data_50HzNotchFilt_bandpass11To15Hz_',subject,'_treatment.mat'];
    end
    cd(foldername)
    if exist(filename)
        load(filename)
    else
        disp(['input doesnt exist, baseline or treatment = ',num2str(baseline_or_treatment),' sub = ',num2str(sub_id)])
        continue
    end 
   
    if baseline_or_treatment == 1
    foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Baseline/After 50Hz Notch Filter/bandpass data SpindleEpochExtendPlusMinus100/'];
    elseif baseline_or_treatment == 2
    foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Treatment/After 50Hz Notch Filter/bandpass data SpindleEpochExtendPlusMinus100/'];
    end
    cd(foldername)
    if exist(subject)==0
        mkdir(subject)
    end
        
    N2_data_1spindle_1subaccu = [];
    N2_data_phase_1spindle_1subaccu = [];
    for iepoch = 1:size(spindle_epochs_start_end_duration_iepoch_1sub,1)
        
        if baseline_or_treatment == 1
        foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Baseline/After 50Hz Notch Filter/bandpass data SpindleEpochExtendPlusMinus100/',subject,'/'];
        filename = ['N2_data_bandpass11To15Hz_SpindleEpochExtendPlusMinus100_NoArousal_subject',subject,'_SpindleEpoch',num2str(iepoch),'_baseline.mat'];
        elseif baseline_or_treatment == 2
        foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Treatment/After 50Hz Notch Filter/bandpass data SpindleEpochExtendPlusMinus100/',subject,'/'];
        filename = ['N2_data_bandpass11To15Hz_SpindleEpochExtendPlusMinus100_NoArousal_subject',subject,'_SpindleEpoch',num2str(iepoch),'_treatment.mat'];
        end  
        cd(foldername)
        corrupt = 0;
        if exist(filename)~=0
            try 
                load(filename)
            catch
                disp(['output corrupt, baseline or treatment = ',num2str(baseline_or_treatment),' sub = ',num2str(sub_id),' spindle=',num2str(iepoch)])
                corrupt = 1;
            end
        end
        if exist(filename)~=0 && corrupt == 0
            disp(['output already exists, baseline or treatment = ',num2str(baseline_or_treatment),' sub = ',num2str(sub_id),' spindle=',num2str(iepoch)]) 
            continue
        end  
        
        spindle_epochs_start_end_duration_iepoch_1spindle = spindle_epochs_start_end_duration_iepoch_1sub(iepoch,:);
        N2_epoch_id = spindle_epochs_start_end_duration_iepoch_1spindle(4);
        N2_epoch_t_start = spindle_epochs_start_end_duration_iepoch_1spindle(1);
        N2_epoch_t_end = spindle_epochs_start_end_duration_iepoch_1spindle(2);
        N2_data_1epoch = data_bandpa_epochs{N2_epoch_id};
        
        % filter by arousal
        arousal_N2_1epoch = arousal_N2_epochs{N2_epoch_id};
        idx_arousal = find(arousal_N2_1epoch==1);
        if nansum(idx_arousal(:))~=0
           N2_data_1epoch(:,idx_arousal) = NaN;
        end       
        % extended spindle epoch*****
        N2_epoch_t_start_ext = N2_epoch_t_start-100;
        N2_epoch_t_end_ext = N2_epoch_t_end+100;
        if N2_epoch_t_start_ext<1 
            continue;
        end
        if N2_epoch_t_end_ext>size(N2_data_1epoch,2)
            continue;
        end
        N2_data_1spindle = N2_data_1epoch(:,N2_epoch_t_start_ext:N2_epoch_t_end_ext);

        if baseline_or_treatment == 1
        foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Baseline/After 50Hz Notch Filter/bandpass data SpindleEpochExtendPlusMinus100/',subject,'/'];
        filename = ['N2_data_bandpass11To15Hz_SpindleEpochExtendPlusMinus100_NoArousal_subject',subject,'_SpindleEpoch',num2str(iepoch),'_baseline.mat'];
        elseif baseline_or_treatment == 2
        foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Treatment/After 50Hz Notch Filter/bandpass data SpindleEpochExtendPlusMinus100/',subject,'/'];
        filename = ['N2_data_bandpass11To15Hz_SpindleEpochExtendPlusMinus100_NoArousal_subject',subject,'_SpindleEpoch',num2str(iepoch),'_treatment.mat'];
        end  
        
        save([foldername,filename],'N2_data_1spindle','-v7.3');
        disp(['file saved, baseline or treatment = ',num2str(baseline_or_treatment),' sub = ',num2str(sub_id),' epoch=',num2str(iepoch)])

    end
end
end

%% phase velocity field over 20x21 2D grid

for baseline_or_treatment = 1:2
for sub_id = 1:9
    
subject_list = [{'01SC'} {'03RA'} {'06AC'} {'08RG'} {'09AS'} {'11NR'} {'12SM'} {'13PR'} {'15DC'} ];
subject = cell2mat(subject_list(sub_id));  
if baseline_or_treatment == 1
    foldername = [main_folder, '/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Baseline'];
    cd(foldername)
    load(['Channel_Location_2D_EEGLAB_20x21_',subject,'_baseline.mat'])
    load(['N2_spindle_epochs_accu_Molle2011_',subject,'_baseline.mat'])
elseif baseline_or_treatment == 2
    foldername = [main_folder, '/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Treatment'];    
    cd(foldername)
    load(['Channel_Location_2D_EEGLAB_20x21_',subject,'_treatment.mat'])
    load(['N2_spindle_epochs_accu_Molle2011_',subject,'_treatment.mat'])
end

for iepoch = 1:size(spindle_epochs_start_end_duration_iepoch_1sub,1)

if baseline_or_treatment == 1
foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Baseline/After 50Hz Notch Filter/bandpass data SpindleEpochExtendPlusMinus100/',subject,'/'];
filename = ['N2_data_bandpass11To15Hz_SpindleEpochExtendPlusMinus100_NoArousal_subject',subject,'_SpindleEpoch',num2str(iepoch),'_baseline.mat'];
elseif baseline_or_treatment == 2
foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Treatment/After 50Hz Notch Filter/bandpass data SpindleEpochExtendPlusMinus100/',subject,'/'];
filename = ['N2_data_bandpass11To15Hz_SpindleEpochExtendPlusMinus100_NoArousal_subject',subject,'_SpindleEpoch',num2str(iepoch),'_treatment.mat'];
end  
cd(foldername)
if exist(filename)==0
    disp(['input missing, baseline or treatment = ',num2str(baseline_or_treatment),' sub = ',num2str(sub_id),' epoch = ',num2str(iepoch)])
    continue
end
load(filename)

N2_data_1spindle_1t = nansum(N2_data_1spindle,1);
N2_data_1spindle_1t(N2_data_1spindle_1t==0) = nan;
idx = find(isnan(N2_data_1spindle_1t));
if nansum(N2_data_1spindle(:))==0 || nansum(idx(:))~=0
    disp(['input empty, baseline or treatment = ',num2str(baseline_or_treatment),' sub = ',num2str(sub_id),' epoch = ',num2str(iepoch)])
    continue
end 

% griddata******************
N2_data_1spindle_valid = N2_data_1spindle(working_channel_id_valid_idx,:);
xCord = [1:21];
yCord = [1:20];
[X,Y] = meshgrid(xCord, yCord);
X_lim = nanmax(xCord);
Y_lim = nanmax(yCord);
data_bandpa_2D_1sample = [];
for iTime= 1:size(N2_data_1spindle_valid,2)
    temp1_x = double(channel_locations_2d_rotate(:,1));
    temp1_y = double(channel_locations_2d_rotate(:,2));
    data_bandpa_2D_1sample(:,:,iTime) = griddata(temp1_x,temp1_y,double(N2_data_1spindle_valid(:,iTime)),double(X),double(Y),'linear') ;
end
% griddata******************
clearvars data_bandpa_1epoch_2D N2_data_1spindle 

% calculate phase
data_bandpa_phase_2D_1sample = nan(size(data_bandpa_2D_1sample));
for irow = 1:size(data_bandpa_2D_1sample,1)
    for icol = 1:size(data_bandpa_2D_1sample,2)
        temp1 = data_bandpa_2D_1sample(irow,icol,:);
        if nansum(temp1(:))== 0
            continue
        end
        data_bandpa_phase_2D_1sample(irow,icol,:) = angle(hilbert(temp1(:)));
    end
end

% calculate phase velocity field
cd([main_folder,'/NeuroPattToolbox-master'])
addpath(genpath(pwd))
cd(main_folder)

sampleFs = 0.002;
max_tt = size(data_bandpa_phase_2D_1sample,3);
max_row = size(data_bandpa_phase_2D_1sample,1);
max_col = size(data_bandpa_phase_2D_1sample,2);

% % % % % % % % % % % % % % % %
% % find velocity field: optical flow
params = setNeuroPattParams(1/sampleFs);
disp('Computing optical flow fields...'); tic

sigInterp_ds2 = data_bandpa_phase_2D_1sample ;
sigInterp_ds2(isnan(data_bandpa_phase_2D_1sample)) = -100 ; % AM this was -10,
vfs = zeros(size(data_bandpa_phase_2D_1sample));
vfs = vfs(:,:,1:end-1);

[vx, vy, csteps] = opticalFlow2(sigInterp_ds2(:,:,:), [], ...
params.opAlpha, params.opBeta, 1);%~params.useAmplitude);   %AM set last value to 1 for phase, 0 for amplitude

vx_n = vx./sqrt(vx.^2+vy.^2);
vy_n = vy./sqrt(vx.^2+vy.^2);
vfs(:,:,:) = vx_n + 1i*vy_n;
clearvars vy_n vx_n sigInterp_ds2

if baseline_or_treatment == 1
    foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Baseline/After 50Hz Notch Filter/20x21/NoArousal/Phase_PhaseVelcity_2Dgriddata/'];
elseif baseline_or_treatment == 2
    foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Treatment/After 50Hz Notch Filter/20x21/NoArousal/Phase_PhaseVelcity_2Dgriddata/'];
end
cd(foldername)
if exist(subject)==0
    mkdir(subject)
end
    
if baseline_or_treatment == 1
    foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Baseline/After 50Hz Notch Filter/20x21/NoArousal/Phase_PhaseVelcity_2Dgriddata/',subject,'/'];
    filename = ['20x21_preprocessed_data_SpindleEpochExtendPlusMinus100_filmis_50HzNotchFilt',subject,'_epoch',num2str(iepoch),'_baseline.mat'];    
elseif baseline_or_treatment == 2
    foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Treatment/After 50Hz Notch Filter/20x21/NoArousal/Phase_PhaseVelcity_2Dgriddata/',subject,'/'];
    filename = ['20x21_preprocessed_data_SpindleEpochExtendPlusMinus100_filmis_50HzNotchFilt',subject,'_epoch',num2str(iepoch),'_treatment.mat'];
end

save([foldername,filename],'N2_data_1spindle_valid','data_bandpa_2D_1sample','data_bandpa_phase_2D_1sample','vfs','-v7.3');    
disp(['bandpass data saved, baseline or treatment = ',num2str(baseline_or_treatment),' sub = ',num2str(sub_id),' epoch = ',num2str(iepoch)])


clearvars data_bandpa_2D_1sample data_bandpa_phase_2D_1sample data_bandpa_1sample N2_data_1spindle_valid N2_data_1spindle

end
end
end
%% spiral detection based on 20x21 phase velocity field

for baseline_or_treatment = 1:2
for sub_id = 1:9
    
subject_list = [{'01SC'} {'03RA'} {'06AC'} {'08RG'} {'09AS'} {'11NR'} {'12SM'} {'13PR'} {'15DC'} ];
subject = cell2mat(subject_list(sub_id));  
if baseline_or_treatment == 1
    foldername = [main_folder, '/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Baseline'];
    cd(foldername)
    load(['N2_spindle_epochs_accu_Molle2011_',subject,'_baseline.mat'])
elseif baseline_or_treatment == 2
    foldername = [main_folder, '/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Treatment'];    
    cd(foldername)
    load(['N2_spindle_epochs_accu_Molle2011_',subject,'_treatment.mat'])
end

for iepoch = 1:size(spindle_epochs_start_end_duration_iepoch_1sub,1)

if baseline_or_treatment == 1
    foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Baseline/After 50Hz Notch Filter/20x21/NoArousal/Phase_PhaseVelcity_2Dgriddata/',subject,'/'];
    filename = ['20x21_preprocessed_data_SpindleEpochExtendPlusMinus100_filmis_50HzNotchFilt',subject,'_epoch',num2str(iepoch),'_baseline.mat'];    
elseif baseline_or_treatment == 2
    foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Treatment/After 50Hz Notch Filter/20x21/NoArousal/Phase_PhaseVelcity_2Dgriddata/',subject,'/'];
    filename = ['20x21_preprocessed_data_SpindleEpochExtendPlusMinus100_filmis_50HzNotchFilt',subject,'_epoch',num2str(iepoch),'_treatment.mat'];
end

cd(foldername)
if exist(filename)==0
    disp(['input missing, baseline or treatment = ',num2str(baseline_or_treatment),' sub = ',num2str(sub_id),' epoch = ',num2str(iepoch)])
    continue
end
load(filename)

% spiral detection
% positive (anticlockwise) spirals *****************
curlz = [];
cav = [];

for time = 1:size(vfs,3)      
    temp1_vx = real(vfs(:,:,time));
    temp1_vy = imag(vfs(:,:,time));
    [x,y] = meshgrid(1:size(temp1_vx,2),1:size(temp1_vx,1));
    [curlz(:,:,time) cav(:,:,time)]  = curl(x,y,temp1_vx,temp1_vy);  % curl value of phase vector field    
end      
cav_filt_pos = curlz;
cav_filt_pos(cav_filt_pos<1) = 0;     % detect only clockwise vorties
cav_filt_pos(isnan(cav_filt_pos)) = 0; % filter out voxels points with curl value > -1

% vortex core detection based on filtered curl values (requires function
% pattDetection_V4.m)

params.minPattTime = 1; % minimum duration of vortex (10 time steps)
params.minPattSize = 3;  % minimum size of vortex (3x3 matrix)

[WCentroids,absoluteTime,instantTotalPower,pattSize,patternIdx] = pattDetection_v4(cav_filt_pos,cav_filt_pos,params,0,'CenterOfMass_Amp');
clearvars instantTotalPower patternIdx pattSize temp1_vx temp1_vy  curlz  cav cav_filt_pos

vortex_size_pos = [];
for ipatt = 1:size(absoluteTime,2) % select 1 vortex pattern detected  
    vortex_size_pos_1patt = [];
      for time_withinpatt = 1:size(absoluteTime{ipatt},2);  % select 1 time point within the duration of the detect vortex pattern
        clearvars i
        y_tempalte = [];
        x_tempalte = [];
        centerofmass_1patt_vxy_angle_dif = [];
        centerofmass_1patt_vxy_angle_dif_fil = [];
        temp1_time = absoluteTime{ipatt};
        time = temp1_time(time_withinpatt); 

        % find the angles of all the phase vectors within the defined 
        % vector at the defined time 
        temp1_vx = real(vfs(:,:,time)); 
        temp1_vy = imag(vfs(:,:,time));
        temp1_vxy_angle = angle(temp1_vx + i.*temp1_vy); 

        % centre of mass position of the defined vector at the definded
        % time
        temp1_centerofmass_1patt = WCentroids{ipatt};
        centerofmass_1patt = temp1_centerofmass_1patt(time_withinpatt,:);

        % calculate the angles of vortex-centre-originated vectors toward all 
        % voxel positions across the flattened cortex
        for irow = 1:20
            for icol = 1:21
                y_tempalte(irow,:) = irow;
                x_tempalte(:,icol) = icol;   
            end
        end
        centerofmass_1patt_vx = round(x_tempalte - centerofmass_1patt(1));
        centerofmass_1patt_vy = round(y_tempalte - centerofmass_1patt(2));
        centerofmass_1patt_vxy_angle = angle(centerofmass_1patt_vx + i.*centerofmass_1patt_vy);

        % lengths of vortex-centre-originated vectors/ distances
        centerofmass_1patt_vxy_abs = sqrt(centerofmass_1patt_vx.*centerofmass_1patt_vx+centerofmass_1patt_vy.*centerofmass_1patt_vy);
        centerofmass_1patt_vxy_abs(centerofmass_1patt_vxy_abs==0) = 1;

        % angle differences between ortex-centre-originated vectors and
        % phase vectors at the same voxel
        centerofmass_1patt_vxy_angle_dif = centerofmass_1patt_vxy_angle - temp1_vxy_angle;
        for irow = 1:20
            for icol = 1:21
                temp1 = centerofmass_1patt_vxy_angle_dif(irow,icol);
                if temp1> pi        % limit the angle differences between -pi and pi
                    centerofmass_1patt_vxy_angle_dif(irow,icol) = centerofmass_1patt_vxy_angle_dif(irow,icol) - 2*pi;
                elseif temp1< -pi   % limit the angle differences between -pi and pi
                    centerofmass_1patt_vxy_angle_dif(irow,icol) = centerofmass_1patt_vxy_angle_dif(irow,icol) + 2*pi;
                end
            end
        end
        centerofmass_1patt_vxy_angle_dif = abs(centerofmass_1patt_vxy_angle_dif); 

        % angle difference threshold: phase vectors with angle differences 
        % >135 or <45 degrees from the centre-orginated-vector are excluded 
        % from the vortex, ideally angle difference should be 90 degrees
        % for a perfect vortex
        angle_threshold_low = pi./2 - 0.7854; %***********%0.7854; %45degrees- %0.5236; % 30degrees below %0.2308; % 15 degrees below
        angle_threshold_high = pi./2 + 0.7854; %*********%0.7854; %45degrees+ %0.2536; % 30degrees above %0.2308; % 15 degrees above
        centerofmass_1patt_vxy_angle_dif_fil = centerofmass_1patt_vxy_angle_dif;
        centerofmass_1patt_vxy_angle_dif_fil(centerofmass_1patt_vxy_angle_dif_fil>angle_threshold_high) = nan;
        centerofmass_1patt_vxy_angle_dif_fil(centerofmass_1patt_vxy_angle_dif_fil<angle_threshold_low) = nan;

        % due to the ambiguity of phase values (and phase vector angles) 
        % close to the vortex centres, we assume 3x3 voxels centred at 
        % the vortex centre have a angle difference of 90 degrees, which
        % means they would always be included in the full-sized vortex
        centerofmass_1patt_round = round(centerofmass_1patt);
        for icol = centerofmass_1patt_round(1)-1:centerofmass_1patt_round(1)+1;
            for irow = centerofmass_1patt_round(2)-1:centerofmass_1patt_round(2)+1;
                if irow <= 20 && icol <=21 && irow > 0 && icol > 0
                centerofmass_1patt_vxy_angle_dif_fil(irow,icol) = 0.5*pi; 
                end
            end
        end

        % incrementally increase the radias from the vortex centre, until 
        % the total count of voxels excluded from the full-sized vortex 
        % exceeds a threshold 
        % threshold = expansion radius * expansion threshold
        % exapansion threshold: stop increment expansion of vortex radius if exceeded
        expansion_threshold = 1;   
        for d = 2:0.5:40    % gradually growing expansion radius, with increments of 0.5 
            % radius filter, remove voxels with distances larger than expansion radius    
            centerofmass_1patt_vxy_abs_fil = centerofmass_1patt_vxy_abs; 
            centerofmass_1patt_vxy_abs_fil(centerofmass_1patt_vxy_abs_fil>d) = nan; 
            temp1 = centerofmass_1patt_vxy_abs_fil./centerofmass_1patt_vxy_abs_fil;
            temp1_count_notnan = nansum(temp1(:)); % count voxels left 

            % combine radius filters and angle difference filter
            temp2 = centerofmass_1patt_vxy_abs_fil .* centerofmass_1patt_vxy_angle_dif_fil; 
            temp2_count_notnan = nansum(temp2(:)./temp2(:)); % count voxels left

            if d == 2
              temp1_count_notnan_PreviousIteration = 0;
            end
            if temp1_count_notnan - temp2_count_notnan > expansion_threshold.*d || temp1_count_notnan_PreviousIteration ==temp1_count_notnan;
                temp2_01 = temp2./temp2;
                temp2_01(isnan(temp2_01)) = 0;
                if centerofmass_1patt_round(1)>size(vfs,2)  % remove spiral at the edge of the cortical map
                    break
                end
                if centerofmass_1patt_round(2)>size(vfs,1)  % remove spiral at the edge of the cortical map
                    break
                end
                break
            end
            temp1_count_notnan_PreviousIteration = temp1_count_notnan;           
        end
        vortex_size_pos_1patt = [vortex_size_pos_1patt;d];

      end
      vortex_size_pos{ipatt} = vortex_size_pos_1patt;
end
vortex_maxsize_pos = [];
for ipatt = 1:size(vortex_size_pos,2)
    temp1 = vortex_size_pos{ipatt};
    max2 = nanmax(temp1);
    vortex_maxsize_pos{ipatt} = max2.*ones(size(temp1,1),1);
end
vortex_filt_pos_ipatt_t_centreXY_size_maxsize = [];
for ipatt = 1:size(WCentroids,2)
    temp1_centerofmass_1patt = WCentroids{ipatt};
    time = absoluteTime{ipatt}';
    radius = vortex_size_pos{ipatt};
    max_radius = vortex_maxsize_pos{ipatt};
    temp2 = [ipatt.*ones(size(time,1),1) time temp1_centerofmass_1patt radius max_radius];
    vortex_filt_pos_ipatt_t_centreXY_size_maxsize = [vortex_filt_pos_ipatt_t_centreXY_size_maxsize;temp2 ];
end

% negative (clockwise) spirals********************
curlz = [];
cav = [];
for time = 1:size(vfs,3)      
    temp1_vx = real(vfs(:,:,time));
    temp1_vy = imag(vfs(:,:,time));
    [x,y] = meshgrid(1:size(temp1_vx,2),1:size(temp1_vx,1));
    [curlz(:,:,time) cav(:,:,time)]  = curl(x,y,temp1_vx,temp1_vy);  % curl value of phase vector field    
end      
cav_filt_nega = curlz;
cav_filt_nega(cav_filt_nega>-1) = 0;     % detect only clockwise vorties
cav_filt_nega(isnan(cav_filt_nega)) = 0; % filter out voxels points with curl value > -1

% vortex core detection based on filtered curl values (requires function
% pattDetection_V4.m)
params.minPattTime = 1; % minimum duration of vortex (10 time steps)
params.minPattSize = 3;  % minimum size of vortex (3x3 matrix)

[WCentroids,absoluteTime,instantTotalPower,pattSize,patternIdx] = pattDetection_v4(cav_filt_nega,cav_filt_nega,params,0,'CenterOfMass_Amp');
clearvars instantTotalPower patternIdx pattSize temp1_vx temp1_vy cav_filt_nega   curlz cav

vortex_size_nega = [];
for ipatt = 1:size(absoluteTime,2) % select 1 vortex pattern detected  
    vortex_size_nega_1patt = [];
      for time_withinpatt = 1:size(absoluteTime{ipatt},2);  % select 1 time point within the duration of the detect vortex pattern
        clearvars i
        y_tempalte = [];
        x_tempalte = [];
        centerofmass_1patt_vxy_angle_dif = [];
        centerofmass_1patt_vxy_angle_dif_fil = [];
        temp1_time = absoluteTime{ipatt};
        time = temp1_time(time_withinpatt); 

        % find the angles of all the phase vectors within the defined 
        % vector at the defined time 
        temp1_vx =  real(vfs(:,:,time)); 
        temp1_vy =  imag(vfs(:,:,time));
        temp1_vxy_angle = angle(temp1_vx + i.*temp1_vy); 

        % centre of mass negaition of the defined vector at the definded
        % time
        temp1_centerofmass_1patt = WCentroids{ipatt};
        centerofmass_1patt = temp1_centerofmass_1patt(time_withinpatt,:);

        % calculate the angles of vortex-centre-originated vectors toward all 
        % voxel negaitions across the flattened cortex
        for irow = 1:20
            for icol = 1:21
                y_tempalte(irow,:) = irow;
                x_tempalte(:,icol) = icol;   
            end
        end
        centerofmass_1patt_vx = round(x_tempalte - centerofmass_1patt(1));
        centerofmass_1patt_vy = round(y_tempalte - centerofmass_1patt(2));
        centerofmass_1patt_vxy_angle = angle(centerofmass_1patt_vx + i.*centerofmass_1patt_vy);

        % lengths of vortex-centre-originated vectors/ distances
        centerofmass_1patt_vxy_abs = sqrt(centerofmass_1patt_vx.*centerofmass_1patt_vx+centerofmass_1patt_vy.*centerofmass_1patt_vy);
        centerofmass_1patt_vxy_abs(centerofmass_1patt_vxy_abs==0) = 1;

        % angle differences between ortex-centre-originated vectors and
        % phase vectors at the same voxel
        centerofmass_1patt_vxy_angle_dif = centerofmass_1patt_vxy_angle - temp1_vxy_angle;
        for irow = 1:20
            for icol = 1:21
                temp1 = centerofmass_1patt_vxy_angle_dif(irow,icol);
                if temp1> pi        % limit the angle differences between -pi and pi
                    centerofmass_1patt_vxy_angle_dif(irow,icol) = centerofmass_1patt_vxy_angle_dif(irow,icol) - 2*pi;
                elseif temp1< -pi   % limit the angle differences between -pi and pi
                    centerofmass_1patt_vxy_angle_dif(irow,icol) = centerofmass_1patt_vxy_angle_dif(irow,icol) + 2*pi;
                end
            end
        end
        centerofmass_1patt_vxy_angle_dif = abs(centerofmass_1patt_vxy_angle_dif); 

        % angle difference threshold: phase vectors with angle differences 
        % >135 or <45 degrees from the centre-orginated-vector are excluded 
        % from the vortex, ideally angle difference should be 90 degrees
        % for a perfect vortex
        angle_threshold_low = pi./2 - 0.7854; %***********%0.7854; %45degrees- %0.5236; % 30degrees below %0.2308; % 15 degrees below
        angle_threshold_high = pi./2 + 0.7854; %*********%0.7854; %45degrees+ %0.2536; % 30degrees above %0.2308; % 15 degrees above
        centerofmass_1patt_vxy_angle_dif_fil = centerofmass_1patt_vxy_angle_dif;
        centerofmass_1patt_vxy_angle_dif_fil(centerofmass_1patt_vxy_angle_dif_fil>angle_threshold_high) = nan;
        centerofmass_1patt_vxy_angle_dif_fil(centerofmass_1patt_vxy_angle_dif_fil<angle_threshold_low) = nan;

        % due to the ambiguity of phase values (and phase vector angles) 
        % close to the vortex centres, we assume 3x3 voxels centred at 
        % the vortex centre have a angle difference of 90 degrees, which
        % means they would always be included in the full-sized vortex
        centerofmass_1patt_round = round(centerofmass_1patt);
        for icol = centerofmass_1patt_round(1)-1:centerofmass_1patt_round(1)+1;
            for irow = centerofmass_1patt_round(2)-1:centerofmass_1patt_round(2)+1;
                if irow <= 20 && icol <=21 && irow > 0 && icol > 0
                centerofmass_1patt_vxy_angle_dif_fil(irow,icol) = 0.5*pi; 
                end
            end
        end

        % incrementally increase the radias from the vortex centre, until 
        % the total count of voxels excluded from the full-sized vortex 
        % exceeds a threshold 
        % threshold = expansion radius * expansion threshold
        % exapansion threshold: stop increment expansion of vortex radius if exceeded
        expansion_threshold = 1;   
        for d = 2:0.5:40    % gradually growing expansion radius, with increments of 0.5 
            % radius filter, remove voxels with distances larger than expansion radius    
            centerofmass_1patt_vxy_abs_fil = centerofmass_1patt_vxy_abs; 
            centerofmass_1patt_vxy_abs_fil(centerofmass_1patt_vxy_abs_fil>d) = nan; 
            temp1 = centerofmass_1patt_vxy_abs_fil./centerofmass_1patt_vxy_abs_fil;
            temp1_count_notnan = nansum(temp1(:)); % count voxels left 

            % combine radius filters and angle difference filter
            temp2 = centerofmass_1patt_vxy_abs_fil .* centerofmass_1patt_vxy_angle_dif_fil; 
            temp2_count_notnan = nansum(temp2(:)./temp2(:)); % count voxels left

            if d == 2
              temp1_count_notnan_PreviousIteration = 0;
            end
            if temp1_count_notnan - temp2_count_notnan > expansion_threshold.*d %|| temp1_count_notnan_PreviousIteration ==temp1_count_notnan;
                temp2_01 = temp2./temp2;
                temp2_01(isnan(temp2_01)) = 0;
                if centerofmass_1patt_round(1)>size(vfs,2)  % remove spiral at the edge of the cortical map
                    break
                end
                if centerofmass_1patt_round(2)>size(vfs,1)  % remove spiral at the edge of the cortical map
                    break
                end
                break
            end
            temp1_count_notnan_PreviousIteration = temp1_count_notnan;               
        end
        vortex_size_nega_1patt = [vortex_size_nega_1patt;d];
      end
      vortex_size_nega{ipatt} = vortex_size_nega_1patt;
end

vortex_maxsize_nega = [];
for ipatt = 1:size(vortex_size_nega,2)
    temp1 = vortex_size_nega{ipatt};
    max2 = nanmax(temp1);
    vortex_maxsize_nega{ipatt} = max2.*ones(size(temp1,1),1);
end
vortex_filt_nega_ipatt_t_centreXY_size_maxsize = [];
for ipatt = 1:size(WCentroids,2)
    temp1_centerofmass_1patt = WCentroids{ipatt};
    time = absoluteTime{ipatt}';
    radius = vortex_size_nega{ipatt};
    max_radius = vortex_maxsize_nega{ipatt};
    temp2 = [ipatt.*ones(size(time,1),1) time temp1_centerofmass_1patt radius max_radius];
    vortex_filt_nega_ipatt_t_centreXY_size_maxsize = [vortex_filt_nega_ipatt_t_centreXY_size_maxsize;temp2 ];
end

if baseline_or_treatment == 1
    foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Baseline/After 50Hz Notch Filter/20x21/NoArousal/spiral detected 2Dgriddata/'];
elseif baseline_or_treatment == 2
    foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Treatment/After 50Hz Notch Filter/20x21/NoArousal/spiral detected 2Dgriddata/'];
end
cd(foldername)
if exist(subject)==0
    mkdir(subject)
end

if baseline_or_treatment == 1
    foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Baseline/After 50Hz Notch Filter/20x21/NoArousal/spiral detected 2Dgriddata/',subject,'/'];
    filename = ['20x21_spiral_detected_SpindleEpochExtendPlusMinus100_filmis_50HzNotchFilt',subject,'_epoch',num2str(iepoch),'_baseline.mat'];    
elseif baseline_or_treatment == 2
    foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Treatment/After 50Hz Notch Filter/20x21/NoArousal/spiral detected 2Dgriddata/',subject,'/'];
    filename = ['20x21_spiral_detected_SpindleEpochExtendPlusMinus100_filmis_50HzNotchFilt',subject,'_epoch',num2str(iepoch),'_treatment.mat'];    
end
    save([foldername,filename],'vortex_filt_nega_ipatt_t_centreXY_size_maxsize','vortex_filt_pos_ipatt_t_centreXY_size_maxsize','-v7.3');    
    disp(['spiral detected saved, baseline or treatment = ',num2str(baseline_or_treatment),' sub = ',num2str(sub_id),' epoch = ',num2str(iepoch)])

    clearvars vfs
    
end
end
end  

%% spiral distrbution template

for baseline_or_treatment = 1:2
for sub_id = 1:9
    
subject_list = [{'01SC'} {'03RA'} {'06AC'} {'08RG'} {'09AS'} {'11NR'} {'12SM'} {'13PR'} {'15DC'} ];
subject = cell2mat(subject_list(sub_id));  
if baseline_or_treatment == 1
    foldername = [main_folder, '/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Baseline'];
    cd(foldername)
    load(['N2_spindle_epochs_accu_Molle2011_',subject,'_baseline.mat'])
elseif baseline_or_treatment == 2
    foldername = [main_folder, '/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Treatment'];    
    cd(foldername)
    load(['N2_spindle_epochs_accu_Molle2011_',subject,'_treatment.mat'])
end

for iepoch = 1:size(spindle_epochs_start_end_duration_iepoch_1sub,1)

if baseline_or_treatment == 1
    foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Baseline/After 50Hz Notch Filter/20x21/NoArousal/spiral detected 2Dgriddata/',subject,'/'];
    filename = ['20x21_spiral_detected_SpindleEpochExtendPlusMinus100_filmis_50HzNotchFilt',subject,'_epoch',num2str(iepoch),'_baseline.mat'];    
elseif baseline_or_treatment == 2
    foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Treatment/After 50Hz Notch Filter/20x21/NoArousal/spiral detected 2Dgriddata/',subject,'/'];
    filename = ['20x21_spiral_detected_SpindleEpochExtendPlusMinus100_filmis_50HzNotchFilt',subject,'_epoch',num2str(iepoch),'_treatment.mat'];    
end

cd(foldername)
if exist(filename)==0
    disp(['input missing, baseline or treatment = ',num2str(baseline_or_treatment),' sub = ',num2str(sub_id),' epoch = ',num2str(iepoch)])
    continue
end
load(filename)

radius_threshold = 6.5; % 95th percentile of the null model
duration_threshold = 114; % 3+ sigma cycles
t_duration_1spindle = spindle_epochs_start_end_duration_iepoch_1sub(iepoch,3) +200;
t_start = 101;
t_end = t_duration_1spindle - 100;
% add duration data
ipatt_nega = vortex_filt_nega_ipatt_t_centreXY_size_maxsize(:,1);
ipatt_nega_unique = unique(ipatt_nega);
duration_nega = nan(size(ipatt_nega));
for ipatt = 1:size(ipatt_nega_unique,1)
    temp1_ipatt =  ipatt_nega_unique(ipatt);
    idx = find(ipatt_nega==temp1_ipatt);
    duration_ipatt = nansum(idx./idx);
    duration_nega(idx) = duration_ipatt;
end
nega_ipatt_t_centreXY_size_maxsize_duration_spindleONLY_1sample = [vortex_filt_nega_ipatt_t_centreXY_size_maxsize duration_nega];

spiral_template_1sample = zeros(20,21);
for i2 = 1:size(vortex_filt_nega_ipatt_t_centreXY_size_maxsize,1)
    x = round(vortex_filt_nega_ipatt_t_centreXY_size_maxsize(i2,3));
    y = round(vortex_filt_nega_ipatt_t_centreXY_size_maxsize(i2,4));
    t = round(vortex_filt_nega_ipatt_t_centreXY_size_maxsize(i2,2));
%                     
    if t<t_start continue; end % skip if spiral not inside spindle epoch
    if t>t_end continue; end

    if x>21 x = 21; end 
    if x<1 x=1; end
    if y>20 y = 20; end 
    if y<1 y=1; end
    % remove spirals smaller than 95% of surrogate data
    max_radius = vortex_filt_nega_ipatt_t_centreXY_size_maxsize(i2,6);
    if max_radius < radius_threshold 
        continue 
    end 
    % remove spirals with duration smaller than x ms
    duration = nega_ipatt_t_centreXY_size_maxsize_duration_spindleONLY_1sample(i2,7);
    if duration < duration_threshold % if spiral rotate less than 1 cycle/ 76ms, remove 
        continue 
    end 

    spiral_template_1sample(y,x) = spiral_template_1sample(y,x) + 1;
end

% add duration data
ipatt_pos = vortex_filt_pos_ipatt_t_centreXY_size_maxsize(:,1);
ipatt_pos_unique = unique(ipatt_pos);
duration_pos = nan(size(ipatt_pos));
for ipatt = 1:size(ipatt_pos_unique,1)
    temp1_ipatt =  ipatt_pos_unique(ipatt);
    idx = find(ipatt_pos==temp1_ipatt);
    duration_ipatt = nansum(idx./idx);
    duration_pos(idx) = duration_ipatt;
end
pos_ipatt_t_centreXY_size_maxsize_duration_spindleONLY_1sample = [vortex_filt_pos_ipatt_t_centreXY_size_maxsize duration_pos];

for i2 = 1:size(vortex_filt_pos_ipatt_t_centreXY_size_maxsize,1)
    x = round(vortex_filt_pos_ipatt_t_centreXY_size_maxsize(i2,3));
    y = round(vortex_filt_pos_ipatt_t_centreXY_size_maxsize(i2,4));
    t = round(vortex_filt_pos_ipatt_t_centreXY_size_maxsize(i2,2));

    if t<t_start continue; end % skip if spiral not inside spindle epoch
    if t>t_end continue; end

    if x>21 x = 21; end 
    if x<1 x=1; end
    if y>20 y = 20; end 
    if y<1 y=1; end
    % remove spirals smaller than 95% of surrogate data
    max_radius = vortex_filt_pos_ipatt_t_centreXY_size_maxsize(i2,6);
    if max_radius < radius_threshold 
        continue 
    end 
    % remove spirals with duration smaller than x ms
    duration = pos_ipatt_t_centreXY_size_maxsize_duration_spindleONLY_1sample(i2,7);
    if duration < duration_threshold % if spiral rotate less than 1 cycle/ 76ms, remove 
        continue 
    end  
    spiral_template_1sample(y,x) = spiral_template_1sample(y,x) + 1;
    ipatt_last = ipatt;
end

if baseline_or_treatment == 1
    foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Baseline/After 50Hz Notch Filter/20x21/NoArousal/spiral distribution 2Dgriddata/'];
elseif baseline_or_treatment == 2
    foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Treatment/After 50Hz Notch Filter/20x21/NoArousal/spiral distribution 2Dgriddata/'];
end
cd(foldername)
if exist(subject)==0
    mkdir(subject)
end
   
if baseline_or_treatment == 1
    foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Baseline/After 50Hz Notch Filter/20x21/NoArousal/spiral distribution 2Dgriddata/',subject,'/'];
    filename = ['20x21_spiral_distribution_SpindleEpochExtendPlusMinus100_',num2str(radius_threshold),'plusRadius_',num2str(duration_threshold),'plusDuration_filmis_50HzNotchFilt',subject,'_epoch',num2str(iepoch),'_baseline.mat'];    
elseif baseline_or_treatment == 2
    foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Treatment/After 50Hz Notch Filter/20x21/NoArousal/spiral distribution 2Dgriddata/',subject,'/'];
    filename = ['20x21_spiral_distribution_SpindleEpochExtendPlusMinus100_',num2str(radius_threshold),'plusRadius_',num2str(duration_threshold),'plusDuration_filmis_50HzNotchFilt',subject,'_epoch',num2str(iepoch),'_treatment.mat'];    
end

save([foldername,filename],'spiral_template_1sample','-v7.3');    
disp(['spiral distribution saved, baseline or treatment = ',num2str(baseline_or_treatment),' sub = ',num2str(sub_id),' epoch = ',num2str(iepoch)])

end
end
end

%% spindle-epoch-avg & subject-avg spiral distribution

% spindle-epoch-avg spiral distribution 
for baseline_or_treatment = 1:2
    for sub_id = 1:9 
        spiral_template_SpindleEpochs_N2epochAccu = [];
        disp(['baseline or treatment = ',num2str(baseline_or_treatment),' sub = ',num2str(sub_id)])
        subject_list = [{'01SC'} {'03RA'} {'06AC'} {'08RG'} {'09AS'} {'11NR'} {'12SM'} {'13PR'} {'15DC'} ];
        subject = cell2mat(subject_list(sub_id));
        % load spindle epoch time
        if baseline_or_treatment == 1
        foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Baseline/After 50Hz Notch Filter/'];
        filename = ['N2_spindle_epochs_accu_Molle2011_',subject,'_baseline.mat'];
        elseif baseline_or_treatment == 2
        foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Treatment/After 50Hz Notch Filter/'];
        filename = ['N2_spindle_epochs_accu_Molle2011_',subject,'_treatment.mat'];
        end
        cd(foldername)
        if exist(filename)
            load(filename)
        else
            disp(['input doesnt exist, baseline or treatment = ',num2str(baseline_or_treatment),' sub = ',num2str(sub_id)])    
            continue
        end      
        % load arousal signal
        if baseline_or_treatment == 1
            foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Baseline/'];
            filename = ['N2_arousal_data_',subject,'_baseline.mat'];
        elseif baseline_or_treatment == 2
            foldername = [main_folder,'/HD_EEG/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Treatment/'];
            filename = ['N2_arousal_data_',subject,'_treatment.mat'];
        end
        cd(foldername)
        if exist(filename)
            load(filename)
        else
            disp(['input doesnt exist, baseline or treatment = ',num2str(baseline_or_treatment),' sub = ',num2str(sub_id)])     
            continue
        end 
        
        radius_threshold = 6.5;
        duration_threshold = 114;
        for iepoch = 1:size(spindle_epochs_start_end_duration_iepoch_1sub,1)
            % remove spindles that are too close to the start and end of N2
            % epoch: <10000 from start and end of 1 N2 epoch************
            N2_epoch_id_1spindle = spindle_epochs_start_end_duration_iepoch_1sub(iepoch,4);
            N2_epoch_t_start_1spindle = spindle_epochs_start_end_duration_iepoch_1sub(iepoch,1);
            N2_epoch_t_end_1spindle = spindle_epochs_start_end_duration_iepoch_1sub(iepoch,2);
            N2_epoch_end = epoch_end(N2_epoch_id_1spindle);
            distance_to_N2epoch_end = N2_epoch_end-N2_epoch_t_end_1spindle;
            distance_to_N2epoch_boundary_threshold = 10000;
            if N2_epoch_t_start_1spindle<distance_to_N2epoch_boundary_threshold || distance_to_N2epoch_end<distance_to_N2epoch_boundary_threshold %70000
                continue
                disp(['spindle too close to boundary, baseline or treatment = ',num2str(baseline_or_treatment),' sub = ',num2str(sub_id), ' epoch=',num2str(iepoch)])    
            end
            % remove spindles that are too close to the start and end of N2
            % epoch: <10000 from start and end of 1 N2 epoch************
            if baseline_or_treatment == 1
                foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Baseline/After 50Hz Notch Filter/20x21/NoArousal/spiral distribution 2Dgriddata/',subject,'/'];
                filename = ['20x21_spiral_distribution_SpindleEpochExtendPlusMinus100_',num2str(radius_threshold),'plusRadius_',num2str(duration_threshold),'plusDuration_filmis_50HzNotchFilt',subject,'_epoch',num2str(iepoch),'_baseline.mat'];    
            elseif baseline_or_treatment == 2
                foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Treatment/After 50Hz Notch Filter/20x21/NoArousal/spiral distribution 2Dgriddata/',subject,'/'];
                filename = ['20x21_spiral_distribution_SpindleEpochExtendPlusMinus100_',num2str(radius_threshold),'plusRadius_',num2str(duration_threshold),'plusDuration_filmis_50HzNotchFilt',subject,'_epoch',num2str(iepoch),'_treatment.mat'];    
            end
            cd(foldername)
            if exist(filename)
                load(filename)
            else
                continue
            end
            spiral_template_SpindleEpochs_N2epochAccu{iepoch,N2_epoch_id_1spindle} = spiral_template_1sample;
        end
        
        if baseline_or_treatment == 1
            foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Baseline/After 50Hz Notch Filter/20x21/NoArousal/spiral distribution analysis/'];
            filename = ['20x21_spiral_distribution_N2epochAccu_SpindleEpochExtendPlusMinus100_Radius',num2str(radius_threshold),'Plus_Duration',num2str(duration_threshold),'PLus_NoArousal_RemoveBoundary',num2str(distance_to_N2epoch_boundary_threshold),'_subject',num2str(subject),'_baseline.mat'];
        elseif baseline_or_treatment == 2
            foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Treatment/After 50Hz Notch Filter/20x21/NoArousal/spiral distribution analysis/']; 
            filename = ['20x21_spiral_distribution_N2epochAccu_SpindleEpochExtendPlusMinus100_Radius',num2str(radius_threshold),'Plus_Duration',num2str(duration_threshold),'PLus_NoArousal_RemoveBoundary',num2str(distance_to_N2epoch_boundary_threshold),'_subject',num2str(subject),'_treatment.mat'];
        end
        save([foldername,filename],'spiral_template_SpindleEpochs_N2epochAccu','-v7.3');  
 
        % subject-avg spiral distribution
        spiral_template_spindles_1N2epochavg_1sub = [];
        for N2_epoch_id_1spindle = 1:size(spiral_template_SpindleEpochs_N2epochAccu,2)
            count = 0;
            spiral_template_spindles_1N2epoch = [];
            for iepoch = 1:size(spiral_template_SpindleEpochs_N2epochAccu,1)
                spiral_template_1sample = spiral_template_SpindleEpochs_N2epochAccu{iepoch,N2_epoch_id_1spindle};
                if nansum(spiral_template_1sample(:))~=0
                   count = count + 1;
                   spiral_template_spindles_1N2epoch(:,:,count) = spiral_template_1sample;
                end
            end
            if nansum(spiral_template_spindles_1N2epoch(:))==0
                continue
            end
            spiral_template_spindles_1N2epochavg = nansum(spiral_template_spindles_1N2epoch,3);
            spiral_template_spindles_1N2epochavg_1sub(:,:,N2_epoch_id_1spindle) = spiral_template_spindles_1N2epochavg;
        end
        spiral_template_spindles_1N2epochavg_1subavg = nansum(spiral_template_spindles_1N2epochavg_1sub,3);
         
        if baseline_or_treatment == 1
            foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Baseline/After 50Hz Notch Filter/20x21/NoArousal/spiral distribution analysis/'];
            filename = ['20x21_spiral_distribution_N2epochAvg_SpindleEpochExtendPlusMinus100_Radius',num2str(radius_threshold),'Plus_Duration',num2str(duration_threshold),'PLus_NoArousal_RemoveBoundary',num2str(distance_to_N2epoch_boundary_threshold),'_subject',num2str(subject),'_baseline.mat'];
        elseif baseline_or_treatment == 2
            foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Treatment/After 50Hz Notch Filter/20x21/NoArousal/spiral distribution analysis/']; 
            filename = ['20x21_spiral_distribution_N2epochAvg_SpindleEpochExtendPlusMinus100_Radius',num2str(radius_threshold),'Plus_Duration',num2str(duration_threshold),'PLus_NoArousal_RemoveBoundary',num2str(distance_to_N2epoch_boundary_threshold),'_subject',num2str(subject),'_treatment.mat'];
        end
        save([foldername,filename],'spiral_template_spindles_1N2epochavg_1sub','spiral_template_spindles_1N2epochavg_1subavg','-v7.3');  
  
end
end


%% short-term (overnight)/long-term (3-month) consistency based on subject/N2-epoch-averaged distribution

for sub_id = 1:9
    for baseline_or_treatment = 1:2
        disp(['baseline or treatment = ',num2str(baseline_or_treatment),' sub = ',num2str(sub_id)])
        subject_list = [{'01SC'} {'03RA'} {'06AC'} {'08RG'} {'09AS'} {'11NR'} {'12SM'} {'13PR'} {'15DC'} ];
        subject = cell2mat(subject_list(sub_id));
        if baseline_or_treatment == 1
            foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Treatment/After 50Hz Notch Filter/20x21/NoArousal/spiral distribution analysis/'];
        elseif baseline_or_treatment == 2
            foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Treatment/After 50Hz Notch Filter/20x21/NoArousal/spiral distribution analysis/'];
        end
        cd(foldername)
        if exist(subject)==0
            mkdir(subject)
        end    
    end
end

radius_threshold = 6.5;
duration_threshold = 114;
distance_to_N2epoch_boundary_threshold = 10000;
        
R2_baseline = [];
R2_treatment = [];
R2_baseline_x_treatment = [];
R2_baseline_x_treatment_1subavg = [];
for sub_id = 1:9 

    for baseline_or_treatment = 1:2     
        disp(['baseline or treatment = ',num2str(baseline_or_treatment),' sub = ',num2str(sub_id)])
        subject_list = [{'01SC'} {'03RA'} {'06AC'} {'08RG'} {'09AS'} {'11NR'} {'12SM'} {'13PR'} {'15DC'} ];
        subject = cell2mat(subject_list(sub_id));
        
        if baseline_or_treatment == 1
            foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Baseline/After 50Hz Notch Filter/20x21/NoArousal/spiral distribution analysis/'];
            filename = ['20x21_spiral_distribution_N2epochAvg_SpindleEpochExtendPlusMinus100_Radius',num2str(radius_threshold),'Plus_Duration',num2str(duration_threshold),'PLus_NoArousal_RemoveBoundary',num2str(distance_to_N2epoch_boundary_threshold),'_subject',num2str(subject),'_baseline.mat'];
        elseif baseline_or_treatment == 2
            foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Treatment/After 50Hz Notch Filter/20x21/NoArousal/spiral distribution analysis/']; 
            filename = ['20x21_spiral_distribution_N2epochAvg_SpindleEpochExtendPlusMinus100_Radius',num2str(radius_threshold),'Plus_Duration',num2str(duration_threshold),'PLus_NoArousal_RemoveBoundary',num2str(distance_to_N2epoch_boundary_threshold),'_subject',num2str(subject),'_treatment.mat'];
        end
        cd(foldername)
        if exist(filename)==0
            continue
        end
        load(filename)

        % spiral short-term consistency 
        if baseline_or_treatment == 1
            spiral_template_N2epoch_accu_baseline = spiral_template_spindles_1N2epochavg_1sub;
            for i2 = 1:size(spiral_template_N2epoch_accu_baseline,3)
                temp1 = spiral_template_N2epoch_accu_baseline(:,:,i2);
                if nansum(temp1(:))==0
                    spiral_template_N2epoch_accu_baseline(:,:,i2) = nan;
                end
            end
            R2 = nan(size(spiral_template_N2epoch_accu_baseline,3),size(spiral_template_N2epoch_accu_baseline,3));
            for iepoch1 = 1:size(spiral_template_N2epoch_accu_baseline,3)
                for iepoch2 = 1:size(spiral_template_N2epoch_accu_baseline,3)
                    temp1 = spiral_template_N2epoch_accu_baseline(:,:,iepoch1);
                    temp2 = spiral_template_N2epoch_accu_baseline(:,:,iepoch2);
                    R = corrcoef(temp1(:),temp2(:));
                    R2(iepoch1,iepoch2) = R(2);
                end
            end
            R2(R2==0) = nan; 
            R2(R2>0.999) = nan;  
            R2_baseline(sub_id) = nanmean(R2(:));
        elseif baseline_or_treatment == 2
            spiral_template_N2epoch_accu_treatment = spiral_template_spindles_1N2epochavg_1sub;
            for i2 = 1:size(spiral_template_N2epoch_accu_treatment,3)
                temp1 = spiral_template_N2epoch_accu_treatment(:,:,i2);
                if nansum(temp1(:))==0
                    spiral_template_N2epoch_accu_treatment(:,:,i2) = nan;
                end
            end
            R2 = nan(size(spiral_template_N2epoch_accu_treatment,3),size(spiral_template_N2epoch_accu_treatment,3));
            for iepoch1 = 1:size(spiral_template_N2epoch_accu_treatment,3)
                for iepoch2 = 1:size(spiral_template_N2epoch_accu_treatment,3)
                    temp1 = spiral_template_N2epoch_accu_treatment(:,:,iepoch1);
                    temp2 = spiral_template_N2epoch_accu_treatment(:,:,iepoch2);
                    R = corrcoef(temp1(:),temp2(:));
                    R2(iepoch1,iepoch2) = R(2);
                end
            end
            R2(R2==0) = nan;
            R2(R2>0.999) = nan;  
            R2_treatment(sub_id) = nanmean(R2(:));
        end
        
        % long-term consistency: baseline x treatment
        if baseline_or_treatment == 2
            exist spiral_template_N2epoch_accu_baseline;
            temp1 = ans;
            if temp1==0
                spiral_template_N2epoch_accu_baseline = [];
            end
            exist spiral_template_N2epoch_accu_treatment;
            temp1 = ans;
            if temp1==0
                spiral_template_N2epoch_accu_treatment = [];
            end    

            if nansum(spiral_template_N2epoch_accu_baseline(:))==0 || nansum(spiral_template_N2epoch_accu_treatment(:))==0
                continue
            end
            R2 = nan(size(spiral_template_N2epoch_accu_baseline,3),size(spiral_template_N2epoch_accu_treatment,3)); 
            for iepoch1 = 1:size(spiral_template_N2epoch_accu_baseline,3)
                for iepoch2 = 1:size(spiral_template_N2epoch_accu_treatment,3)
                    temp1_baseline = spiral_template_N2epoch_accu_baseline(:,:,iepoch1);
                    temp1_treatment = spiral_template_N2epoch_accu_treatment(:,:,iepoch2);
                    if nansum(temp1_baseline(:))==0 || nansum(temp1_treatment(:))==0
                        continue
                    end
                    R = corrcoef(temp1_baseline(:),temp1_treatment(:));
                    R2(iepoch1,iepoch2) = R(2);
                end
            end
            R2(R2==0) = nan;
            R2(R2>0.999) = nan;
            R2_baseline_x_treatment(sub_id) = nanmean(R2(:));
            
            spiral_template_N2epoch_accu_baseline_1subavg = nansum(spiral_template_N2epoch_accu_baseline,3);
            spiral_template_N2epoch_accu_treatment_1subavg = nansum(spiral_template_N2epoch_accu_treatment,3);
           
            R = corrcoef(spiral_template_N2epoch_accu_baseline_1subavg(:),spiral_template_N2epoch_accu_treatment_1subavg(:));
            R2_baseline_x_treatment_1subavg(sub_id) = R(2);
        end
    end
end

filename = ['short_long_term_stability_N2epochAvg_20x21_noarousal_Radius',num2str(radius_threshold),'plus_Duration',num2str(duration_threshold),'plus_RemoveBoundary',num2str(distance_to_N2epoch_boundary_threshold),'_9sub_baseline_treatment.mat'];
foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Baseline/After 50Hz Notch Filter/20x21/NoArousal/spiral distribution analysis/'];
save([foldername,filename],'R2_baseline_x_treatment','R2_treatment','R2_baseline','R2_baseline_x_treatment_1subavg')
foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Treatment/After 50Hz Notch Filter/20x21/NoArousal/spiral distribution analysis/'];
save([foldername,filename],'R2_baseline_x_treatment','R2_treatment','R2_baseline','R2_baseline_x_treatment_1subavg')

%% spiral long-term consistency predicts age/memory

filename = ['short_long_term_stability_N2epochAvg_20x21_noarousal_Radius',num2str(radius_threshold),'plus_Duration',num2str(duration_threshold),'plus_RemoveBoundary',num2str(distance_to_N2epoch_boundary_threshold),'_9sub_baseline_treatment.mat'];
foldername = [main_folder,'/PRJ-HdEEG_OSA_CPAP/UWM/N2 data clean/Baseline/After 50Hz Notch Filter/20x21/NoArousal/spiral distribution analysis/'];
cd(foldername)
load(filename)
R2_baseline_x_treatment_subavg = R2_baseline_x_treatment_1subavg;

age_9sub = [46.8528   53.5667   64.8444   49.5417   44.8889   53.6028   43.6861   53.3694   41.6667];
temp1_baseline_pre_score = [20    16    19    20    28    13    31    15    25];
temp1_baseline_post_score = [19    11    17    20    27     8    30    16    26];
temp1_treatment_pre_score = [20    18    19    30    32    17    31    22    26];
temp1_treatment_post_score = [22    13    15    27    32    14    31    16    24];
temp1_baseline_postpre_score = temp1_baseline_post_score./temp1_baseline_pre_score;
temp1_treatment_postpre_score = temp1_treatment_post_score./temp1_treatment_pre_score;

% baseline-treatment avg %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp1_baseline_treatment_pre_score_avg = (temp1_baseline_pre_score + temp1_treatment_pre_score)./2;
temp1_baseline_treatment_post_score_avg = (temp1_baseline_post_score + temp1_treatment_post_score)./2;
temp1_baseline_treatment_postpre_score_avg = temp1_baseline_treatment_post_score_avg ./ temp1_baseline_treatment_pre_score_avg;
mdl_BaselineTreatmentAvg_presleep_recall_R2_baseline_x_treatment = fitlm(R2_baseline_x_treatment_subavg(:),temp1_baseline_treatment_pre_score_avg(:));
mdl_BaselineTreatmentAvg_postsleep_recall_R2_baseline_x_treatment = fitlm(R2_baseline_x_treatment_subavg(:),temp1_baseline_treatment_post_score_avg(:));
mdl_BaselineTreatmentAvg_postpresleep_recall_R2_baseline_x_treatment = fitlm(R2_baseline_x_treatment_subavg(:),temp1_baseline_treatment_postpre_score_avg(:));
correlation_BaselineTreatmentAvg_presleep_recall_R2_baseline_x_treatment = corrcoef(R2_baseline_x_treatment_subavg(:),temp1_baseline_treatment_pre_score_avg(:));
correlation_BaselineTreatmentAvg_postsleep_recall_R2_baseline_x_treatment = corrcoef(R2_baseline_x_treatment_subavg(:),temp1_baseline_treatment_post_score_avg(:));
correlation_BaselineTreatmentAvg_postpresleep_recall_R2_baseline_x_treatment = corrcoef(R2_baseline_x_treatment_subavg(:),temp1_baseline_treatment_postpre_score_avg(:));
mdl_age_memory_postprechange_baseline_treatment_avg = fitlm(R2_baseline_x_treatment_subavg,age_9sub(:));
correlation_age_memory_postprechange_baseline_treatment_avg = corrcoef(R2_baseline_x_treatment_subavg,age_9sub(:));

mdl_Baseline_postpresleep_recall_R2_baseline_x_treatment = fitlm(R2_baseline_x_treatment_subavg(:),temp1_baseline_postpre_score(:));
R2_Baseline_postpresleep_recall_R2_baseline_x_treatment = corrcoef(R2_baseline_x_treatment_subavg(:),temp1_baseline_postpre_score(:));
mdl_Treatment_postpresleep_recall_R2_baseline_x_treatment = fitlm(R2_baseline_x_treatment_subavg(:),temp1_treatment_postpre_score(:));
R2_Treatment_postpresleep_recall_R2_baseline_x_treatment = corrcoef(R2_baseline_x_treatment_subavg(:),temp1_treatment_postpre_score(:));

figure(1)
scatter(R2_baseline_x_treatment_subavg(:),temp1_baseline_treatment_post_score_avg(:)./temp1_baseline_treatment_pre_score_avg(:),[],'r','filled')
hold on
plot(mdl_BaselineTreatmentAvg_postpresleep_recall_R2_baseline_x_treatment)
hold off
xlabel(['spiral long-term consistency'])
ylabel(['memory retention (baseline + treatment)'])
title(['spiral long-term consistency predicts memory retention'])

% spiral short-term consistency: baseline vs treatment
R2_baseline_subavg_avg = nanmean(R2_baseline(:));
R2_baseline_subavg_std = nanstd(R2_baseline(:));
R2_treatment_subavg_avg = nanmean(R2_treatment(:));
R2_treatment_subavg_std = nanstd(R2_treatment(:));
figure(2)
err_std = [R2_baseline_subavg_std R2_treatment_subavg_std];
avg = abs([R2_baseline_subavg_avg R2_treatment_subavg_avg]);
data = [R2_baseline_subavg(:) R2_treatment_subavg(:)];
x = [1 2];
XTicklabel = {'baseline','treatment'};
bar(x,avg)
hold on
er = errorbar(x,avg,err_std,err_std)
er.Color = [0 0 0];
er.LineStyle = 'none';
hold on
scatter(ones(size(R2_baseline_subavg(:))).*(1+(rand(size(R2_baseline_subavg(:)))-0.5)/3),R2_baseline_subavg(:),'r','filled')
scatter(ones(size(R2_treatment_subavg(:))).*(1+(rand(size(R2_treatment_subavg(:)))+2.5)/3),R2_treatment_subavg(:),'r','filled')
hold off
title(['spiral short-term consistency, baseline vs treatment'])

% spiral long-term consistency predicts age
figure(3)
scatter(R2_baseline_x_treatment_subavg(:),age_9sub(:),[],'r','filled')
hold on
plot(mdl_age_memory_postprechange_baseline_treatment_avg)
hold off
xlabel(['spiral long-term consistency'])
ylabel(['age'])
title(['spiral long-term consistency predicts age'])

% spiral long-term consistency is the dominant predictor of memory retention 
% multiple linear regression: spiral long-term consistency + age predict memory: 
X = [ones(9,1) R2_baseline_x_treatment_subavg(:) age_9sub(:)];
y = temp1_baseline_treatment_post_score_avg(:)./temp1_baseline_treatment_pre_score_avg(:);
[b,bint,r,rint,stats] = regress(y,X);
R_square_R2_plus_age_predict_memory = stats(1);
% comparing to linear regression with one independent variable (spiral
% long-term consistency)
X = [ones(9,1) R2_baseline_x_treatment_subavg(:)];
y = temp1_baseline_treatment_post_score_avg(:)./temp1_baseline_treatment_pre_score_avg(:);
[b,bint,r,rint,stats] = regress(y,X);
R_square_R2_predict_memory = stats(1);
% unique contribution by age alone
R_square_age_predict_memory_unique = R_square_R2_plus_age_predict_memory - R_square_R2_predict_memory;

% Model with only age as sole predictor
X = [ones(9,1) age_9sub(:)];
[b, bint, r, rint, stats] = regress(y, X);
R_square_age_only = stats(1);
% Unique contribution of R2 (spiral long-term consistency)
R_square_R2_predict_memory_unique = R_square_R2_plus_age_predict_memory - R_square_age_only;

