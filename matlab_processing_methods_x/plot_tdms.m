

clear all
close all

%% give initial thresholds
Fs=1000; % Sampling frequency of the signal
FOI=0.5; % cut off frequency for filtering (in Hz)

%Which plot you would like to see?
plot_raw_volume=0; %1=yes/0=no
plot_filterd_volume=1; %1=yes/0=no
plot_speed=1; %1=yes/0=no

out_put_format='.pdf'; %give the format of the output figure (e.g.: .png). Deafult is .fig


%% read the files from folder

dir_name = uigetdir;
Path = dir_name;
dir_file = dir(dir_name);
num_file = size(dir_file);


Suck_all_time=[];
addpath(dir_name)
cntr_subj=0;

for i=1:num_file(1)
    File = dir_file(i).name;
    [~,name,ext] = fileparts([dir_name,'/',File]);
    
    xls_e = strcmp(ext, '.tdms');
    
    if (xls_e > 0)
              cntr_subj=cntr_subj+1;
        
        FileName{cntr_subj}=File;
        TDMS = TDMS_readTDMSFile(fullfile([dir_name,'/',File]));
        
        
            nofSamples=size(TDMS.data{3},2);
            Data1=zeros(size(TDMS.data,2)-3,nofSamples);
            for idx=4:size(TDMS.data,2)
            Data1(idx-3,:)=TDMS.data{idx};
            end

            t=(1:nofSamples)./Fs;

            fidInfo = fopen(fullfile(dir_name,[name,'.info']), 'r');
            rawLines = textscan(fidInfo,'%s','Delimiter','\n','Whitespace','');
            fclose(fidInfo);
            
            rawLines = rawLines{1};
            params = containers.Map;
            
            for i = 1:numel(rawLines)
                line = strtrim(rawLines{i});
                if isempty(line) || startsWith(line,'%')
                    continue
                end
                % Match: "Key: Value"
                tokens = regexp(line,'^([^:]+):\s*([-\d\.eE]+)$','tokens');
                if ~isempty(tokens)
                    key = strtrim(tokens{1}{1});
                    val = str2double(tokens{1}{2});
                    params(key) = val;
                end
            end
            
            % Debug: list all available keys
            disp(params.keys)

            weight=zeros(4,nofSamples);
            height=zeros(4,nofSamples);
            temp=zeros(4,nofSamples);

            for IDX_unit=1:4
                %apply the parameters form the output linearization
                weight(IDX_unit,:) = (Data1(1+(IDX_unit-1)*3,:)).*C{1,2}...
                (26+(IDX_unit-1)*2)+C{1,2}(27+(IDX_unit-1)*2)-C{1,2}...
                (46+(IDX_unit-1)*2);
                height(IDX_unit,:) = (Data1(2+(IDX_unit-1)*3,:)).*C{1,2}...
                (18+(IDX_unit-1)*2)+C{1,2}(19+(IDX_unit-1)*2)-C{1,2}...
                (47+(IDX_unit-1)*2);
                temp(IDX_unit,:) = (Data1(3+(IDX_unit-1)*3,:)).*C{1,2}...
                (34+(IDX_unit-1)*2)+C{1,2}(35+(IDX_unit-1)*2);
            end
            
            for idx=4:size(TDMS.data,2)
                Data1(idx-3,:)=TDMS.data{idx};
            end

            % --- Derive available units from Data1 rows instead of assuming 4 ---
            nRows  = size(Data1,1);
            nUnits = floor(nRows/3);
            if nUnits < 1
                error('TDMS has no complete unit triplets (rows=%d).', nRows);
            end
            if nUnits < 4
                warning('Only %d unit(s) present in TDMS (rows=%d). Looping over %d units.', nUnits, nRows, nUnits);
            end

            rate=1;
            weight_ds=zeros(4,ceil(nofSamples/rate));
            height_ds=zeros(4,ceil(nofSamples/rate));
            temp_ds=zeros(4,ceil(nofSamples/rate));
            
            for IDX_unit = 1:nUnits
            weight_ds(IDX_unit,:) = weight(IDX_unit,:);
            height_ds(IDX_unit,:) = height(IDX_unit,:);
            temp_ds(IDX_unit,:) = temp(IDX_unit,:);
            end

            weight_comp=zeros(4,ceil(nofSamples/rate));
            A_tube=2*pi*(C{1,2}(10)/2)^2-2*pi*(C{1,2}(11)/2)^2;
            A_hose=2*pi*(C{1,2}(12)/2)^2-2*pi*(C{1,2}(13)/2)^2;
            A=A_tube+A_hose;
            for IDX_unit=1:4
            Density(IDX_unit)=C{1,2}(42+(IDX_unit-1));
            weight_comp(IDX_unit,:)=A.*height_ds(IDX_unit,:).*Density(IDX_unit);
            end
            % Appliy correction for weight difference by buonancy
            Volume=weight_ds-weight_comp;
            Volume=ones(size(Volume)).*max(Volume)-Volume;
            t1=rate/Fs:rate/Fs:nofSamples/(Fs/rate)/rate;
            t2=0:rate/Fs:nofSamples/(Fs/rate)/rate;
            
            Volume=Volume(1,:);
            V=Volume(1,:);
            data_filt = bst_bandpass_fft(Volume(1,:)*1000, Fs, 0,FOI, 1, 1);
            Speed=gradient(data_filt)' .*(Fs/rate);
            
            if plot_raw_volume==1
                f=figure;
                plot(t,data_filt)
                xlabel('time (s)')
                ylabel('volume (ml)')
                title (['Volume_raw_',File],'Interpreter', 'none')
                saveas(f,['Volume_raw_',File,out_put_format])
                close(f)
            end
            
            if plot_filterd_volume==1
                f=figure;
                plot(t,data_filt)
                xlabel('time (s)')
                ylabel('volume (ml)')
                title (['Volume_filtered_',File],'Interpreter', 'none')
                saveas(f,['Volume_filtered_',File,out_put_format])
                close(f)
            end
            
            if plot_speed==1
                f=figure;
                plot(t,Speed)
                xlabel('time (s)')
                ylabel('speed (ml/s)')
                title (['Speed_',File],'Interpreter', 'none')
                saveas(f,['Speed_',File,out_put_format])
                close(f)
            end
    end
end
