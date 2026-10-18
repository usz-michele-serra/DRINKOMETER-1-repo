
clear all
close all


%% give initial thresholds - here you can modify the parameters
Fs=1000; %Hz; Sampling frequency of the signal
FOI=0.5; %Hz; cut off frequency for filtering (in Hz)
artefact_speed_threshold=50; %ml/sec; Speed signal is considered as an artefact above the threshold (absolute value). The detected interval is deleted and interpolated linearly.
min_amplitude_of_suck=0.5; %ml/sec; speed is set to zero below this threshold
min_peak_distance=0.5; %sec; the minimal distance between peaks sucks (prevent to detect noine as suck) 
min_suck_amount=2; %ml; below this threshold the sucks considered as noise and not included in the analysis.
name_of_the_refill_table='refill_table.xlsx';
refill_table=0; %1=yes, 0=no
want_figures=1; %1=yes, 0=no
histogram_bins=100; %number of columns (bins) in the histogramm for the GM

%% read the files from folder

dir_name = uigetdir;
Path = dir_name;
dir_file = dir(dir_name);
num_file = size(dir_file);

%% read refill data
if refill_table
    [interval_refill,filenames_refill] = xlsread(name_of_the_refill_table);
    filenames_refill=filenames_refill(2:end,1);
end

did_debug_plot = false;

Suck_all_time=[];
addpath(dir_name)
cntr_subj=0;
 cntr_no_subj=0;
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

            N = numel(t);
            win = 2000;  % samples around artefact

            % --- Read .info into key-value pairs (params) ---
            fidInfo = fopen(fullfile(dir_name,[name,'.info']), 'r');
            rawLines = textscan(fidInfo,'%s','Delimiter','\n','Whitespace','');
            fclose(fidInfo);
            rawLines = rawLines{1};
            
            params = containers.Map;
            for j = 1:numel(rawLines)
                line = strtrim(rawLines{j});
                if isempty(line) || startsWith(line,'%'), continue; end
                tok = regexp(line,'^([^:]+):\s*([-\d\.eE]+)$','tokens','once');
                if ~isempty(tok)
                    params(strtrim(tok{1})) = str2double(tok{2});
                end
            end
            % disp(params.keys)   % <- uncomment once to see available keys

            weight=zeros(4,nofSamples);
            height=zeros(4,nofSamples);
            temp=zeros(4,nofSamples);

            % Use exact key names from your .info sample:
            %   UnitX_Scale_Amp, UnitX_Scale_Offset
            %   UnitX_US_Amp,    UnitX_US_Offset
            %   UnitX_Temp_Amp,  UnitX_Temp_Offset
            %   UnitX_Start_Weight, UnitX_Start_Height
            
            for IDX_unit=1:4
                scaleAmp  = params(sprintf('Unit%d_Scale_Amp',    IDX_unit));
                scaleOff  = params(sprintf('Unit%d_Scale_Offset', IDX_unit));
                startW    = params(sprintf('Unit%d_Start_Weight', IDX_unit));
            
                usAmp     = params(sprintf('Unit%d_US_Amp',       IDX_unit));
                usOff     = params(sprintf('Unit%d_US_Offset',    IDX_unit));
                startH    = params(sprintf('Unit%d_Start_Height', IDX_unit));
            
                tAmp      = params(sprintf('Unit%d_Temp_Amp',     IDX_unit));
                tOff      = params(sprintf('Unit%d_Temp_Offset',  IDX_unit));
            
                weight(IDX_unit,:) = Data1(1+(IDX_unit-1)*3,:).*scaleAmp + scaleOff - startW;
                height(IDX_unit,:) = Data1(2+(IDX_unit-1)*3,:).*usAmp    + usOff   - startH;
                temp(IDX_unit,:)   = Data1(3+(IDX_unit-1)*3,:).*tAmp     + tOff;
            end

            rate=1;
            weight_ds=zeros(4,ceil(nofSamples/rate));
            height_ds=zeros(4,ceil(nofSamples/rate));
            temp_ds=zeros(4,ceil(nofSamples/rate));
            for IDX_unit=1:4
            weight_ds(IDX_unit,:) = weight(IDX_unit,:);
            height_ds(IDX_unit,:) = height(IDX_unit,:);
            temp_ds(IDX_unit,:) = temp(IDX_unit,:);
            end

            weight_comp=zeros(4,ceil(nofSamples/rate));
            A_tube = pi*(params('Metal_tube_diam_out(m)')/2)^2 ...
                   - pi*(params('Metal_tube_diam_in(m)')/2)^2;
            
            A_hose = pi*(params('Rubber_tube_diam_out(m)')/2)^2 ...
                   - pi*(params('Rubber_tube_diam_in(m)')/2)^2;
            
            A = A_tube + A_hose;
            
            for IDX_unit=1:4
                Density(IDX_unit) = params(sprintf('Density_Unit%d',IDX_unit));
                weight_comp(IDX_unit,:) = A .* height_ds(IDX_unit,:) .* Density(IDX_unit);
            end

            % --- Volume ---
            Volume = weight_ds - weight_comp;
            Volume = ones(size(Volume)).*max(Volume) - Volume;  % your invert
            
            % NEW: baseline to zero per channel and clip tiny negatives
            Volume = Volume - Volume(:,1);     % start at 0 ml at t=0
            Volume(Volume < 0) = 0;
            
            fprintf('\n[%s]\nV0(min,max)=[%.4g, %.4g]  Vol(min,max)=[%.4g, %.4g]\n', ...
                File, min((weight_ds - weight_comp),[],'all'), max((weight_ds - weight_comp),[],'all'), ...
                      min(Volume,[],'all'),                 max(Volume,[],'all'));


%% refill
if refill_table
            L=zeros(1,length(filenames_refill));
            for b=1:length(filenames_refill)
                
                tf = regexpi(File,filenames_refill{b});
                if ~isempty(tf)
                    L(b)=1;
                end
                
            end
            t_r=[];
            if sum(L)>0
                hol=find(L);
                for b=1:sum(L)
                    t_r(b)=interval_refill(hol(b),1);
                    
                    % deal with the refill
                    Volume((interval_refill(hol(b),2)*Fs):end) = Volume((interval_refill(hol(b),2)*Fs):end) + ...
                        (Volume(interval_refill(hol(b),1)*Fs)-Volume(interval_refill(hol(b),2)*Fs)); %
                    
                    % in small interval around the refill make the signal
                    % equal to prevent spurious sucks
                     h1=(interval_refill(hol(b),1)*Fs)-500:(interval_refill(hol(b),1)*Fs);
                     h2=(interval_refill(hol(b),2)*Fs):(interval_refill(hol(b),2)*Fs)+500;
                    Volume([h1,h2])=mean(Volume([h1,h2]));
                        Volume((interval_refill(hol(b),1)*Fs):(interval_refill(hol(b),2)*Fs)) = mean(Volume([h1,h2]));
                end
                
            end
end
%% end refill

%% filter
       
            data_filt=bst_bandpass_fft(Volume(1,:)*1000, 1000, 0,FOI, 1, 1);
            speed=gradient(data_filt)' .*(Fs/rate);
            amount_consumed=data_filt;
%% interpolation
            [ind, ~] = find(abs(speed) > artefact_speed_threshold);
            vege = ind(:);  % column vector of sample indices
            
            interp_speed = speed;

            if cntr_subj == 1 && want_figures==1



                if want_figures==1 && ~did_debug_plot
                    figure;
                    subplot(2,1,1); plot(t, Volume(1,:)*1000); ylabel('ml');
                    title(['Volume (preview) — ', File],'Interpreter','none');
                    subplot(2,1,2);
                    sp = gradient(bst_bandpass_fft(Volume(1,:)*1000, Fs, 0, FOI, 1, 1))' .* Fs;
                    plot(t, sp); ylabel('ml/s'); xlabel('s'); title('Speed (preview)');
                    did_debug_plot = true;
                end
                
                sp = gradient(bst_bandpass_fft(Volume(1,:)*1000, Fs, 0, FOI, 1, 1))' .* Fs;
                plot(t, sp); ylabel('ml/s'); xlabel('s'); title('Speed (pre-art)');
            end

            interp_amount_consumed = amount_consumed;
            
            while ~isempty(vege)
                i0 = vege(1);
            
                % Clamp initial window to valid indices
                ind_start = max(1, i0 - win);
                ind_end   = min(N, i0 + win);
                ind_interval = ind_start:ind_end;
            
                % Grow window while it captures more artefacts, but always clamp
                last_len = 0;
                while true
                    [~, ia] = intersect(vege, ind_interval);
                    if numel(ia) > last_len
                        last_len = numel(ia);
                        ind_end = min(N, ind_end + win);
                        ind_interval = ind_start:ind_end;
                    else
                        break
                    end
                end
            
                % Remove those artefacts from the list
                vege(ia) = [];
            
                % Build time vectors safely
                time_interp = t(ind_interval);
                ok_idx = setdiff(1:N, ind_interval);
                time_ok = t(ok_idx);
            
                % Interpolate speed
                good_speed = speed(ok_idx);
                signal_interp = interp1(time_ok, good_speed, time_interp, 'linear', 'extrap');
                interp_speed(ind_interval) = signal_interp;
            
                % Interpolate amount_consumed
                good_amount = amount_consumed(ok_idx);
                signal_interp = interp1(time_ok, good_amount, time_interp, 'linear', 'extrap');
                interp_amount_consumed(ind_interval) = signal_interp;
            end
            
            interp_speed(interp_speed < min_amplitude_of_suck) = 0;
interp_speed(end) = 0;
            
            
            %% find zeros, make zeros
if refill_table          
           if sum(L)>0
                hol=find(L);
                for b=1:sum(L) %here the accidently arising sucks (due to noise) were deleted
                     h1=(interval_refill(hol(b),1)*Fs)-500:(interval_refill(hol(b),1)*Fs);
                     h2=(interval_refill(hol(b),2)*Fs):(interval_refill(hol(b),2)*Fs)+500;
                     interp_speed(h1(1):h2(end))=0;
                     ind=0;
                     cntr=1;
                     while ind==0
                         if interp_speed(h1(1)-cntr)==0
                             ind=1;
                         else
                             interp_speed(h1(1)-cntr)=0;
                             cntr=cntr+1;
                         end
                     end
                     
                     ind=0;
                     cntr=1;
                     while ind==0 
                         if interp_speed(h2(end)+cntr)==0 | h2(end)+cntr==length(interp_speed)-1
                             ind=1;
                         else
                             interp_speed(h2(end)+cntr)=0;
                             cntr=cntr+1;
                         end
                     end
                     
                end
                %%
               interp_speed=interp_speed';
                hh=0;
                hol2=[];
                for b=1:sum(L) %find closest nonzero speeds around refills
                    
                     h1=(interval_refill(hol(b),1)*Fs-hh)-500:(interval_refill(hol(b),1)*Fs-hh);
                     h2=(interval_refill(hol(b),2)*Fs-hh):(interval_refill(hol(b),2)*Fs-hh)+500;
                     ind=0;
                     cntr=1;
                     while ind==0
                         if interp_speed(h1(1)-cntr)~=0 | h1(1)-cntr==1
                             ind=1;
                             hol2(b,1)=h1(1)-cntr;
                         else
                             cntr=cntr+1;
                         end
                     end
                     
                     ind=0;
                     cntr=1;
                     while ind==0
                         if interp_speed(h2(end)+cntr)~=0 | h2(end)+cntr==length(interp_speed)-1
                             ind=1;
                             hol2(b,2)=h2(end)+cntr;
                         else
                             cntr=cntr+1;
                         end
                     end
                     
                      % delete the pauses around refill sucks and add 3000
                      % datapoints
                     interp_speed=[(interp_speed(1:hol2(b,1))),(interp_speed(hol2(b,2):end))];
                     interp_amount_consumed=[(interp_amount_consumed(1:hol2(b,1))),(interp_amount_consumed(hol2(b,2):end))];
                     hh=hh+(hol2(b,2)-hol2(b,1));
                end
                
           end
end
            t=(1:length(interp_amount_consumed))/Fs;
            
            
%%

        [pks,locs,w,proms] = findpeaks(interp_speed,t,'MinPeakProminence',...
             min_amplitude_of_suck,'MinPeakDistance',min_peak_distance,'Annotate','peaks'); %,'MinPeakWidth',min_peak_width
        
        [C,ai,ib] = intersect(t,locs);
        
        
        
        Basis_ind=zeros(2,length(ai));
        
        for a=1:length(ai)
            
            ind=0;
            cntr=1;
            while ind==0
                
                if interp_speed(ai(a)+cntr)>0
                    cntr=cntr+1;
                else
                    Basis_ind(1,a)=ai(a)+cntr;
                    ind=1;
                end
                
                if ismember(ai(a)+cntr,ai)
                    
                    [M,I]=min(interp_speed(ai(a):ai(a)+cntr));
                    Basis_ind(1,a)=ai(a)+I-1;
                    ind=1;
                    
                end
                
            end
            
            
            ind=0;
            cntr=1;
            while ind==0
                
                if interp_speed(ai(a)-cntr)>0
                    cntr=cntr+1;
                else
                    Basis_ind(2,a)=ai(a)-cntr;
                    ind=1;
                end
                
                if ismember(ai(a)-cntr,ai)

                    [M,I]=min(interp_speed(ai(a)-cntr:ai(a)));
                    Basis_ind(2,a)=ai(a)-length(interp_speed(ai(a)-cntr:ai(a)))+I;
                    ind=1;
                    
                end
                
            end
        end
        
              %% consumed for sucks
        suck_consumed=(interp_amount_consumed(Basis_ind(1,:)-1)-interp_amount_consumed(Basis_ind(2,:)+1))';
        

        
        %%
        T_suck=(t(Basis_ind(1,:)-1)-t(Basis_ind(2,:)+1))';
        small_suck=find(suck_consumed<min_suck_amount);
        suck_consumed(small_suck)=[];
        T_suck(small_suck)=[];
        Vmax=pks;
        Vmax(small_suck)=[];
        TimeVamx=locs;
        locs(small_suck)=[];
        TimeVamx(small_suck)=[];
        Basis_ind(:,small_suck)=[];
        
        % ---- SAFETY GUARD ----
        no_peaks = isempty(locs);
        no_basis = isempty(Basis_ind) || size(Basis_ind,2) == 0;
        if no_peaks || no_basis
            warning('No valid sucks detected in %s (peaks=%d, basis cols=%d). Skipping.', ...
                    File, numel(locs), size(Basis_ind,2));
            continue  % skip this file safely
        end
        % ----------------------

        start_end_point = [min(Basis_ind(2,:)), max(Basis_ind(1,:))];
        
        
        
        %% intersuck interval
        suck_interval=zeros(1,size(Basis_ind,2)-1);
        
        for a=1:size(Basis_ind,2)-1
            
            suck_interval(a)=t(Basis_ind(2,a+1))-t(Basis_ind(1,a));
            
        end
        O(cntr_subj,1)=sum(suck_consumed);
        O(cntr_subj,2)=interp_amount_consumed(end);
        O(cntr_subj,3)=sum(suck_interval)./max(t);


        timestamp=t;
        %%
        if want_figures==1   

                f=figure;

                % Compute a pre-artifact speed for scale estimation (use your signal name if different)
                speed_for_scale = interp_speed;                 % if you already have it
                pos = speed_for_scale(speed_for_scale > 0);
                if isempty(pos)
                    min_prom = 0.2;                             % fallback
                else
                    p95 = prctile(pos, 95);
                    % Adaptive prominence: 20% of p95 but clamped to [0.1, 0.8] ml/s
                    min_prom = max(0.1, min(0.8, 0.20 * p95));
                end
                % Adaptive minimal distance: keep your default but allow down to 0.3 s
                min_dist = max(0.3, min_peak_distance);
                
                [pks,locs,w,proms] = findpeaks( ...
                    interp_speed, t, ...
                    'MinPeakProminence', min_prom, ...
                    'MinPeakDistance',   min_dist, ...
                    'Annotate','peaks');

                findpeaks(interp_speed,timestamp,'MinPeakProminence',...
                min_amplitude_of_suck,'MinPeakDistance',min_peak_distance,'Annotate','peaks'); %'MinPeakWidth',min_peak_width,

                hold on

                plot(timestamp(Basis_ind(1,:)),interp_speed(Basis_ind(1,:)),'r*')

                hold on

                plot(timestamp(Basis_ind(2,:)),interp_speed(Basis_ind(2,:)),'g*')

                xlabel('time (s)')
                ylabel('speed (ml/s)')
                title (File)

                saveas(f,[File,'.pdf'])
                close(f)

        end

        Suck_all_time=[Suck_all_time suck_interval];
    end
end

%% IBI


%% make the loge
Suck_all_time2=Suck_all_time;
Suck_all_time2(Suck_all_time2==0)=[];
LOG_suck=log(Suck_all_time2)';


%% find the number of GM with the 
for index_j=1:10
options = statset('MaxIter', 1000);  % Increase maximum iterations

Q1 = prctile(LOG_suck, 25);
Q3 = prctile(LOG_suck, 75);
IQR_value = Q3 - Q1;

% Remove outliers that are beyond 1.5 times the IQR
LOG_suck_filtered = LOG_suck(all(LOG_suck > (Q1 - 1.5 * IQR_value) & LOG_suck < (Q3 + 1.5 * IQR_value), 2), :);

number_of_comp = 2;  % Set the desired number of components

gm = fitgmdist(LOG_suck_filtered, number_of_comp, 'Options', options, 'RegularizationValue', 1e-2);


%options = statset('MaxIter',100000,'Display','final');
%gm = fitgmdist(LOG_suck,index_j,'Options',options);
cost(index_j)=gm.AIC;
end

[AIC, number_of_comp]=min(cost);

gm = fitgmdist(LOG_suck,number_of_comp,'Options',options);
mu = gm.mu;
sd = squeeze(gm.Sigma);
heigh=1;%[10*gm.ComponentProportion(1),10*gm.ComponentProportion(2),10*gm.ComponentProportion(3),60,60];

% --- Histogram (PDF) with GMM components overlaid ---
% Make histogram centers and component pdfs for later use
[counts, edges] = histcounts(LOG_suck, histogram_bins, 'Normalization','pdf');
centers = edges(1:end-1) + diff(edges)/2;   % replaces old "x"

L = zeros(gm.NumComponents, numel(centers));
for k = 1:gm.NumComponents
    mu_k  = gm.mu(k);
    var_k = gm.Sigma(:,:,k);          % 1x1xK in 1-D
    sd_k  = sqrt(var_k);
    L(k,:) = gm.ComponentProportion(k) * normpdf(centers, mu_k, sd_k);
end

figure
bar(centers, counts, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', [0.5 0.5 0.5]); hold on

% Overlay each Gaussian component correctly
for k = 1:gm.NumComponents
    mu_k = gm.mu(k);                       % 1-by-1
    var_k = gm.Sigma(:,:,k);               % 1-by-1-by-K -> 1-by-1
    sd_k = sqrt(var_k);                    % scalar
    yk = gm.ComponentProportion(k) * normpdf(centers, mu_k, sd_k);
    plot(centers, yk, 'k', 'LineWidth', 3);
end

xlabel('Log_e(inter-suck interval length)','FontSize',25)
ylabel('Probability density','FontSize',25)
set(gca,'fontsize',15)
xlim([min(centers) max(centers)])
title('GMM fit')

        
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
title('GM')

figure
plot (cost)
xlabel ('Number of components','FontSize',25)
ylabel ('AIC','FontSize',25)
set(gca,'fontsize',15)
title('AIC')

disp('AIC and NUMBER OF COMPONENTS:')
disp(['The minimal AIC is ', num2str(min(cost))])
disp(['The number of components are ', num2str(number_of_comp)])

disp('CURVE INTERSECTIONS:')

for c1 = 1:number_of_comp
    for c2 = c1+1:number_of_comp
        d = L(c1,:) - L(c2,:);                     % difference of the two PDFs
        signChange = d(1:end-1) .* d(2:end) <= 0;  % look for zero crossings
        idx = find(signChange);

        for m = idx
            x1 = x(m);   x2 = x(m+1);
            y1 = d(m);   y2 = d(m+1);
            if y2 ~= y1
                xstar = x1 + (0 - y1) * (x2 - x1) / (y2 - y1);  % linear interp
            else
                xstar = x1;  % degenerate case
            end
            fprintf('Curve_%d intersect Curve_%d at %.6f second\n', c1, c2, exp(xstar));
        end
    end
end





