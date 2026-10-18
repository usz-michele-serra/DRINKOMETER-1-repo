clear all
close all

%% give initial thresholds - here you can modify the parameters
Fs=1000; %Hz; Sampling frequency of the signal
FOI=0.5; %Hz; cut off frequency for filtering (in Hz)
artefact_speed_threshold=20.0; %ml/sec; Speed signal is considered as an artefact above the threshold (absolute value). The detected interval is deleted and interpolated linearly.
min_amplitude_of_suck=0.5; %ml/sec; speed is set to zero below this threshold
min_peak_distance=1.0; %sec; the minimal distance between peaks sucks (prevent to detect noine as suck) 
min_suck_amount=1.0; %ml; below this threshold the sucks considered as noise and not included in the analysis.


min_burst_interval=2.05; %sec; set IBI


refill_pause=min_burst_interval; % the refill intervall is replace by the IBI
name_of_the_refill_table='refill_table.xlsx';
refill_table=0; %1=yes, 0=no
want_figures=0; %1=yes, 0=no



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

            fidInfo=fopen([dir_name,'/', name,'.info']);
            frewind(fidInfo)
            C=textscan(fidInfo,'%s %f' , 'CommentStyle',{'%'},'Delimiter',{':', '\n'});
            fclose(fidInfo);

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
            t=rate/Fs:rate/Fs:nofSamples/(Fs/rate)/rate;
            Volume=Volume(1,:);
            V=Volume(1,:);

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
            [ind,val]=find(abs(speed)>artefact_speed_threshold);
            vege=ind;
            interp_speed=speed;
            interp_amount_consumed=amount_consumed;
            while ~isempty(vege)
                %% find close artefacts together
                ind_interval=vege(1)-2000:vege(1)+2000; %mekkora intervallum legyen, ahol a közöen elõforduló mûtermékeket interpolálja?
                indi=0;
                ia_e=1;
                while indi==0
                [C,ia,ib] = intersect(vege,ind_interval); %ia=vege; ib=ind_interval
                if (length(ia)-length(ia_e))>0
                    ia_e=ia;
                    ind_interval=(ind_interval(1)) : (ind_interval(end)+2000);
                    
                else
                    indi=1;
                    
                end
                end
                vege(ia)=[];

                time_interp=t(ind_interval);
                time_ok=setdiff(t,time_interp);

                %speed
                good_speed=speed;
                good_speed(ind_interval)=[];
                signal_interp = interp1(time_ok,good_speed,time_interp,'linear');
                interp_speed(ind_interval)=signal_interp;

                %amount_consumed
                good_amount_consumed=amount_consumed;
                good_amount_consumed(ind_interval)=[];
                signal_interp = interp1(time_ok,good_amount_consumed,time_interp,'linear');
                interp_amount_consumed(ind_interval)=signal_interp;

            end
            interp_speed(interp_speed<min_amplitude_of_suck)=0;
            interp_speed(end)=0;
            
            
            %% find zeros, make zeros
%             error('here')
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
                     interp_speed=[(interp_speed(1:hol2(b,1))),(zeros(1,refill_pause*Fs)),(interp_speed(hol2(b,2):end))];
                     interp_amount_consumed=[(interp_amount_consumed(1:hol2(b,1))),(ones(1,refill_pause*Fs)*interp_amount_consumed(hol2(b,1))),(interp_amount_consumed(hol2(b,2):end))];
                     hh=hh+(hol2(b,2)-hol2(b,1))-refill_pause*Fs;
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
        
        start_end_point(1)=min(Basis_ind(2,:));
        start_end_point(2)=max(Basis_ind(1,:));
        
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





%%  all table
        
%total volume consumed
OUT_all.total_volume_consumed(cntr_subj,1)=interp_amount_consumed(start_end_point(2)) - interp_amount_consumed(start_end_point(1));
        
%total time
OUT_all.total_time(cntr_subj,1)=timestamp(start_end_point(2)) - timestamp(start_end_point(1));
        
%avg speed
OUT_all.avg_speed(cntr_subj,1)=(interp_amount_consumed(start_end_point(2)) - interp_amount_consumed(start_end_point(1)))/...
(timestamp(start_end_point(2)) - timestamp(start_end_point(1)));
        
%Vmax
OUT_all.Vmax_all(cntr_subj,:)= mean_sd_se_median_iqr( Vmax );
        
%number of suck
OUT_all.suck_number(cntr_subj,1)= length(suck_consumed);
                
%suck consumed
OUT_all.suck_consumed(cntr_subj,:)= mean_sd_se_median_iqr(suck_consumed);
     
% inter_sucking interval
OUT_all.inter_suck_interval(cntr_subj,:)= mean_sd_se_median_iqr(suck_interval(suck_interval<min_burst_interval));
        
% interburst interval
OUT_all.inter_burst_interval(cntr_subj,:)= mean_sd_se_median_iqr(suck_interval(suck_interval>min_burst_interval));
        
where_burst=find(suck_interval>min_burst_interval);
start_burst=Basis_ind(2,[1,where_burst+1]);
end_burst=Basis_ind(1,[where_burst, length(Basis_ind)]);
        
% number of suck per burst
sucks_in_bursts=[where_burst,length(suck_interval)+1]-[0,where_burst];
OUT_all.number_of_suck_per_burst(cntr_subj,:)= mean_sd_se_median_iqr(sucks_in_bursts);
        
% suck rate
OUT_all.suck_rate(cntr_subj,:)= mean_sd_se_median_iqr(suck_consumed./T_suck);
        
% burst intervals
burst_intervals=timestamp(end_burst)-timestamp(start_burst);
OUT_all.burst_intervals(cntr_subj,:)= mean_sd_se_median_iqr(burst_intervals);
        
% burst consumed
burst_consumed=interp_amount_consumed(end_burst)-interp_amount_consumed(start_burst);
OUT_all.burst_consumed(cntr_subj,:)= mean_sd_se_median_iqr(burst_consumed);
        
% burst speed
OUT_all.burst_speed(cntr_subj,:)= mean_sd_se_median_iqr(burst_consumed./burst_intervals);
        
%number of burst
OUT_all.burst_number(cntr_subj,1)= length(where_burst)+1;

%suck duration
OUT_all.suck_duration(cntr_subj,:)= mean_sd_se_median_iqr(T_suck);

        
%% first 15 sec
[M,index_end_15]=min(abs(timestamp-(15+timestamp(start_end_point(1)))));
        
%total volume consumed - 15
OUT_15.total_volume_consumed(cntr_subj,1)=interp_amount_consumed(index_end_15) - interp_amount_consumed(start_end_point(1));
        
%avg speed -15
OUT_15.avg_speed(cntr_subj,1)=(interp_amount_consumed(index_end_15) - interp_amount_consumed(start_end_point(1)))/...
(timestamp(index_end_15) - timestamp(start_end_point(1)));
        
%Vmax -15
Vmax_15=Vmax(TimeVamx<=timestamp(index_end_15));
OUT_15.Vmax_all(cntr_subj,:)= mean_sd_se_median_iqr( Vmax_15 );
        
%number of suck - 15
OUT_15.suck_number(cntr_subj,1)= length(Vmax_15);
        
%suck consumed - 15
OUT_15.suck_consumed(cntr_subj,:)= mean_sd_se_median_iqr( suck_consumed(1:length(Vmax_15)) );
        
%suck rate - 15
OUT_15.suck_rate(cntr_subj,:)= mean_sd_se_median_iqr( suck_consumed(1:length(Vmax_15))./ T_suck(1:length(Vmax_15)));
     
% inter_sucking interval - 15
suck_interval_15=suck_interval(1:length(Vmax_15)-1);
OUT_15.inter_suck_interval(cntr_subj,:)= mean_sd_se_median_iqr(suck_interval_15(suck_interval_15<min_burst_interval));

%suck duration
OUT_15.suck_duration(cntr_subj,:)= mean_sd_se_median_iqr(T_suck(1:length(Vmax_15)));


%% first 30 sec
[M,index_end_30]=min(abs(timestamp-(30+timestamp(start_end_point(1)))));
        
%total volume consumed - 30
OUT_30.total_volume_consumed(cntr_subj,1)=interp_amount_consumed(index_end_30) - interp_amount_consumed(start_end_point(1));
        
%avg speed -30
OUT_30.avg_speed(cntr_subj,1)=(interp_amount_consumed(index_end_30) - interp_amount_consumed(start_end_point(1)))/...
(timestamp(index_end_30) - timestamp(start_end_point(1)));
        
%Vmax -30
Vmax_30=Vmax(TimeVamx<=timestamp(index_end_30));
OUT_30.Vmax_all(cntr_subj,:)= mean_sd_se_median_iqr( Vmax_30 );
        
%number of suck - 30
OUT_30.suck_number(cntr_subj,1)= length(Vmax_30);
        
%suck consumed - 30
OUT_30.suck_consumed(cntr_subj,:)= mean_sd_se_median_iqr( suck_consumed(1:length(Vmax_30)) );
        
%suck rate - 30
OUT_30.suck_rate(cntr_subj,:)= mean_sd_se_median_iqr( suck_consumed(1:length(Vmax_30))./ T_suck(1:length(Vmax_30)));
     
% inter_sucking interval - 30
suck_interval_30=suck_interval(1:length(Vmax_30)-1);
OUT_30.inter_suck_interval(cntr_subj,:)= mean_sd_se_median_iqr(suck_interval_30(suck_interval_30<min_burst_interval));

%suck duration - 30
OUT_30.suck_duration(cntr_subj,:)= mean_sd_se_median_iqr(T_suck(1:length(Vmax_30)));

%% first 60 sec
[M,index_end_60]=min(abs(timestamp-(60+timestamp(start_end_point(1)))));
        
%total volume consumed - 60
OUT_60.total_volume_consumed(cntr_subj,1)=interp_amount_consumed(index_end_60) - interp_amount_consumed(start_end_point(1));
        
%avg speed -60
OUT_60.avg_speed(cntr_subj,1)=(interp_amount_consumed(index_end_60) - interp_amount_consumed(start_end_point(1)))/...
(timestamp(index_end_60) - timestamp(start_end_point(1)));
        
%Vmax -60
Vmax_60=Vmax(TimeVamx<=timestamp(index_end_60));
OUT_60.Vmax_all(cntr_subj,:)= mean_sd_se_median_iqr( Vmax_60 );
        
%number of suck - 60
OUT_60.suck_number(cntr_subj,1)= length(Vmax_60);
        
%suck consumed - 60
OUT_60.suck_consumed(cntr_subj,:)= mean_sd_se_median_iqr( suck_consumed(1:length(Vmax_60)) );
        
%suck rate - 60
OUT_60.suck_rate(cntr_subj,:)= mean_sd_se_median_iqr( suck_consumed(1:length(Vmax_60))./ T_suck(1:length(Vmax_60)));
     
% inter_sucking interval - 60
suck_interval_60=suck_interval(1:length(Vmax_60)-1);
OUT_60.inter_suck_interval(cntr_subj,:)= mean_sd_se_median_iqr(suck_interval_60(suck_interval_60<min_burst_interval));

%suck duration - 60
OUT_60.suck_duration(cntr_subj,:)= mean_sd_se_median_iqr(T_suck(1:length(Vmax_60)));
        
%% first burst
if ~isempty(where_burst)
    end_of_first_burst=Basis_ind(1,where_burst(1));
else
    end_of_first_burst=start_end_point(2);
end
        
%total volume consumed - 1st burst
OUT_firstB.total_volume_consumed(cntr_subj,1)=interp_amount_consumed(end_of_first_burst) - interp_amount_consumed(start_end_point(1));
        
%total time - 1st burst
OUT_firstB.total_time(cntr_subj,1)=timestamp(end_of_first_burst) - timestamp(start_end_point(1));
        
%avg speed - 1st burst
OUT_firstB.avg_speed(cntr_subj,1)=(interp_amount_consumed(end_of_first_burst) - interp_amount_consumed(start_end_point(1)))/...
(timestamp(end_of_first_burst) - timestamp(start_end_point(1)));
    
%Vmax - 1st burst
Vmax_firstB=Vmax(TimeVamx<=timestamp(end_of_first_burst));
OUT_firstB.Vmax_all(cntr_subj,:)= mean_sd_se_median_iqr( Vmax_firstB );

%number of suck - 1st burst
OUT_firstB.suck_number(cntr_subj,1)= length(Vmax_firstB);
        
%suck consumed - 1st burst
OUT_firstB.suck_consumed(cntr_subj,:)= mean_sd_se_median_iqr( suck_consumed(1:length(Vmax_firstB)) );
        
%suck rate - 1st burst
OUT_firstB.suck_rate(cntr_subj,:)= mean_sd_se_median_iqr( suck_consumed(1:length(Vmax_firstB)) ./T_suck(1:length(Vmax_firstB)));
%inter_sucking interval - 1st burst
suck_interval_firstB=suck_interval(1:length(Vmax_firstB)-1);
OUT_firstB.inter_suck_interval(cntr_subj,:)= mean_sd_se_median_iqr(suck_interval_firstB);

%suck duration - 1st burst
OUT_firstB.suck_duration(cntr_subj,:)= mean_sd_se_median_iqr(T_suck(1:length(Vmax_firstB)));


%% all burst table

        where_burst=find(suck_interval>min_burst_interval);
        
        start_burst=Basis_ind(2,[1,where_burst+1]);
        end_burst=Basis_ind(1,[where_burst, length(Basis_ind)]);
        
        for a=1:20
            
            if a > size(start_burst)
                OUT_Bursts.total_volume_consumed(cntr_subj,1,a)=NaN;
                OUT_Bursts.total_time(cntr_subj,1,a)=NaN;
                OUT_Bursts.avg_speed(cntr_subj,1,a)=NaN;
                OUT_Bursts.Vmax(cntr_subj,:,a)= [NaN NaN NaN NaN NaN];
                OUT_Bursts.suck_number(cntr_subj,1,a)= NaN;
                OUT_Bursts.suck_consumed(cntr_subj,:,a)= [NaN NaN NaN NaN NaN];
                OUT_Bursts.suck_rate(cntr_subj,:,a)= [NaN NaN NaN NaN NaN];
                OUT_Bursts.inter_suck_interval(cntr_subj,:,a)= [NaN NaN NaN NaN NaN];
                OUT_Bursts.suck_duration(cntr_subj,:,a)= [NaN NaN NaN NaN NaN];
            else
            

            %% it counts 20 burst
            
            OUT_Bursts.total_volume_consumed(cntr_subj,1,a)=interp_amount_consumed(end_burst(a)) - interp_amount_consumed(start_burst(a));
        
            %total time
            OUT_Bursts.total_time(cntr_subj,1,a)=timestamp(end_burst(a)) - timestamp(start_burst(a));
        
            %avg speed
            OUT_Bursts.avg_speed(cntr_subj,1,a)=OUT_Bursts.total_volume_consumed(cntr_subj,1,a)/OUT_Bursts.total_time(cntr_subj,1,a);
    
            %Vmax
            Vmax_Burst=Vmax(TimeVamx>=timestamp(start_burst(a)) & TimeVamx<=timestamp(end_burst(a)));
            OUT_Bursts.Vmax(cntr_subj,:,a)= mean_sd_se_median_iqr( Vmax_Burst );
        
            %number of suck
            OUT_Bursts.suck_number(cntr_subj,1,a)= length(Vmax_Burst);
        
            %suck consumed 
            OUT_Bursts.suck_consumed(cntr_subj,:,a)= mean_sd_se_median_iqr(...
                suck_consumed(TimeVamx>=timestamp(start_burst(a)) & TimeVamx<=timestamp(end_burst(a))));
            
            %suck rate 
            OUT_Bursts.suck_rate(cntr_subj,:,a)= mean_sd_se_median_iqr(...
                suck_consumed(TimeVamx>=timestamp(start_burst(a)) & TimeVamx<=timestamp(end_burst(a)))./...
                T_suck(TimeVamx>=timestamp(start_burst(a)) & TimeVamx<=timestamp(end_burst(a))));
        
            %inter_sucking interva
            hol_b=find(TimeVamx>=timestamp(start_burst(a)) & TimeVamx<=timestamp(end_burst(a)));
            hol_b(end)=[];
            suck_interval_Burst=suck_interval(hol_b);
            OUT_Bursts.inter_suck_interval(cntr_subj,:,a)= mean_sd_se_median_iqr(suck_interval_Burst);
            
            %suck duration
            OUT_Bursts.suck_duration(cntr_subj,:,a)= mean_sd_se_median_iqr(T_suck(TimeVamx>=timestamp(start_burst(a)) & TimeVamx<=timestamp(end_burst(a))));
            
            end
            
        end
        
        
%% first half - second half
        where_burst=find(suck_interval>min_burst_interval);
        
        start_burst=Basis_ind(2,[1,where_burst+1]);
        end_burst=Basis_ind(1,[where_burst, length(Basis_ind)]);

        half_of_bursts=floor(length(start_burst)/2);
        if half_of_bursts>=1
%             cntr_subj=cntr_subj+1;
%             FileName{cntr_subj}=File;
            First_half(1)=1;
            First_half(2)=half_of_bursts;
            Second_half(1)=length(start_burst)-(half_of_bursts-1);
            Second_half(2)=length(start_burst);
%             error('iit')
            
            t_volume=[];
            t_time=[];
            t_speed=[];
            t_vmax=[];
            t_suck_number=[];
            t_suck_consumed=[];
            t_suck_interval=[];
            t_burst=[];
            t_suck_duration=[];
            
            for a=1:length(start_burst)

                t_volume(a)=interp_amount_consumed(end_burst(a)) - interp_amount_consumed(start_burst(a));
                t_time(a)=timestamp(end_burst(a)) - timestamp(start_burst(a));
                t_speed(a)=t_volume(a)/t_time(a);
                t_vmax(a)=mean(Vmax(TimeVamx>=timestamp(start_burst(a)) & TimeVamx<=timestamp(end_burst(a))));
                t_suck_number(a)=length(Vmax(TimeVamx>=timestamp(start_burst(a)) & TimeVamx<=timestamp(end_burst(a))));
                t_suck_consumed(a)=mean(suck_consumed(TimeVamx>=timestamp(start_burst(a)) & TimeVamx<=timestamp(end_burst(a))));
                t_suck_rate(a)=t_suck_consumed(a)./mean(T_suck(TimeVamx>=timestamp(start_burst(a)) & TimeVamx<=timestamp(end_burst(a))));
                hol_b=find(TimeVamx>=timestamp(start_burst(a)) & TimeVamx<=timestamp(end_burst(a)));
                hol_b(end)=[];
                t_suck_interval(a)=mean(suck_interval(hol_b));
                t_suck_duration(a)=mean(T_suck(TimeVamx>=timestamp(start_burst(a)) & TimeVamx<=timestamp(end_burst(a))));
                

            end
            if length(First_half)>1
                for a=1:length(start_burst)-1
                    t_burst(a)=timestamp(start_burst(a+1))-timestamp(end_burst(a));
                end
                half_of_ibi=floor(length(t_burst)/2);
                First_ibi(1)=1;
                First_ibi(2)=half_of_ibi;
                Second_ibi(1)=length(t_burst)-(half_of_ibi-1);
                Second_ibi(2)=length(t_burst);
                B_t_1=t_burst(First_ibi(1):First_ibi(2));
                B_t_2=t_burst(Second_ibi(1):Second_ibi(2));
            else
                B_t_1=NaN;
                B_t_2=NaN;
            end
                
            

            %
            FistLast_Bursts.volume(cntr_subj,1,1)=interp_amount_consumed(end_burst(First_half(2))) - interp_amount_consumed(start_burst(First_half(1)));
            FistLast_Bursts.volume(cntr_subj,1,2)=interp_amount_consumed(end_burst(Second_half(2))) - interp_amount_consumed(start_burst(Second_half(1)));

            FistLast_Bursts.time(cntr_subj,1,1)=timestamp(end_burst(First_half(2))) - timestamp(start_burst(First_half(1)));
            FistLast_Bursts.time(cntr_subj,1,2)=timestamp(end_burst(Second_half(2))) - timestamp(start_burst(Second_half(1)));

            FistLast_Bursts.speed(cntr_subj,1,1)=FistLast_Bursts.volume(cntr_subj,1,1)/FistLast_Bursts.time(cntr_subj,1,1);
            FistLast_Bursts.speed(cntr_subj,1,2)=FistLast_Bursts.volume(cntr_subj,1,2)/FistLast_Bursts.time(cntr_subj,1,2);
%             error('itt')
            FistLast_Bursts.total_volume_consumed(cntr_subj,:,1)=mean_sd_se_median_iqr(t_volume(First_half(1):First_half(2)));
            FistLast_Bursts.total_volume_consumed(cntr_subj,:,2)=mean_sd_se_median_iqr(t_volume(Second_half(1):Second_half(2)));
        
            %total time
            FistLast_Bursts.total_time(cntr_subj,:,1)=mean_sd_se_median_iqr(t_time(First_half(1):First_half(2)));
            FistLast_Bursts.total_time(cntr_subj,:,2)=mean_sd_se_median_iqr(t_time(Second_half(1):Second_half(2)));

            %avg speed
            FistLast_Bursts.avg_speed(cntr_subj,:,1)=mean_sd_se_median_iqr(t_speed(First_half(1):First_half(2)));
            FistLast_Bursts.avg_speed(cntr_subj,:,2)=mean_sd_se_median_iqr(t_speed(Second_half(1):Second_half(2)));

            %Vmax
            FistLast_Bursts.Vmax(cntr_subj,:,1)=mean_sd_se_median_iqr(t_vmax(First_half(1):First_half(2)));
            FistLast_Bursts.Vmax(cntr_subj,:,2)=mean_sd_se_median_iqr(t_vmax(Second_half(1):Second_half(2)));
        
            %number of suck
            FistLast_Bursts.suck_number(cntr_subj,:,1)=mean_sd_se_median_iqr(t_suck_number(First_half(1):First_half(2)));
            FistLast_Bursts.suck_number(cntr_subj,:,2)=mean_sd_se_median_iqr(t_suck_number(Second_half(1):Second_half(2)));
        
            %suck consumed
            FistLast_Bursts.suck_consumed(cntr_subj,:,1)= mean_sd_se_median_iqr(t_suck_consumed(First_half(1):First_half(2)));
            FistLast_Bursts.suck_consumed(cntr_subj,:,2)=mean_sd_se_median_iqr(t_suck_consumed(Second_half(1):Second_half(2)));
            
            %suck rate
            FistLast_Bursts.suck_rate(cntr_subj,:,1)= mean_sd_se_median_iqr(t_suck_rate(First_half(1):First_half(2)));
            FistLast_Bursts.suck_rate(cntr_subj,:,2)=mean_sd_se_median_iqr(t_suck_rate(Second_half(1):Second_half(2)));
        
            %inter_sucking interval - 1st burst
            F=t_suck_interval(First_half(1):First_half(2));
            F(isnan(F))=[];
            FistLast_Bursts.inter_suck_interval(cntr_subj,:,1)=mean_sd_se_median_iqr(F);
            
            F=t_suck_interval(Second_half(1):Second_half(2));
            F(isnan(F))=[];
            FistLast_Bursts.inter_suck_interval(cntr_subj,:,2)=mean_sd_se_median_iqr(F);
            
            %inter_burst interval
            FistLast_Bursts.inter_burst_interval(cntr_subj,:,1)=mean_sd_se_median_iqr(B_t_1);
            FistLast_Bursts.inter_burst_interval(cntr_subj,:,2)=mean_sd_se_median_iqr(B_t_2);
            
            %suck-duration
            FistLast_Bursts.suck_duration(cntr_subj,:,1)=mean_sd_se_median_iqr(t_suck_duration(First_half(1):First_half(2)));
            FistLast_Bursts.suck_duration(cntr_subj,:,2)=mean_sd_se_median_iqr(t_suck_duration(Second_half(1):Second_half(2)));
            
           %number of suck all
            FistLast_Bursts.suck_number_all(cntr_subj,1,1)=sum(t_suck_number(First_half(1):First_half(2)));
            FistLast_Bursts.suck_number_all(cntr_subj,1,2)=sum(t_suck_number(Second_half(1):Second_half(2)));
        
            
        else
            FistLast_Bursts.volume(cntr_subj,1,1)=0;
            FistLast_Bursts.volume(cntr_subj,1,2)=0;

            FistLast_Bursts.time(cntr_subj,1,1)=0;
            FistLast_Bursts.time(cntr_subj,1,2)=0;

            FistLast_Bursts.speed(cntr_subj,1,1)=0;
            FistLast_Bursts.speed(cntr_subj,1,2)=0;
%             error('itt')
            FistLast_Bursts.total_volume_consumed(cntr_subj,:,1)=zeros(1,5);
            FistLast_Bursts.total_volume_consumed(cntr_subj,:,2)=zeros(1,5);
        
            %total time
            FistLast_Bursts.total_time(cntr_subj,:,1)=zeros(1,5);
            FistLast_Bursts.total_time(cntr_subj,:,2)=zeros(1,5);

            %avg speed
            FistLast_Bursts.avg_speed(cntr_subj,:,1)=zeros(1,5);
            FistLast_Bursts.avg_speed(cntr_subj,:,2)=zeros(1,5);

            %Vmax
            FistLast_Bursts.Vmax(cntr_subj,:,1)=zeros(1,5);
            FistLast_Bursts.Vmax(cntr_subj,:,2)=zeros(1,5);
        
            %number of suck
            FistLast_Bursts.suck_number(cntr_subj,:,1)=zeros(1,5);
            FistLast_Bursts.suck_number(cntr_subj,:,2)=zeros(1,5);
        
            %suck consumed
            FistLast_Bursts.suck_consumed(cntr_subj,:,1)=zeros(1,5);
            FistLast_Bursts.suck_consumed(cntr_subj,:,2)=zeros(1,5);
            
            %suck rate
            FistLast_Bursts.suck_rate(cntr_subj,:,1)= zeros(1,5);
            FistLast_Bursts.suck_rate(cntr_subj,:,2)=zeros(1,5);
        
            %inter_sucking interval - 1st burst
            FistLast_Bursts.inter_suck_interval(cntr_subj,:,1)=zeros(1,5);
            
            FistLast_Bursts.inter_suck_interval(cntr_subj,:,2)=zeros(1,5);
            
            %inter_burst interval
            FistLast_Bursts.inter_burst_interval(cntr_subj,:,1)=zeros(1,5);
            FistLast_Bursts.inter_burst_interval(cntr_subj,:,2)=zeros(1,5);
            
            % suck duration
            FistLast_Bursts.suck_duration(cntr_subj,:,1)=zeros(1,5);
            FistLast_Bursts.suck_duration(cntr_subj,:,2)=zeros(1,5);
            
            % number of suck all
            FistLast_Bursts.suck_number_all(cntr_subj,1,1)=0;
            FistLast_Bursts.suck_number_all(cntr_subj,1,2)=0;

        end
        %%


    end
end

%% write all
Sum_out_all=cat(2,OUT_all.total_volume_consumed, OUT_all.total_time, OUT_all.avg_speed, OUT_all.Vmax_all, OUT_all.suck_number, OUT_all.burst_number,...
    OUT_all.burst_intervals,OUT_all.burst_consumed,OUT_all.burst_speed,OUT_all.suck_consumed,OUT_all.suck_rate, OUT_all.inter_suck_interval,...
OUT_all.inter_burst_interval ,OUT_all.number_of_suck_per_burst,OUT_all.suck_duration);

Label={'total_volume_consumed',...
    'total_drinking_time',...
    'avg_speed',...
    'Vmax_mean','Vmax_SD','Vmax_SE','Vmax_Med','Vmax_IQR',...
    'Number_of_suck'...
    'Number_of_burst'...
    'Burst_time_mean','Burst_time_SD','Burst_time_SE','Burst_time_Med','Burst_time_IQR'...
    'Burst_amount_mean','Burst_amount_SD','Burst_amount_SE','Burst_amount_Med','Burst_amount_IQR'...
    'Burst_speed_mean','Burst_speed_SD','Burst_speed_SE','Burst_speed_Med','Burst_speed_IQR'...
    'Suck_consumed_mean','Suck_consumed_SD','Suck_consumed_SE','Suck_consumed_Med','Suck_consumed_IQR',...
    'Suck_rate_mean','Suck_rate_SD','Suck_rate_SE','Suck_rate_Med','Suck_rate_IQR',...
    'inter_sucking_interval_mean','inter_sucking_interval_SD','inter_sucking_interval_SE','inter_sucking_interval_Med','inter_sucking_interval_IQR',...
    'inter_burst_interval_mean','inter_burst_interval_SD','inter_burst_interval_SE','inter_burst_interval_Med','inter_burst_interval_IQR',...
    'Number_of_suck_per_burst_mean','Number_of_suck_per_burst_SD','Number_of_suck_per_burst_SE','Number_of_suck_per_burst_Med','Number_of_suck_per_burst_IQR',...
    'suck_duration_mean','suck_duration_SD','suck_duration_SE','suck_duration_Med','suck_duration_IQR'};

W=cat(1,Label,num2cell(Sum_out_all));
W=cat(2,(['FileNames', FileName])',W);

writecell(W, [erase(['Full_Consumption','_',num2str(min_amplitude_of_suck),'_',num2str(min_peak_distance),'_',num2str(min_suck_amount),'_',...
    num2str(min_burst_interval),'_',date,],'.'),'.xlsx'])
% 
%% write 15

Sum_out_15=cat(2,OUT_15.total_volume_consumed, OUT_15.avg_speed, OUT_15.Vmax_all, OUT_15.suck_number,...
    OUT_15.suck_consumed, OUT_15.suck_rate, OUT_15.inter_suck_interval,OUT_15.suck_duration);

Label_15={'total_volume_consumed',...
    'avg_speed',...
    'Vmax_mean','Vmax_SD','Vmax_SE','Vmax_Med','Vmax_IQR',...
    'Number_of_suck'...
    'Suck_consumed_mean','Suck_consumed_SD','Suck_consumed_SE','Suck_consumed_Med','Suck_consumed_IQR',...
    'Suck_rate_mean','Suck_rate_SD','Suck_rate_SE','Suck_rate_Med','Suck_rate_IQR',...
    'inter_sucking_interval_mean','inter_sucking_interval_SD','inter_sucking_interval_SE','inter_sucking_interval_Med','inter_sucking_interval_IQR',...
    'suck_duration_mean','suck_duration_SD','suck_duration_SE','suck_duration_Med','suck_duration_IQR'};

W_15=cat(1,Label_15,num2cell(Sum_out_15));
W_15=cat(2,(['FileNames', FileName])',W_15);

writecell(W_15,[erase(['First15sec','_',num2str(min_amplitude_of_suck),'_',num2str(min_peak_distance),'_',num2str(min_suck_amount),'_',...
    num2str(min_burst_interval),'_',date,],'.'),'.xlsx'])
%% write 30

Sum_out_30=cat(2,OUT_30.total_volume_consumed, OUT_30.avg_speed, OUT_30.Vmax_all, OUT_30.suck_number,...
    OUT_30.suck_consumed, OUT_30.suck_rate, OUT_30.inter_suck_interval,OUT_30.suck_duration);

Label_30={'total_volume_consumed',...
    'avg_speed',...
    'Vmax_mean','Vmax_SD','Vmax_SE','Vmax_Med','Vmax_IQR',...
    'Number_of_suck'...
    'Suck_consumed_mean','Suck_consumed_SD','Suck_consumed_SE','Suck_consumed_Med','Suck_consumed_IQR',...
    'Suck_rate_mean','Suck_rate_SD','Suck_rate_SE','Suck_rate_Med','Suck_rate_IQR',...
    'inter_sucking_interval_mean','inter_sucking_interval_SD','inter_sucking_interval_SE','inter_sucking_interval_Med','inter_sucking_interval_IQR',...
    'suck_duration_mean','suck_duration_SD','suck_duration_SE','suck_duration_Med','suck_duration_IQR'};

W_30=cat(1,Label_30,num2cell(Sum_out_30));
W_30=cat(2,(['FileNames', FileName])',W_30);

writecell(W_30,[erase(['First30sec','_',num2str(min_amplitude_of_suck),'_',num2str(min_peak_distance),'_',num2str(min_suck_amount),'_',...
    num2str(min_burst_interval),'_',date,],'.'),'.xlsx'])
%% write 60

Sum_out_60=cat(2,OUT_60.total_volume_consumed, OUT_60.avg_speed, OUT_60.Vmax_all...
    , OUT_60.suck_number, OUT_60.suck_consumed, OUT_60.suck_rate, OUT_60.inter_suck_interval,OUT_60.suck_duration);

Label_60={'total_volume_consumed',...
    'avg_speed',...
    'Vmax_mean','Vmax_SD','Vmax_SE','Vmax_Med','Vmax_IQR',...
    'Number_of_suck'...
    'Suck_consumed_mean','Suck_consumed_SD','Suck_consumed_SE','Suck_consumed_Med','Suck_consumed_IQR',...
    'Suck_rate_mean','Suck_rate_SD','Suck_rate_SE','Suck_rate_Med','Suck_rate_IQR',...
    'inter_sucking_interval_mean','inter_sucking_interval_SD','inter_sucking_interval_SE','inter_sucking_interval_Med','inter_sucking_interval_IQR',...
    'suck_duration_mean','suck_duration_SD','suck_duration_SE','suck_duration_Med','suck_duration_IQR'};

W_60=cat(1,Label_60,num2cell(Sum_out_60));
W_60=cat(2,(['FileNames', FileName])',W_60);

writecell(W_60,[erase(['First60sec','_',num2str(min_amplitude_of_suck),'_',num2str(min_peak_distance),'_',num2str(min_suck_amount),'_',...
    num2str(min_burst_interval),'_',date,],'.'),'.xlsx'])
%% write firstB

Sum_out_firstB=cat(2,OUT_firstB.total_volume_consumed,OUT_firstB.total_time, OUT_firstB.avg_speed, OUT_firstB.Vmax_all, OUT_firstB.suck_number,...
    OUT_firstB.suck_consumed, OUT_firstB.suck_rate, OUT_firstB.inter_suck_interval,OUT_firstB.suck_duration);

Label_firstB={'total_volume_consumed',...
    'total_drinking_time',...
    'avg_speed',...
    'Vmax_mean','Vmax_SD','Vmax_SE','Vmax_Med','Vmax_IQR',...
    'Number_of_suck'...
    'Suck_consumed_mean','Suck_consumed_SD','Suck_consumed_SE','Suck_consumed_Med','Suck_consumed_IQR',...
    'Suck_rate_mean','Suck_rate_SD','Suck_rate_SE','Suck_rate_Med','Suck_rate_IQR',...
    'inter_sucking_interval_mean','inter_sucking_interval_SD','inter_sucking_interval_SE','inter_sucking_interval_Med','inter_sucking_interval_IQR',...
    'suck_duration_mean','suck_duration_SD','suck_duration_SE','suck_duration_Med','suck_duration_IQR'};

W_firstB=cat(1,Label_firstB,num2cell(Sum_out_firstB));
W_firstB=cat(2,(['FileNames', FileName])',W_firstB);

writecell(W_firstB,[erase(['FirstBurst','_',num2str(min_amplitude_of_suck),'_',num2str(min_peak_distance),'_',num2str(min_suck_amount),'_',...
    num2str(min_burst_interval),'_',date,],'.'),'.xlsx'])
%% write all the burst table
Sum_out_Bursts=[];
for a=1:20

    Sum_out_Bursts(a,:,:)=cat(2,squeeze(OUT_Bursts.total_volume_consumed(:,1,a))...
        ,squeeze(OUT_Bursts.total_time(:,1,a))...
        ,squeeze(OUT_Bursts.avg_speed(:,1,a))...
        ,squeeze(OUT_Bursts.Vmax(:,:,a))...
        ,squeeze(OUT_Bursts.suck_number(:,1,a))...
        ,squeeze(OUT_Bursts.suck_consumed(:,:,a))...
        ,squeeze(OUT_Bursts.suck_rate(:,:,a))...
        ,squeeze(OUT_Bursts.inter_suck_interval(:,:,a))...
        ,squeeze(OUT_Bursts.suck_duration(:,:,a)));
 
    Label_Burst{a}={[num2str(a),'_total_volume_consumed'],...
    [num2str(a),'_total_drinking_time'],...
    [num2str(a),'_avg_speed'],...
    [num2str(a),'_Vmax_mean'],[num2str(a),'_Vmax_SD'],[num2str(a),'_Vmax_SE'],[num2str(a),'_Vmax_Med'],[num2str(a),'_Vmax_IQR'],...
    [num2str(a),'_Number_of_suck']...
    [num2str(a),'_Suck_consumed_mean'],[num2str(a),'_Suck_consumed_SD'],[num2str(a),'_Suck_consumed_SE'],[num2str(a),'_Suck_consumed_Med'],[num2str(a),'_Suck_consumed_IQR'],...
    [num2str(a),'_Suck_rate_mean'],[num2str(a),'_Suck_rate_SD'],[num2str(a),'_Suck_rate_SE'],[num2str(a),'_Suck_rate_Med'],[num2str(a),'_Suck_rate_IQR'],...
    [num2str(a),'_inter_sucking_interval_mean'],[num2str(a),'_inter_sucking_interval_SD'],[num2str(a),'_inter_sucking_interval_SE'],[num2str(a),'_inter_sucking_interval_Med'],[num2str(a),'_inter_sucking_interval_IQR'],...
    [num2str(a),'_Suck_duration_mean'],[num2str(a),'_Suck_duration_SD'],[num2str(a),'_Suck_duration_SE'],[num2str(a),'_Suck_duration_Med'],[num2str(a),'_Suck_duration_IQR']};
    
end

O=[];
L={};
for a=1:20
    O=[O,squeeze(Sum_out_Bursts(a,:,:))];
    L=[L,Label_Burst{a}];
end
W_B=cat(1,L,num2cell(O));
W_B=cat(2,(['FileNames', FileName])',W_B);


writecell(W_B,[erase(['All_Burst','_',num2str(min_amplitude_of_suck),'_',num2str(min_peak_distance),'_',num2str(min_suck_amount),'_',...
    num2str(min_burst_interval),'_',date,],'.'),'.xlsx'])

%% write first half/second half

Sum_out_Bursts=[];
for a=1:2

    Sum_out_Bursts(a,:,:)=cat(2, squeeze(FistLast_Bursts.volume(:,1,a))...
, squeeze(FistLast_Bursts.time(:,1,a))...
, squeeze(FistLast_Bursts.speed(:,1,a))...
        ,squeeze(FistLast_Bursts.total_volume_consumed(:,:,a))...
        ,squeeze(FistLast_Bursts.total_time(:,:,a))...
        ,squeeze(FistLast_Bursts.avg_speed(:,:,a))...
        ,squeeze(FistLast_Bursts.Vmax(:,:,a))...
        ,squeeze(FistLast_Bursts.suck_number(:,:,a))...
        ,squeeze(FistLast_Bursts.suck_consumed(:,:,a))...
        ,squeeze(FistLast_Bursts.suck_rate(:,:,a))...
        ,squeeze(FistLast_Bursts.inter_suck_interval(:,:,a))...
        ,squeeze(FistLast_Bursts.inter_burst_interval(:,:,a))....
        ,squeeze(FistLast_Bursts.suck_duration(:,:,a))...
        ,squeeze(FistLast_Bursts.suck_number_all(:,:,a)));

end

Label_Burst={['firsthalf','_total_volume_consumed'],...
    ['firsthalf','_total_drinking_time'],...
    ['firsthalf','_avg_speed'],...
    ['firsthalf','_burst_volume_consumed_mean'],['firsthalf','_burst_volume_consumed_SD'],['firsthalf','_burst_volume_consumed_SE'],['firsthalf','_burst_volume_consumed_Med'],['firsthalf','_burst_volume_consumed_IQR'],...
    ['firsthalf','_burst_time_mean'],['firsthalf','_burst_time_SD'],['firsthalf','_burst_time_SE'],['firsthalf','_burst_time_Med'],['firsthalf','_burst_time_IQR'],...
    ['firsthalf','_burst_speed_mean'],['firsthalf','_burst_speed_SD'],['firsthalf','_burst_speed_SE'],['firsthalf','_burst_speed_Med'],['firsthalf','_burst_speed_IQR'],...
    ['firsthalf','_Vmax_mean'],['firsthalf','_Vmax_SD'],['firsthalf','_Vmax_SE'],['firsthalf','_Vmax_Med'],['firsthalf','_Vmax_IQR'],...
    ['firsthalf','_Number_of_suck_mean'],['firsthalf','_Number_of_suck_SD'],['firsthalf','_Number_of_suck_SE'],['firsthalf','_Number_of_suck_Med'],['firsthalf','_Number_of_suck_IQR'],...
    ['firsthalf','_Suck_consumed_mean'],['firsthalf','_Suck_consumed_SD'],['firsthalf','_Suck_consumed_SE'],['firsthalf','_Suck_consumed_Med'],['firsthalf','_Suck_consumed_IQR'],...
    ['firsthalf','_Suck_rate_mean'],['firsthalf','_Suck_rate_SD'],['firsthalf','_Suck_rate_SE'],['firsthalf','_Suck_rate_Med'],['firsthalf','_Suck_rate_IQR'],...
    ['firsthalf','_inter_sucking_interval_mean'],['firsthalf','_inter_sucking_interval_SD'],['firsthalf','_inter_sucking_interval_SE'],['firsthalf','_inter_sucking_interval_Med'],['firsthalf','_inter_sucking_interval_IQR'],...
    ['firsthalf','_inter_burst_interval_mean'],['firsthalf','_inter_burst_interval_SD'],['firsthalf','_inter_burst_interval_SE'],['firsthalf','_inter_burst_interval_Med'],['firsthalf','_inter_burst_interval_IQR'],...
    ['firsthalf','_Suck_duration_mean'],['firsthalf','_Suck_duration_SD'],['firsthalf','_Suck_duration_SE'],['firsthalf','_Suck_duration_Med'],['firsthalf','_Suck_duration_IQR'],...
    ['firsthalf','_suck_number_all'],...
    ['secondhalf','_total_volume_consumed'],...
    ['secondhalf','_total_drinking_time'],...
    ['secondhalf','_avg_speed'],...
    ['secondhalf','_burst_volume_consumed_mean'],['secondhalf','_burst_volume_consumed_SD'],['secondhalf','_burst_volume_consumed_SE'],['secondhalf','_burst_volume_consumed_Med'],['secondhalf','_burst_volume_consumed_IQR'],...
    ['secondhalf','_burst_time_mean'],['secondhalf','_burst_time_SD'],['secondhalf','_burst_time_SE'],['secondhalf','_burst_time_Med'],['secondhalf','_burst_time_IQR'],...
    ['secondhalf','_burst_speed_mean'],['secondhalf','_burst_speed_SD'],['secondhalf','_burst_speed_SE'],['secondhalf','_burst_speed_Med'],['secondhalf','_burst_speed_IQR'],...
    ['secondhalf','_Vmax_mean'],['secondhalf','_Vmax_SD'],['secondhalf','_Vmax_SE'],['secondhalf','_Vmax_Med'],['secondhalf','_Vmax_IQR'],...
    ['secondhalf','_Number_of_suck_mean'],['secondhalf','_Number_of_suck_SD'],['secondhalf','_Number_of_suck_SE'],['secondhalf','_Number_of_suck_Med'],['secondhalf','_Number_of_suck_IQR'],...
    ['secondhalf','_Suck_consumed_mean'],['secondhalf','_Suck_consumed_SD'],['secondhalf','_Suck_consumed_SE'],['secondhalf','_Suck_consumed_Med'],['secondhalf','_Suck_consumed_IQR'],...
    ['secondhalf','_Suck_rate_mean'],['secondhalf','_Suck_rate_SD'],['secondhalf','_Suck_rate_SE'],['secondhalf','_Suck_rate_Med'],['secondhalf','_Suck_rate_IQR'],...
    ['secondhalf','_inter_sucking_interval_mean'],['secondhalf','_inter_sucking_interval_SD'],['secondhalf','_inter_sucking_interval_SE'],['secondhalf','_inter_sucking_interval_Med'],['secondhalf','_inter_sucking_interval_IQR'],...
    ['secondhalf','_inter_burst_interval_mean'],['secondhalf','_inter_burst_interval_SD'],['secondhalf','_inter_burst_interval_SE'],['secondhalf','_inter_burst_interval_Med'],['secondhalf','_inter_burst_interval_IQR'],...
    ['secondhalf','_Suck_duration_mean'],['secondhalf','_Suck_duration_SD'],['secondhalf','_Suck_duration_SE'],['secondhalf','_Suck_duration_Med'],['secondhalf','_Suck_duration_IQR'],...
    ['secondhalf','_suck_number_all']};



O=[];
for a=1:2
    O=[O,squeeze(Sum_out_Bursts(a,:,:))];
end


W_B=cat(1,Label_Burst,num2cell(O));
W_B=cat(2,(['FileNames', FileName])',W_B);

writecell(W_B, [erase(['FirstHalf_SecondHalf_Burst','_',num2str(min_amplitude_of_suck),'_',num2str(min_peak_distance),'_',num2str(min_suck_amount),'_',...
    num2str(min_burst_interval),'_',date,],'.'),'.xlsx'])
