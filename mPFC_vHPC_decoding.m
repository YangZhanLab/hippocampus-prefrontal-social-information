% This function was used for PFC vs. RE decoding analysis
function mPFC_RE_decoding

% listed are example files for illustration only
folder_name_mPFC={
   '20170525#4', [1,5,9,13,17,21,29] []; 
   '20170523#1', [1,5,9,13,17,21,25,29] []; 
    };

folder_name_vHPC={
    '20210712#11', [1,5,13,17,21,29] []; 
    '20210722#11', [1,5,9,13,17,21,29] [];
    };

dir_names_compared={[pwd '\data\mPFC\produce_single_cell_data_all_mPFC\'],...
    [pwd '\data\vHPC\produce_single_cell_data_all_vHPC\']};

folder_name_compared={folder_name_mPFC,folder_name_vHPC};

file_prefix_compared={'mPFCspike_','vHPCspike_'};

dir_cate={'mPFC','vHPC'};

% bin_length_compared=[100,300,500]; % ms
bin_length_compared=[75,100,150]; % ms
pop_cell_dim_compared=[1:10];

pop_cell_dim_max=pop_cell_dim_compared(end);

decode_percent_all_compared=[];
decode_percent_overall_compared=[];
decode_percent_overall_valid_session_compared=[];

%% Section I: decoding
if false % set as TRUE to re-run the decoding using the illustrative files
%if true

for idim=1:length(pop_cell_dim_compared)
pop_cell_dim=pop_cell_dim_compared(idim);

for idir=1:length(dir_names_compared)
dir_name=dir_names_compared{idir};
folder_name=folder_name_compared{idir};
file_prefix=file_prefix_compared{idir};

social_class={ 'socialA','socialB';
              'socialA','novel';
              'socialA','empty';
              'socialB','novel';
              'socialB','empty';
              'novel','empty' };
stimulusID= {'socialA';'empty';'socialB';'novel'};

[folderLength,~]=size(folder_name);    
datalengthMili=3000;
rvAll_pairs=cell(folderLength,6);
rvAll_names=cell(folderLength,2);

vectorAll_during=cell(folderLength,4);
vectorAll_pre=cell(folderLength,4);

vectorAll_during_norm=cell(folderLength,4);
vectorAll_during_normss=cell(folderLength,1);


for ibin=1:length(bin_length_compared)
bin_length=bin_length_compared(ibin);

label_all_test=[];
label_all_judge=[];

kk_used = [];
for kk=1:folderLength
    
    clear raster_data raster_datab raster_labels dataVector2 databVector2 dataVector_pre2 dataVector_during2
    name_now=folder_name{kk,1};
    name_a=[file_prefix name_now(1:8) '_m' name_now(10:end) '*.mat'];     
    
    name_cell_a=dir([dir_name name_a]);
    
    if isempty(name_cell_a) 
        warning(['''' name_a '''not found!']);
        continue; 
    end
    
    load([name_cell_a(1).folder '\' name_cell_a(1).name],'raster_labels');
    
    TrialLength=length(raster_labels.stimulus_ID);
    stimulusID_now = raster_labels.stimulus_ID;
    
    if exist('stimulusIDGroup','var')
        for istimCur=1:length(stimulusID_now)
            for istimType=1:length(stimulusID)
                if ~isempty(strfind(stimulusID{istimType},stimulusID_now{istimCur}))
                    stimulusID_now{istimCur}=stimulusID{istimType};
                end
            end
        end
    end 

    % only use sessions with more than 8 cells
    if length(name_cell_a)>=pop_cell_dim_max
          kk_used =[kk_used kk];
          % construct vector averaged by pre and during time
          dataVector_pre=zeros(length(name_cell_a),TrialLength);
          dataVector_during=zeros(length(name_cell_a),TrialLength);
          databVector_pre=zeros(length(name_cell_a),TrialLength);
          databVector_during=zeros(length(name_cell_a),TrialLength);
          % vector for all cells
           dataVector=zeros(length(name_cell_a),datalengthMili*2/bin_length,TrialLength);
           databVector=zeros(length(name_cell_a),datalengthMili*2/bin_length,TrialLength);
          
           dataVector_mean=zeros(length(name_cell_a),1);
           % convert to counts from 0/1,using bin size
           
           
        for mm=1:length(name_cell_a)
            clear raster_data raster_datab raster_labels raster_during_now raster_pre_now vector_length_now covVector_now covVector_now2 
            cellName_now=name_cell_a(mm).name;
            
            load([name_cell_a(mm).folder '\' cellName_now])
   
            data_count=zeros(TrialLength,datalengthMili/bin_length*2);
            datab_count=zeros(TrialLength,datalengthMili/bin_length*2);
            data_fr=zeros(TrialLength,datalengthMili/bin_length*2);
            datab_fr=zeros(TrialLength,datalengthMili/bin_length*2);
            for jj=1:TrialLength
                 % use 100 ms window to count
                 data_count(jj,:) = histcounts(find(raster_data(jj,:)==1),0:bin_length:datalengthMili*2) *1000/bin_length;
                 datab_count(jj,:) = histcounts(find(raster_datab(jj,:)==1),0:bin_length:datalengthMili*2) *1000/bin_length;

                  data_fr(jj,:)=data_count(jj,:)*1000/bin_length; % firing rate (Hz)
                  datab_fr(jj,:)=datab_count(jj,:)*1000/bin_length; % firing rate (Hz)
            end
            data_pre_now = mean( data_fr(:,1:datalengthMili/bin_length),2);
            data_during_now = mean( data_fr(:,datalengthMili/bin_length+1:datalengthMili/bin_length*2),2);
            datab_pre_now = mean( datab_fr(:,1:datalengthMili/bin_length),2);
            datab_during_now = mean( datab_fr(:,datalengthMili/bin_length+1:datalengthMili/bin_length*2),2);
            
            dataVector_pre(mm,:)=data_pre_now;
            dataVector_during(mm,:)=data_during_now;
            databVector_pre(mm,:)=datab_pre_now;
            databVector_during(mm,:)=datab_during_now;
            
            dataVector(mm,:,:)=data_fr';
            databVector(mm,:,:)=datab_fr';
            
            dataVector_mean(mm,1)=mean(mean(data_fr));
        end

        dataVector2 = dataVector(dataVector_mean>=0.5,:,:);
        databVector2 = databVector(dataVector_mean>=0.5,:,:);
        dataVector_during2 = dataVector_during(dataVector_mean>=0.5,:);
        dataVector_pre2 = dataVector_pre(dataVector_mean>=0.5,:);
        
        % further restriction: only using sessions with certain number
        % of cells all having a mean firing rate above the threshold
        if sum(dataVector_mean>=0.5)>=pop_cell_dim_max
           
            % 选取前面若干维的细胞
            dataVector2 = dataVector2(1:pop_cell_dim,:,:);
            databVector2 = databVector2(1:pop_cell_dim,:,:);
            dataVector_during2 = dataVector_during2(1:pop_cell_dim,:);
            dataVector_pre2 = dataVector_pre2(1:pop_cell_dim,:);
            
            vector_all=[];
            vector_label_all=[];
            
            for ss=1:length(stimulusID)
                
                raster_index=[];
                
                for jj=1:TrialLength
                    if strcmp(stimulusID_now{1,jj},stimulusID{ss})
                        raster_index = [ raster_index jj];
                    end
                end
                Nraster_now = length(raster_index);
                
                raster_during_now=dataVector_during2(:,raster_index);
                raster_pre_now=dataVector_pre2(:,raster_index);
                
                pre_mean = mean(dataVector2(:,1:3,raster_index),2);
                vector_now = dataVector2(:,datalengthMili/bin_length+2:datalengthMili/bin_length*2,raster_index)...
                    -repmat(pre_mean,1,datalengthMili/bin_length-1,1);
                
                if isempty(vector_now)
                    continue;
                end
                
                [c_a,c_b,c_c] = size(vector_now);
                
                vector_combine = reshape(vector_now,c_a,c_b*c_c);
                
                vector_label = ones(c_a,c_b*c_c)*ss;
                
                vector_all=[vector_all vector_combine];
                vector_label_all=[vector_label_all vector_label(1,:)];
            end
            
            
            % decoding
            features=vector_all';
            labels=vector_label_all';
            label_classes=unique(labels);
            
            
            session_test_labels=[];
            session_judge_labels=[];
            session_judge_scores=[];  % 0~1
            
            %### cross validation : leave-one-out
            for ifea=1:length(labels)
                test_index=false(size(labels));
                test_index(ifea)=true;
                
                
                test_features=features(test_index,:);
                test_labels=labels(test_index,:);
                train_features=features(~test_index,:);
                train_labels=labels(~test_index,:);                
                
                %# Decoding
                % p(s|f1,f2,...,fn) ∝ p(s)p(f1|s)p(f2|s)...p(fn|s),
                % s--stimulus, f--feature
                p_sf_numer=[];
                for iclass=1:length(label_classes)
                    fea=train_features(train_labels==label_classes(iclass),:);
                    m=mean(fea,1);
                    v=var(fea,0,1);
                    
                    p_s=mean(train_labels==label_classes(iclass));
                    for itestfea=1:size(test_features,1)
                        p_fs=exp(-(test_features(itestfea,:)-m).^2./v/2)./sqrt(2*pi*v);
                        p_sf_numer(itestfea,iclass)=p_s*prod(p_fs);
                    end
                end
                p_sf=p_sf_numer./(sum(p_sf_numer,2)*ones(1,size(p_sf_numer,2)));                
                [~,s_max]=max(p_sf,[],2);
                judge_labels=label_classes(s_max);
                sc=p_sf;     
                
                session_test_labels=[session_test_labels;test_labels];
                session_judge_labels=[session_judge_labels;judge_labels];
                session_judge_scores=[session_judge_scores;sc];
            end
                        
            
            label_all_test{end+1,1}=session_test_labels;
            label_all_judge{end+1,1}=session_judge_labels;            
            
        end
    end
 end
 
 valid_PFC_fileNo=length(label_all_test);
 
decode_percent_all = zeros(length(stimulusID),length(stimulusID),valid_PFC_fileNo);
decode_class_all = zeros(length(stimulusID),length(stimulusID),valid_PFC_fileNo);
for nn=1:valid_PFC_fileNo
    label_test_now = label_all_test{nn,1};
    label_judge_now = label_all_judge{nn,1};
    [Ntest,Nrandom] = size(label_test_now);
    label_test_now2 = reshape(label_test_now,Ntest*Nrandom,1);
    label_judge_now2 = reshape(label_judge_now,Ntest*Nrandom,1);
    for kk=1:length(stimulusID)
         label_test_now11 = label_test_now2(label_test_now2 == kk);
          label_judge_now11 = label_judge_now2(label_test_now2 == kk);
         for jj=1:length(stimulusID)
             decode_percent_all(kk,jj,nn) = length( find( label_judge_now11 ==jj) ) /length(label_test_now11);
             decode_class_all(kk,jj,nn) = length( find( label_judge_now11 ==jj) );
         end
    end
end

validIdx=true(1,valid_PFC_fileNo);
for nn=1:valid_PFC_fileNo
    if sum(sum(isnan(decode_percent_all(:,:,nn)))) ...
            || mean(decode_percent_all(:,1,nn)==1)==1
        validIdx(nn)=false;
    end
end

decode_percent_mean=mean(decode_percent_all(:,:,validIdx),3);

decode_percent_overall = zeros(1,valid_PFC_fileNo);
for nn=1:valid_PFC_fileNo
    label_test_now = label_all_test{nn,1};
    label_judge_now = label_all_judge{nn,1};    
    decode_percent_overall(nn) = mean(label_test_now==label_judge_now);
end

decode_percent_all_compared{idim,idir,ibin}=decode_percent_all(:,:,validIdx);
decode_percent_overall_compared{idim,idir,ibin}=decode_percent_overall;
decode_percent_overall_valid_session_compared{idim,idir,ibin}=decode_percent_overall(validIdx);

end
end
end
save([pwd '\data\decode_results_Eular_dist_4_class_all.mat']);
end

%% Section II: plot using saved results
load([pwd '\data\decode_results_Eular_dist_4_class_all.mat'])

% Confusion matrix
bin_compared=100; %ms
dim_compared=10;
clims=[0,0.7];
figure; 
for idir=1:length(dir_names_compared)
    subplot(1,length(dir_names_compared),idir);
    conf_mat=decode_percent_all_compared{pop_cell_dim_compared==dim_compared,idir,bin_length_compared==bin_compared};
    conf_mat_mean=mean(conf_mat,3);
    imagesc(conf_mat_mean,clims);
    title(dir_cate{idir});
    xlabel('Predicted Classes');
    ylabel('True Classes');
    set(gca,'xtick',1:length(stimulusID),'xticklabel',{'S1','E','S2','N'});
    set(gca,'ytick',1:length(stimulusID),'yticklabel',{'S1','E','S2','N'});
end
colorbar('eastoutside');

% population dimension effects
bin_compared=100; %ms
acc_compared=[];
color_compared={'b','r'};
figure; 
for idir=1:length(dir_names_compared)
    v=[];
    g=[];
    for idim=1:length(pop_cell_dim_compared)
        conf_mat=decode_percent_all_compared{idim,idir,bin_length_compared==bin_compared};
        acc=[];
        for isample=1:size(conf_mat,3)
            acc(:,isample)=diag(conf_mat(:,:,isample));
        end
        acc_compared{idir,idim}=mean(acc,1);        
        v=[v,acc_compared{idir,idim}];
        g=[g,repmat(idim,1,length(acc_compared{idir,idim}))];
    end    
    boxplot(v*100,g,'colors',color_compared{idir}, 'Widths',0.2,'position',(1:length(pop_cell_dim_compared))+0.4*(idir-1.5));
    hold on;
end
for idim=1:length(pop_cell_dim_compared)
    p=ranksum(acc_compared{1,idim},acc_compared{2,idim});
    hold on;
    if p<0.05        
        plot(idim+[-0.2,-0.2,0.2,0.2],[80,81,81,80],'k-','linewidth',0.5);
        plot(idim,83,'k*');
    end
end
box_cate=findall(gca,'Tag','Box');
legend(box_cate([length(pop_cell_dim_compared)+1,1]),dir_cate)
set(gca,'XLim',[0 length(pop_cell_dim_compared)+1],'xtick',1:length(pop_cell_dim_compared),'xticklabel',cellstr(num2str(pop_cell_dim_compared')));
set(gca,'YLim',[15 90]);
xlabel('Population dimension');
ylabel('Mean predicted accuracy (%)');
title(['Decoding performance under ' num2str(bin_compared) 'ms-bin condition']);

% bin length effects
dim_compared=10;
acc_compared=[];
color_compared={'b','r'};
figure;
for idir=1:length(dir_names_compared)
    v=[];
    g=[];
    for ibin=1:length(bin_length_compared)
        conf_mat=decode_percent_all_compared{pop_cell_dim_compared==dim_compared,idir,ibin};
        acc=[];
        for isample=1:size(conf_mat,3)
            acc(:,isample)=diag(conf_mat(:,:,isample));
        end
        acc_compared{idir,ibin}=mean(acc,1);        
        v=[v,acc_compared{idir,ibin}];
        g=[g,repmat(ibin,1,length(acc_compared{idir,ibin}))];
    end
    boxplot(v*100,g,'colors',color_compared{idir}, 'Widths',0.2,...
        'position',(1:length(bin_length_compared))+0.4*(idir-1.5));
    hold on;
end
for ibin=1:length(bin_length_compared)
    p=ranksum(acc_compared{1,ibin},acc_compared{2,ibin});
    hold on;
    if p<0.05
        plot(ibin+[-0.2,-0.2,0.2,0.2],[80,81,81,80],'k-','linewidth',0.5);
        plot(ibin,83,'k*');
    end
end
box_cate=findall(gca,'Tag','Box');
legend(box_cate([4,1]),dir_cate)
set(gca,'XLim',[0 4],'xtick',1:length(bin_length_compared),'xticklabel',cellstr(strcat(num2str(bin_length_compared'),'ms')));
set(gca,'YLim',[15 90]);
xlabel('Bin length');
ylabel('Mean predicted accuracy (%)');
title(['Decoding performance under ' num2str(dim_compared) '-cell population condition']);


