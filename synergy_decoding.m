clear 
% 
datafiles=[dir('.\synergy\mPFC\*.mat'); ...
    dir('.\synergy\vHPC\*.mat')];

numBalancedCells=8;
datafiles=BalanceCells(datafiles,numBalancedCells,{'c'}); % balance the number of cells

fileinfo=regexp(replace({datafiles(:).name}','_','-'),...
            '(?<brainarea>\w*)-(?<session>\w*)-(?<mouse>\w*)-(?<featurename>\S+).mat','names');
fileinfo=[fileinfo{:}]';
v={};
for fld=fieldnames(fileinfo)'
    v.(fld{1})=replace({fileinfo(:).(fld{1})}','-','_');
end
fileinfo=v;

bandedFeatures={...
      ...'hilbert_amp'...
      ...'hilbert_phase'...
      ...'power'...
      ...'power_relative'...
      'power_z'...
    }';

bandUsed={...
     '2_4hz'...
     '5_12hz'...
     '16_25hz'...
     '31_61hz'...
     '62_122hz'...
     '123_153hz'...
    }';

featureUsed=[
    {'c'}...
   ... reshape(strcat(repmat(bandedFeatures',length(bandUsed),1),...
    ...'_',repmat(bandUsed,1,length(bandedFeatures))),1,[])...
    ]';

stimulusIDUsed={'socialA','socialB','novel','empty'};
 
trainingManner='binwise'; % 'binwise'

brainAreaCompared={'mPFC','vHPC'}; 

bRunShuffle=true;  % whether to run shuffling of stimulus IDs
nSuffleTimes=100;    % the number of shuffle time

predictAccuracy=[];  % brain area x session x mouse 
predictConfusionMatrix=[]; % brain area x session x mouse x Confusion Matrix
predictDescription={}; % brain area x session x mouse
predictAccuracy_shuffle=[]; % brain area x session x mouse x nshuffle
predictConfusionMatrix_shuffle=[];  % brain area x session x mouse x Confusion Matrix x nshuffle

for bshfl=[false, true]
    if bshfl
        if bRunShuffle
            nrloop=nSuffleTimes;            
        else
            nrloop=0;
        end
    else
        nrloop=1;
    end
    for irloop=1:nrloop

        for iba=1:length(brainAreaCompared)  % select brain area
            isel_ba=contains(fileinfo.brainarea,brainAreaCompared{iba},'IgnoreCase',true);
        
            ss=unique(lower(fileinfo.session(isel_ba)));
            % for iss=1:3  % session
            for iss=1:length(ss)  % session
                isel_ss=isel_ba&...
                    matches(fileinfo.session,ss{iss},IgnoreCase=true);
        
                mice=unique(lower(fileinfo.mouse(isel_ss)));
                for im=1:length(mice)  % mouse            
                    isel_im=isel_ss&...
                        matches(fileinfo.mouse,mice{im},IgnoreCase=true);
        
                    if ~bshfl
                        disp(['Processing #' brainAreaCompared{iba} ', ' ss{iss} ', ' mice{im} '...']);
                    else
                        disp(['(shuffle loop ' num2str(irloop) ') Processing #' brainAreaCompared{iba} ', ' ss{iss} ', ' mice{im} '...']);
                    end     
        
                    isel_fea=isel_im&...  % features
                        contains(fileinfo.featurename,featureUsed,IgnoreCase=true);
        
                    iPesudoCell=find(isel_fea);
        
                    binWidth=200; % ms
                    tBaseline=[0,0.6]; %s
        
                    features=[];  % nFeature x tTime x nSample
                    labels={};  % nSample x 1
                    trange=[]; % tTime x 1 (s)
                    feasType=fileinfo.featurename(iPesudoCell);
        
                    for ic=1:length(iPesudoCell)
                        d=load([datafiles(iPesudoCell(ic)).folder '\' datafiles(iPesudoCell(ic)).name]);
                        itrial=matches(d.events_labels.stimulus_ID,stimulusIDUsed,"IgnoreCase",true);                
                        
                        v=d.lfp_data(itrial,:); % nSample x tTime
        
                        % binned values (mean value)
                        nbinWidth=binWidth*d.infos.Fs/1000;
                        v=v(:,1:floor(size(v,2)/nbinWidth)*nbinWidth);
                        v=permute(mean(reshape(v',nbinWidth,[],size(v,1)),1,'omitmissing'),[3,2,1]);
                        t=(0:size(v,2)-1)*binWidth/1000;  % s                
        
                        % remove baseline
                        v=v-mean(v(:,t>tBaseline(1) & t<=tBaseline(2)),2);
                        
                        tDecode=[1.4,4];
                        v=v(:,t>tDecode(1) & t<=tDecode(2));
                        t=t(t>tDecode(1) & t<=tDecode(2));
                        
                        if isempty(trange)
                            trange=t;
                        elseif ~isequal(trange,t)
                            error(['Time labels are inconsistent with the others for file: ' datafiles(iPesudoCell(ic)).folder '\' datafiles(iPesudoCell(ic)).name '!']);
                        end
        
                        % process different kinds of features
                        ftype=fileinfo.featurename(ic);
                        if contains(ftype,'power','IgnoreCase',true)
                            if contains(ftype,'power_relative','IgnoreCase',true)
        
                                % power_relative 
        
                            elseif contains(ftype,'power_z','IgnoreCase',true)
        
                                % power_z 
        
                            else
        
                                % power 
        
                            end
                        elseif contains(ftype,'hilbert_amp','IgnoreCase',true)
        
                            % hilbert_amp
        
                        elseif contains(ftype,'hilbert_phase','IgnoreCase',true)
        
                            % hilbert_phase
        
                        elseif contains(ftype,'c','IgnoreCase',true)
        
                            % spike
                            v=v/binWidth*1000; % Hz
        
                        else
                            error(['Unable to process feature type: ' ftype '!']);
                        end
        
                        % normalize 
                        v=rescale(v,"InputMin",min(v),"InputMax",max(v)); % 0~1 for each time bin
        
                        if ~isempty(v)
                            features=cat(3,features,v);  % nSample x tTime x nFeature
                            if isempty(labels)
                                labels=d.events_labels.stimulus_ID(itrial)';
                            else
                                if ~isequal(labels,d.events_labels.stimulus_ID(itrial)')
                                    error(['Trial labels are inconsistent with the others for file: ' datafiles(iPesudoCell(ic)).folder '\' datafiles(iPesudoCell(ic)).name '!']);
                                end
                            end
                        end
                    end


                    features=permute(features,[3,2,1]);  % nSample x tTime x nFeature format transfer into nFeature x tTime x nSample
      
       
                    % pool bins together
                    labels=repmat(labels,size(features,2),1);
                    features=reshape(permute(features,[1,3,2]),size(features,1),[])'; % nSample x nFeature; average in tTime

                    % convert labels to digital values
                    dlabels=zeros(size(labels));
                    for istimulus=1:length(stimulusIDUsed)
                        idx=strcmpi(labels,stimulusIDUsed{istimulus});
                        dlabels(idx)=istimulus;
                    end
                    labels=dlabels;            
                    
                    confMatix=[];
                    buff_features=features;
                    buff_labels=labels; 

                    % set the number of loops for permutation of trials
                    nloop=1;
                    
                    for iloop=1:nloop

                        labels=buff_labels;
                        features=buff_features;

                        if bshfl  % shuffle the labels of stimulus ID
                            labels=labels(randperm(length(labels)));
                        end
                        
                        trueLabels=[];
                        predictLabels=[];
                        predictScores=[];  % 0~1
                        label_classes=unique(labels);            
                        %### cross validation : leave-one-out
                        for ifea=1:length(labels)
                            test_index=false(size(labels));
                            test_index(ifea)=true;                
                            
                            test_features=features(test_index,:);
                            test_labels=labels(test_index,:);
                            train_features=features(~test_index,:);
                            train_labels=labels(~test_index,:);    

 
                             % % # Decoding 1: Joint Bayes classifier (?? test only)
                            % %p(s|f1,f2,...,fn) 鈭?p(s)p(f1|s)p(f2|s)...p(fn|s),
                            % %s--stimulus, f--feature
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
            
                            trueLabels=[trueLabels;test_labels];
                            predictLabels=[predictLabels;judge_labels];
                            predictScores=[predictScores;sc];
                        end            
            
                        confMatix(:,:,iloop)=confusionmat(trueLabels(:),predictLabels(:));
                    end

                    % evaluate the accuracy
                    if ~isempty(confMatix)
                        confMatix=ceil(mean(confMatix,3));
                        predictDescription{iba,iss,im}=[brainAreaCompared{iba} '-' ss{iss} '-' mice{im}];
                        if ~bshfl
                            predictConfusionMatrix(iba,iss,im,:,:)=confMatix;
                            predictAccuracy(iba,iss,im)=mean(diag(confMatix)./sum(confMatix,2)); % balanced accuracy
                        else
                            predictConfusionMatrix_shuffle(iba,iss,im,:,:,irloop)=confMatix;
                            predictAccuracy_shuffle(iba,iss,im,irloop)=mean(diag(confMatix)./sum(confMatix,2)); % balanced accuracy
                        end
                    end

                end
                end
           end
      end
end

predictValid=~cellfun(@isempty,predictDescription); % brain area x session x mouse
accList=predictAccuracy(:); accList=accList(predictValid); % sample x 1
descripList=predictDescription(:); descripList=descripList(predictValid); % sample x 1

accAll={};
confMatixAll=[]; % confusion matrix x brain area
accAll_shfl={};
confMatixAll_shfl=[]; % confusion matrix x brain area
for iba=1:length(brainAreaCompared)    
    iValid=reshape(predictValid(iba,:,:),1,[]);
    v=reshape(predictAccuracy(iba,:,:),1,[]);
    accAll{1,iba}=v(iValid);
    confMatixAll(:,:,iba)=permute(sum(sum(predictConfusionMatrix(iba,:,:,:,:),2),3),[4,5,1,2,3]);

    if bRunShuffle
        v=reshape(predictAccuracy_shuffle(iba,:,:,:),[],size(predictAccuracy_shuffle,4));
        accAll_shfl{1,iba}=mean(v(iValid,:),2)';
        vm=permute(sum(sum(predictConfusionMatrix_shuffle(iba,:,:,:,:,:),2),3),[4,5,6,1,2,3]);
        confMatixAll_shfl(:,:,iba)=round(mean(vm,3))';
    end
end

figure;
subplot(121);
plotgroups(accAll,{''},brainAreaCompared,'scatter','omitzero');
ylim([0,1]);
ylabel('accuracy');
title('Decoding accuracies');

if bRunShuffle
    subplot(122);
    plotgroups(accAll_shfl,{''},brainAreaCompared,'scatter','omitzero');
    ylim([0,1]);
    ylabel('accuracy');
    title('Decoding accuracies (shuffled)');
end

figure;
for iba=1:length(brainAreaCompared)
    subplot(2,length(brainAreaCompared),iba);
    confusionchart(confMatixAll(:,:,iba));
    title(['Total confusion matrix for ' brainAreaCompared{iba}]);

    if bRunShuffle
        subplot(2,length(brainAreaCompared),iba+2);
        confusionchart(confMatixAll_shfl(:,:,iba));
        title(['Total confusion matrix for ' brainAreaCompared{iba} '(shuffled)']);
    end
end

all_accuracy = [accAll{1,1}',accAll{1,2}',accAll_shfl{1,1}',accAll_shfl{1,2}'];

function rfiles=BalanceCells(dfiles,ncell,celltype)
% Screen files to balance the number of cells in each session
%   Input:   
%            dfiles    the file list structure returned by the dir function.
%            ncell     the number of cell of each session after balance
%            celltype  the cell type to be balanced. E.g., {'spike'}
%
%   Output:
%            rfiles   the screened files with balanced number of cells for
%                     each session.

fileinfo=regexp(replace({dfiles(:).name}','_','-'),...
            '(?<brainarea>\w*)-(?<session>\w*)-(?<mouse>\w*)-(?<featurename>\S+).mat','names');
fileinfo=[fileinfo{:}]';
v={};
for fld=fieldnames(fileinfo)'
    v.(fld{1})=replace({fileinfo(:).(fld{1})}','-','_');
end
fileinfo=v;

datalengthMili=3000;
bin_length = 300;
stimulusID ={'socialA','socialB','novel','empty'};

% balance the number of cells
uBrainArea=unique(fileinfo.brainarea);
isel_cell=false(size(fileinfo.featurename));  

for iba=1:length(uBrainArea)  % select brain area
    isel_ba=contains(fileinfo.brainarea,uBrainArea{iba},'IgnoreCase',true);
    ss=unique(lower(fileinfo.session(isel_ba)));
    for iss=1:length(ss)  % session
        isel_ss=isel_ba&...
            matches(fileinfo.session,ss{iss},IgnoreCase=true);

        disp(['Balancing cells #' uBrainArea{iba} ', ' ss{iss} '...']);

        mice=unique(lower(fileinfo.mouse(isel_ss)));        
        for im=1:length(mice)  % mouse
            isel_im=isel_ss&...
                matches(fileinfo.mouse,mice{im},IgnoreCase=true);            

            % select files with c or lfp
            isel_fea_c = isel_im & contains(fileinfo.featurename,celltype,IgnoreCase=true);
            isel_fea_lfp = isel_im & ~contains(fileinfo.featurename,celltype,IgnoreCase=true);

            idx_cell_c = find(isel_fea_c);
            idx_cell_lfp = find(isel_fea_lfp);

            % select files in c files
            if length(idx_cell_c)>=ncell
                % calculate mean
                dataVector_pre = zeros(length(idx_cell_c), 1);
                dataVector_during = zeros(length(idx_cell_c), 1);
                for mm = 1:length(idx_cell_c)
                    clear lfp_data events_labels raster_during_now raster_pre_now vector_length_now covVector_now covVector_now2 
                    cellName_now = dfiles(idx_cell_c(mm)).name;

                    load([dfiles(idx_cell_c(mm)).folder '\' cellName_now]);

                    data_count = zeros(1, datalengthMili/bin_length*2);
                    data_fr = zeros(1, datalengthMili/bin_length*2);

                    trialLength = size(lfp_data, 1);
                    data_count = zeros(trialLength, datalengthMili/bin_length*2);

                    for jj = 1:trialLength
                         % use 100 ms window to count\
                         data_count(jj,:) = histcounts(find(lfp_data(jj,:)==1),0:bin_length:datalengthMili*2) *1000/bin_length;
                         data_fr(jj,:) = data_count(jj,:) * 1000/bin_length; % firing rate (Hz)
                    end

                    data_pre_now = mean(mean( data_fr(:,1:datalengthMili/bin_length),2));
                    data_during_now = mean(mean( data_fr(:,datalengthMili/bin_length+1:datalengthMili/bin_length*2),2));

                    dataVector_pre(mm) = data_pre_now;
                    dataVector_during(mm) = data_during_now;
                end

                % sort
                [sorted_pre, sort_idx_pre] = sort(dataVector_pre, 'descend');
                [sorted_during, sort_idx_during] = sort(dataVector_during, 'descend');

                % select n cells
                selected_idx = sort_idx_during(1:ncell);   
                isel_cell(idx_cell_c(selected_idx)) = true;
            else
                isel_cell(isel_im)=false; % abandon the mouse if not enough cells were found
                disp(['Less than ' num2str(ncell) ' cells were found for mouse #'  ', ' mice{im} '! Skipped!']); 
            end

            isel_cell(idx_cell_lfp) = true;
        end
    end
end

rfiles = dfiles(isel_cell);
end


