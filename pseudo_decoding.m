clear 

datafiles=[dir('.\data\pseudo\mPFC\*.mat'); ...
    dir('.\data\pseudo\vHPC\*.mat')];


fileinfo=regexp(replace({datafiles(:).name}','_','-'),...
            '(?<brainarea>\w*)-(?<session>\w*)-(?<mouse>\w*)-(?<featurename>\S+).mat','names');
fileinfo=[fileinfo{:}]';
v={};
for fld=fieldnames(fileinfo)'
    v.(fld{1})=replace({fileinfo(:).(fld{1})}','-','_');
end
fileinfo=v;

% LFP features
bandedFeatures={...
    ...'hilbert_amp'...
    ...'hilbert_phase'...
    'power'...
    'power_relative'...
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
lfpFeatures=reshape(strcat(repmat(bandedFeatures',length(bandUsed),1),...
    '_',repmat(bandUsed,1,length(bandedFeatures))),1,[]);

% combine LFP & c features
featureUsed=[
    {'c'}...
    ...lfpFeatures
    ]';


stimulusIDUsed={'socialA','novel'};

brainAreaCompared={'mPFC','vHPC'}; 

binWidth=200; % ms
tBaseline=[0,1]; %s
tTraining=[1.5,4]; %s. used for training decoder

trainingManner='binwise'; % 'binwise'

bSpatialFiltering=false;  % whether to do spatial filtering for two classes
npcs=10; % the number of projected components retrieved by spatial filtering
pcAcc={};
pcAcc_shuffle={};
pcInds={1:npcs}; % use the first several PCs
for ipc=1:length(pcInds)


% pseudo-population trials across sessions
ncCellUsed=10;
nLFPUnitUsed=0;
ntrialUsed=6;
nloop=1000;
predictAccuracy=[];  % brain area x time x nloop
predictConfusionMatrix=[]; % brain area x Confusion Matrix x time x nloop
predictAccuracy_shuffle=[]; % brain area x time x nshuffle
predictConfusionMatrix_shuffle=[];  % brain area x Confusion Matrix x time x nshuffle
trange=[]; % tTime x 1 (s)
datafiles_buff=datafiles;
fileinfo_buff=fileinfo;
for iloop=1:nloop
    for bshfl=[false, true]

        datafiles=datafiles_buff;
        fileinfo=fileinfo_buff;

        % select files with enough trials for each stimulus
        isel_file=false(size(datafiles));
        for iba=1:length(brainAreaCompared)
            isel_ba=contains(fileinfo.brainarea,brainAreaCompared{iba},'IgnoreCase',true);
            ss=unique(lower(fileinfo.session(isel_ba)));
            for iss=1:length(ss)  % session
                isel_ss=isel_ba&...
                    matches(fileinfo.session,ss{iss},IgnoreCase=true);

                mice=unique(lower(fileinfo.mouse(isel_ss)));
                for im=1:length(mice)  % mouse
                    isel_im=isel_ss&...
                        matches(fileinfo.mouse,mice{im},IgnoreCase=true);

                    isel_fea=isel_im&...  % features
                        contains(fileinfo.featurename,featureUsed,'IgnoreCase',true);

                    iPesudoCell=find(isel_fea);
                    if ~isempty(iPesudoCell)
                        iTypicalFile=iPesudoCell(1);  % assume all pseudo cells recorded in the same session have the same number of trials for each stimulus
                        d=load([datafiles(iTypicalFile).folder '\' datafiles(iTypicalFile).name]);

                        bEnoughTrials=true;
                        for istim=1:length(stimulusIDUsed)
                            itrial=find(matches(d.events_labels.stimulus_ID,stimulusIDUsed{istim},"IgnoreCase",true));
                            if length(itrial)<ntrialUsed
                                bEnoughTrials=false;
                                disp(['Not enough trials for #' brainAreaCompared{iba} ', ' ss{iss} ', mouse #'  mice{im} '! Skipped!']);
                                break;
                            end                            
                        end
                        if bEnoughTrials
                            isel_file(isel_im)=true;                            
                        end
                    end

                end
            end
        end
        datafiles=datafiles(isel_file);
        for fld=fieldnames(fileinfo)'
            fileinfo.(fld{1})=fileinfo.(fld{1})(isel_file);
        end

        % select files with balanced cells & LFP units among conditions
        isel_file=false(size(datafiles));
        for iba=1:length(brainAreaCompared)
            isel_ba=contains(fileinfo.brainarea,brainAreaCompared{iba},'IgnoreCase',true);

            rng('shuffle');

            % c: randomly select cells across sessions
            idx_c_cells=find(isel_ba & contains(fileinfo.featurename,'c','IgnoreCase',true));            
            if length(idx_c_cells)<ncCellUsed
                error(['Totally ' num2str(length(idx_c_cells)) ' c cells were found for #' brainAreaCompared{iba} ...
                    ', less than the required number of cells, ' num2str(ncCellUsed) '!']);
            end
            if ~isempty(idx_c_cells) && sum(contains(lower(featureUsed),'c'))~=0
                isel_c=idx_c_cells(randperm(length(idx_c_cells),ncCellUsed));
            else
                isel_c=[];
            end

            % LFP: randomly select LFP units (with all associated features) across sessions.
            idx_lfp_unit=find(isel_ba & contains(fileinfo.featurename,lfpFeatures{1},'IgnoreCase',true));
            if length(idx_lfp_unit)<nLFPUnitUsed
                error(['Totally ' num2str(length(idx_lfp_unit)) ' LFP units were found for #' brainAreaCompared{iba} ...
                    ', less than the required number of LFP units, ' num2str(nLFPUnitUsed) '!']);
            end
            lfpFeatureUsed=featureUsed(~contains(featureUsed,'c','IgnoreCase',true));            
            if ~isempty(lfpFeatureUsed)
                isel_lfp_unit=idx_lfp_unit(randperm(length(idx_lfp_unit),nLFPUnitUsed)); 
                lfpfiles={datafiles(isel_lfp_unit).name}';
                for ilfp=1:length(lfpFeatureUsed)
                    lfpfiles=cat(1,lfpfiles,regexprep({datafiles(isel_lfp_unit).name}',lfpFeatures{1},lfpFeatureUsed{ilfp},'ignorecase'));
                end
                isel_lfp=find(matches({datafiles(:).name}',lfpfiles,"IgnoreCase",true));
            else
                isel_lfp=[];
            end
            
            isel_file(isel_c)=true;
            isel_file(isel_lfp)=true;

            disp([brainAreaCompared{iba} ': ' num2str(length(isel_c)) '/' num2str(length(idx_c_cells)) ' cells, '...
                num2str(length(isel_lfp)/length(lfpFeatureUsed)) '/' num2str(length(idx_lfp_unit)) ' LFP uints (' num2str(length(isel_lfp)) ' LFP features) were randomly selected!']);
        end
        datafiles=datafiles(isel_file);
        for fld=fieldnames(fileinfo)'
            fileinfo.(fld{1})=fileinfo.(fld{1})(isel_file);
        end

        % pseudo-population analysis
        for iba=1:length(brainAreaCompared)  % select brain area
            isel_ba=contains(fileinfo.brainarea,brainAreaCompared{iba},'IgnoreCase',true);

            rng('shuffle');

            % prepare pseudo-population trials  
            disp('Preparing pseudo-population trials ...');
            
            features=[];  % nSample x nFeature x tTime (or nTrial x nPseudo-population x nBin)
            labels={};  % nSample x 1            
            feasType={}; % nFeature x 1

            ss=unique(lower(fileinfo.session(isel_ba)));
            for iss=1:length(ss)  % session
                isel_ss=isel_ba&...
                    matches(fileinfo.session,ss{iss},IgnoreCase=true);

                mice=unique(lower(fileinfo.mouse(isel_ss)));
                for im=1:length(mice)  % mouse
                    isel_im=isel_ss&...
                        matches(fileinfo.mouse,mice{im},IgnoreCase=true);

                    if bshfl
                        disp(['pc No.# ' num2str(ipc) '(shuffle loop ' num2str(iloop) '/' num2str(nloop)...
                            ') Processing #' brainAreaCompared{iba} ', ' ss{iss} ', ' mice{im} '...']);
                    else
                        disp(['pc No.# ' num2str(ipc) '(loop ' num2str(iloop) '/' num2str(nloop)...
                            ') Processing #' brainAreaCompared{iba} ', ' ss{iss} ', ' mice{im} '...']);
                    end

                    isel_fea=isel_im&...  % features
                        contains(fileinfo.featurename,featureUsed,'IgnoreCase',true);

                    iPesudoCell=find(isel_fea);

                    itrialSel=[];
                    for ic=1:length(iPesudoCell) % simultaneously recorded cells
                        ftype=datafiles(iPesudoCell(ic)).name;
                        d=load([datafiles(iPesudoCell(ic)).folder '\' datafiles(iPesudoCell(ic)).name]);

                        % randomly select trials for each stimulus.                        
                        if isempty(itrialSel) ...
                                || contains(ftype,'c','IgnoreCase',true) % destroy noise correlation between cells by shuffling the trials for each cell (globally & LFP unaffected)
                            itrialSel=[];
                            for istim=1:length(stimulusIDUsed)
                                itrial=find(matches(d.events_labels.stimulus_ID,stimulusIDUsed{istim},"IgnoreCase",true));
                                if length(itrial)>=ntrialUsed
                                    itrial=itrial(randperm(length(itrial),ntrialUsed));
                                    itrialSel=[itrialSel,itrial];
                                else
                                    itrialSel=[];
                                    break;
                                end
                            end
                        end
                        if isempty(itrialSel)
                            disp(['Not enough trials for #' brainAreaCompared{iba} ', ' ss{iss} ', mouse #'  mice{im} '! Skipped!']);
                            break;
                        end

                        v=d.lfp_data(itrialSel,:); % nSample x tTime

                        % binned values (mean value)
                        nbinWidth=binWidth*d.infos.Fs/1000;
                        v=v(:,1:floor(size(v,2)/nbinWidth)*nbinWidth);
                        v=permute(mean(reshape(v',nbinWidth,[],size(v,1)),1,'omitmissing'),[3,2,1]);
                        t=(0:size(v,2)-1)'*binWidth/1000;  % s

                        % remove baseline
                        v=v-mean(v(:,t>tBaseline(1) & t<=tBaseline(2)),2);

                        if isempty(trange)
                            trange=t;
                        elseif ~isequal(trange,t)
                            error(['Time labels are inconsistent with the others for file: ' datafiles(iPesudoCell(ic)).folder '\' datafiles(iPesudoCell(ic)).name '!']);
                        end

                        % process different kinds of features
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

                            % c
                            v=v/binWidth*1000; % Hz

                        else
                            error(['Unable to process feature type: ' ftype '!']);
                        end

                        % normalize
                        v=rescale(v,"InputMin",min(v),"InputMax",max(v)); % 0~1 for each time bin

                        if ~isempty(v)
                            features=cat(3,features,v);
                            if isempty(labels)
                                labels=d.events_labels.stimulus_ID(itrialSel)';
                            else
                                if ~isequal(labels,d.events_labels.stimulus_ID(itrialSel)')
                                    error(['Trial labels are inconsistent with the others for file: ' datafiles(iPesudoCell(ic)).folder '\' datafiles(iPesudoCell(ic)).name '!']);
                                end
                            end
                            feasType=cat(1,feasType,ftype);
                        end
                    end
                end
            end
            features=permute(features,[1,3,2]); % nSample x nFeature x tTime
            features_spatial_filtered=[]; % features of each sample after spatial filtering . nSample x nFeature x tTime

            % decoding analysis            
            disp('Decoding ...');

            % convert labels to digital values
            dlabels=zeros(size(labels));
            for istimulus=1:length(stimulusIDUsed)
                idx=strcmpi(labels,stimulusIDUsed{istimulus});
                dlabels(idx)=istimulus;
            end
            labels=dlabels;

            confMatix=[];

            trueLabels=[];      % nSample x tTime
            predictLabels=[];   % nSample x tTime
            predictScores=[];  % 0~1. nSample x tTime x dScore
            label_classes=unique(labels);
            [nsmpl, nfea, ntime]=size(features);

            %### cross validation #2: leave-one-sample-pair-out
            % (auto-balance the sample sizes of each classes)
            nclassSampleIdx={};
            for icl=1:length(label_classes)
                nclassSampleIdx{icl}=find(labels==label_classes(icl));
            end
            nsampleCombs=min(cellfun(@length,nclassSampleIdx));
            for ismpl=1:nsampleCombs
                test_index=false(size(labels));
                for icl=1:length(label_classes)
                    test_index(nclassSampleIdx{icl}(ismpl))=true;
                end

                % divide into train & test sets
                test_features=features(test_index,:,:);  % nSample x nFeature x tTime
                test_labels=labels(test_index);
                train_features=features(~test_index,:,:); % nSample x nFeature x tTime
                train_labels=labels(~test_index);

                % specify the training bins
                train_features=train_features(:,:,trange>=tTraining(1) & trange<=tTraining(2));
 
                if bshfl  % shuffle the labels of stimulus ID
                    train_labels=train_labels(randperm(length(train_labels)));
                end

                % spatial filtering
                if bSpatialFiltering
                    pcs_all=[];
                    ulabels=unique(train_labels);
                    if length(ulabels)~=2
                        error(['Spatial filtering is only supported for binary classification! ' num2str(length(ulabels)) ' classes are found!']);
                    end
                    if ndims(features)~=3
                        error('A 3-dimensional features set (nSample x nFeature x tTime) is required for running the spatial filtering!');
                    end
                    splabels=train_labels; splabels(splabels==ulabels(1))=1; splabels(splabels~=ulabels(1))=-1;
                    sratio=sum(splabels==1)/sum(splabels==-1);
                    if sratio>=1/5 && sratio<=5 ...
                            && sum(splabels==1)>=2 && sum(splabels==-1)>=2
                        [mdl,pcs,ffea]=SpRayleigh(permute(train_features,[2,3,1]),splabels,npcs,'fisher','nozscore');                        
                        train_features=real(permute(pcs(pcInds{ipc},:,:),[3,1,2])); % use PCs                     
                        [test_features,ffea]=TransformFea(permute(test_features,[2,3,1]),mdl);
                        test_features=real(permute(test_features(pcInds{ipc},:,:),[3,1,2])); % use PCs
                        features_spatial_filtered(test_index,:,:)=permute(ffea,[3,1,2]);
                    else
                        continue;
                    end
                end

                % # Decoding
                switch lower(trainingManner)
                    case 'binwise' %  bin-wise training

                        train_labels=repmat(train_labels,size(train_features,3),1); % nSample x 1
                        train_features=reshape(permute(train_features,[2,1,3]),size(train_features,2),[])'; % nSample x nFeature
                        test_labels=repmat(test_labels,size(test_features,3),1); % nSample x 1
                        test_features=reshape(permute(test_features,[2,1,3]),size(test_features,2),[])'; % nSample x nFeature

                        % # Decoding (Naive Bayesian)
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
                                p_fs(v<1e-5)=1; % exclude insignificant features with extremely low variances
                                p_sf_numer(itestfea,iclass)=p_s*prod(p_fs);
                            end
                        end
                        p_sf=p_sf_numer./(sum(p_sf_numer,2)*ones(1,size(p_sf_numer,2)));
                        rnd_idx=randperm(size(p_sf,2))';
                        [~,s_max]=max(p_sf(:,rnd_idx),[],2);  % use randomized indices to avoid bias induced by possible selection tendency of the max function. 
                        judge_labels=label_classes(rnd_idx(s_max));
                        sc=p_sf;

                        % reserve the temporal order
                        ntesttrials=sum(test_index);
                        test_labels=reshape(test_labels,ntesttrials,ntime);   % nSample x tTime
                        judge_labels=reshape(judge_labels,ntesttrials,ntime); % nSample x tTime
                        sc=reshape(sc,ntesttrials,ntime,[]);                  % nSample x tTime x dScore

                    otherwise
                        error(['Unsupported training manner: ''' trainingManner '''!']);
                end

                trueLabels=cat(1,trueLabels,test_labels);
                predictLabels=cat(1,predictLabels,judge_labels);
                predictScores=cat(1,predictScores,sc);
            end

            for it=1:ntime
                confMatix=confusionmat(trueLabels(:,it),predictLabels(:,it));
                if bshfl
                    predictConfusionMatrix_shuffle(iba,:,:,it,iloop)=confMatix;
                    predictAccuracy_shuffle(iba,it,iloop)=mean(diag(confMatix)./sum(confMatix,2)); % balanced accuracy
                else
                    predictConfusionMatrix(iba,:,:,it,iloop)=confMatix;
                    predictAccuracy(iba,it,iloop)=mean(diag(confMatix)./sum(confMatix,2)); % balanced accuracy
                end
            end


            % visualization of features (before & after spatial filtering)
            if 0 ...
                    && ~isempty(features_spatial_filtered) ...  % nSample x nFeature x tTime
                    && ~bshfl
                
                % # Vis 1:  show the cells' distribution in the categorical firing space
                if length(stimulusIDUsed)>2
                    error('Visualization of all cells''s distribution is only supported for binary categories!');
                end

                if iloop==1
                    if iba==1
                        hfig_cell_dist=figure;
                    end

                    feas_all_cells{iba}=[];  % before spatial filtering. {brain area}[nFeature x nStimulus]
                    feas_all_cells_spatial_filtered{iba}=[];% after spatial filtering. {brain area}[nFeature x nStimulus]

                    feastype_all_cells{iba}=[];    % {brain area}[nFeature x 1]
                    feas_repeat_all_cells{iba}=[]; % repeat times. {brain area}[nFeature x 1]                    
                    
                else
                    figure(hfig_cell_dist);
                end
                feastype_all_cells{iba}=cat(1,feastype_all_cells{iba},...
                    setdiff(feasType,feastype_all_cells{iba}));
                [~,icell]=intersect(feastype_all_cells{iba},feasType);
                feas_all_cells{iba}(setdiff(icell,1:size(feas_all_cells{iba},1)),1:length(stimulusIDUsed))=0;
                feas_all_cells_spatial_filtered{iba}(setdiff(icell,1:size(feas_all_cells_spatial_filtered{iba},1)),1:length(stimulusIDUsed))=0;
                feas_repeat_all_cells{iba}(setdiff(icell,1:size(feas_repeat_all_cells{iba},1)),1)=0;
                feas=[];  % nFeature x stimulus
                feassp=[];
                for istimulus=1:length(stimulusIDUsed)
                    feas(:,istimulus)=mean(mean(features(...
                        labels==istimulus,:,:),1),3)';
                    feassp(:,istimulus)=real(mean(mean(features_spatial_filtered(...
                        labels==istimulus,:,:),1),3))';
                end
                feas_all_cells{iba}(icell,:)=feas_all_cells{iba}(icell,:)+feas;
                feas_all_cells_spatial_filtered{iba}(icell,:)=feas_all_cells_spatial_filtered{iba}(icell,:)+feassp;
                feas_repeat_all_cells{iba}(icell)=feas_repeat_all_cells{iba}(icell)+1;
                for bspatialfiltered=[false,true]
                    if bspatialfiltered
                        feas=feas_all_cells_spatial_filtered{iba}./feas_repeat_all_cells{iba};
                    else
                        feas=feas_all_cells{iba}./feas_repeat_all_cells{iba};
                    end
                    subplot(2,length(brainAreaCompared),iba+bspatialfiltered*2);
                    hs=[];hf=[];
                    for iftype=1:length(featureUsed)
                        ifea=contains(feastype_all_cells{iba},featureUsed{iftype},'IgnoreCase',true);
                        x=feas(ifea,1);
                        y=feas(ifea,2);
                        hs(iftype)=scatter(x,y);
                        f=fit(x,y,'poly1');
                        hold on;
                        hf(iftype)=plot(x,f(x),'linewidth',1.5);
                    end
                    xl=xlim;
                    hh=plot(xl,xl,'--','Color',[0.7,0.7,0.7],'linewidth',1.5);
                    xlabel(stimulusIDUsed{1});
                    ylabel(stimulusIDUsed{2});
                    hold off;                    
                    axis equal;
                    if bspatialfiltered
                        title([brainAreaCompared{iba} ' (after spatial filtering)']);
                    else
                        if iba==1
                            legend([hs,hf,hh],[featureUsed; strcat('fit-',featureUsed); {'x=y'}],'Interpreter','none','Location','bestoutside');
                        end
                        title([brainAreaCompared{iba} ' (before spatial filtering)']);                        
                    end
                end

            end

        end
    end
end


pcAcc{ipc}=predictAccuracy;
pcAcc_shuffle{ipc}=predictAccuracy_shuffle;
end
if length(pcInds)>1  % compare & plot decoding accuracies of each PCs
    figure;
    hl=[];
    hlshf=[];
    for iba=1:length(brainAreaCompared)
        m=[]; m_shf=[];
        s=[]; s_shf=[];
        for ipc=1:length(pcInds)
            m(ipc)=mean(mean(pcAcc{ipc}(iba,:,:),2),3);  % brain area x time x nloop
            s(ipc)=std(mean(pcAcc{ipc}(iba,:,:),2),[],3)/sqrt(nloop);

            m_shf(ipc)=mean(mean(pcAcc_shuffle{ipc}(iba,:,:),2),3);  % brain area x time x nloop
            s_shf(ipc)=std(mean(pcAcc_shuffle{ipc}(iba,:,:),2),[],3)/sqrt(nloop);
        end
        hl(iba)=plot(m);
        hold on;
        plotshaded(1:length(pcInds),[m+s;m-s],get(hl(iba),'Color'));
        hlshf(iba)=plot(m_shf,'--');
        plotshaded(1:length(pcInds),[m_shf+s_shf;m_shf-s_shf],get(hlshf(iba),'Color'));
    end
    legend([hl,hlshf],[brainAreaCompared,strcat(brainAreaCompared,'-shuffle')]);
    xticks(1:length(pcInds)); xticklabels(strcat('PC',num2str((1:length(pcInds))')));
    ylabel('Accuracy');
    title(['Decoding accuracies of spatial filtering PCs (' strjoin(stimulusIDUsed,' vs. ') ')']);
end


% plot time course
figure;
hl=[];
hlshf=[];
iselt=trange>0 & trange<=4;  
t=trange(iselt);
for iba=1:length(brainAreaCompared)  
    m=mean(predictAccuracy(iba,iselt,:),3);
    s=std(predictAccuracy(iba,iselt,:),[],3)/sqrt(nloop);
    hl(iba)=plot(t,m);
    hold on;
    plotshaded(t',[m+s;m-s],get(hl(iba),'Color'));    

    m=mean(predictAccuracy_shuffle(iba,iselt,:),3);
    s=std(predictAccuracy_shuffle(iba,iselt,:),[],3)/sqrt(nloop);
    hlshf(iba)=plot(t,m,'--');
    plotshaded(t',[m+s;m-s],get(hlshf(iba),'Color'));    
end
legend([hl,hlshf],[brainAreaCompared,strcat(brainAreaCompared,'-shuffle')]);
xlabel('Time(s)');
ylabel('Accuracy');
title(['Time course of decoding accuracies (' strjoin(stimulusIDUsed,' vs. ') ')']);


% plot mean values
accAll={};
confMatixAll=[]; % confusion matrix x brain area
confMatixAll_shfl=[];
for iba=1:length(brainAreaCompared)
    accAll{1,iba}=permute(mean(predictAccuracy(iba,:,:),2),[2,3,1]);
    accAll{2,iba}=permute(mean(predictAccuracy_shuffle(iba,:,:),2),[2,3,1]);
    confMatixAll(:,:,iba)=uint32(permute(mean(sum(predictConfusionMatrix(iba,:,:,:,:),4),5),[2,3,1,4,5]));
    confMatixAll_shfl(:,:,iba)=uint32(permute(mean(sum(predictConfusionMatrix_shuffle(iba,:,:,:,:),4),5),[2,3,1,4,5]));
end

figure;
plotgroups(accAll,{'normal','shuffled'},brainAreaCompared,'scatter','omitzero');
ylim([0,1]);
ylabel('accuracy');
title('Decoding accuracies');

figure;
for iba=1:length(brainAreaCompared)
    subplot(2,length(brainAreaCompared),iba);
    confusionchart(confMatixAll(:,:,iba));
    title(['Total confusion matrix for ' brainAreaCompared{iba}]);

    subplot(2,length(brainAreaCompared),iba+2);
    confusionchart(confMatixAll_shfl(:,:,iba));
    title(['Total confusion matrix for ' brainAreaCompared{iba} '(shuffled)']);
end

