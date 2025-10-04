function h=plotgroups(v,glabel,tlabel,varargin)
% h=plotgroups(V,GLABEL,TLABEL,OPTIONS) plots bar or box graphs arranged in
% a GROUP x TYPE manner, where each TYPE of data within a GROUP will be
% plotted together. 
%
% V: a cell array of data vectors, {GROUP x TYPE}[1 x SAMPLES].
% GLABEL: a cell array of characters, labels of each GROUP, shown as x-axis labels;
% TLABEL: a cell array of characters, labels of each TYPE, shown as legends.
%
% options:
%    'bar'     --  bar plot with error bars (default)
%    'box'     --  box plot
%    'scatter' --  show scatters
%    'omitzero'    -- omit all zero values when plotting and performing analysis.
%    ('err','std'|'sem'(default)) -- (for bar plot) type of the error bar. 
%                  'std' -- standard error; 
%                  'sem' -- standard error of mean (default).
%    ('sigtest',true(default)|false) -- test for significance between types within each group.
%    ('testmethod','ranksum'(default)|'ttest'|'anova1'|'kruskalwallis') test method. 'ttest', paired-sample t-test.
%    ('tail','both'(default)|'right'|'left') -- 
%                   sig. test type: two-sided, right-sided, or left-sided
%    ('offset',voffset) -- offset of the plots relative to the default positions
%    ('alpha',valpha)   --  (for bar plot) alpha value of face transparency
%    ('linestyle',lstyle)  -- (for box plot) linestyple, compatible with BOXPLOT
%    ('linewidth',lwidth)  -- (for box plot) linewidth, compatible with BOXPLOT
%    ('colorlist',clist) -- cull array of 3-element color vectors,
%                           customized color list specifying color of each TYPE.
%
%  E.g.,
%          v=repmat({rand(1,100)},2,3);
%          plotgroups(v,{'grp1','grp2'},{'type1','type2','type3'});
%

% by mz.
% 2021/4/20

% Update: add comments and an option for customized color list.
% by mz.
% 2022/5/27

% Update: for constant values, i.e., there is only one sample in some of GROUP
% x TYPE bins, no error bar will be presented for bar plots.
% by mz.
% 2022/7/24

% Update: add support for significance test between each type-pair within
% each group.
% by mz.
% 2022/8/19

% Update: add support for one-way ANOVA test, non-parametric Kruskal-Wallis test,
% and multiple comparison test among types within each group. 
% by mz.
% 2023/5/17

% Update: add an option, 'omitzero', for omitting zero values when plotting and performing analysis.
% by mz.
% 2024/12/2

plottype='bar';
bshowScatter=false;
bomitZero=false;
voffset=0;
valpha=0.3;
lstyle='-';
lwidth=0.5;
errtype='sem';
sigtest=true;
testmethod='ranksum';
tail='both';
colorList={[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],...
    [0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],...
    [0.6350 0.0780 0.1840]};
narg=1;
while narg<=length(varargin)
    switch(lower(varargin{narg}))
        case 'bar'
        case 'box'
            plottype='box';
        case 'scatter'
            bshowScatter=true;
        case 'omitzero'
            bomitZero=true;
        case 'offset'
            narg=narg+1;
            voffset=varargin{narg};
        case 'alpha'
            narg=narg+1;
            valpha=varargin{narg};
        case 'linestyle'
            narg=narg+1;
            lstyle=varargin{narg};         
        case 'linewidth'
            narg=narg+1;
            lwidth=varargin{narg};            
        case 'colorlist'
            narg=narg+1;
            colorList=varargin{narg};
        case 'err'
            narg=narg+1;
            errtype=varargin{narg};
        case 'sigtest'
            narg=narg+1;
            sigtest=varargin{narg};
        case 'testmethod'
            narg=narg+1;
            testmethod=varargin{narg};
        case 'tail'
            narg=narg+1;
            tail=varargin{narg};
        otherwise
            error('Unsupported option!');
    end
    narg=narg+1;
end

if length(glabel)~=size(v,1) ||...
        length(tlabel)~=size(v,2)
    error('The number of groups or types dose match the corresponding data dimension!'); 
end

if bomitZero  % omit all zero values
     for iv=1:numel(v)
         v{iv}=v{iv}(v{iv}~=0);
     end
end

vplotpos=[];  % GROUP x TYPE
vplotwidth=0;
switch(plottype)
    case 'bar'   % bar plot
        x=(1:length(glabel))+voffset;
        y=[];
        err=[];
        for itp=1:length(tlabel)
            for igrp=1:length(glabel)
                y(igrp,itp)=mean(v{igrp,itp},'omitnan');
                switch lower(errtype)
                    case 'std'
                        err(igrp,itp)=std(v{igrp,itp},'omitnan');
                    case 'sem' 
                        err(igrp,itp)=std(v{igrp,itp},'omitnan')/sqrt(sum(~isnan(v{igrp,itp}))); % S.E.M.
                    otherwise
                        error('Unsupported error bar type!');
                end
            end
        end
        h=bar(x,y,0.7,'FaceAlpha',valpha);
        for ih=1:length(h) % match the colorlist
            h(ih).FaceColor=colorList{mod(ih-1,length(colorList))+1};
            h(ih).EdgeColor=colorList{mod(ih-1,length(colorList))+1};
        end
        hold on;
        for itp=1:length(tlabel)
            barpos=h(itp).XEndPoints;
            fcolor=h(itp).FaceColor;
            ishowbar=err(:,itp)>1e-10;            
            if bshowScatter
                for igrp=1:length(glabel)
                    d=v{igrp,itp};                    
                    scatter(barpos(igrp)*ones(size(d)),d,'filled',...
                        'MarkerEdgeColor',fcolor,...
                        'MarkerFaceColor',fcolor,...
                        'MarkerFaceAlpha',1,'SizeData',3,...
                        'jitter','on','jitterAmount',h(itp).BarWidth/20);
                end
            end
            errorbar(barpos(ishowbar),y(ishowbar,itp),err(ishowbar,itp),'.','linewidth',1,...
                'Capsize',6,'Color',[0,0,0]); 
            vplotpos(:,itp)=barpos';
            vplotwidth=h(itp).BarWidth;
        end
        set(gca,'fontsize',10);
        
    case 'box'  % box plot
        boxwidth=0.05;
        interblank=0.3;
        intvl=(1-interblank)/length(tlabel);
        startpos=-(length(tlabel)-1)/2*intvl+voffset;
        for itp=1:length(tlabel) 
            d=[]; g=[];
            for igrp=1:length(glabel)
                d=[d,v{igrp,itp}];
                g=[g,repmat(igrp,1,length(v{igrp,itp}))];
            end
            boxplot(d,g,'colors',colorList{mod(itp-1,length(colorList))+1}, 'Widths',boxwidth,...
                'position',(1:length(glabel))+startpos+(itp-1)*intvl);    
            hold on;
            if bshowScatter
                for igrp=1:length(glabel)
                    d=v{igrp,itp};
                    scatter(igrp*ones(size(d))+startpos+(itp-1)*intvl,d,'filled',...
                        'MarkerEdgeColor',colorList{mod(itp-1,length(colorList))+1},...
                        'MarkerFaceColor',colorList{mod(itp-1,length(colorList))+1},...
                        'MarkerFaceAlpha',0.3,'SizeData',20,...
                        'jitter','on','jitterAmount',0.02);
                end
            end
            vplotpos(:,itp)=[1:length(glabel)]'+startpos+(itp-1)*intvl;
            vplotwidth=boxwidth;
        end
        box_cate=findall(gca,'Tag','Box');
        set(box_cate(1:length(tlabel)*length(glabel)),'linestyle',lstyle,'linewidth',lwidth);
        h=box_cate(length(tlabel)*length(glabel):-length(glabel):1)';
end

% # sig. test
if sigtest
    ymax=max([v{:}],[],'all'); ymin=min([v{:}],[],'all');
    ystep=(ymax-ymin)/30;
    yoff=ymax+ystep;
    for igrp=1:length(glabel)
        ppairs=[];
        if strcmpi(testmethod,'anova1')|| ... % one-way ANOVA among types within each group
                strcmpi(testmethod,'kruskalwallis')
            g={};
            for itp=1:length(tlabel)
                g{itp}=repmat(itp,size(v{igrp,itp}));
            end
            vv=cat(2,v{igrp,:});
            gg=cat(2,g{:});
            eval(['[p,~,stats]=' testmethod '(vv,gg,''off'');']);
            if p<0.1
                sigmark='*';if p<0.05;sigmark='**';if p<0.01;sigmark='***';end;end
                glabel{igrp}=[glabel{igrp} sigmark];
            end
            c = multcompare(stats,'Display','off');
            for ic=1:size(c,1)
                ppairs(c(ic,1),c(ic,2))=c(ic,6);
            end
        else
            for itp1=1:length(tlabel)
                for itp2=itp1+1:length(tlabel)
                    if length(v{igrp,itp1})==1  % single point vs. group points
                        [~,I]=sort([v{igrp,itp1},v{igrp,itp2}],2,'descend');
                        switch tail
                            case 'both'
                                p=1-abs(find(I==1)-length(I)/2)/length(I)/2;
                            case 'right'
                                p=find(I==1)/length(I);
                            case 'left'
                                p=1-find(I==1)/length(I);
                        end
                    elseif length(v{igrp,itp2})==1 % group points vs. single point
                        [~,I]=sort([v{igrp,itp2},v{igrp,itp1}],2,'descend');
                        switch tail
                            case 'both'
                                p=1-abs(find(I==1)-length(I)/2)/length(I)/2;
                            case 'right'
                                p=1-find(I==1)/length(I);
                            case 'left'
                                p=find(I==1)/length(I);
                        end
                    else  % group vs. group points
                        switch lower(testmethod)
                            case 'ranksum'
                                [p,ht,stats]=ranksum(v{igrp,itp1},v{igrp,itp2},'tail',tail);
                            case 'ttest'
                                [ht,p,ci,stats]=ttest(v{igrp,itp1},v{igrp,itp2},'tail',tail);
                            otherwise
                                error('Unsupported significance test method!');
                        end
                    end
                    ppairs(itp1,itp2)=p;
                end
            end
        end
        for itp1=1:length(tlabel)
            for itp2=itp1+1:length(tlabel)
                p=ppairs(itp1,itp2);
                if p<0.1
                    sigmark='*';if p<0.05;sigmark='**';if p<0.01;sigmark='***';end;end
                    b=vplotpos(igrp,itp1);
                    e=vplotpos(igrp,itp2);
                    plot([b,b,e,e],[yoff,yoff+ystep,yoff+ystep,yoff],'k-');
                    yoff=yoff+1.5*ystep;
                    hold on;
                    text((b+e)/2,yoff,sigmark,'FontSize',8,'HorizontalAlignment','center');
                    yoff=yoff+ystep;
                end
            end
        end
    end
    ylim([ymin,yoff+0.01]);
end

legend(h,tlabel);
set(gca,'xtick',1:length(glabel),'xticklabel',glabel);



