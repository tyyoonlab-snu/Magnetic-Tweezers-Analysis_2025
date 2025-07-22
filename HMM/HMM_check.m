%%%%%%% Baum-Welch algorithm  %%%%%%%%
%                                    %
%Modified by "Hyun-Kyu Choi" from SNU%
%                                    %
%------------------------------------%
%                                    %
%              |/ / / /              %
%          //  -    -   //           %
%             (  @ @  )              %
% +---------oOOo- & -oOOo--------+   %
% |          Best Regards        |   %
% +-------------------Oooo------+    %
%           oooO\    /(    )         %
%          (    )      )  /          %
%            )  (     (__/           %
%            (__/                    %
%------------------------------------%
if 1
    clear all
    close all
    clc
end
%% Loading data

% data(:,t,i) = observation vector at time t in sequence 'i-th'
% Q = num. hidden states
% M = num. mixture components
% numex = num. sequence (num. independent traces)
% O = Dimensionality of variable (e.g,(x,y,z) is O = 3)
% T = data length per each sequence

% Every loading data must be synchronized! == all seqeunce must have same scaled factor (ex> nm,um etc) and same positions of final/initial line

%%%%%%%%%%%%%%
species = 3;
%%%%%%%%%%%%%%%
if 1
    dpath = 'C:\Users\Owner\Desktop\Beta2-adrenergic receptor_GPCR\A Data and analysis\HMM\Reshaped\Apo with 2 mM TCEP\5pN';
    epath = [dpath,'\State position'];
    cd(dpath);
    
    SRread = dir('Stepresults*.mat');
    SR = struct(SRread);
    A = load(SR(1).name);
    B =A.HMMtrace_total;
    HMMtrace_total = B;
    
    dpmention = find(dpath == '\');
    
    disp(['Loading --',dpath(dpmention(end-1)+1:end),'_for HMM analysis'])
    
    rawdata = dir('*.txt');
    S = struct(rawdata);
    
    numdata= size(S,1);
    partitionposi = zeros(numdata,1);
    dataraw={};
    
    for i=1:numdata
        partitionposi(i) = find(S(i).name =='B');
        datalabel(i) = S(i).name(partitionposi(i)+1);
        dataraw{i} = load(S(i).name)';
    end
    
    for i = 1:numdata-1
        if datalabel(i) ~= datalabel(i+1)
            partition(i) = i;
        else
            partition(i) = 0;
        end
    end
    %% Setting conditions for HMM
    cpath = 'C:\Users\Owner\Desktop\GlpG\HMM_MJ';
    cd(cpath);
    saving = 1;
    SNRcorrection = 1;
    fps = 1200 %Hz;
    Precondition = 5 % str2num(input('Filtering process ? please, input median window size (Hz) [raw = 1200Hz] : ','s'));
    
    for i=1:numdata
        datafiltered{i} = medfilt1(dataraw{i},round(fps./Precondition));
    end
    Precondition1 = 1;%str2num(input('Do you know number of states in the system ? (yes = 1, no = otherwise): ','s'));
    Precondition2 = 1 ;%str2num(input('Do you know positions of states in the system ? (yes = 1, no = otherwise): ','s'));
    Precondition3 = 1 ;%str2num(input('Which one you want to analyze ? (Trace by Trace = 1, Bead by Bead = otherwise): ','s'));
    
    M = 1 ;  % Single Gaussian component at each state (1D information means M = 1)
    
    
    normalSD = zeros(numdata,1);
    smoothSD = zeros(numdata,1);
    SNR = zeros(numdata,1);
    SNR2 = zeros(numdata,1);
    %%
    cd(epath);
    
    rawposi= dir('*posi.txt');
    Sposi = struct(rawposi);
    mu_input = load(Sposi.name);
    if 1
        rawsig= dir('*sig.txt');
        Ssig = struct(rawsig);
        sig_input = load(Ssig.name);
    end
    cd(dpath);
    %% data arrangement
    
    for k = 1:numel(dataraw)
        Precondition0 =  Precondition;
        normalSD(k) = std(dataraw{k}(end-(2*fps/Precondition0):end));
        smoothSD(k) = std(medfilt1(dataraw{k}(end-(2*fps/Precondition0):end),round(fps/Precondition0)));
        SNR(k) = 4.5/normalSD(k);
        SNR2(k) = 4.5/smoothSD(k);
        
        if SNRcorrection
            if SNR2(k)>= 2
                Precondition0 =  Precondition0;
            else
                Precondition0 = floor(Precondition0.*(SNR2(k)./2));
            end
        end
        filterwindow(k)=Precondition0
    end
end
%% Vis
Ns = 4;

for k = 1:numel(dataraw)
    i=1;
    figure(k);
    subplot(3,Ns,1:Ns)
    hold on; box on; grid on;
    plot((1:length(dataraw{k}(fps/(3*filterwindow(k)):end)))./fps,dataraw{k}(fps/(3*filterwindow(k)):end),'color',[0.9 0.9 0.9],'LineWidth',0.2)
    plot((1:length(dataraw{k}(fps/(3*filterwindow(k)):end)))./fps,medfilt1(dataraw{k}(fps/(3*filterwindow(k)):end),round(fps./Precondition)),'k','LineWidth',0.5)
    plot((1:length(HMMtrace_total{k}{1}))./fps,HMMtrace_total{k}{1},'r','LineWidth',0.5)
    set(gca,'fontsize', 12, 'linewidth', 1, 'fontweight', 'bold','TickLength',[0.02 0.02])
    %set(gcf,'Renderer','painters')
    xlabel('Time (s)');     ylabel('Extension (nm)')
    %% find state posi
    zstate{k} = mu_input(:,k);
    for j = 1:numel(HMMtrace_total{k}{i})-1
        trans_point{k}(j) = HMMtrace_total{k}{i}(j+1) - HMMtrace_total{k}{i}(j);
    end
    
    disp('Total number of transitions')
    ind{k}= find(trans_point{k}~=0);
    numel(ind{k})
    
    for q = 1:numel(ind{k})
        if q == 1
            % present state check
            delzlist = zstate{k} - HMMtrace_total{k}{i}(ind{k}(q));
            pre_state = find(abs(delzlist)==min(abs(delzlist))); % 1 == Uz, 2== I1, 3 == I2 ...
            
            % check un/re folding directionality
            trans_direc = trans_point{k}(ind{k}(q))./abs(trans_point{k}(ind{k}(q))); % +1 = unfolding, -1=folding
            
            % collect dwell time data
            dwelltime{k}{pre_state}{1+0.5*(1+trans_direc)} = (ind{k}(q))/fps; % s;
            
        else
            % present state check
            delzlist = zstate{k} - HMMtrace_total{k}{i}(ind{k}(q));
            pre_state = find(abs(delzlist)==min(abs(delzlist))); % 1 == Uz, 2== I1, 3 == I2 ...
            % check un/re folding directionality
            trans_direc = trans_point{k}(ind{k}(q))./abs(trans_point{k}(ind{k}(q))); % +1 = unfolding, -1=folding
            
            % collect dwell time data
            temptdwell =(ind{k}(q)-ind{k}(q-1))/fps;
            dwelltime_temp{k}{pre_state}{1+0.5*(1+trans_direc)}(q) = temptdwell;
        end
    end
end

for k = 1:numel(dataraw)
    for i = 1:numel(dwelltime_temp{k})
        for j = 1:numel(dwelltime_temp{k}{i})
            ind_tmp = dwelltime_temp{k}{i}{j} ~= 0
            dwelltime_indi{k}{i}{j} = dwelltime_temp{k}{i}{j}(ind_tmp)
            figure(k);
            subplot(3,Ns,4*j+i)
            hist(dwelltime_indi{k}{i}{j})
        end
    end
end

%% put indi_data set into single one
dwellt_UtoI1= {};
dwellt_I1toI2 = {};
dwellt_I2toN= {};
dwellt_NtoI2= {};
dwellt_I2toI1= {};
dwellt_I1toU= {};

for k = 1:numel(dataraw)
    for i = 1:numel(dwelltime_indi{k})
        if i == 1
            if isempty(dwelltime_indi{k}{i}) == 1
                dwellt_UtoI1{k} = 0;
            else
                dwellt_UtoI1{k} = dwelltime_indi{k}{1}{1};
            end
        elseif i == 2
            if size(dwelltime_indi{k}{i},2) == 1
                dwellt_I1toI2{k} = dwelltime_indi{k}{2}{1};
                dwellt_I1toU{k} = 0;
            else
                dwellt_I1toI2{k} = dwelltime_indi{k}{2}{1};
                dwellt_I1toU{k} = dwelltime_indi{k}{2}{2};
            end
        elseif i == 3
            if size(dwelltime_indi{k}{i},2) == 1
                dwellt_I2toN{k} = dwelltime_indi{k}{3}{1};
                dwellt_I2toI1{k} = 0;
            else
                dwellt_I2toN{k} = dwelltime_indi{k}{3}{1};
                dwellt_I2toI1{k} = dwelltime_indi{k}{3}{2};
            end
        elseif i == 4
            dwellt_NtoI2{k} = dwelltime_indi{k}{4}{2};
        end
    end
end

t_all{1} = cell2mat(dwellt_UtoI1)';
t_all{2}= cell2mat(dwellt_I1toI2)';
t_all{3} = cell2mat(dwellt_I2toN)';
t_all{4} = cell2mat(dwellt_NtoI2)';
t_all{5}= cell2mat(dwellt_I2toI1)';
t_all{6}= cell2mat(dwellt_I1toU)';
%%%%%%%%%%%%%%
cutoff = 0.05;%
%%%%%%%%%%%%%%

figure(111); hold all; box on;
for i = 1:numel(t_all)
    if i == 4
        ind_cutoff = t_all{i} > cutoff & t_all{i} < 10;
        
        t_cut =t_all{i}(ind_cutoff);
    else
        ind_cutoff = t_all{i} > cutoff;
        
        t_cut =t_all{i}(ind_cutoff);
    end
    
    if species == 1 % WT
        if floor(round(max(t_cut) - min(t_cut))) < 8
            [count,timebin] = hist(t_cut,12);
        else
            if i == 3
                [count,timebin] = hist(t_cut,floor(round(max(t_cut) - min(t_cut))*1.8));
            else
                [count,timebin] = hist(t_cut,floor(round(max(t_cut) - min(t_cut))));
            end
        end
        scaling = [1.4 3.5 1.7 1 2 1.05]% Uz 5[1.5 2.8 3.7 3 3 1.3]  %Uz 6 [1.8 1.8 3.75 3 3 1] % WT [1.4 3.5 1.7 1 2 1.05]
        
        count(1) = count(1)./scaling(i);
        timebin = timebin.*scaling(i)
    elseif species == 2 % A206G
        if floor(round(max(t_cut) - min(t_cut))) < 8
            [count,timebin] = hist(t_cut,10);
        else
            if i == 3
                [count,timebin] = hist(t_cut,floor(round(max(t_cut) - min(t_cut)))*0.8);
            else
                [count,timebin] = hist(t_cut,floor(round(max(t_cut) - min(t_cut))));
            end
        end
        scaling = [1.3 3 1.5 1 1.75 1.14];
        count(1) = count(1)./scaling(i);
        if i == 2
            count(2) = count(2)*0.75;
            count(3) = count(3)*1.5;
            timebin = timebin.*scaling(i)
        else
            timebin = timebin.*scaling(i)
        end
    elseif species == 3 % L155A
        if floor(round(max(t_cut) - min(t_cut))) < 8
            if i == 6
                [count,timebin] = hist(t_cut,floor(round(max(t_cut) - min(t_cut)))*2);
            else
                [count,timebin] = hist(t_cut,15);
            end
        else
            if i == 3
                [count,timebin] = hist(t_cut,floor(round(max(t_cut) - min(t_cut)))*1.8);
            elseif i == 5
                [count,timebin] = hist(t_cut,floor(round(max(t_cut) - min(t_cut)))*2);
            elseif i == 4
                [count,timebin] = hist(t_cut,floor(round(max(t_cut) - min(t_cut))));
            else
                [count,timebin] = hist(t_cut,floor(round(max(t_cut) - min(t_cut))));
            end
        end
        scaling = [1.4 3 2.8 2.8 1.5 1];
        count(1) = count(1)./scaling(i);
        timebin = timebin.*scaling(i)
        
    elseif species == 4 %F121EL133E
        if floor(round(max(t_cut) - min(t_cut))) < 8
            [count,timebin] = hist(t_cut,10);
        else
            if i == 3
                [count,timebin] = hist(t_cut./1.8,floor(round(max(t_cut) - min(t_cut))*0.3));
            else
                [count,timebin] = hist(t_cut,floor(round(max(t_cut) - min(t_cut))));
            end
        end
        
        scaling = [1.55 1.6 1 1 3 0.7]; % WT scaling = [1.4 3.5 1.7 1 2 1.05];
        count(1) = count(1)./scaling(i);
        
        timebin = timebin.*scaling(i)
    elseif species == 5 %F121EL133E 5pN
        if floor(round(max(t_cut) - min(t_cut))) < 8
            if i == 5
                [count,timebin] = hist(t_cut,10);
                count(3) = 1.8*count(3);
            else
                [count,timebin] = hist(t_cut,10);
            end
        else
            if i == 3
                [count,timebin] = hist(t_cut,floor(round(max(t_cut) - min(t_cut)))*1.8);
            else
                [count,timebin] = hist(t_cut,floor(round(max(t_cut) - min(t_cut))));
            end
        end
        
        scaling = [1.7 2.2 3.4 1 3.5 1.05]; % WT scaling = [1.4 3.5 1.7 1 2 1.05];
        count(1) = count(1)./scaling(i);
        timebin = timebin.*scaling(i)
        
    end
    
    count_tot = sum(count);
    
    normcount = count./count_tot;
    
    
    param = {'k','Nconst'};
    startpoint = [0.1,0.1];
    syms tdistb k Nconst t;
    tdistb = Nconst*exp(-k*t);
    dwellsexp = fittype(char(tdistb),'independent','t','coefficient',param);
    fitres_se = fit(timebin',normcount',dwellsexp,'startpoint',startpoint,'lower',[0.01,0.01],'upper',[100,100])
    
    kfit(i,1) = fitres_se.k;
    Nconstfit(i) = fitres_se.Nconst;
    fit_error = confint(fitres_se);
    kfit_error = fit_error(:,1);
    kfit(i,2) = abs(kfit_error(2,1) - kfit_error(1,1))/2 ;
    
    %t_se_result=[kfit' Afit'];
    tbint = round(max(timebin)-min(timebin));
    
    
    tb = 0.05:0.05:tbint+1.05;
    tmpz1 = feval(fitres_se,tb);
    
    
    
    if i < numel(t_all)/2+1
        subplot(2,numel(t_all)/2,i); hold on; box on;
        bar(timebin, normcount);
        plot(tb,tmpz1,'r')
        
        title(['State',num2str(i),' to ','state',num2str(i+1)])
        text(timebin(end-1),0.5,['n = ', num2str(round(count_tot))])
        text(timebin(end-1),0.75,[num2str(kfit(i,1)),' +- ',num2str(kfit(i,2)),' 1/s'])
        
        xlim([0 round(max(timebin))+1])
        ylim([0 ceil(max(normcount))])
    else
        subplot(2,numel(t_all)/2,10-i); hold on; box on;
        bar(timebin, normcount);
        plot(tb,tmpz1,'r');
        
        title(['State',num2str(8-i),' to ','state',num2str(7-i)])
        text(timebin(end-1),0.5,['n = ', num2str(round(count_tot))])
        text(timebin(end-1),0.75,[num2str(kfit(i,1)),' +- ',num2str(kfit(i,2)),' 1/s'])
        
        xlim([0 round(max(timebin))+1])
        ylim([0 ceil(max(normcount))])
    end
end
kfit(:,1)
set(gcf,'Renderer','painters')
