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
species = 1;  % 1 == WT
%%%%%%%%%%%%%%%


if 1
    dpath = 'C:\Users\Owner\Desktop\Beta2-AR_GPCR\A Data and analysis\HMM\Reshaped\Following-up paper\Apo with 2 mM TCEP\4pN';
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
trans_point_sp = {};
for k = 1:numel(dataraw)
    i=1;
    figure(k);
    subplot(3,Ns,1:Ns)
    hold on; box on; grid on;
    plot((1:length(dataraw{k}(fps/(3*filterwindow(k)):end)))./fps,dataraw{k}(fps/(3*filterwindow(k)):end),'color',[0.9 0.9 0.9],'LineWidth',0.2)
    plot((1:length(dataraw{k}(fps/(3*filterwindow(k)):end)))./fps,medfilt1(dataraw{k}(fps/(3*filterwindow(k)):end),round(fps./Precondition)),'k','LineWidth',0.5)
    plot((1:length(HMMtrace_total{k}{i}))./fps,HMMtrace_total{k}{1},'r','LineWidth',0.5)
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
    %%%%%%%%%%%%%%%%
    cutoff = 1/10;% second %
    %%%%%%%%%%%%%%%%
    
    %% Reformatting HMM trace after cutting-off 
    if 1
        for q = 1:numel(ind{k})
            if  q == 1
                tempdt{k}(q) =  (ind{k}(q))/fps;
            else
                tempdt{k}(q) = (ind{k}(q)-ind{k}(q-1))/fps;
            end
        end
        indcutoff = 0;
        indcutoff = find(tempdt{k} < 0.05);
        HMMtrace_del{k}{i} = zeros(1,numel(HMMtrace_total{k}{i}));
        for qq = 1:numel(indcutoff)
            if indcutoff(qq) == 1
                HMMtrace_del{k}{i}(ind{k}(indcutoff(qq))) = 0;
            else
                HMMtrace_del{k}{i}(ind{k}(indcutoff(qq)-1)+1:ind{k}(indcutoff(qq))) = HMMtrace_total{k}{i}(ind{k}(indcutoff(qq)))-HMMtrace_total{k}{i}(ind{k}(indcutoff(qq)-1));
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for j = 1:numel(HMMtrace_total{k}{i})
            HMMtrace_total_cut{k}{i}(j) = HMMtrace_total{k}{i}(j) -  HMMtrace_del{k}{i}(j);
        end
        figure(k);
        subplot(3,Ns,1:Ns)
        plot((1:length(HMMtrace_total_cut{k}{i}))./fps,HMMtrace_total_cut{k}{i},'color',[1,0.5,0.1],'LineWidth',0.5);
    end
    % Find dwell time at every transition
    
    for j = 1:numel(HMMtrace_total_cut{k}{i})-1
        trans_point_cut{k}(j) = HMMtrace_total_cut{k}{i}(j+1) - HMMtrace_total_cut{k}{i}(j);
    end
    
    ind_cut{k}= find(trans_point_cut{k}~=0);
    numel(ind_cut{k})
    
    for q = 1:numel(ind_cut{k})
        if q == 1
            % present state check
            delzlist = zstate{k} - HMMtrace_total_cut{k}{i}(ind_cut{k}(q));
            pre_state = find(abs(delzlist)==min(abs(delzlist))); % 1 == Uz, 2== I1, 3 == I2 ...
            
            % check un/re folding directionality
            trans_direc = trans_point_cut{k}(ind_cut{k}(q))./abs(trans_point_cut{k}(ind_cut{k}(q))); % +1 = unfolding, -1=folding
            
            % collect dwell time data
            dwelltime{k}{pre_state}{1+0.5*(1+trans_direc)} = (ind_cut{k}(q))/fps; % s;
            
        else
            % present state check
            delzlist = zstate{k} - HMMtrace_total_cut{k}{i}(ind_cut{k}(q));
            pre_state = find(abs(delzlist)==min(abs(delzlist))); % 1 == Uz, 2== I1, 3 == I2 ...
            % check un/re folding directionality
            trans_direc = trans_point_cut{k}(ind_cut{k}(q))./abs(trans_point_cut{k}(ind_cut{k}(q))); % +1 = unfolding, -1=folding
            
            % collect dwell time data
            temptdwell =(ind_cut{k}(q)-ind_cut{k}(q-1))/fps;
            dwelltime_temp{k}{pre_state}{1+0.5*(1+trans_direc)}(q) = temptdwell;
        end
    end
    %pause;
end

for k = 1:numel(dataraw)
    for i = 1:numel(dwelltime_temp{k})
        for j = 1:numel(dwelltime_temp{k}{i})
            ind_tmp = dwelltime_temp{k}{i}{j} ~= 0
            dwelltime_indi{k}{i}{j} = dwelltime_temp{k}{i}{j}(ind_tmp)
            figure(k);
            subplot(3,Ns,Ns*j+i)
            hist(dwelltime_indi{k}{i}{j})
            ylabel('Count'); xlabel('Dwell time (s)')
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
            dwellt_NtoI2{k} = dwelltime_indi{k}{6}{2};
        end
    end
end

t_all{1} = cell2mat(dwellt_UtoI1)';
t_all{2}= cell2mat(dwellt_I1toI2)';
t_all{3}= cell2mat(dwellt_I2toN)';

t_all{4} = cell2mat(dwellt_NtoI2)';
t_all{5}= cell2mat(dwellt_I2toI1)';
t_all{6}= cell2mat(dwellt_I1toU)';

figure(111); hold all; box on;

for i = 1:numel(t_all)
    
    t_cut =t_all{i};
    
    if species == 1 % WT
        if floor(round(max(t_cut) - min(t_cut))) < 9
            [count,timebin] = hist(t_cut,12);
        else
            [count,timebin] = hist(t_cut,floor(round(max(t_cut) - min(t_cut))));
        end
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
        text(timebin(2),0.5,['n = ', num2str(round(count_tot))])
        text(timebin(2),0.75,[num2str(kfit(i,1)),' +- ',num2str(kfit(i,2)),' 1/s'])
        
        xlim([0 round(max(timebin))+1])
        ylim([0 ceil(max(normcount))])
    else
        subplot(2,numel(t_all)/2,(3*Ns-2-i)); hold on; box on;
        bar(timebin, normcount);
        plot(tb,tmpz1,'r');
        
        title(['State',num2str(2*Ns-i),' to ','state',num2str(2*Ns-i)])
        text(timebin(2),0.5,['n = ', num2str(round(count_tot))])
        text(timebin(2),0.75,[num2str(kfit(i,1)),' +- ',num2str(kfit(i,2)),' 1/s'])
        
        xlim([0 round(max(timebin))+1])
        ylim([0 ceil(max(normcount))])
    end
end
kfit(:,1)
kfit(:,1)
sc = 45;
figure(112); hold on; grid on;
scatter(1:Ns-1,kfit(1:Ns-1,1),sc,'ro')
scatter(Ns:2*Ns-2,kfit(Ns:2*Ns-2,1),sc,'bo')
errorbar(1:Ns-1,kfit(1:Ns-1,1),kfit(1:Ns-1,2))
errorbar(Ns:2*Ns-2,kfit(Ns:2*Ns-2,1),kfit(Ns:2*Ns-2,2))
set(gca,'yscale','log')
ylim([0.1 10])
set(gcf,'Renderer','painters')

