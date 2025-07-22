clear all
close all
clc

cpath = 'C:\Users\tylab\Desktop\HK_R\HMM_MJ';
dpath = 'C:\Users\tylab\Desktop\HK_R\HMM_MJ\30% PG\WT\2018.06.12 _ Refolding trace (1.2kHz)\5pN';
fpath = 'C:\Users\tylab\Desktop\HK_R\HMM_MJ\Synchronized\WT'; % 5 Hz
%%%% Intermediate --> scaling = 50nm
%%%% Refolding --> scaling = 24nm at 5pN
%%
saving = 0;
%%
cd(dpath);
RTraw=dir(['*.txt']);
S = struct(RTraw);
nfile = size(RTraw,1);

la = 30000;
display(['Set acquisition rate (Hz)']);
fps = str2num(input('Set acquisition rate (Hz): ','s'));
fc = 60;
fs2 = 20;
fs3 = 10;
fs = 5;
fs3 = 30;

RTz = {};
smoothRTz={};
Transz ={};
Norz = {};
Cutz = {};

for i = 1:nfile
    RTz{i}=load(S(i).name);
end

cd(cpath)
a = 22;  %Normailized facfor 26 for 6pN, 22 for 5pN

Minz = {};

for j = 1:nfile
    smoothRTz{j} = medfilt1(RTz{j},fps/fs);
    Minz{j} = min(smoothRTz{j});
    Transz{j} = RTz{j} - Minz{j} + a/4;
    Norz{j} = (Transz{j})/a;
end

tis = {};

refoldingz = 0.15;
notrefoldingz = 0.85;

if 1
    for i = 1:nfile
        m = 0.01;
        x = 0.05;
        figure(100+i);
        f=figure(100+i);
        Screensize = get( groot, 'Screensize' );
        set(f, 'Position', [600, 100, Screensize(3)*0.65, Screensize(4)*0.85 ]);
        figp = get(gcf,'Position');
        set(0,'DefaultFigurePosition', figp);
        display(['Shift traces of cycle' num2str(i)])
        
        while 1
            f;
            plot(Norz{i},'Color',[0.7 0.7 0.7],'LineWidth',1.5);
            hold on;
            axis([0 length(Norz{i}) -0.4 1.4]);
            plot(get(gca,'xlim'), [refoldingz refoldingz],'r','LineWidth',4);
            hold on;
            %plot(get(gca,'xlim'),[refoldingI1 refoldingI1],'Color',[1 0.5 0],'LineWidth',2);
            %plot(get(gca,'xlim'),[refoldingI2 refoldingI2],'g','LineWidth',2);
            plot(get(gca,'xlim'),[notrefoldingz notrefoldingz],'b','LineWidth',4);
            hold on;
            plot(medfilt1(Norz{i},fps/fs3),'y','LineWidth',1.5);
            plot(medfilt1(Norz{i},fps/fs),'k','LineWidth',1.5);
            
            axis([0 length(Norz{i}) -0.4 1.4]);
            grid on;
            set(gca, 'fontsize', 16, 'linewidth', 1, 'fontweight', 'bold','TickLength',[0.02 0.02]);
            xlabel('Time (s)'); ylabel('Normalized extension (nm)');
            
            ans = getkey
            if ans == 28
                Norz{i} = Norz{i} - m;
            elseif ans == 29
                Norz{i} =Norz{i} + m;
            elseif ans == 30
                Norz{i} = Norz{i} + x;
            elseif ans == 31
                Norz{i} = Norz{i} - x;
            else
                break
            end
            clf(f);
        end
        pause;
        getpts
        tis{i} = ans;
        ao= tis{i}(1);
        co = tis{i}(end);
        Mnorz(i) = mean(Norz{i}(ao:co));
        Norz{i} = Norz{i}- (Mnorz(i) - refoldingz);
        
        if length(Norz{i}) < la
            N(i) = 1;
            Cutz{i}{1} = Norz{i};
        else
            N(i) = ceil(length(Norz{i})./la);
            for j = 1:N(i)-1
                Cutz{i}{j} = Norz{i}(1+la*(j-1):la*(j));
            end
            
            for j = N(i)
                Cutz{i}{j} = Norz{i}(1+la*(j-1):end);
            end
        end
        close(100+i);
    end
    
      for i = 1:nfile
        m = 0.25;
        x = 0.5;
        figure(10+i);
        f=figure(10+i);
        Screensize = get( groot, 'Screensize' );
        set(f, 'Position', [600, 100, Screensize(3)*0.65, Screensize(4)*0.85 ]);
        figp = get(gcf,'Position');
        set(0,'DefaultFigurePosition', figp);
        display(['Shift traces of cycle' num2str(i)])
        
        while 1
            f;
            plot(stt_a206g{j},'k','LineWidth',1.5);
            hold on;
            axis([0 length(stt_a206g{j}) -5 85]);
            plot(get(gca,'xlim'), [zup zup],'m','LineWidth',3);
            hold on;
            plot(medfilt1(stt_a206g{j},20),'y','LineWidth',1.5);
            axis([0 length(stt_a206g{j}) -5 85]);
            grid on;
            set(gca, 'fontsize', 16, 'linewidth', 1, 'fontweight', 'bold','TickLength',[0.02 0.02]);
            xlabel('Time (s)'); ylabel('Normalized extension (nm)');
            
            ans = getkey
            if ans == 28
                stt_a206g{j} = stt_a206g{j} - m;
            elseif ans == 29
                stt_a206g{j} =stt_a206g{j} + m;
            elseif ans == 30
                stt_a206g{j} = stt_a206g{j} + x;
            elseif ans == 31
                stt_a206g{j} = stt_a206g{j} - x;
            else
                break
            end
            clf(f);
        end
        pause;
        getpts
        tis{i} = ans;
        ao= tis{i}(1);
        co = tis{i}(end);
        Mnorz(i) = mean(stt_a206g{j}(ao:co));
        stt_a206g{j} = stt_a206g{j}- (Mnorz(i) - zup);
        
        if length(stt_a206g{j}) < la
            N(i) = 1;
            Cutz{i}{1} = stt_a206g{j};
        else
            N(i) = ceil(length(stt_a206g{j})./la);
            for j = 1:N(i)-1
                Cutz{i}{j} = stt_a206g{j}(1+la*(j-1):la*(j));
            end
            
            for j = N(i)
                Cutz{i}{j} = stt_a206g{j}(1+la*(j-1):end);
            end
        end
        close(100+j);
    end
end

C =[0:0.05:1; 0:0.05:1; 0:0.05:1]';

if 1
    for i= 1:nfile
        for j = 1:N(i)
            figure(100*i+j);
            f=figure(100*i+j);
            Screensize = get( groot, 'Screensize' );
            set(f, 'Position', [600, 150, Screensize(3)*0.65, Screensize(4)*0.85 ]);
            figp = get(gcf,'Position');
            set(0,'DefaultFigurePosition', figp);
            
            plot(Cutz{i}{j},'color',C(j,:),'LineWidth',1.5);
            hold on;
            plot(medfilt1(Cutz{i}{j},floor(fps/fc)),'y','LineWidth',2);
            hold on;
            plot(medfilt1(Cutz{i}{j},fps/fs3),'g','LineWidth',1.5);
            plot(medfilt1(Cutz{i}{j},fps/fs),'b','LineWidth',2);
            grid on;
            set(gca, 'fontsize', 16, 'linewidth', 1, 'fontweight', 'bold','TickLength',[0.02 0.02]);
            axis([0 length(Cutz{i}{j})+2 -0.5 1.5]);
            xlabel('Frame (#)'); ylabel('Normailized extension');
            plot(get(gca,'xlim'), [refoldingz refoldingz],'r','LineWidth',4);
            hold on;
            plot(get(gca,'xlim'),[notrefoldingz notrefoldingz],'b','LineWidth',4);
        end
    end
    figure(777);
    f=figure(777);
    Screensize = get( groot, 'Screensize' );
    set(f, 'Position', [600, 150, Screensize(3)*0.65, Screensize(4)*0.85 ]);
    figp = get(gcf,'Position');
    set(0,'DefaultFigurePosition', figp);
    
    plot(cell2mat(Norz'),'color',C(j,:),'LineWidth',1.5);
    hold on;
    plot(medfilt1(cell2mat(Norz'),floor(fps/fc)),'y','LineWidth',2);
    hold on;
    plot(medfilt1(cell2mat(Norz'),fps/fs),'b','LineWidth',2);
    grid on;
    set(gca, 'fontsize', 16, 'linewidth', 1, 'fontweight', 'bold','TickLength',[0.02 0.02]);
    axis([0 length(cell2mat(Norz'))+2 -0.5 1.5]);
    plot(get(gca,'xlim'), [refoldingz refoldingz],'r','LineWidth',4);
    hold on;
    plot(get(gca,'xlim'),[notrefoldingz notrefoldingz],'b','LineWidth',4);
    xlabel('Frame (#)'); ylabel('Normailized extension');
end

sss = find(dpath == '8');

if saving
    for i = 1:nfile
        for j =1:N(i)
            cd(fpath)
            Norsmothz{i}= medfilt1(Cutz{i}{j},floor(fps/fc));
            temp = Norsmothz{i}(floor(fps/(fc*3)):end-floor(fps/(fc*3)));
            save(['60Hz window-',dpath(sss-1:sss+6),'_',S(i).name(1:7),'_',num2str(i),'_',num2str(j),'.txt'],'temp','-ascii')
        end
    end
    
    
    for i = 1:nfile
        for j =1:N(i)
            cd(fpath)
            Norsmothz{i}= medfilt1(Cutz{i}{j},fps/fs);
            temp = Norsmothz{i}(floor(fps/(fs*3)):end-floor(fps/(fs*3)));
            save(['5Hz window-',dpath(sss-1:sss+6),'_',S(i).name(1:7),'_',num2str(i),'_',num2str(j),'.txt'],'temp','-ascii')
        end
    end
    
    for i = 1:nfile
        for j =1:N(i)
            cd(fpath)
            Norsmothz{i}= Cutz{i}{j};
            temp = Norsmothz{i};
            save([num2str(fps),'Hz window-',dpath(sss-1:sss+6),'_',S(i).name(1:7),'_',num2str(i),'_',num2str(j),'.txt'],'temp','-ascii')
        end
    end
    
     for i = 1:nfile
        for j =1:N(i)
            cd(fpath)
            Norsmothz{i}= medfilt1(Cutz{i}{j},fps/fs2);
            temp = Norsmothz{i}(floor(fps/(fs2*3)):end-floor(fps/(fs2*3)));
            save(['20Hz window-',dpath(sss-1:sss+6),'_',S(i).name(1:7),'_',num2str(i),'_',num2str(j),'.txt'],'temp','-ascii')
        end
    end
    
    
    for i = 1:nfile
        for j =1:N(i)
            cd(fpath)
            Norsmothz{i}= medfilt1(Cutz{i}{j},fps/fs3);
            temp = Norsmothz{i}(floor(fps/(fs3*3)):end-floor(fps/(fs3*3)));
            save(['10Hz window-',dpath(sss-1:sss+6),'_',S(i).name(1:7),'_',num2str(i),'_',num2str(j),'.txt'],'temp','-ascii')
        end
    end
    cd(cpath)
end

if 0
    
    clear all
    close all
    clc
    
    cpath = 'C:\Users\tylab\Desktop\HK_R';
    dpath = 'C:\Users\tylab\Desktop\HK_R\HMM_HS\Reshaping_WT (HS + 30% PG)\5pN\60Hz';
    epath = 'C:\Users\tylab\Desktop\HK_R\HMM_HS\Reshaping_WT (HS)\5pN\16.67-ms window';
    
    %%%% Intermediate --> scaling = 50nm
    %%%% Refolding --> scaling = 24nm at 5pN
    %%
    saving = 1;
    %%
    
    refoldingz = 0.15;
    notrefoldingz = 0.85;
    
    cd(dpath);
    RTrawhs=dir(['*.txt']);
    Shs = struct(RTrawhs);
    nfilehs = size(RTrawhs,1);
    
    RThs = {};
    
    
    for i = 1:nfilehs
        RThs{i}=load(Shs(i).name);
    end
    
    
    cd(epath);
    RTraw=dir(['*.txt']);
    S = struct(RTraw);
    nfile = size(RTraw,1);
    
    RT = {};
    
    for i = 1:nfile
        RT{i}=load(S(i).name);
    end
    
    cd(cpath);
    
    
    RT_total = cell2mat(RT');
    RT_totalhs = cell2mat(RThs');
    
    figure(1)
    
    subplot(2,1,1)
    plot((1:length(RT_total))./ 1200,RT_total,'k','LineWidth',1.5);
    hold on;
    plot((1:length(RT_total))./1200,medfilt1(RT_total,10),'y','LineWidth',1.5);
    grid on;
    plot(get(gca,'xlim'), [refoldingz refoldingz],'r','LineWidth',4);
    hold on;
    plot(get(gca,'xlim'),[notrefoldingz notrefoldingz],'b','LineWidth',4);
    set(gca, 'fontsize', 16, 'linewidth', 1, 'fontweight', 'bold','TickLength',[0.02 0.02]);
    xlabel('Time (s)'); ylabel('Normailized extension');
    %title('Normal MT-60Hz tracking')
    legend('60Hz raw trace','167-ms median-filtered','Location','Best')
    legend(gca,'boxoff')
    subplot(2,1,2)
    plot((1:length(RT_totalhs))./1200,RT_totalhs,'k','LineWidth',1.5);
    hold on;
    plot((1:length(RT_totalhs))./1200,medfilt1(RT_totalhs,120),'y','LineWidth',1.5);
    grid on;
    plot(get(gca,'xlim'), [refoldingz refoldingz],'r','LineWidth',4);
    hold on;
    plot(get(gca,'xlim'),[notrefoldingz notrefoldingz],'b','LineWidth',4);
    set(gca, 'fontsize', 16, 'linewidth', 1, 'fontweight', 'bold','TickLength',[0.02 0.02]);
    title('High-speed MT-1200Hz tracking')
    xlabel('Time (s)'); ylabel('Normailized extension');
    legend('60Hz-smoothed trace from 1200Hz raw trace','100-ms median-filtered (median-filter twice)','Location','Best')
    legend(gca,'boxoff')
    
    figure(2)
    plot((1:length(RT_total))./1200,RT_total,'g','LineWidth',1.5);
    hold on;
    plot((1:length(RT_totalhs))./1200,RT_totalhs,'b','LineWidth',1.5);
    hold on;
    grid on;
    set(gca, 'fontsize', 16, 'linewidth', 1, 'fontweight', 'bold','TickLength',[0.02 0.02]);
    
    xlabel('Time (s)'); ylabel('Normailized extension');
end
