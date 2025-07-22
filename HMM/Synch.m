clear all
close all
clc

cpath = 'C:\Users\Owner\Desktop\GlpG\HMM_MJ';
dpath = 'C:\Users\Owner\Desktop\GlpG\HMM_MJ\30% PG\F121EL133E\Raw\6pN';
fpath = 'C:\Users\Owner\Desktop\GlpG\HMM_MJ\Synchronized\F121EL133E\6pN'; % 5 Hz
%%%% Intermediate --> scaling = 50nm
%%%% Refolding --> scaling = 24nm at 5pN
%%
saving = 1;
correction = 1;
%%
cd(dpath);
RTraw=dir(['*.txt']);
S = struct(RTraw);
nfile = size(RTraw,1);

fps = 1200;
fs = 10;

RTz = {};
smoothRTz={};
Transz ={};
Norz = {};
Cutz = {};

for i = 1:nfile
    RTz{i}=load(S(i).name);
end

cd(cpath)


for j = 1:nfile
    Transz{j} = RTz{j}
    Norz{j} = (Transz{j});
end

tis = {};
refoldingz = 20;
notrefoldingz = 38;

if 1
    for i = 1:nfile
        if correction
            figure(100+i);
            plot(Norz{i},'Color',[0.7 0.7 0.7],'LineWidth',1);
            hold on;
            axis([0 length(Norz{i}) min(Norz{i})-10 max(Norz{i})+10]);
            plot(get(gca,'xlim'), [refoldingz refoldingz],'r','LineWidth',4);
            hold on;
            plot(get(gca,'xlim'),[notrefoldingz notrefoldingz],'b','LineWidth',4);
            hold on;
            plot(medfilt1(Norz{i},fps/fs),'y','LineWidth',1.25);
            
            axis([0 length(Norz{i}) min(Norz{i})-10 max(Norz{i})+10]);
            grid on;
            set(gca, 'fontsize', 16, 'linewidth', 1, 'fontweight', 'bold','TickLength',[0.02 0.02]);
            xlabel('Time (s)'); ylabel('Normalized extension (nm)');
            pause;
            getpts
            tis3{i} = ans;
            
            zc1 = round(tis3{i}(1));
            zc2 = round(tis3{i}(2));
            zc3 = round(tis3{i}(end-1));
            zc4 = round(tis3{i}(end));
            zdrift = (mean(Norz{i}(zc3:zc4))-mean(Norz{i}(zc1:zc2)))./(zc4-zc1) ;% nm/frame(#)
            for loc = 1:length(Norz{i})
                Norz{i}(loc) = Norz{i}(loc) - zdrift*loc ;
            end
        end
        pause;
        getpts
        tis{i} = ans;
        ao= tis{i}(1);
        co = tis{i}(end);
        Mnorz(i) = mean(Norz{i}(ao:co));
        Norz{i} = Norz{i}- (Mnorz(i) - refoldingz);
        
        Cutz{i} = Norz{i};
        
        close(100+i);
        figure(200+i);
        plot(Norz{i},'Color',[0.7 0.7 0.7],'LineWidth',1);
        hold on;
        axis([0 length(Norz{i}) min(Norz{i})-10 max(Norz{i})+10]);
        plot(get(gca,'xlim'), [refoldingz refoldingz],'r','LineWidth',4);
        hold on;
        plot(get(gca,'xlim'),[notrefoldingz notrefoldingz],'b','LineWidth',4);
        hold on;
        plot(medfilt1(Norz{i},fps/fs),'y','LineWidth',1.25);
        
        axis([0 length(Norz{i}) min(Norz{i})-10 max(Norz{i})+10]);
        grid on;
        set(gca, 'fontsize', 16, 'linewidth', 1, 'fontweight', 'bold','TickLength',[0.02 0.02]);
        xlabel('Time (s)'); ylabel('Normalized extension (nm)');
    end     
end

C =[0:0.05:1; 0:0.05:1; 0:0.05:1]';

if 1
    for j = 1:nfile
        tdomain(j) = length(Cutz{j});
    end
    for i= 1:nfile
        if i == 1
            figure(100);
            plot((1:length(Cutz{i})),Cutz{i},'LineWidth',1);
            hold on;
            plot((1:length(Cutz{i})),medfilt1(Cutz{i},fps/fs),'y','LineWidth',1.25);
            grid on;
            set(gca, 'fontsize', 16, 'linewidth', 1, 'fontweight', 'bold','TickLength',[0.02 0.02]);
        else
            figure(100);
            tdomain2 = sum(tdomain(1:i-1));
            plot((tdomain2+1:tdomain2+length(Cutz{i})),Cutz{i},'LineWidth',1);
            hold on;
            plot((tdomain2+1:tdomain2+length(Cutz{i})),medfilt1(Cutz{i},fps/fs),'y','LineWidth',1.25);
            grid on;
            set(gca, 'fontsize', 16, 'linewidth', 1, 'fontweight', 'bold','TickLength',[0.02 0.02]);
        end
    end
    axis([0 sum(tdomain) 0 60]);
    xlabel('Frame (#)'); ylabel('Normailized extension');
    plot(get(gca,'xlim'), [refoldingz refoldingz],'r','LineWidth',4);
    hold on;
    plot(get(gca,'xlim'),[notrefoldingz notrefoldingz],'b','LineWidth',4);
end



sss = find(dpath == '8');

if 0
    for i = 1:nfile
        cd(fpath)
        Norsmothz{i}= Cutz{i};
        temp = Norsmothz{i};
        save([S(i).name,'.txt'],'temp','-ascii')
    end
end
pause;
cd(cpath)
