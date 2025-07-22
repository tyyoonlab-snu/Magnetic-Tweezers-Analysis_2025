
%%%%%%% Baum-Welch algorithm  %%%%%%%%
%
if 1
    clear all
    close all
    clc
end
%% Explain code
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initiation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[prior0, transmat0, mixmat0, mu0, sigma0] =  init_mhmm(data, Q, M,'full',0)

% init_mhmm(data, Q, M, cov_type, left_right)
% INIT_MHMM Compute initial param. estimates for an HMM with mixture of Gaussian outputs.
% [init_state_prob, transmat, obsmat, mixmat, mu, Sigma] = init_mhmm(data, Q, M, cov_type, left_right)
%
% Inputs:
% data{1}(:,t) = observation vector at time t in sequence l
% Q = num. hidden states
% M = num. mixture components
% cov_type = 'full', 'diag' or 'spherical'
% left_right = 1 if the model is a left-to-right HMM, 0 otherwise

%%%%%%%%%%%%%%%%%%%%%%%%% Generate input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Learning Hmm for specifying initial input parameters

% [LL1, prior1, transmat1, mu1, sigma1, mixmat1, muhis] = learn_mhmm(data, prior0, transmat0, mu0, sigma0, mixmat0, 'max_iter', 'thresh', 'verbose', 'cov_type', 'static')

% [LL, PRIOR, TRANSMAT, MU, SIGMA, MIXMAT, MU_HISTORY] = LEARN_MHMM(DATA, PRIOR0, TRANSMAT0,
% MU0, SIGMA0, MIXMAT0) computes the ML estimates of the following parameters,
% where, for each time t, Q(t) is the hidden state, M(t) is the mixture component, and Y(t) is the observation.
%   prior(i) = Pr(Q(1) = i),
%   transmat(i,j) = Pr(Q(t+1)=j | Q(t)=i)
%   mixmat(j,k) = Pr(M(t)=k | Q(t)=j)
%   mu(:,j,k) = E[Y(t) | Q(t)=j, M(t)=k ]
%   Sigma(:,:,j,k) = Cov[Y(t) | Q(t)=j, M(t)=k]
% PRIOR0 is the initial estimate of PRIOR, etc.
% To learn an HMM with a single Gaussian output, just set mixmat = ones(Q,1).

%%%%%%%%%%%%%%%%%%%  Analysis using HMM in detail  %%%%%%%%%%%%%%%%%%%%%%%%

%[LL, prior, transmat, mu, sigma, mixmat] = mhmm_em(data, prior1, transmat1, mu1, sigma1,[], 'max_iter', 400, 'thresh', 1e-4, 'verbose', 1,'cov_type', 'full', 'adj_prior', 1, 'adj_trans', 1, 'adj_mix', 1,'adj_mu', 1, 'adj_Sigma', 1);
% , 'verbose', 1,'cov_type','full','max_iter', 500,'adj_prior',1,'adj_mix',1,'adj_mu',1,'adj_Sigma',1, 'adj_trans', 1
% LEARN_MHMM Compute the ML parameters of an HMM with (mixtures of) Gaussians output using EM.
% [ll_trace, prior, transmat, mu, sigma, mixmat] = learn_mhmm(data, ...
%   prior0, transmat0, mu0, sigma0, mixmat0, ...)
%
% Notation: Q(t) = hidden state, Y(t) = observation, M(t) = mixture variable
%
% INPUTS:
% data{ex}(:,t) or data(:,t,ex) if all sequences have the same length
% prior(i) = Pr(Q(1) = i),
% transmat(i,j) = Pr(Q(t+1)=j | Q(t)=i)
% mu(:,j,k) = E[Y(t) | Q(t)=j, M(t)=k  ]
% Sigma(:,:,j,k) = Cov[Y(t) | Q(t)=j, M(t)=k]
% mixmat(j,k) = Pr(M(t)=k | Q(t)=j) : set to [] or ones(Q,1) if only one mixture component
%
% Optional parameters may be passed as 'param_name', param_value pairs.
% Parameter names are shown below; default values in [] - if none, argument is mandatory.
%
% 'max_iter' - max number of EM iterations [10]
% 'thresh' - convergence threshold [1e-3]
% 'verbose' - if 1, print out loglik at every iteration [1]
% 'cov_type' - 'full', 'diag' or 'spherical' ['full']
%
% To clamp some of the parameters, so learning does not change them:
% 'adj_prior' - if 0, do not change prior [1]
% 'adj_trans' - if 0, do not change transmat [1]
% 'adj_mix' - if 0, do not change mixmat [1]
% 'adj_mu' - if 0, do not change mu [1]
% 'adj_Sigma' - if 0, do not change Sigma [1]
%
% If the number of mixture components differs depending on Q, just set  the trailing
% entries of mixmat to 0, e.g., 2 components if Q=1, 3 components if Q=2,
% then set mixmat(1,3)=0. In this case, B2(1,3,:)=1.0.


%% Loading data

% data(:,t,i) = observation vector at time t in sequence 'i-th'
% Q = num. hidden states
% M = num. mixture components
% numex = num. sequence (num. independent traces)
% O = Dimensionality of variable (e.g,(x,y,z) is O = 3)
% T = data length per each sequence


newpath = ''; % set path, where the traces are
cpath = ''; % set path, where the subcodes are
 
Precondition = str2num(input('Filtering process ? please, input median window size (Hz) [raw = 1200Hz] : ','s'));

cd(newpath);
disp(['Loading --',newpath,'_for HMM analysis'])

rawdata = dir('*.txt');
S = struct(rawdata);

numdata= size(S,1);
dataraw={};

for i=1:numdata
    dataraw{i} = load(S(i).name)';
end

for i = 1:numdata-1
    partition(i) = i;
end

%% Setting conditions for HMM

cd(cpath);
saving = 1;
SNRcorrection = 2;
fps = 1200; %Hz;
Qmax = 13;
for i=1:numdata
    datafiltered{i} = medfilt1(dataraw{i},round(fps./Precondition));
end

% change these preconditions!!!!!!!!!!!!!!!!!
Precondition1 = str2num(input('Do you know number of states in the system ? (yes = 1, no = otherwise): ','s'));
Precondition2 = str2num(input('Do you know the positions of states in the system ? (yes = 1, no = otherwise): ','s'));
Precondition3 = str2num(input('Which one you want to analyze ? (Trace by Trace = 1, Bead by Bead = otherwise): ','s'));

%% Partitioning by number of beads or trace
if Precondition3 ~=1
    dividingposi = find(partition~=0);
    for i = 1:length(dividingposi)+1
        if i == length(dividingposi)+1
            for j = 1:numdata-dividingposi(end)
                dividing{i}(j) = dividingposi(end)+j;
            end
        elseif i == 1
            for j = 1:dividingposi(1)
                dividing{i}(j) = j;
            end
        else
            for j = 1:dividingposi(i)-dividingposi(i-1)
                dividing{i}(j) = dividingposi(i-1)+j;
            end
        end
    end
end


%%%%%%%%%%%%%%  Model parameters    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = 1 ;  % Single Gaussian component at each state (1D information means M = 1)

mu_total = {};
mixmat_total = {};
prior_total = {};
sigma_total = {};
transmat_total = {};
Gaussprob_total = {};
HMMtrace_total = {};
AIC = {};
normalSD = zeros(numdata,1);
smoothSD = zeros(numdata,1);close
SNR = zeros(numdata,1);
SNR2 = zeros(numdata,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start HMM analysis !
if Precondition1 == 1 % Known number of state
            Q = str2num(input('Number of state: ','s'));
    %% ATP 10Hz, 7 state standard

        mu=input('What is the location of the positions? '); % e.g [0;10;40]

    kr_total = zeros(Q-1,numdata);
    ku_total = zeros(Q-1,numdata);
    Qcheck = Q;
    if Precondition2 == 1 % Known position of state
        if Precondition3 == 1
            filterwindow = zeros(numdata,1);


            mu_input = repmat(mu, 1, numdata);

            if 1
                defined_sig=[1];
                sig_input = repmat(defined_sig, numel(mu), numdata);

            end
            cd(cpath);
            stdDeviations = zeros(1, numdata); % Preallocate array for standard deviations
            for k = 1:numdata
                disp(['%%% Known state number and positions %%%']);
                disp(['%%%%% Trace',num2str(k),'- Starting analysis %%%%%%']);
                disp(['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%']);
                data ={};
                Precondition0 = Precondition;
                normalSD(k) = std(dataraw{k}(end-(2*fps/Precondition0):end)); % standard deviation
                smoothSD(k) = std(medfilt1(dataraw{k}(end-(2*fps/Precondition0):end),round(fps/Precondition)));
                SNR(k) = 4.5/normalSD(k);
                SNR2(k) = 4.5/smoothSD(k);
                if SNRcorrection
                    if SNR2(k) >= 2
                        Precondition0 =  Precondition;
                    else
                        Precondition0 = Precondition;
                        %Precondition0 = floor(Precondition0.*(SNR2(k)./2));
                    end

                    disp (['window size: ', num2str(Precondition)])

                end
                filterwindow(k)=Precondition0;
                %Precondition0 = 2;
                filtdata = medfilt1(dataraw{k},round(fps/Precondition0));
                data{1} = filtdata((fps/(3*Precondition0):end));

                [prior0, transmat0, mixmat0, mu0, sigma0] =  init_mhmm(data, Q, M,'full',0);
                disp(['Assume ',num2str(Q),'-state system_','Finish Initialization'])

                mu0 = mu_input(:,k)';
                sigma0(1,1,:) = sig_input(:,k);
                [LL1, prior1, transmat1, mu1, sigma1, mixmat1, muhis] = learn_mhmm(data, prior0, transmat0, mu0, sigma0, mixmat0, 'max_iter', 'thresh', 'verbose', 'cov_type', 'static');
                disp(['Assume ',num2str(Q),'-state system_','Finish Setting input parameters'])
                mu1 = mu0;
                sigma1(1,1,:) = sig_input(:,k);
                [LL, prior, transmat, mu, sigma, mixmat] = mhmm_em(data, prior1, transmat1, mu1, sigma1,[], 'max_iter', 400, 'thresh',1e-4, 'verbose', 1,'cov_type', 'full', 'adj_prior', 1, 'adj_trans', 1, 'adj_mix', 1,'adj_mu',0, 'adj_Sigma',0);
                disp(['Assume ',num2str(Q),'-state system_','Finish HMM running'])

                HMM_find; % Extract model parameters

                mu_total{k} = mu_result;
                mixmat_total{k}= mixmat_result;
                prior_total{k}= prior_result;
                sigma_total{k}= sigma_result;
                transmat_total{k}= transmat_result;
                Gaussprob_total{k}= Gaussprob;
                HMMtrace_total{k}= HMMtrace;
                AIC{k} = LL; % Already log scale
                BIC(k) = -2*AIC{k}(end) + Q*(Q+2)*log(length(cell2mat(data)));  % number of model parameters*log(number of data points)

                for jj = 1:Q-1
                    kr_total(jj,k) = transmat_total{k}(jj,jj+1);
                    ku_total(jj,k) = transmat_total{k}(Q-jj+1,Q-jj);
                end
                for jj = 1:Q
                    z_track(jj,k) = mu_total{k}(jj) ;
                end

            end

            rgbscale = linspace(0,1,numdata);
            rgbflat = zeros(numdata,1);
            rgbscaleinv = sort(rgbscale,'descend');
            c = [rgbscaleinv' rgbflat rgbscaleinv']; %c =[ 1 0 0 ; 1 0.6 0.2 ; 0 0.6 0.2];
            cm = flipud(c);
            for k = 1:numdata
                if 0
                    figure(1111*k);
                    box on; grid on; hold on;
                    plot(AIC{k},'r','LineWidth',3);
                    set(gca,'fontsize', 12, 'linewidth', 1, 'fontweight', 'bold','TickLength',[0.02 0.02])
                    xlabel('Number of iteration (#)');
                end


                I1 = (max(HMMtrace_total{k}{1})-min(HMMtrace_total{k}{1}))*0.7 + min(HMMtrace_total{k}{1});
                I2 = (max(HMMtrace_total{k}{1})-min(HMMtrace_total{k}{1}))*0.41 + min(HMMtrace_total{k}{1});
                expt1 = linspace(I1,I1,length(HMMtrace_total{k}{1}));
                expt2 = linspace(I2,I2,length(HMMtrace_total{k}{1}));

                figure(100+k);hold on; box on; grid on;
                plot((1:length(dataraw{k}(fps/(3*filterwindow(k)):end)))./fps,dataraw{k}(fps/(3*filterwindow(k)):end),'Color', [0.9 0.9 0.9],'LineWidth',0.2)
                plot((1:length(dataraw{k}(fps/(3*filterwindow(k)):end)))./fps,medfilt1(dataraw{k}(fps/(3    *filterwindow(k)):end),round(fps./Precondition)),'k','LineWidth',0.5)
                %plot((1:length(HMMtrace_total{k}{1}))./fps,expt1,'color',[0 1 0],'LineWidth',1)
                %plot((1:length(HMMtrace_total{k}{1}))./fps,expt2,'color',[1 0.3 0],'LineWidth',1)
                plot((1:length(HMMtrace_total{k}{1}))./fps,HMMtrace_total{k}{1},'r','LineWidth',0.5)
                set(gca,'fontsize', 12, 'linewidth', 1, 'fontweight', 'bold','TickLength',[0.02 0.02])
                set(gcf,'Renderer','painters')
                xlabel('Time (s)')
                ylabel('Extension (nm)')
                stdDeviations(k) = std(dataraw{k}(fps/(3*filterwindow(k)):end) - HMMtrace_total{k}{1});
                sc = 40;
                for i =1:Q
                    figure(1); hold on; box on; grid on;
                    scatter(k,mu_total{k}(i),sc,'ko','filled');hold on; box on; grid on;
                end
                figure(1);
                set(gca,'xtick',0:1:numdata+1)
                xlim([0 numdata+2])
                xlabel('Trace to Trace variance')
                ylabel('Extension (nm)')

                figure(2); hold on; box on; grid on;
                for jj = 1:Q-1
                    scatter(jj, kr_total(jj,k),sc,cm(k),'o','filled')
                    scatter(jj+Q-1,ku_total(jj,k),sc,cm(k),'o','filled')
                    Name=[num2str(k),'th:',S(k).name];
                    t=text (1, kr_total(1,k),Name)  ;
                end
                set(gca,'xtick',[]);
                set(gca, 'YScale', 'log')
                xlim([0 2*(Q-1)+1])
                ylabel('Kinetic rate (1/s)')
            end
            figure(1);hold on;
            for jj = 1:Q
                scatter(numdata+1,mean(z_track(jj,:)),3*sc,'ro','filled');
                errorbar(numdata+1,mean(z_track(jj,:)),std(z_track(jj,:)));
                % text(x_position, y_position, ['S(' num2str(k) ').name'], 'Color', cm(k));
            end
            xlim([0 numdata+2])

            figure(2); hold on;
            for jj = 1:Q-1
                scatter(jj, mean(kr_total(jj,:)),3*sc,'d','filled')
                scatter(jj+Q-1, mean(ku_total(jj,:)),3*sc,'d','filled')
            end
            set(gca,'xtick',[]);
            set(gca, 'YScale', 'log')
            xlim([0 2*(Q-1)+1])
            ylim([0.01 100])
            ylabel('Kinetic rate (1/s)')

            figure(3);
            scatter(1:length(filterwindow),filterwindow,5*sc,'ko')

        else

            mu=input("what is the state defined?  eg. [0; 17; 54 ;78 ;111] " );
            mu_input = repmat(mu, 1, numdata);
            cd(newpath);
            filterwindow = zeros(length(dividing),1);
            for k = 1:length(dividing)
                disp(['%%% Known state number and positions %%%']);
                disp(['%%%%% Bead',num2str(k),'- Starting analysis %%%%%%%']);
                disp(['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%']);

                Precondition0 =  Precondition;
                normalSD(k) = std(dataraw{k}(end-(2*fps/Precondition0):end));
                smoothSD(k) = std(medfilt1(dataraw{k}(end-(2*fps/Precondition0):end),round(fps/Precondition0)));
                SNR(k) = 4.5/normalSD(k);
                SNR2(k) = 4.5/smoothSD(k);

                if SNRcorrection
                    if SNR2(k)>= 2
                        Precondition0 = 4;
                    else
                        Precondition0 = 4;
                    end
                end
                filterwindow(k)=Precondition0;

                data ={};
                for i = 1:length(dividing{k})
                    datarawarr{k}{i} = dataraw{dividing{k}(i)}(fps/(3*Precondition0):end);
                    data{i} = medfilt1(dataraw{dividing{k}(i)}(fps/(3*Precondition0):end),round(fps/Precondition0));
                end

                [prior0, transmat0, mixmat0, mu0, sigma0] =  init_mhmm(data, Q, M,'full',0);
                disp(['Assume ',num2str(Q),'-state system_','Finish Initialization'])

                mu_input2 = 0;
                mu_input2 = mu_input(dividing{k},:)
                mu0 = [mean(mu_input2(:,1)) mean(mu_input2(:,2)) mean(mu_input2(:,3)) mean(mu_input2(:,4))];


                [LL1, prior1, transmat1, mu1, sigma1, mixmat1, muhis] = learn_mhmm(data, prior0, transmat0, mu0, sigma0, mixmat0, 'max_iter', 'thresh', 'verbose', 'cov_type', 'static');
                disp(['Assume ',num2str(Q),'-state system_','Finish Setting input parameters'])
                mu1 = mu0;

                [LL, prior, transmat, mu, sigma, mixmat] = mhmm_em(data, prior1, transmat1, mu1, sigma1,[], 'max_iter', 400, 'thresh',1e-4, 'verbose', 1,'cov_type', 'full', 'adj_prior', 1, 'adj_trans', 1, 'adj_mix', 1,'adj_mu',0, 'adj_Sigma', 1);
                disp(['Assume ',num2str(Q),'-state system_','Finish HMM running'])

                HMM_find; % Extract model parameters

                mu_total{k} = mu_result;
                mixmat_total{k}= mixmat_result;
                prior_total{k}= prior_result;
                sigma_total{k}= sigma_result;
                transmat_total{k}= transmat_result;
                Gaussprob_total{k}= Gaussprob;
                HMMtrace_total{k}= HMMtrace;
                AIC{k} = LL; % Already log scale
                BIC(k) = -2*AIC{k}(end) + Q*(Q+2)*log(length(cell2mat(data)));  % number of model parameters*log(number of data points)
                for jj = 1:Q-1
                    kr_total(jj,k) = transmat_total{k}(jj,jj+1);
                    ku_total(jj,k) = transmat_total{k}(Q-jj+1,Q-jj);
                end
            end

            for k = 1:length(dividing)
                if 0
                    figure(1111*k);
                    box on; grid on; hold on;
                    plot(AIC{k},'r','LineWidth',3);
                    set(gca,'fontsize', 12, 'linewidth', 1, 'fontweight', 'bold','TickLength',[0.02 0.02])
                    xlabel('Number of iteration (#)');
                end

                for i = 1:size(HMMtrace_total{k},2)
                    figure(100*k+i);hold on; box on; grid on;
                    plot((1:length(datarawarr{k}{i}))./fps,datarawarr{k}{i},'Color', [17 17 17]/255,'LineWidth',0.25)
                    plot((1:length(datarawarr{k}{i}))./fps,medfilt1(datarawarr{k}{i},round(fps./Precondition)),'k','LineWidth',0.5)
                    plot((1:length(HMMtrace_total{k}{i}))./fps,HMMtrace_total{k}{i},'r','LineWidth',0.5)
                    set(gca,'fontsize', 12, 'linewidth', 1, 'fontweight', 'bold','TickLength',[0.02 0.02])
                end

                sc = 40;
                for i =1:Q
                    figure(1); hold on; box on; grid on;
                    scatter(k,mu_total{k}(i),sc,'ko','filled');hold on; box on; grid on;
                end
                figure(1);
                set(gca,'xtick',0:1:length(dividing)+1)
                xlim([0 length(dividing)+1])
                xlabel('Bead to Bead variance')
                ylabel('Extension (nm)')

                figure(2); hold on; box on; grid on;
                for jj = 1:Q-1
                    scatter(jj, kr_total(jj,k),sc,'o','filled')
                    scatter(jj+Q-1,ku_total(jj,k),sc,'o','filled')
                end
                set(gca,'xtick',[]);
                set(gca, 'YScale', 'log')
                xlim([0 7])
                ylabel('Kinetic rate (1/s)')
            end
            figure(2); hold on;
            for jj = 1:Q-1
                scatter(jj, mean(kr_total(jj,:)),3*sc,'d','filled')
                scatter(jj+Q-1, mean(ku_total(jj,:)),3*sc,'d','filled')
            end
            set(gca,'xtick',[]);
            set(gca, 'YScale', 'log')
            xlim([0 7])
            ylabel('Kinetic rate (1/s)')

            figure(3);
            scatter(1:length(filterwindow),filterwindow,5*sc,'ko')
        end

    else % Unknown position of state
        if Precondition3 == 1
            filterwindow = zeros(numdata,1);
            for k = 1:numdata
                disp(['%%% Known state number and but unknown positions %%%']);
                disp(['%%%%% Trace',num2str(k),'- Starting analysis %%%%%%']);
                disp(['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%']);
                data ={};
                Precondition0 = Precondition;
                normalSD(k) = std(dataraw{k}(end-(2*fps/Precondition0):end));
                smoothSD(k) = std(medfilt1(dataraw{k}(end-(2*fps/Precondition0):end),round(fps/Precondition0)));
                SNR(k) = 4.5/normalSD(k);
                SNR2(k) = 4.5/smoothSD(k);

                if SNRcorrection
                    if SNR2(k)>= 2
                        Precondition0 =  Precondition0 ;
                    else
                        Precondition0 =  Precondition0 ;
                        % Precondition0 = floor(Precondition0.*(SNR2(k)./2))
                    end
                end
                filterwindow(k)=Precondition0;
                data{1} = medfilt1(dataraw{k}(fps/(3*Precondition0):end),round(fps./Precondition0));
                [prior0, transmat0, mixmat0, mu0, sigma0] =  init_mhmm(data, Q, M,'full',0);
                disp(['Assume ',num2str(Q),'-state system_','Finish Initialization'])

                [LL1, prior1, transmat1, mu1, sigma1, mixmat1, muhis] = learn_mhmm(data, prior0, transmat0, mu0, sigma0, mixmat0, 'max_iter', 'thresh', 'verbose', 'cov_type', 'static');
                disp(['Assume ',num2str(Q),'-state system_','Finish Setting input parameters'])

                [LL, prior, transmat, mu, sigma, mixmat] = mhmm_em(data, prior1, transmat1, mu1, sigma1,[], 'max_iter', 400, 'thresh',1e-4, 'verbose', 1,'cov_type', 'full', 'adj_prior', 1, 'adj_trans', 1, 'adj_mix', 1,'adj_mu',1, 'adj_Sigma', 1);
                disp(['Assume ',num2str(Q),'-state system_','Finish HMM running'])

                HMM_find; % Extract model parameters

                mu_total{k} = mu_result;
                mixmat_total{k}= mixmat_result;
                prior_total{k}= prior_result;
                sigma_total{k}= sigma_result;
                transmat_total{k}= transmat_result;
                Gaussprob_total{k}= Gaussprob;
                HMMtrace_total{k}= HMMtrace;
                AIC{k} = LL; % Already log scale
                BIC(k) = -2*AIC{k}(end) + Q*(Q+2)*log(length(cell2mat(data)));  % number of model parameters*log(number of data points)
                for jj = 1:Q-1
                    kr_total(jj,k) = transmat_total{k}(jj,jj+1);
                    ku_total(jj,k) = transmat_total{k}(Q-jj+1,Q-jj);
                end

                for jj = 1:Q
                    z_track(jj,k) = mu_total{k}(jj) ;
                end
            end

            for k = 1:numdata
                if 0
                    figure(1111*k);
                    box on; grid on; hold on;
                    plot(AIC{k},'r','LineWidth',3);
                    set(gca,'fontsize', 12, 'linewidth', 1, 'fontweight', 'bold','TickLength',[0.02 0.02])
                    xlabel('Number of iteration (#)');
                end
                figure(100+k);hold on; box on; grid on;
                plot((1:length(dataraw{k}(fps/(3*filterwindow(k)):end)))./fps,dataraw{k}(fps/(3*filterwindow(k)):end),'Color', [0.80,0.80,0.80],'LineWidth',0.25)
                plot((1:length(dataraw{k}(fps/(3*filterwindow(k)):end)))./fps,medfilt1(dataraw{k}(fps/(3*filterwindow(k)):end),round(fps./Precondition)),'k','LineWidth',0.5)
                plot((1:length(HMMtrace_total{k}{1}))./fps,HMMtrace_total{k}{1},'r','LineWidth',0.5)
                set(gca,'fontsize', 12, 'linewidth', 1, 'fontweight', 'bold','TickLength',[0.02 0.02])
                sc = 40;
                for i =1:Q
                    figure(1); hold on; box on; grid on;
                    scatter(k,mu_total{k}(i),sc,'ko','filled');hold on; box on; grid on;
                end
                figure(1);
                set(gca,'xtick',0:1:numdata+1)
                xlim([0 numdata+2])
                xlabel('Trace to Trace variance')
                ylabel('Extension (nm)')

                figure(2); hold on; box on; grid on;
                for jj = 1:Q-1
                    scatter(jj, kr_total(jj,k),sc,'o','filled')
                    scatter(jj+Q-1,ku_total(jj,k),sc,'o','filled')
                end
                set(gca,'xtick',[]);
                set(gca, 'YScale', 'log')
                xlim([0 2*(Q-1)+1])
                ylabel('Kinetic rate (1/s)')
            end
            figure(1);hold on;
            for jj = 1:Q
                scatter(numdata+1,mean(z_track(jj,:)),3*sc,'ro','filled');
                errorbar(numdata+1,mean(z_track(jj,:)),std(z_track(jj,:)));
            end
            xlim([0 numdata+2])

            figure(2); hold on;
            for jj = 1:Q-1
                scatter(jj, mean(kr_total(jj,:)),3*sc,'d','filled')
                scatter(jj+Q-1, mean(ku_total(jj,:)),3*sc,'d','filled')
            end
            set(gca,'xtick',[]);
            set(gca, 'YScale', 'log')
            xlim([0 2*(Q-1)+1])
            ylabel('Kinetic rate (1/s)')

            figure(3);
            scatter(1:length(filterwindow),filterwindow,5*sc,'ko')

        else
            filterwindow = zeros(length(dividing),1);
            for k = 1:length(dividing)
                disp(['%%% Known state number and but unknown positions %%%']);
                disp(['%%%%% Bead',num2str(k),'- Starting analysis %%%%%%%']);
                disp(['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%']);

                Precondition0 =  Precondition
                normalSD(k) = std(dataraw{k}(end-(2*fps/Precondition0):end));
                smoothSD(k) = std(medfilt1(dataraw{k}(end-(2*fps/Precondition0):end),round(fps/Precondition0)));
                SNR(k) = 4.5/normalSD(k);
                SNR2(k) = 4.5/smoothSD(k);

                if SNRcorrection
                    if SNR2(k)>= 2
                        Precondition0 =  Precondition0;
                    else
                        Precondition0 =  Precondition0;
                        % Precondition0 = floor(Precondition0.*(SNR2(k)./2))
                    end
                end
                filterwindow(k)=Precondition0;

                data ={};
                for i = 1:length(dividing{k})
                    datarawarr{k}{i} = dataraw{dividing{k}(i)}(fps/(3*Precondition0):end);
                    data{i} = medfilt1(dataraw{dividing{k}(i)}(fps/(3*Precondition0):end),round(fps/Precondition0));
                end

                [prior0, transmat0, mixmat0, mu0, sigma0] =  init_mhmm(data, Q, M,'full',0);
                disp(['Assume ',num2str(Q),'-state system_','Finish Initialization'])

                [LL1, prior1, transmat1, mu1, sigma1, mixmat1, muhis] = learn_mhmm(data, prior0, transmat0, mu0, sigma0, mixmat0, 'max_iter', 'thresh', 'verbose', 'cov_type', 'static');
                disp(['Assume ',num2str(Q),'-state system_','Finish Setting input parameters'])

                [LL, prior, transmat, mu, sigma, mixmat] = mhmm_em(data, prior1, transmat1, mu1, sigma1,[], 'max_iter', 400, 'thresh', 1e-4, 'verbose', 1,'cov_type', 'full', 'adj_prior', 1, 'adj_trans', 1, 'adj_mix', 1,'adj_mu',1, 'adj_Sigma', 1);
                disp(['Assume ',num2str(Q),'-state system_','Finish HMM running'])

                HMM_find; % Extract model parameters

                mu_total{k} = mu_result;
                mixmat_total{k}= mixmat_result;
                prior_total{k}= prior_result;
                sigma_total{k}= sigma_result;
                transmat_total{k}= transmat_result;
                Gaussprob_total{k}= Gaussprob;
                HMMtrace_total{k}= HMMtrace;
                AIC{k} = LL; % Already log scale
                BIC(k) = -2*AIC{k}(end) + Q*(Q+2)*log(length(cell2mat(data)));  % number of model parameters*log(number of data points)

                for jj = 1:Q-1
                    kr_total(jj,k) = transmat_total{k}(jj,jj+1);
                    ku_total(jj,k) = transmat_total{k}(Q-jj+1,Q-jj);

                    for jj = 1:Q
                        z_track(jj,k) = mu_total{k}(jj) ;
                    end
                end

                for k = 1:length(dividing)
                    figure(1111*k);
                    box on; grid on; hold on;
                    plot(AIC{k},'r','LineWidth',3);
                    set(gca,'fontsize', 12, 'linewidth', 1, 'fontweight', 'bold','TickLength',[0.02 0.02])
                    xlabel('Number of iteration (#)');

                    for i = 1:size(HMMtrace_total{k},2)
                        figure(100*k+i);hold on; box on; grid on;
                        plot((1:length(datarawarr{k}{i}))./fps,datarawarr{k}{i},'Color', [0.80,0.80,0.80],'LineWidth',0.25)
                        plot((1:length(datarawarr{k}{i}))./fps,medfilt1(datarawarr{k}{i},round(fps./Precondition)),'k','LineWidth',0.5)
                        plot((1:length(HMMtrace_total{k}{i}))./fps,HMMtrace_total{k}{i},'r','LineWidth',0.5)
                        set(gca,'fontsize', 12, 'linewidth', 1, 'fontweight', 'bold','TickLength',[0.02 0.02])
                    end
                    sc = 40;
                    for i =1:Q
                        figure(1); hold on; box on; grid on;
                        scatter(k,mu_total{k}(i),sc,'ko','filled');hold on; box on; grid on;
                    end
                    figure(1);
                    set(gca,'xtick',0:1:length(dividing)+1)
                    xlim([0 length(dividing)+1])
                    xlabel('Bead to Bead variance')
                    ylabel('Extension (nm)')

                    figure(2); hold on; box on; grid on;
                    for jj = 1:Q-1
                        scatter(jj, kr_total(jj,k),sc,'o','filled')
                        scatter(jj+Q-1,ku_total(jj,k),sc,'o','filled')
                    end
                end

                figure(2); hold on;
                for jj = 1:Q-1
                    scatter(jj, mean(kr_total(jj,:)),3*sc,'d','filled')
                    scatter(jj+Q-1, mean(ku_total(jj,:)),3*sc,'d','filled')
                end
                set(gca,'xtick',[]);
                set(gca, 'YScale', 'log')
                xlim([0 2*(Q-1)+1])
                ylabel('Kinetic rate (1/s)')

                figure(3);
                scatter(1:length(filterwindow),filterwindow,5*sc,'ko')
            end
        end
    end
else
    % Unknown number of state
    data = datafiltered;
    for Q = 1:Qmax
        disp(['%% Nothing known; Need a general Inspection %%']);
        disp(['%%%%% If the system has ',num2str(Q),' states - Starting analysis %%%']);
        disp(['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%']);

        Qcheck(Q) = Q;
        [prior0, transmat0, mixmat0, mu0, sigma0] =  init_mhmm(data, Q, M,'full',0);
        disp(['Assume ',num2str(Q),'-state system_','Finish Initialization'])

        [LL1, prior1, transmat1, mu1, sigma1, mixmat1, muhis] = learn_mhmm(data, prior0, transmat0, mu0, sigma0, mixmat0, 'max_iter', 'thresh', 'verbose', 'cov_type', 'static');
        disp(['Assume ',num2str(Q),'-state system_','Finish Setting input parameters'])

        [LL, prior, transmat, mu, sigma, mixmat] = mhmm_em(data, prior1, transmat1, mu1, sigma1,[], 'max_iter', 400, 'thresh',1e-4, 'verbose', 1,'cov_type', 'full', 'adj_prior', 1, 'adj_trans', 1, 'adj_mix', 1,'adj_mu',1, 'adj_Sigma', 1);
        disp(['Assume ',num2str(Q),'-state system_','Finish HMM running'])

        HMM_find % Extract model parameters

        mu_total{Q} = mu_result;
        mixmat_total{Q}= mixmat_result;
        prior_total{Q}= prior_result;
        sigma_total{Q}= sigma_result;
        transmat_total{Q}= transmat_result;
        Gaussprob_total{Q}= Gaussprob;
        HMMtrace_total{Q}= HMMtrace;
        AIC{Q} = LL; % Already log scale
        BIC(Q) = -2*AIC{Q}(end) + Q*(Q+2)*log(length(cell2mat(data)));  % number of model parameters*log(number of data points)
    end

    sc=250;
    figure(7777);
    for p =1:Qmax
        % for p =1:8
        box on; grid on; hold on;
        scatter(p,BIC(p),sc,'ro','filled');hold on;
    end
    plot(1:Qmax,BIC(1:Qmax),'Color', [0.80,0.80,0.80],'LineWidth',1)
    set(gca,'fontsize', 12,'xtick',1:1:Qmax+1, 'linewidth', 1, 'fontweight', 'bold','TickLength',[0.02 0.02])
    xlabel('Number of states (#)');
    xlim([1 Qmax+1])
    %ylim([6e+6 7e+6])
    if size(dataraw,2) > ceil(sqrt(size(dataraw,2)))*floor(sqrt(size(dataraw,2)))
        psx = ceil(sqrt(size(dataraw,2)))
        psy = ceil(sqrt(size(dataraw,2)))
    elseif  size(dataraw,2) < ceil(sqrt(size(dataraw,2)))*floor(sqrt(size(dataraw,2)))
        psx = ceil(sqrt(size(dataraw,2)))
        psy = floor(sqrt(size(dataraw,2)))
    elseif size(dataraw,2) == ceil(sqrt(size(dataraw,2)))*floor(sqrt(size(dataraw,2)))
        psx = ceil(sqrt(size(dataraw,2)))
        psy = floor(sqrt(size(dataraw,2)))
    end

    for p = 1:Qmax
        figure(p)
        for i = 1:size(HMMtrace,2)
            subplot(psx,psy,i); hold on; box on; grid on;
            plot((1:length(dataraw{i}))./fps,dataraw{i},'Color', [0.80,0.80,0.80],'LineWidth',0.25)
            plot((1:length(data{i}))./fps,data{i},'k','LineWidth',0.5)
            plot((1:length(HMMtrace_total{p}{i}))./fps,HMMtrace_total{p}{i},'r','LineWidth',0.5)
            set(gca,'fontsize', 12, 'linewidth', 1, 'fontweight', 'bold','TickLength',[0.02 0.02])
            minz_set{p}(i) = min(dataraw{i});
            maxz_set{p}(i) = max(dataraw{i});
        end
        figure(p);
        for i = 1:size(HMMtrace,2)
            subplot(psx,psy,i)
            ylim([floor(min(minz_set{p}))-2 ceil(max(maxz_set{p}))+2])
        end
        figure(6666);hold on;
        for j = 1:p
            scatter(p,mu_total{p}(j),sc,'bd','filled');hold on; box on; grid on;
        end
        figure(6666); axis([0 Qmax+1 0 80])
    end
end



new_folder_name = ['k=',num2str(Q),'_W=', num2str(Precondition),'_endreal_HMM_',date];
disp (new_folder_name)
mkdir(fullfile(newpath, new_folder_name));
% Update dpath to point to the new folder
d_in_path = fullfile(newpath, new_folder_name);


%% Extract model parameters
if saving
    cd(d_in_path);
    save(['Mean_by ',num2str(Precondition),'Hz window','.mat'],'mu_total','-mat')

    save(['Mixmat_by ',num2str(Precondition),'Hz window','.mat'],'mixmat_total','-mat')

    save(['Prior_by ',num2str(Precondition),'Hz window','.mat'],'prior_total','-mat')

    save(['Sigma_by ',num2str(Precondition),'Hz window','.mat'],'sigma_total','-mat')

    save(['Kinetic matrix_by ',num2str(Precondition),'Hz window','.mat'],'transmat_total','-mat')

    save(['Gaussprob_by',num2str(Precondition),'Hz window','.mat'],'Gaussprob_total','-mat')

    save(['Stepresults_by ',num2str(Precondition),'Hz window','.mat'],'HMMtrace_total','-mat')

    save(['AIC_by ',num2str(Precondition),'Hz window','.mat'],'AIC','-mat')

    save(['BIC_by ',num2str(Precondition),'Hz window','.mat'],'BIC','-mat')

end



