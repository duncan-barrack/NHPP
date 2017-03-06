clear
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Example 3 - Sinusoidal and linear rate function%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %set parameters
t1 = 0; %LH bound of time interval
t2 = 4 * pi; %RH bound of time interva
nb = 100; %number of cubic basis splines
t = linspace(0, t2, 1000)'; %time vector (used for plotting NHPP rate functions)
ns = 20; %no. of samples for each class
labels=[ones(1, ns), 2*ones(1, ns)];% class lables
opts = optimoptions(@fmincon, 'Display', 'iter', 'GradObj', 'on','MaxFunEvals', 1e7, 'TolFun', 1e-20, 'MaxIter', 1000);  %options for fmincon solver


%rate function for class 1
lambdat = @(x) 100*sin(x/2).^2; %define rate function for NHPP
rng(1)
a1 = (lambdat(linspace(t1, t2, 100)) + 10*randn(1, 100))'; 
lambdaTrue1=lambdat(t); 

%generate some data
data1 = arrayfun(@(x) NHPP(lambdat, t1, t2), 1:ns, 'UniformOutput', false);

%rate function for class 2
lambdat = @(x) 100*sin(x).^2; %define rate function for NHPP
rng(1)
a2 = (lambdat(linspace(t1, t2, 100)) + 10*randn(1, 100))';
lambdaTrue2=lambdat(t);

%generate some data
data2 = arrayfun(@(x) NHPP(lambdat, t1, t2), 1:ns, 'UniformOutput', false);

%rate function for class 3
lambdat = @(x) ((x<t2/4)* 20 + (x<t2/2 & x>= t2/4)* 40 + (x<3*t2/4 & x>= t2/2)* 60 + (x>=3*t2/4)*80); %define rate function for NHPP
rng(1)
a3 = (lambdat(linspace(t1, t2, 100)) + 10*randn(1, 100))';  
lambdaTrue3=lambdat(t); 

%generate some data
data3 = arrayfun(@(x) NHPP(lambdat, t1, t2), 1:ns, 'UniformOutput', false);

%rate function for class 4
lambdat = @(x) ((x<t2/4)* 80 + (x<t2/2 & x>= t2/4)* 60 + (x<3*t2/4 & x>= t2/2)* 40 + (x>=3*t2/4)*20); 
rng(1)
a4 = (lambdat(linspace(t1, t2, 100)) + 10*randn(1, 100))';
lambdaTrue4=lambdat(t);

%generate some data
data4 = arrayfun(@(x) NHPP(lambdat, t1, t2), 1:ns, 'UniformOutput', false);

% combine data
data = [data1, data2, data3, data4];

% cluster the data
[sp_fn, mp, q] = NHPP_cluster(data, 4, t1, t2, opts, nb, [a1, a2, a3, a4], [0.25, 0.25, 0.25, 0.25], 10, []);

lambdaHat1 = fnval(sp_fn(1), t); %NHPP estimate for class 1
lambdaHat2 = fnval(sp_fn(2), t); %NHPP estimate for class 2
lambdaHat3 = fnval(sp_fn(3), t); %NHPP estimate for class 3
lambdaHat4 = fnval(sp_fn(4), t); %NHPP estimate for class 4

%%% Plot Results

%%% Raster plots
%%% Class 1
figure
for i=1:ns
    samp = datasample(data1{i}, 200);
    plot([samp;samp],[i*ones(size(samp));(0+(i-1)).*ones(size(samp))] +0.5,'k-')
    hold on    
end
xlim([t1, t2])
ylim([0.5,(ns +0.5)])
xlabel('time', 'FontSize',18,'Interpreter','latex')
ylabel('Sample no.', 'FontSize',18,'Interpreter','latex')
set(gca,'fontsize',16)
title('Class 1 event times','FontSize',18,'Interpreter','latex')

%%% Class 2
figure
for i=1:ns
    samp = datasample(data2{i}, 200);
    plot([samp;samp],[i*ones(size(samp));(0+(i-1)).*ones(size(samp))] +0.5,'k-')
    hold on    
end
xlim([t1, t2])
ylim([0.5,(ns +0.5)])
xlabel('time', 'FontSize',18,'Interpreter','latex')
ylabel('Sample no.', 'FontSize',18,'Interpreter','latex')
set(gca,'fontsize',16)
title('Class 2 event times','FontSize',18,'Interpreter','latex')

%%% Class 3
figure
for i=1:ns
    samp = datasample(data3{i}, 200);
    plot([samp;samp],[i*ones(size(samp));(0+(i-1)).*ones(size(samp))] +0.5,'k-')
    hold on    
end
xlim([t1, t2])
ylim([0.5,(ns +0.5)])
xlabel('time', 'FontSize',18,'Interpreter','latex')
ylabel('Sample no.', 'FontSize',18,'Interpreter','latex')
set(gca,'fontsize',16)
title('Class 3 event times','FontSize',18,'Interpreter','latex')

%%% Class 4
figure
for i=1:ns
    samp = datasample(data4{i}, 200);
    plot([samp;samp],[i*ones(size(samp));(0+(i-1)).*ones(size(samp))] +0.5,'k-')
    hold on    
end
xlim([t1, t2])
ylim([0.5,(ns +0.5)])
xlabel('time', 'FontSize',18,'Interpreter','latex')
ylabel('Sample no.', 'FontSize',18,'Interpreter','latex')
set(gca,'fontsize',16)
title('Class 4 event times','FontSize',18,'Interpreter','latex')

%%% NHPP estimates
figure
plot(t,lambdaTrue1,'r','linewidth',2) %plot actual lambda
hold on
plot(t,lambdaHat1,'-.','linewidth',2); %plot estimate
fig2_leg=legend('$\lambda_{1}(t)$','$\widehat{\lambda_{1}}(t)$');
set(fig2_leg,'FontSize',18,'Interpreter','latex')
xlim([t1,t2])
ylim([0,1.4*round(max(lambdaTrue1))])
xlabel('$t$','FontSize',18,'Interpreter','latex')
ylabel('$\lambda(t)$','FontSize',18,'Interpreter','latex')
plot_title=sprintf('class 1 rate function');
set(gca,'fontsize',16)
title(plot_title,'FontSize',18,'Interpreter','latex')

figure
plot(t,lambdaTrue2,'r','linewidth',2) %plot actual lambda
hold on
plot(t,lambdaHat2,'-.','linewidth',2); %plot estimate
fig2_leg=legend('$\lambda_{2}(t)$','$\widehat{\lambda_{2}}(t)$');
set(fig2_leg,'FontSize',18,'Interpreter','latex')
xlim([t1,t2])
ylim([0,1.4*round(max(lambdaTrue2))])
xlabel('$t$','FontSize',18,'Interpreter','latex')
ylabel('$\lambda(t)$','FontSize',18,'Interpreter','latex')
plot_title=sprintf('class 2 rate function');
set(gca,'fontsize',16)
title(plot_title,'FontSize',18,'Interpreter','latex')

figure
plot(t,lambdaTrue3,'r','linewidth',2) %plot actual lambda
hold on
plot(t,lambdaHat3,'-.','linewidth',2); %plot estimate
fig2_leg=legend('$\lambda_{3}(t)$','$\widehat{\lambda_{3}}(t)$');
set(fig2_leg,'FontSize',18,'Interpreter','latex')
xlim([t1,t2])
ylim([0,1.4*round(max(lambdaTrue3))])
xlabel('$t$','FontSize',18,'Interpreter','latex')
ylabel('$\lambda(t)$','FontSize',18,'Interpreter','latex')
plot_title=sprintf('class 3 rate function');
set(gca,'fontsize',16)
title(plot_title,'FontSize',18,'Interpreter','latex')

figure
plot(t,lambdaTrue4,'r','linewidth',2) %plot actual lambda
hold on
plot(t,lambdaHat4,'-.','linewidth',2); %plot estimate
fig2_leg=legend('$\lambda_{4}(t)$','$\widehat{\lambda_{4}}(t)$');
set(fig2_leg,'FontSize',18,'Interpreter','latex')
xlim([t1,t2])
ylim([0,1.4*round(max(lambdaTrue4))])
xlabel('$t$','FontSize',18,'Interpreter','latex')
ylabel('$\lambda(t)$','FontSize',18,'Interpreter','latex')
plot_title=sprintf('class 4 rate function');
set(gca,'fontsize',16)
title(plot_title,'FontSize',18,'Interpreter','latex')

%%% Obtain posterior probabilites and classification rate for test data
labelsH = round([mp(1:ns, 1); 2*mp((1 + nts):2*ns, 2); 3*mp((1 + 2*ns):3*ns, 3); 4*mp((1 + 3*nts):4*ns, 4)]); %predictions for class labels
ac = classperf(labels, labelsH);
disp(sprintf('Piecewise linear and sinusoidal rate functions example - Classification accuracy %d %%', ac.CorrectRate*100))

toc





