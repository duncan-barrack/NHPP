clear
tic %start timer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%Example 1 - Sinusoidal rate function%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %set parameters
t1 = 0; %LH bound of time interval
t2 = 4 * pi; %RH bound of time interva
nb = 100; %number of cubic basis splines
t = linspace(0, t2, 1000)'; %time vector (used for plotting NHPP rate functions)
ntr = 10; %no. of training samples for each class
nts = 10; %no. of test samples for each class
labels=[ones(1, ntr), 2*ones(1, ntr)];% class lables
opts = optimoptions(@fmincon, 'Display', 'iter', 'GradObj', 'on','MaxFunEvals',1e7,'TolFun',1e-5);  %options for fmincon solver
opts = optimoptions(@fmincon,'Display','iter','GradObj','on','MaxFunEvals',1e7,'TolFun',1e-20,'TolX',1e-20);

%rate function for class 1
lambdat = @(x) 100*sin(x/2).^2; %define rate function for NHPP
lambdaTrue1=lambdat(t); 

%generate some test and training data
train1 = arrayfun(@(x) NHPP(lambdat, t1, t2), 1:ntr, 'UniformOutput', false);
test1 = arrayfun(@(x) NHPP(lambdat, t1, t2), 1:nts, 'UniformOutput', false);

%rate function for class 2
lambdat = @(x) 100*sin(x).^2; %define rate function for NHPP
lambdaTrue2=lambdat(t);

%generate some data
train2 = arrayfun(@(x) NHPP(lambdat, t1, t2), 1:ntr, 'UniformOutput', false);
test2 = arrayfun(@(x) NHPP(lambdat, t1, t2), 1:nts, 'UniformOutput', false);  

% combine data
train = [train1, train2];
test = [test1, test2];

%%% Obtain NHPP estimates for training data
[spn_fn] = NHPP_train(train, labels, t1, t2, nb, opts);

lambdaHat1 = fnval(spn_fn(1), t); %NHPP estimate for class 1
lambdaHat2 = fnval(spn_fn(2), t); %NHPP estimate for class 2

%%% Plot Results

%%% Raster plots
%%% Training set - Class 1
figure
for i=1:ntr
    samp = datasample(train1{i}, 200);
    plot([samp;samp],[i*ones(size(samp));(0+(i-1)).*ones(size(samp))] +0.5,'k-')
    hold on    
end
xlim([t1, t2])
ylim([0.5,(ntr +0.5)])
xlabel('time', 'FontSize',18,'Interpreter','latex')
ylabel('Sample no.', 'FontSize',18,'Interpreter','latex')
set(gca,'fontsize',16)
title('Class 1 training set event times','FontSize',18,'Interpreter','latex')

%%% Training set - Class 2
figure
for i=1:ntr
    samp = datasample(train2{i}, 200);
    plot([samp;samp],[i*ones(size(samp));(0+(i-1)).*ones(size(samp))] +0.5,'k-')
    hold on    
end
xlim([t1, t2])
ylim([0.5,(ntr +0.5)])
xlabel('time', 'FontSize',18,'Interpreter','latex')
ylabel('Sample no.', 'FontSize',18,'Interpreter','latex')
set(gca,'fontsize',16)
title('Class 2 training set event times','FontSize',18,'Interpreter','latex')

%%% Test set - Class 1
figure
for i=1:ntr
    samp = datasample(test1{i}, 200);
    plot([samp;samp],[i*ones(size(samp));(0+(i-1)).*ones(size(samp))] +0.5,'k-')
    hold on    
end
xlim([t1, t2])
ylim([0.5,(ntr +0.5)])
xlabel('time', 'FontSize',18,'Interpreter','latex')
ylabel('Sample no.', 'FontSize',18,'Interpreter','latex')
set(gca,'fontsize',16)
title('Class 1 test set event times','FontSize',18,'Interpreter','latex')

%%% Test set - Class 2
figure
for i=1:ntr
    samp = datasample(test2{i}, 200);
    plot([samp;samp],[i*ones(size(samp));(0+(i-1)).*ones(size(samp))] +0.5,'k-')
    hold on    
end
xlim([t1, t2])
ylim([0.5,(ntr +0.5)])
xlabel('time', 'FontSize',18,'Interpreter','latex')
ylabel('Sample no.', 'FontSize',18,'Interpreter','latex')
set(gca,'fontsize',16)
title('Class 2 test set event times','FontSize',18,'Interpreter','latex')


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

%%% Obtain posterior probabilites and classification rate for test data
mp = NHPP_test(test, spn_fn);
labelsH = round([mp(1:nts, 1); 2*mp((1 + nts):2*nts, 2)]); %predictions for class labels
ac = classperf(labels, labelsH);
disp(sprintf('Sinusoid rate functions example - Classification accuracy %d %%', ac.CorrectRate*100))

toc %stop timer


