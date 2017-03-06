clear
tic %start timer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%Example 3 - Sinusoidal and linear rate functions%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %set parameters
t1 = 0; %LH bound of time interval
t2 = 4 * pi; %RH bound of time interva
nb = 100; %number of cubic basis splines
t = linspace(0, t2, 1000)'; %time vector (used for plotting NHPP rate functions)
ntr = 10; %no. of training samples for each class
nts = 10; %no. of test samples for each class
labels=[ones(1,10),2*ones(1,10),3*ones(1,10),4*ones(1,10)];% class labels
opts = optimoptions(@fmincon, 'Display', 'iter', 'GradObj', 'on','MaxFunEvals',1e7,'TolFun',1e-5);  %options for fmincon solver

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

%rate function for class 3
lambdat = @(x) ((x<t2/4)* 20 + (x<t2/2 & x>= t2/4)* 40 + (x<3*t2/4 & x>= t2/2)* 60 + (x>=3*t2/4)*80); %define rate function for NHPP
lambdaTrue3=lambdat(t);

%generate some test and training data
train3 = arrayfun(@(x) NHPP(lambdat, t1, t2), 1:ntr, 'UniformOutput', false);
test3 = arrayfun(@(x) NHPP(lambdat, t1, t2), 1:nts, 'UniformOutput', false);

%rate function for class 4
lambdat = @(x) ((x<t2/4)* 80 + (x<t2/2 & x>= t2/4)* 60 + (x<3*t2/4 & x>= t2/2)* 40 + (x>=3*t2/4)*20); %define rate function for NHPP
lambdaTrue4=lambdat(t);

%generate some data
train4 = arrayfun(@(x) NHPP(lambdat, t1, t2), 1:ntr, 'UniformOutput', false);
test4 = arrayfun(@(x) NHPP(lambdat, t1, t2), 1:nts, 'UniformOutput', false);  

% combine data
train = [train1, train2, train3, train4];
test = [test1, test2, test3, test4];

% Obtain NHPP estimates
[spn_fn] = NHPP_train(train, labels, t1, t2, nb, opts);

lambdaHat1 = fnval(spn_fn(1), t); %NHPP estimate for class 1
lambdaHat2 = fnval(spn_fn(2), t); %NHPP estimate for class 2
lambdaHat3 = fnval(spn_fn(3), t); %NHPP estimate for class 3
lambdaHat4 = fnval(spn_fn(4), t); %NHPP estimate for class 4
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

%%% Training set - Class 3
figure
for i=1:ntr
    samp = datasample(train3{i}, 200);
    plot([samp;samp],[i*ones(size(samp));(0+(i-1)).*ones(size(samp))] +0.5,'k-')
    hold on    
end
xlim([t1, t2])
ylim([0.5,(ntr +0.5)])
xlabel('time', 'FontSize',18,'Interpreter','latex')
ylabel('Sample no.', 'FontSize',18,'Interpreter','latex')
set(gca,'fontsize',16)
title('Class 3 training set event times','FontSize',18,'Interpreter','latex')

%%% Training set - Class 4
figure
for i=1:ntr
    samp = datasample(train4{i}, 200);
    plot([samp;samp],[i*ones(size(samp));(0+(i-1)).*ones(size(samp))] +0.5,'k-')
    hold on    
end
xlim([t1, t2])
ylim([0.5,(ntr +0.5)])
xlabel('time', 'FontSize',18,'Interpreter','latex')
ylabel('Sample no.', 'FontSize',18,'Interpreter','latex')
set(gca,'fontsize',16)
title('Class 4 training set event times','FontSize',18,'Interpreter','latex')

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

%%% Test set - Class 3
figure
for i=1:ntr
    samp = datasample(test3{i}, 200);
    plot([samp;samp],[i*ones(size(samp));(0+(i-1)).*ones(size(samp))] +0.5,'k-')
    hold on    
end
xlim([t1, t2])
ylim([0.5,(ntr +0.5)])
xlabel('time', 'FontSize',18,'Interpreter','latex')
ylabel('Sample no.', 'FontSize',18,'Interpreter','latex')
set(gca,'fontsize',16)
title('Class 3 test set event times','FontSize',18,'Interpreter','latex')

%%% Test set - Class 4
figure
for i=1:ntr
    samp = datasample(test4{i}, 200);
    plot([samp;samp],[i*ones(size(samp));(0+(i-1)).*ones(size(samp))] +0.5,'k-')
    hold on    
end
xlim([t1, t2])
ylim([0.5,(ntr +0.5)])
xlabel('time', 'FontSize',18,'Interpreter','latex')
ylabel('Sample no.', 'FontSize',18,'Interpreter','latex')
set(gca,'fontsize',16)
title('Class 4 test set event times','FontSize',18,'Interpreter','latex')


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
plot_title=sprintf('Class 1 rate function');
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
plot_title=sprintf('Class 2 rate function');
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
plot_title=sprintf('Class 3 rate function');
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
plot_title=sprintf('Class 4 rate function');
set(gca,'fontsize',16)
title(plot_title,'FontSize',18,'Interpreter','latex')

%%% Obtain posterior probabilites and classification rate for test data
mp = NHPP_test(test, spn_fn);
labelsH = round([mp(1:nts, 1); 2*mp((1 + nts):2*nts, 2); 3*mp((1 + 2*nts):3*nts, 3); 4*mp((1 + 3*nts):4*nts, 4)]); %predictions for class labels
ac = classperf(labels, labelsH);
disp(sprintf('Piecewise linear and sinusoidal rate functions example - Classification accuracy %d %%', ac.CorrectRate*100))


toc %stop timer


