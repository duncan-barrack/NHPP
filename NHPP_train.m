function [sp_fn] = NHPP_train(train, labels, t1, t2, nb, opts)
%Performs maximum likelihood classification on NHPP data

%%%%%%Input arguments%%%%%%%
% train           - cell array of training set of all event times
% labels          - vector contain class labels of training data (integer startting from 1 ,e.g. [1,1,2,3,1])
% t1              - scaler for the left hand bound of the time interval over which the NHPP data is to be generated
% t2              - scaler for the right hand bound of the time interval over which the NHPP data is to be generated
% nb              - number of cubic basis spline functions to use 
% opts            - options for the fmincon solver

%%%%%%Output arguments%%%%%%%
%sp_fn            - NHPP B-spline functions for each training class

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Training data Checks %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %ensure train is an array
    if iscell(train)==0
        error('train argument must be a cell array!')
    else
    end

    %ensure train is not empty
    count_emp=0;
    for i=1:size(train,2)
        if isempty(train{i})==0
            count_emp=count_emp+1;
        else
        end
    
    end
    if count_emp==0
        error('train array is empty')
    else
    end

    %ensure train contains column vectors
    for i=1:size(train,2)
        if isrow(train{i})==1
        train{i}=train{i}';
        else
        end    
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%get no. of training classes%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    nClass=length(unique(labels)); %number of classses
    count=0;
    classes=1:nClass;
    %check labels vector is in correct format
    for i=1:nClass
        count=count+length(find(labels==classes(i)));
    end
    if count==size(train,2);
    else
        error('Error with labels vector. Check dimension or if in correct format (integers from 1, e.g. labels=[1,1,2,3,1])')
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%get estimate of rate functions for training data%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    sp_fn = arrayfun(@(x) estimate_coeffs(train, t1, t2, nb, x, labels, opts), unique(labels));
     
   
end
 


function [sp_fn] = estimate_coeffs(train, t1, t2, nb, class, labels, opts)
disp(sprintf('Fitting training set for Class %d', class))
train_idx = labels == class;
train = train(train_idx);
knots = linspace(t1, t2, nb); % define knots
hk = knots(2) - knots(1); %get distance between knots 
knots = [knots(1)-2*hk, knots(1)-hk, knots, knots(end)+hk, knots(end)+2*hk];
sp = spmak(knots, eye(nb)); % define b-spline basis 
c_at_t = cellfun(@(x) fnval(sp, x), train, 'UniformOutput', false);

sp_basis = spmak(knots, eye(nb));
sp_con = spmak(knots, ones(1, nb));
t = linspace(t1, t2, nb*100); %time vector
tcon = linspace(t1, t2, nb*10);
C_at_t2 = fnval(sp_basis, t)';
C_at_t2 = trapz(C_at_t2); %\int c_i(t2) dt for each basis function i
C_at_t2=(t2/length(t))*C_at_t2'; %take transpose 

% Optmisiation to obtain NHPP function
disp(sprintf('Running optimisation to obtain B-spline coefficients.'))
a0 = ones(nb,1); % initial guess for B- spline coefficients. 
a = fmincon(@(a) NHPP_of(a, c_at_t, C_at_t2), a0, [], [], [], [], [], [], @(a) NHPP_con(a, sp_con, tcon), opts);
sp_fn = spmak(knots, a');


end


