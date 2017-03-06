function [sp_fn, mp, q] = NHPP_cluster(data, k, t1, t2, varargin)
%Clusters NHPP data into k classes using the EM algorithm

%%%%%%%%%%%%Mandatory Input arguments%%%%%%% 
% data          - cell array containing data
% k             - number of clusters
% t1            - scaler for the left hand bound of the time interval over which the NHPP data is to be generated
% t2            - scaler for the right hand bound of the time interval over which the NHPP data is to be generated

%%%%%%%%%%%%Optional Input arguments (specified in varargin)%%%%%%% 
% fmincon_opts  - options for the fmincon solver (default is none)
% nb            - number of cubic bspline basis functions to use (default 50)
% a0            - initial guesses for coefficients of the basis functions for the k classes (default is ones(1,k*nb).*rand(1,k*nb);)
% pk0           - initial guess for mixing probabilites (default is 1/k)
% MaxIter       - maximum number of iterations for the EM algorimth (default is 12)
% AbsTol        - Tolerance for the absolute difference in auxillary function for termination of alogorithm (default is 0.0001)  

%%%%%%Output arguments%%%%%%%
% sp_fn         - NHPP B-spline functions for each k
% mp            - posterior probabilities for eaxh sampe
% q             - value of auxiliary function
  

    %check if data argument is a cell array
    if iscell(data)==0
        error('data argument must be a cell array')
    end
    
    %check if any samples contain no events
    for i = 1:size(data, 2)
        if (length(data{i})) == 0
            error('Sample %d in data cell array contains no events!', i)
        end
    end
    
    % Check ther are at most 6 optional inputs arguments
    numvarargs = length(varargin);
    if numvarargs ~= 6
         error('NHPP_cluster requires 10 inputs. Optional arguments must be empty if you wish to use default settings, e.g. nb =[]');
    end
    
    % set defaults for optional inputs
    optargs{1} = optimoptions(@fmincon,'HessFcn', []); %n.b. HessFcn is not set by default
    optargs{2} = 50;
    optargs{4} = ones(1,k)/k;
    optargs{5} = 12;
    optargs{6} = 0.0001;
    
    %Obtain initial guesses of B-spline coefficients if not provided by
    %user
    if length(varargin{2}) > 0 && length(varargin{3}) ==0
        varargin{3} = ones(optargs{2},k).*rand(optargs{2},k);
    else
    end
    
    %overwrite optional arguments with user specified arguments
    for i=1:numvarargs
        if length(varargin{i}) > 0
            optargs{i} = varargin{i};
        end
    end
              
    % Place optional args in variable names
    [opts, nb, a0, pk0, MaxIter, AbsTol] = optargs{:};
    
    % Error checks for user supplied arguments    
    if size(a0, 2) ~= k || size(a0, 1) ~= nb
        error('Coefficients matrix a0 in wrong format or no. of coefficients provided inconsistent values for k and nb\n k is %d \n nb is %d \n a0 must be %dx%d', k, nb, nb, k)
    end
    
    if size(pk0, 2) ~= k
        error('Number of values in pk0 must equal k')
     end
    
    if sum(pk0) ~= 1
        error('Values in pk0 must sum to 1')
    end
    
    
    %%% initialisation
    
    % set B-spline knots
    knots = linspace(t1, t2, nb); 
    hk = knots(2) - knots(1); %get distance between knots 
    knots = [knots(1)-2*hk, knots(1)-hk, knots, knots(end)+hk, knots(end)+2*hk];
    
    for i=1:size(a0, 2)
        sp_fn(i) = spmak(knots, a0(:, i)'); %spline function    
        sp_con(i) = spmak(knots, ones(1, nb)); %spline function for constraint
    end
    
    sp = spmak(knots, eye(nb)); % define b-spline basis 
    c_at_t = cellfun(@(x) fnval(sp, x), data, 'UniformOutput', false); 
    sp_basis = spmak(knots, eye(nb));
    t = linspace(t1, t2, nb*100); %time vector
    tcon = linspace(t1, t2, nb*10);
    C_at_t2 = fnval(sp_basis, t)';
    C_at_t2 = trapz(C_at_t2); %\int c_i(t2) dt for each basis function i
    C_at_t2=(t2/length(t))*C_at_t2'; %take transpose 
    
    disp(sprintf('Running EM algorithm with %d B-spline functions for each NHPP in mixture model of size %d', nb, k))
    for iter=1:MaxIter
        disp(sprintf('EM algorithm Iteration %d ......', iter))
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%Expectation Step%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %set initial conditions
        if iter==1
            pk=pk0;
            a=a0;            
        else
        end
        mp = mps(data, sp_fn, pk); %membership probabilities
  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%Maximisation Step%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
        %flatten array containing coefficient
        a=reshape(a,numel(a),1);
        disp(sprintf('Running optimisation to obtain B-spline estimates.')) 
        
        % last function auxilliary function value
        if iter==1
            ql = NHPP_of_EM(a, c_at_t, C_at_t2, mp, pk, k);
        else
            ql = q;
        end
        
        % get B-spline coefficients
        [a, q] = fmincon(@(a) NHPP_of_EM(a, c_at_t, C_at_t2, mp, pk, k), a, [], [], [], [], [], [], @(a) NHPP_con_EM(a, sp_con, tcon), opts);
        
        %  get mixing probabilities 
        pk = mix(mp, k); 
        
        % reshape B-spline coefficient vector
        a= reshape(a,length(a)/k,k);             
    
        %update spline functions with new coefficient values
        for i=1:size(a, 2)
            sp_fn(i) = spmak(knots, a(:, i)'); %spline function        
        end
    
        %check stopping criteria       
        if abs(-(q - ql)) < AbsTol
            disp('EM stopping criteria satisfied')
            mp = mps(data, sp_fn, pk); % final membership probabilities
            q = -q;
            break;
        end
    
    end
    disp('Maximum number of EM iterations reached')
    q = -q;
        
end

function [out] = mps(data, sp_fn, pk)
 %%% Membership probabilities
    
     m_ts = arrayfun(@(x) m_t(x), sp_fn);
    for i = 1:length(sp_fn)
        liks{i} = cellfun(@(x) lik(m_ts(i), x, sp_fn(i)), data, 'UniformOutput', false);
    end
    mpn = sym(zeros(length(data), length(sp_fn)));
    mpd = sym(zeros(length(data), 1));
    for i = 1:length(data)
        for j = 1:length(sp_fn)
            mpn(i, j) = pk(j)*liks{j}{i};
            mpd(i) = mpd(i) + pk(j)*liks{j}{i};
        end
    end
    
    mp = zeros(length(data), length(sp_fn));
    for i = 1:length(sp_fn)
        mp(:, i) = double(mpn(:, i)./mpd);
    end
    out = mp;
end

function [out] = lik(C_at_t2, data, NHPP_fn)
    out = arrayfun(@(x) fnval(NHPP_fn, x), data);
    out = prod(sym(out));
    out = exp(sym(-C_at_t2)) * out;
end

function [out] = m_t(NHPP_fn)
    t1 = NHPP_fn.knots(3);
    t2 = NHPP_fn.knots(end - 2);
    t = linspace(t1, t2, NHPP_fn.number*100); %time vector
    C_at_t2 = fnval(NHPP_fn, t)';
    C_at_t2 = trapz(C_at_t2); %\int c_i(t2) dt for each basis function i
    out = (t2/length(t))*C_at_t2';
end

function [out] = mix(mp, k)
    for i=1:k
            out(i)=0;
            for j=1:size(mp, 1)
                out(i)=out(i)+mp(j, i);
            end
            out(i)=out(i)/size(mp, 1);        
    end
end

