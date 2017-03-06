function [mp] = NHPP_test(test, sp_fn)
%Provides log_liklehood and membership probabilities for test data 

%%%%%%Input arguments%%%%%%%
% test            - cell array of training set of all event times
% sp_fn           - NHPP B-spline functions (1 per class)

%%%%%%Output arguments%%%%%%%
%mp               - memberhsip probabilites of each class
%log_liks         - Log liklehoods for test data for each class

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Test data Checks %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %ensure test is an array
    if iscell(test)==0
        error('test argument must be a cell array!')
    else
    end

    %ensure test is not empty
    count_emp=0;
    for i=1:size(test,2)
        if isempty(test{i})==0
            count_emp=count_emp+1;
        else
        end
    
    end
    if count_emp==0
        error('test array is empty')
    else
    end

    %ensure test contains column vectors
    for i=1:size(test,2)
        if isrow(test{i})==1
        test{i}=test{i}';
        else
        end    
    end     


    m_ts = arrayfun(@(x) m_t(x), sp_fn);
    for i = 1:length(sp_fn)
        liks{i} = cellfun(@(x) lik(m_ts(i), x, sp_fn(i)), test, 'UniformOutput', false);
    end
    mpn = sym(zeros(length(test), length(sp_fn)));
    mpd = sym(zeros(length(test), 1));
    for i = 1:length(test)
        for j = 1:length(sp_fn)
            mpn(i, j) = liks{j}{i};
            mpd(i) = mpd(i) + liks{j}{i};
        end
    end
    
    mp = zeros(length(test), length(sp_fn));
    for i = 1:length(sp_fn)
        mp(:, i) = double(mpn(:, i)./mpd);
    end
    
end

function [out] = lik(C_at_t2, data, NHPP_fn)
    out = arrayfun(@(x) fnval(NHPP_fn, x), data);
    out = prod(sym(out));
    out = exp(-C_at_t2) * out;
end

function [out] = m_t(NHPP_fn)
    t1 = NHPP_fn.knots(3);
    t2 = NHPP_fn.knots(end - 2);
    t = linspace(t1, t2, NHPP_fn.number*100); %time vector
    C_at_t2 = fnval(NHPP_fn, t)';
    C_at_t2 = trapz(C_at_t2); %\int c_i(t2) dt for each basis function i
    out = (t2/length(t))*C_at_t2';
end

