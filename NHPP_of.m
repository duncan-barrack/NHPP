function [f, gradf] = NHPP_of(a, g_at_t, g_at_Tau)

% Objective function
c = a' * g_at_Tau;
Sum = cellfun(@(x) NHPPes(a, x) - c, g_at_t);    
f =  - sum(Sum);

% Gradiant
gradf=[];
if nargout>1
    bar = cellfun(@(x) NHPPesGrad(a, x), g_at_t, 'UniformOutput', false);  
    bar = cell2mat(bar);
    gradf = size(g_at_t, 2) * g_at_Tau - sum(bar, 2);   
end

function [out] = NHPPes(a, g_at_t)
C = num2cell(g_at_t,1);
out = cellfun(@(x) log(a'*x), C); 
out = sum(real(out));

function [out] = NHPPesGrad(a, g_at_t)
C = num2cell(g_at_t,1);
out = cellfun(@(x) (x / (a'*x)), C, 'UniformOutput', false); 
out = cell2mat(out);
out = sum(out, 2);


