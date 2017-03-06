function [f, gradf] = NHPP_of_EM(a, g_at_t, g_at_tau, mp, pk, k)

nb = size(g_at_tau, 1); %no. of basis

%%% Objective function
Sumc = [];
for i = 1:k
    c = a(i*nb-nb+1:i*nb)' * g_at_tau;
    Sum = cellfun(@(x) NHPPes(a(i*nb-nb+1:i*nb), x) - c + log(pk(i)), g_at_t); 
    Sum = mp(:, i)' .* Sum;
    Sumc = [Sumc, Sum];
end
f =  - sum(Sumc);

%Gradient (
if nargout>1
    bar_c = [];
    for j=1:k        
        for i=1:size(g_at_t,2)
            bar=0;
            for T_ = 1:size(g_at_t{i},2)
                bar = bar + (g_at_t{i}(:,T_)/(a(j*nb-nb+1:j*nb)'*g_at_t{i}(:,T_)));
            end
            bar_c=[bar_c, mp(i,j)*(-bar+g_at_tau)];
            
        end
    end
    
    gradf=[];
    
    for j=1:k
        gradf=[gradf; bar_c(:,j*size(g_at_t,2)-size(g_at_t,2)+1:j*size(g_at_t,2))];
    end
    gradf = sum(gradf');
   
    
end
end

function [out] = NHPPes(a, g_at_t)
    C = num2cell(g_at_t,1);
    out = cellfun(@(x) log(a'*x), C); 
    out = sum(real(out));
end

