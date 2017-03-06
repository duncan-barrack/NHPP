function [f, gradf] = NHPP_con_EM(a, sp, t)

% Inequality constraint
f = [];
nb = sp(1).number;
k = size(sp, 2);
for i=1:k
    sp(i).coefs = a(i*nb-nb+1:i*nb)';
    f = [f, -fnval(sp(i), t)];
end
gradf=[];

%%n.b. I have included the gradient funtion for the constraint below but
%%have commented it out because it appears to slow down the optimisation. Most likley
%%because of the fnval call.
%if nargout>1
%    ll = size(t, 2);
%    gradf = zeros(nb * k, k * ll);
%    for i = 1:k
%        sp_basis(i) = spmak(sp(i).knots, eye(sp(i).number));
%        gradf(i*nb-nb+1:i*nb, i*ll-ll+1:i*ll) = fnval(sp_basis(i), t);        
%    end
%end

