function [f, gradf] = NHPP_con(a, sp, t)

% Inequality constraint
sp.coefs = a';
f = -fnval(sp, t);
gradf=[];

%%n.b. I have included the gradient funtion for the constraint below but
%%have commented it out because it appears to slow down the optimisation. Most likley
%%because of the fnval call.
%if nargout>1
%sp_basis = spmak(sp.knots, eye(sp.number));
%gradf = fnval(sp_basis, t);
%end

