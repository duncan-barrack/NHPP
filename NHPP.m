function [ti] = NHPP(lambdat,t0,tend)
%generates event times for a non-homogenius Poisson process (NHPP)

%%%%%%Input arguments%%%%%%%
% lambdat - rate function for the NHPP
% t0      - scaler for the left hand bound of the time interval over which the NHPP data is to be generated
% tend    - scaler for the right hand bound of the time interval over which the NHPP data is to be generated

%%%%%%Output arguments%%%%%%%
% ti    - vector containing the event times  

%%%%%%Examples%%%%%%%
%lambdat = @(x) 10*sin(x).^2; %define rate function for NHPP
%events=NHPP(lambdat,0,10);
%
%lambdat = @(x) 2*exp(x); %define rate function for NHPP
%events=NHPP(lambdat,0,10);

%check if rate function is positive
if lambdat(fminbnd(lambdat,t0,tend))<0
    error('NHPP rate function provided is not positive everywhere')
else
end    

tMax= fminbnd(@(x) -lambdat(x),t0,tend); %find the value of t which maximises the rate function lambdat
MaxLambda=lambdat(tMax); %max value of lambdat

%initialise time vector and count vector
ti=[];
t=t0;
i=0;
while t < tend
    U=rand(1,1);
    t=t-(1/MaxLambda)*(log(U));   
    lambdat_val=lambdat(t);
    u=rand(1,1);
        if u<=lambdat_val/MaxLambda
            i=i+1;
            ti(i)=t;
        else
        end    
end


ti=ti(find(ti<tend));%remove any values greater than tend



end

