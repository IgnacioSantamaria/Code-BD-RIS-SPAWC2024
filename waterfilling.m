function [p, mu] = waterfilling(b,Pt)

% Waterfilling algorithm
% Adapted from Palomar & Fonollosa TSP paper "Practical algorithms for a
% family of waterfilling problems"

% Pt = Total power (linear scale)
% b = row vector with the inverse SNRs sigma^2/lambda_i^2, i= 1,..,L

% p = optimal power vector allocated to the L channels
% mu = waterlevel


L = length(b);
a = ones(1,L);

%Reorder substreams with increasing bi
[~,index] = sort(b);
b = b(index);  

%Loop to obtain the optimal waterlevel
b(L+1) = +Inf; 
L_ = L;
while L_>=1
   mu = b(L_);
   if mu < b(L_+1)  &&  mu < ( Pt + sum(b(1:L_)) )/( sum(a(1:L_)) )
       break; %accept hypothesis
   else
      L_= L_-1; %reject hypothesis and form a new one
   end
end

%Compute the definite waterlevel and the power allocation
mu = ( Pt + sum(b(1:L_)) ) / sum(a(1:L_));
p = max(0, mu*a(1:L) - b(1:L));
      
%Recover the original ordering
p(index)=p;

% Water-level
mu=1/mu;