function [c,ceq] = glb_nonlin_constraint(x,N,J,mu)  
     %N = 2; J = 4; mu = 2*ones(N,1);     
%      C = zeros(N,N*(J+1)); 
%      for idx =1:N
%          C(idx,idx:N:N*J) = ones(J,1)';
%      end
%      c = C*x(1:N*J) - mu'*x(:,J+1);
    c = x(:,1:J)*ones(J,1) - x(:,J+1).*mu; 
    ceq = [];
end

