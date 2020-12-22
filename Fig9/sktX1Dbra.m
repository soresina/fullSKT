function out=sktX1Dbra(p,u)
% out
    out=[u(1);...      % u(-lx)
         u(p.np+1);... % v(-lx)
         (p.Om)^(-1/2)*sum(p.mat.M(1:p.np,1:p.np)*abs(u(1:p.np).^2))^(1/2);...% ||u||_L^2; 
         (p.Om)^(-1/2)*sum(p.mat.M(1:p.np,1:p.np)*abs(u(p.np+1:2*p.np).^2))^(1/2)];% ||v||_L^2];
end
