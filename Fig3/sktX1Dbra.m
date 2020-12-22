function out=sktX1Dbra(p,u)
% out 
    out=[u(1);...      % 1st component in -lx
         u(p.np+1);... % 2nd component in -lx
         (p.Om)^(-1/2)*sum(p.mat.M(1:p.np,1:p.np)*abs(u(1:p.np).^2))^(1/2);...% ||u_1||_L^2; 
         (p.Om)^(-1/2)*sum(p.mat.M(1:p.np,1:p.np)*abs(u(p.np+1:2*p.np).^2))^(1/2)];% ||u_1||_L^2]; 
end
