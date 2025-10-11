function [x,y] = Random_Gauss_generator(N,sigma,mu)
arguments
    N (1,1) 
    sigma (1,1) = 1
    mu (1,1) = 0
end
% Random_Gauss_generator: 产生满足二维Gauss分布的随机数
% 概率密度函数:p(x,y) = 1/(2*pi) * exp(-(x.^2+y.^2)/(2*sigma))
% 也可以输出一维Gauss分布随机数
% p(x) = 1/sqrt(2*pi*sigma^2) * exp(-(x-mu)^2/(2*sigma))

x0 = Random_generator_16807(2*N);
u_1 = x0(1:N).*2-1;
v_1 = x0(N+1:2*N).*2-1;
r = sqrt(u_1.^2+v_1.^2);
Index = r < 1;
x = u_1(Index)./r(Index).*sqrt(-2.*log(r(Index).^2)).*sigma + mu;
y = v_1(Index)./r(Index).*sqrt(-2.*log(r(Index).^2)).*sigma;
end