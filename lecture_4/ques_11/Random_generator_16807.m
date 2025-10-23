function [x] = Random_generator_16807(N)
% 用16807产生器 返回一个含有N个数的随机数序列

% 初始参数定义
a = 16807;m = 2147483647;q = 127773;r = 2836;
I = zeros(1,N);
% 用当前时间构造初始值
cl = fix(clock());
I(1) = cl(1) + 70*(cl(2) + 12 * (cl(3) +...
    31 * (cl(4) + 23 * (cl(5) + 59 * cl(6)))));

% 产生随机数
for i = 2 : N
    I(i) = a * mod(I(i-1), q) - r * floor(I(i-1)/q);
    if I(i) < 0
       I(i) = I(i) + m;
    end
end
x = I/m;
end