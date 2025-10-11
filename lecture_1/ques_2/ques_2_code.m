tic
% 初始参数定义
a = 16807;m = 2147483647;q = 127773;r = 2836;
N = 10000000;l = 55;
I = zeros(1,N+l+2);
X = zeros(size(I));m1 = 2^32;
p1 = zeros(1,10);
p2 = zeros(1,10);
% 循环10次,每次计算平均值

for k = 1:10
    % 用当前时间构造初始值
    cl = fix(clock());
    I(1) = cl(1) + 70*(cl(2) + 12 * (cl(3) +...
        31 * (cl(4) + 23 * (cl(5) + 59 * cl(6)))));

    % 16807产生器
    for i = 2 : N+2
        I(i) = a * mod(I(i-1), q) - r * floor(I(i-1)/q);
        if I(i) < 0
            I(i) = I(i) + m;
        end
    end

    % Fibonacci 延迟产生器
    X(1:l) = I(1:l);
    for i = l+1:N+l
        X(i) = mod(X(i-24)+X(i-55),m1);
    end

    p1(k) = sum(I(l+1:N+l) > I(l+3:N+l+2) &...
        I(l+3:N+2+l) > I(l+2:N+1+l))./N;
    p2(k) = sum(X(l+1:N+l) > X(l+3:N+l+2) &...
        X(l+3:N+2+l) > X(l+2:N+1+l))./N;
    disp(['seeding = ',num2str(I(1)),' p1 = ',num2str(p1(k),'%.5f'),...
        ' p2 = ',num2str(p2(k),'%.5f')])

    pause(2)
end
disp('average:')
disp(['16807:' ,num2str(sum(p1)/10,'%.5f')])
disp(['Fibonacci: ',num2str(sum(p2)/10,'%.5f')])

toc

