tic
% 初始参数定义
a = 16807;m = 2147483647;q = 127773;r = 2836;
N = 20000000;l = 10;
I = zeros(1,N+l);
% 用当前时间构造初始值
cl = fix(clock());
I(1) = cl(1) + 70*(cl(2) + 12 * (cl(3) +...
    31 * (cl(4) + 23 * (cl(5) + 59 * cl(6)))));

% 产生随机数
for i = 2 : N+10
    I(i) = a * mod(I(i-1), q) - r * floor(I(i-1)/q);
    if I(i) < 0
       I(i) = I(i) + m;
    end
end
x = I/m;

% <x^{k}>均匀性检测
N_1 = [100,1000,10000,100000,1000000,10000000];
xk = zeros(size(N_1,2),5);
err = zeros(size(xk));

for i = 1:size(N_1,2)
    disp(['N_1 = ',num2str(N_1(i))])
    for k = 1:5
        xk(i,k) = sum(x(1:N_1(i)).^k)/N_1(i);
        err(i,k) = abs(xk(i,k) - 1/(1+k));
        disp(['<x^',num2str(k),'> = ',num2str(xk(i,k)),' expectation = ',num2str(1/(1+k))])
    end
end

% 以k = 5为例,拟合err 与 lg(N)
figure;
hold on;
lg_n = log10(N_1);
lg_err = log10(err(:,5));
plot(lg_n,lg_err,'-o','LineWidth',2,'Color','b')
xlabel('lg(N)','FontSize',20);
ylabel('lg(error)','FontSize',20);
p = polyfit(lg_n,lg_err,1);
f1 = polyval(p,lg_n);
plot(lg_n,f1,'LineWidth',2,'Color','r')
legend('lg(error)','拟合')
text(lg_n(3),f1(3),['\leftarrow y = ',num2str(p(1)) , ' x + ',num2str(p(2))],...
    'Fontsize',16);
hold off;

% C(l)测量二维独立性
% C(l) = (<x_{n+l}*x_{n}> - <x_n>^2)/(<x_n^{2}> - <x_n>^2)
cl = zeros(1,9);
for l1 = 2:l
    xnl_xn = x(1:N)*x(1+l1:N+l1)'/N;
    x_n = sum(x(1:N))/N;
    x_n_2 = x(1:N)*x(1:N)'/N;
    cl(l1-1) = (xnl_xn - x_n^2)/(x_n_2 - x_n^2);
    disp(['C(',num2str(l1),') = ',num2str(cl(l1-1))])
end

% 卡方检测
for K = 2:10
    m_k = N/K;mu = K-1;% 自由度
    x_k = floor(x*K);
    [counts,centers] = hist(x_k,0:K-1);
    chi_2 = sum((counts - m_k).^2/m_k);
    fun = @(t) t.^(mu/2-1).*exp(-t/2);
    P = integral(fun,0,chi_2)*2.^(-mu/2)./gamma(mu/2);
    disp(['K = ',num2str(K),'     chi_2 = ',num2str(chi_2),...
        '     P = ',num2str(P)])
end

% 采点，绘制前200个点
% figure;
% hold on;
% N1 = 5000;
% l = 10;
% scatter(x(1:N1),x(1+l:N1+l),'filled','b')
% xlabel('X','FontSize',20);
% ylabel('Y','FontSize',20);
% axis equal;
% axis([-0.1,1.1,-0.1,1.1])
% title(['Plane distribution of random numbers(N = ',num2str(N1),...
%     ' ,interval = ',num2str(l),')' ],'FontSize',20)
toc;