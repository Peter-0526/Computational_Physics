tic
% 初始参量定义
N = 10000000;
b = 0.1;

a = 0.258.*sqrt(b);
% 选取合适的alpha使其满足条件
alpha = 1.6;
% \xi_1满足Gauss分布
sigma = 1/sqrt(2*a);
xi_1 = Random_Gauss_generator(N,sigma);
pause(2);
xi_2 = Random_generator_16807(length(xi_1));
fun1 = @(x) 1./(sqrt(2*pi).*sigma).*exp(-x.^2/(2.*sigma.^2)).*alpha;
fun2 = @(x) sqrt(2)./((1 + b.*x.^4)*pi.*b^(-1/4));
xi_2 = xi_2.*fun1(xi_1);
x0 = floor(-5.*sigma):0.01:ceil(5.*sigma);
y0_gauss = fun1(x0);
y0_lorenze = fun2(x0);

plot(x0,y0_gauss,x0,y0_lorenze,'LineWidth',2);
hold on;
scatter(xi_1(1:1000),xi_2(1:1000),2,'filled','b');
legend('Gauss','lorenze','FontSize',14);
figure;

% 假定绘制区间为[-5\sigma,5\sigma]

plot(x0,y0_gauss,x0,y0_lorenze,'LineWidth',3);

xlabel('X','FontSize',18);
ylabel('Y','FontSize',18);
hold on;

Index = find(xi_2 <= fun2(xi_1));
xi_x = xi_1(Index);
xi_y = xi_2(Index);
scatter(xi_x(1:1000),xi_y(1:1000),2,'filled','r')
histogram(xi_x,'Normalization','pdf')
legend('Gauss','lorenze','hist','FontSize',14);
title(['取样点概率密度直方图','(b = ',num2str(b),' a = ',num2str(a),')'],...
    'FontSize',18)

disp(['Sample_first : ',num2str(length(xi_1)/N)]);
disp(['Sample_second : ',num2str(length(xi_x)/length(xi_1))]);
disp(['Sample_total : ',num2str(length(xi_x)/N)]);
toc