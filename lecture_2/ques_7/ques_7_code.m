% 读取数据,得到log(N)关于能量(E)分布的概率密度函数
A = readmatrix("data_7.txt");
A(:,2) = log10(A(:,2));
Sum = sum(A(:,2));
A(:,2) = A(:,2)./Sum;

% 产生随机数xi_1,xi_2
N = 10000000;
xi_1 = Random_generator_16807(N);
pause(2);
xi_2 = Random_generator_16807(N);

% 直接抽样法
tic
C = cumsum(A(:,2));
xi_x1 = find_index(C,xi_1) + A(1,1); % CDF的反函数
toc    
% 绘图
figure;
plot(A(:,1),A(:,2),'LineWidth',3)
hold on;
histogram(xi_x1,length(A),'FaceColor','y','Normalization','pdf');
title('直接抽样直方图','FontSize',18);
legend('data','hist1','FontSize',14);
xlabel('E(eV)','FontSize',14);
ylabel('log N(proability)','FontSize',14);

% 舍选法
tic
M = max(A(:,2));
xi_y2 = M * xi_2;
Index = zeros(1,N);
for i = 1:N
    if xi_y2(i) <= A(ceil(xi_1(i).*length(A)),2)
        Index(i) = 1;
    end
end
Index1 = find(Index == 1);
xi_y2 = xi_y2(Index1);
xi_x2 = ceil(xi_1(Index1).*length(A)) + A(1,1) - 1;
toc
% 绘图
figure
plot(A(:,1),A(:,2),'LineWidth',3)
hold on
histogram(xi_x2,length(A),'FaceColor','r','Normalization','pdf');
title('舍选法抽样直方图','FontSize',18);
legend('data','hist2','FontSize',14);
xlabel('E(eV)','FontSize',14);
ylabel('log N(proability)','FontSize',14);