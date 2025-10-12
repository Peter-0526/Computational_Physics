tic
N = 100000;
X = Random_generator_16807(50*N);
gauss = @(x) exp(-x.^2/2)/(sqrt(2*pi));
x_0 = -4:0.02:4;

Y1 = zeros(size(X(1:N)));   % 均匀分布
Y2 = zeros(size(X(1:N)));   % 指数分布
Y3 = zeros(size(X(1:N)));   % 泊松分布
C3 = zeros(size(1,20));
Y4 = zeros(size(X(1:N)));   % 二项分布
C4 = zeros(size(1:21));

% 生成[0,1]区间的均匀分布
t1 = tiledlayout(4,4);

for i = 1:50
    Y1 = Y1 + X((i-1)*N+1:i*N);
    if i == 2 || i == 5 || i == 10 || i == 50
        mu = sum(Y1/i)./(N);
        sigma = sqrt((sum((Y1/i).^2)/N-mu^2));
        nexttile
        histogram((Y1/i-mu)/(sigma),100,'Normalization','pdf');hold on;
        plot(x_0,gauss(x_0),'LineWidth',2);
        title(['均匀分布','N = ',num2str(i)])
        legend('Experiment','Gauss')
    end
end

% 指数分布X ~ Exp(\lambda) ,\lamaba = 2
fun2 = @(x,lamaba) -log(1-x)./lamaba;
lamaba1 = 2;
for i = 1:50
    Y2 = Y2 + fun2(X((i-1)*N+1:i*N),lamaba1);
    if i == 2 || i == 5 || i == 10  || i == 50
        mu = sum(Y2/i)./(N);
        sigma = sqrt((sum((Y2/i).^2)/N-mu^2));
        nexttile
        histogram((Y2/i-mu)/(sigma),100,'Normalization','pdf');hold on;
        plot(x_0,gauss(x_0),'LineWidth',2);
        title(['指数分布','N = ',num2str(i)])
        legend('Experiment','Gauss')        
    end
end

% 泊松分布 X ~ P(\lamaba) lamaba = 1
fun3 = @(lamaba,k) exp(-lamaba).*lamaba.^k./factorial(k);
lamaba2 = 4.5;
% 累计分布函数
C3(1) = fun3(lamaba2,0);
for i = 2:20
    C3(i) = C3(i-1)+fun3(lamaba2,i-1);
end
Indexs = find_index(C3,X);
for i = 1:50
    Y3 = Y3 + Indexs((i-1)*N+1:i*N);
    if i == 2 || i == 5 || i==10 || i == 50
        mu = sum(Y3/i)./(N);
        sigma = sqrt((sum((Y3/i).^2)/N-mu^2));
        nexttile
        histogram((Y3/i-mu)/(sigma),50,'Normalization','pdf');hold on;
        plot(x_0,gauss(x_0),'LineWidth',2);
        title(['泊松分布','N = ',num2str(i)])
        legend('Experiment','Gauss')        
    end
end

% 二项分布 B(20,0.3)
fun4 = @(n,p,k) factorial(n)./(factorial(n-k).*factorial(k)).*p.^k.*(1-p).^(n-k);
n = 20;
p = 0.3;
% 累计分布函数
C4(1) = fun4(n,p,0);
for i = 2:21
    C4(i) = C4(i-1)+fun4(n,p,i-1);
end
Indexs_2 = find_index(C4,X);
for i = 1:50
    Y4 = Y4 + Indexs_2((i-1)*N+1:i*N);
    if i == 2 || i == 5 || i == 10 || i == 50
        mu = sum(Y4/i)./(N);
        sigma = sqrt((sum((Y4/i).^2)/N-mu^2));
        nexttile
        histogram((Y4/i-mu)/(sigma),50,'Normalization','pdf');hold on;
        plot(x_0,gauss(x_0),'LineWidth',2);
        title(['二项分布','N = ',num2str(i)])
        legend('Experiment','Gauss')        
    end
end

toc