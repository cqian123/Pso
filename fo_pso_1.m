%%完成一个分数阶的粒子群算法，在速度项加入gl离散差分
clc;clear;close all;
%% 初始化种群
f= @(x)x .* sin(x) .* cos(2 * x) - 2 * x .* sin(3 * x); % 函数表达式
figure(1);
fplot(f,[0,20]);

%gl离散化粒子群的阶次alpha
alpha = 0.9; %阶次取0.9试试看


N = 50;                         % 初始种群个数
d = 1;                          % 空间维数
ger = 100;                      % 最大迭代次数     
limit = [0, 20];                % 设置位置参数限制,这里的你看f函数的横坐标0-20为限制
%？？
vlimit = [-1, 1];               % 设置速度限制(不是系数，就是速度大小)
%w = alpha;                        % 惯性权重 在分数阶这边是没用了，靠分数阶阶次来构造相应的系数
c1 = 0.5;                       % 自我学习因子
c2 = 0.5;                       % 群体学习因子 

%因为这里是单变量，故粒子群的维度只有1，循环只进行1次
for i = 1:d
    %针对每个维度，在位置限制之内线性随机分布
    x = limit(i, 1) + (limit(i, 2) - limit(i, 1)) * rand(N, d);%初始种群的位置，初始化50个粒子的位置
end
v = rand(N, ger);                  % 初始种群的速度 ,记录种群每一列代表当前每个粒子达到的速度(第一列表示N个点的初始速度，但初始化矩阵做100列的，后面的在迭代过程中替换掉，是无用数据)
xm = x;                          % 每个个体的历史最佳位置(先记录为初始位置)
ym = zeros(1, d);                % 种群的历史最佳位置，就是所有粒子在当前所有迭代过程中能够找到的最佳位置，是坐标(只有一个点，故是1*d维的)

fxm = f(x);               % 每个个体的历史最佳适应度，是函数值
fym = -inf;                      % 种群历史最佳适应度，是函数值，代表负无穷，程序是为了找max值
hold on
plot(xm, f(xm), 'ro');title('初始状态图');

%到这里为止做好了所有随机点的初始化
figure(2)
%% 群体更新
iter = 1;%记录当前迭代次数
k = 1;
%创建数组用来记录系数
coff(k) = alpha;

record = zeros(ger, 1);          % 记录器
while iter <= ger
     %分数阶的系数G-L定义下的
     coff(k+1) = - coff(k)*(alpha - k)/(k+1);
     tempcoff = fliplr(coff(1:k)); %可能不需要(将数组从左向右翻转)
     
     fx = f(x) ; % 个体当前适应度   
     for i = 1:N      
        if fxm(i) < fx(i)
            fxm(i) = fx(i);     % 更新个体历史最佳适应度
            xm(i,:) = x(i,:);   % 更新个体历史最佳位置
        end 
     end
    if fym < max(fxm)
        [fym, nmax] = max(fxm);   % 更新群体历史最佳适应度
        ym = xm(nmax, :);      % 更新群体历史最佳位置
    end
    
    temp = zeros(N,1);%%临时变量
    for i=1:k %用前面的速度信息去更新下一次的速度，要用一个循环的
    temp = temp + v(:,i)*tempcoff(i);
    end   
    v(:,k+1) = temp + c1 * rand * (xm - x) + c2 * rand * (repmat(ym, N, 1) - x);% 速度更新
    
    % 边界速度处理
    v(v > vlimit(2)) = vlimit(2);
    v(v < vlimit(1)) = vlimit(1);
    
    x = x + v(:,k+1);% 位置更新
    % 边界位置处理
    x(x > limit(2)) = limit(2);
    x(x < limit(1)) = limit(1);
    
    record(iter) = fym;%最大值记录
    
    %显示出50个粒子每一次迭代时的对应的搜索过程
    x0 = 0 : 0.01 : 20;
    plot(x0, f(x0), 'b-', x, f(x), 'ro');title('状态位置变化')
    pause(0.1)  %pause(a)表示程序暂停a秒后继续执行
    
    iter = iter+1;
    k = k + 1;
end
figure(3);plot(record);title('收敛过程')

x0 = 0 : 0.01 : 20;
figure(4);plot(x0, f(x0), 'b-', x, f(x), 'ro');title('最终状态位置')

disp(['最大值：',num2str(fym)]);
disp(['变量取值：',num2str(ym)]);