
clc;clear;close all;
f= @(x) - (x - 10) .^ 2 + x .* sin(x) .* cos(2 * x) - 5 * x .* sin(3 * x) ; % 适应度函数表达式(求这个函数的最大值)  
figure(1);
fplot(f, [0 20], 'b-');                 % 画出初始图像 

d = 1;                           % 空间维数(该例子是1维)  
N = 30;                          % 初始种群个数         
x_limit = [0, 20];               % 设置位置限制  
v_limit = [-1, 1];               % 设置速度限制    

x = x_limit(1) + (x_limit(2) - x_limit(1)) * rand(N, d); %初始每个粒子的位置  
v = rand(N, d);                  % 初始每个粒子的速度  

pbest = x;                       % 初始化每个个体的历史最佳位置 
gbest = zeros(1, d);             % 初始化种群的历史最佳位置  

fp_best = zeros(N, 1);           % 初始化每个个体的历史最佳适应度 为 0  
fg_best = -inf;                  % 初始化种群历史最佳适应度 为 负无穷  

iter = 50;                       % 最大迭代次数
w = 0.8;                         % 惯性权重  
c1 = 0.5;                        % 自我学习因子  
c2 = 0.2;                        % 群体学习因子 

hold on;
plot(x, f(x), 'ro');
title('初始状态图');  

  
record = zeros(iter, 1);        % 记录器(用于记录 fg_best 的变化过程)
figure(2);
i = 1; 
while i <= iter  
    fx = f(x) ;                 % 计算个体当前适应度     
    for j = 1:N        
        if fp_best(j) < fx(j) 
            fp_best(j) = fx(j); % 更新个体历史最佳适应度
            pbest(j) = x(j);    % 更新个体历史最佳位置 
        end   
    end
    if fg_best < max(fp_best) 
        [fg_best,ind_max] = max(fp_best);     % 更新群体历史最佳适应度
        gbest = pbest(ind_max);               % 更新群体历史最佳位置，其中 ind_max 是最大值所在的下标
        % 注：将上面的两个式子换成下面两个是不可以的
        % [gbest, ind_max] = max(pbest);      % 更新群体历史最佳位置，其中 ind_max 是最大值所在的下标
        % fg_best = fp_best(ind_max);         % 更新群体历史最佳适应度   
    end  
    v = v * w + c1 * rand() * (pbest - x) + c2 * rand() * (repmat(gbest, N, 1) - x); % 速度更新
    % 注： repmat(A,r1,r2):可将矩阵 扩充 为每个单位为A的r1*r2的矩阵
    
    % 边界速度处理  
    v(v > v_limit(2)) = v_limit(2);
    v(v < v_limit(1)) = v_limit(1);  
    
    x = x + v;% 位置更新  
    
    % 边界位置处理  
    x(x > x_limit(2)) = x_limit(2);  
    x(x < x_limit(1)) = x_limit(1);  
    
    record(i) = fg_best;%最大值记录  
    
    % 画动态展示图
    zuo_x = 0 : 0.01 : 20;  
    plot(zuo_x, f(zuo_x), 'b-', x, f(x), 'ro');
    title('状态位置变化')  
    pause(0.1)  
    i = i + 1;
%     if mod(i,10) == 0   % 显示进度
%         i
%     end
end  
figure(3);
plot(record);
title('收敛过程'); 
zuo_x = 0 : 0.01 : 20;  

figure(4);
plot(zuo_x, f(zuo_x), 'b-', x, f(x), 'ro');
title('最终状态图')

disp(['最佳适应度：',num2str(fg_best)]);  
disp(['最佳粒子的位置x：',num2str(gbest)]);  

