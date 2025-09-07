function [gbest,gbestval,convergence,fitcount]= AMTSA_func(fhd,Dimension,Particle_Number,Max_Gen,VRmin,VRmax,varargin)
%[gbest,gbestval,fitcount]= PSO_func('f8',3500,200000,30,30,-5.12,5.12)

N=Particle_Number; % initial number of trees in an area.
low=ceil(0.1*N);          % is used by determining the number of seeds produced for a tree
high=ceil(0.25*N);        % is used by determining the number of seeds produced for a tree
D= Dimension; % Dimensionality of the problem
ST=0.1;           % Search tendency parameter used for equation selection 搜索趋势参数
% ST_initial = 0.9;  % 初始值，高探索
% ST_final = 0.1;    % 结束值，高开发
% ST_initial = 0.9;  % 初始值，高探索
% ST_final = 0.1;    % 结束值，高开发

dmin=VRmin;        % The lower bound of search space
dmax=VRmax;         % The upper bound of search space
fes=N;              % is the fes counter.
% max_fes=500000;    % is the termination condition, 与PSO的函数评估次数相等
trees=zeros(N,D);    % tree population
obj=zeros(1,N);      % objective function value for each tree
t=1;
convergence = ones(1,Max_Gen);

% 设置 Weibull 分布的参数
%     lambda = (dmax - dmin) / 2; % 尺度参数
%     k = 1 + sqrt(D / 10); % 形状参数

% prime_number_min = D * 2 + 3; %%
% while 1
%     if isprime(prime_number_min) == 1
%         break;
%     else
%         prime_number_min = prime_number_min + 1;
%     end
% end
% for i=1:N
%     for j=1:D
%        r = mod(2 * cos(2 * pi * j / prime_number_min) * i, 1); % 对应维度的 r
%        trees(i, j) = dmin + r * (dmax - dmin); % 计算树的位置
%     end
% end
% beta=1.0;
% sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);


for i=1:N
    for j=1:D
        trees(i,j)=dmin+(dmax-dmin)*rand;   % The trees location on the search space
%           trees(i, j) = dmin + wblrnd(lambda, k);
    end
end
obj=feval(fhd,trees',varargin{:});

%%% The determination of best tree location in the stand
%必须注意到，在取obj时，该算法与STASA不一样
[values,indis]=min(obj);
indis=indis(end);
best=obj(indis);
best_params=trees(indis,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  while  (t <= Max_Gen) %(fes<max_fes)  % Max_Gen 最大迭代次数未使用
% fitness_diff = max(obj) - min(obj);
% ST_iter = ST_initial - (ST_initial - ST_final) * (t / Max_Gen);
% ST_fitness = 0.1 + 0.8 * (1 - exp(-fitness_diff / mean(obj)));
% ST = (ST_iter + ST_fitness) / 2;

%%%%%%%%%%%%%%%%%%   SEED PRODUCTION   %%%%%%%%%%%%%%%%%%
      beta = 0.2+(1.2-0.2)*(1-(t/Max_Gen)^3)^2;                        % Eq.(14.2)
      alpha = abs(beta.*sin((3*pi/2+sin(3*pi/2*beta))));  
      mean_obj = mean(obj);
      std_obj = std(obj);
      lambda = (dmax - dmin) * (1 + std_obj / mean_obj);  % Adjust lambda based on fitness variance
      k = 1 + sqrt(D / 10) * (1 + (mean_obj - best) / mean_obj); % Adjust k based on fitness improvement
      
%       step_size = alpha;
 
%一种
%      beta = 0.2+(1.2-0.2)*(1-(t/Max_Gen)^3)^2;                        % Eq.(14.2)
%      alpha = abs(beta.*sin((3*pi/2+sin(3*pi/2*beta))));              % Eq.(14.1)

%二种
% 初始步长 alpha 设为 1.0
% alpha = 1.0;
% alpha = 1.0 + (1.0 - 1 / (1 + exp(-Dimension / 10))); 


% 调整参数 beta 的取值
% beta = 0.2 + (1.2 - 0.2) * (1 - (t / Max_Gen)^3)^2;

% % 调整初始步长 alpha
% initial_alpha = 1.0; % 初始步长
% % 调整参数 beta
% initial_beta = 1.0; % 初始 beta 值
% beta = initial_beta;
% % 动态调整步长参数 alpha
% alpha = abs(beta * sin((3 * pi / 2 + sin(3 * pi / 2 * beta)))); % 根据公式(14.1)计算 alpha
% alpha = max(alpha, initial_alpha); % 初始步长应足够大，有助于跳出局部最优解


% 改进后的 alpha 和 beta 计算
% 使用指数函数动态调整 alpha 和 beta
% z=5 * ( Dimension / 50);
% exp_factor = 5; % 指数函数的系数，控制调整速率
% beta = 0.2 + (1.2 - 0.2) * (1 - (t / Max_Gen)^3)^2;
% alpha = abs(beta .* sin((3 * pi / 2 + sin(3 * pi / 2 * beta)))) * exp(-exp_factor * t / Max_Gen);


% % 引入自适应因子，根据搜索进展动态调整参数变化速率
% adapt_factor = 1 + exp(-exp_factor * t / Max_Gen); % 自适应因子，控制参数变化速率
% beta = beta * adapt_factor;
% alpha = alpha * adapt_factor;


%      ro = alpha.*(2*rand-1);
%         mean_obj = mean(obj);
%         std_obj = std(obj);
%         lambda = lambda * (1 + std_obj / mean_obj);  % Adjust lambda based on fitness variance
%         k = k * (1 + (mean_obj - best) / mean_obj); % Adjust k based on fitness improvement

        for i=1:N
            ns = ceil(wblrnd(lambda, k)); %weibull分布
            ns=ns(1);
            seeds=zeros(ns,D);
           
            for j=1:ns %生成种子
                % Determination of the neighbour for ith tree.    
                   r=fix(rand*N)+1;
                   while(i==r)
                        r=fix(rand*N)+1;
                   end
                %%% SEEDS ARE CREATING
                for d=1:D
                    

%                      u=rand*sigma;
%                      v=rand;
%                      step=u./abs(v).^(1/beta);
%                      stepsize=step*((Max_Gen-t)/Max_Gen);
                    if(rand<ST)
                        % Calculate step2              
                    %这一串中的seeds(j,d)需要修改，特别重要
%                      step1 = 0.5 * sign(rand(1,D) - 0.5) .* norm(best_params -trees(r,d));
%                      u = rand; 
%                      wbld = B * (-log(1 - u))^(1/A); 
%                      seeds(j,d) = trees(i,d) + step1* wbld ;


%                        step_length = alpha / (t^beta);
% 改进后的 step_length 计算
% 动态调整 step_length，并引入自适应因子和随机性
% step_exp_factor = 3; % 指数函数的系数，控制调整速率
% adapt_factor = 1 + exp(-step_exp_factor * t / Max_Gen); % 自适应因子，控制参数变化速率

% % 计算动态调整后的 step_length
%                        step_length = alpha / (t^beta)+ randn * 0.1; % 添加随机扰动项，增加随机性
%                        seeds(j,d) = trees(i,d) + step_length * (best_params(d) - trees(r,d));


%                        step_length = norm(best_params(d)-trees(r,d)) * (rand * 0.5 + 0.5); % 使用一定范围内的随机步长
%                         seeds(j,d)=trees(i,d)+(best_params(d)-trees(r,d))*ro; 
%                        exponential_value = exprnd(1); % 生成指数分布的随机数
%                          seeds(j, d) = trees(i, d) + step_length * (2 * rand - 1) * exponential_value;
%                         seeds(j,d)=trees(i,d)+(best_params(d)-trees(r,d))*(rand-0.5)*2; 


%                        adaptive_step_length = step_length / (1 + 0.1 * abs(trees(i, d) - trees(r, d))); % 根据当前位置动态调整步长
%                        seeds(j, d) = trees(i, d) + adaptive_step_length * (2 * rand - 1)*(best_params(d)-trees(r,d));
                       
                        

                        step_length = alpha / (t^beta);
                        seeds(j,d) = trees(i,d) + step_length * (best_params(d) - trees(r,d)) ;
                       
                       
                       
%                       seeds(j,d)=trees(i,d)+(best_params(d)-trees(r,d))*step_size;     
%                          competition_strength = rand * k; % 根据 Weibull 分布的形状参数调整竞争强度
%                         seeds(j, d) = trees(i, d) + competition_strength * (best_params(d) - trees(r, d)) * (rand - 0.5) * 2;



                        if(seeds(j,d)>dmax || seeds(j,d)<dmin)
                            seeds(j,d)=dmin+(dmax-dmin)*rand;
                        end
                    else                   
%                         step2 = 0.1 * sign(rand(1,D) - 0.5) ;
%                         u = rand; 
%                         wbld = B * (-log(1 - u))^(1/A); 
% %                         seeds(j,d) = trees(i,d) + (step2 * wbld * (trees(i,d)-trees(r,d)))/t;
%                         seeds(j,d) = trees(i,d) + step2.*wbld.*(trees(i,d)-trees(r,d));
                      



                       seeds(j,d)=trees(i,d)+(trees(i,d)-trees(r,d))*(rand-0.5)*2 ;
                      
                       
                       
%                         seeds(j, d) = trees(i, d) + (trees(i, d) - trees(r, d)) * (rand - 0.5) * 2;
                       
                       
                       if(seeds(j,d)>dmax || seeds(j,d)<dmin)
                        seeds(j,d)=dmin+(dmax-dmin)*rand;
                        end
                    end
                end
%                 objs(j)=feval(objfun,seeds(j,:));
%                 objs(i)=feval(fhd,seeds(j,:)',varargin{:});
            end
            
            objs=feval(fhd,seeds',varargin{:});
            
%              [mintohum,indis]=min(objs);
%             indis=indis(end);
%             best = obj(indis);
%             if(objs(indis)<obj(i))
%                 trees(i,:)=seeds(indis,:);
%                 obj(i)=objs(indis);
%             end
       

        [mintohum,indis]=min(objs);
        if(mintohum<obj(i))
                trees(i,:)=seeds(indis,:);
                obj(i)=mintohum;
                best_params=seeds(indis,:);
        end
     %%%%%修改
        
            fes=fes+ns;            
        end
%      [~, indis] = min(obj);
%     indis = indis(end);
%     temp_best = obj(indis);
%     if (temp_best < best)
%         best = temp_best;
%         best_params = trees(indis, :);
%     end
     
     %   改进尝试
%         [fVal,SortIndex]=sort(objs);     %min和sort也要考虑
       
        %又一新的尝试添加
%          newtrees=zeros(N,D);    % tree population
        
      
%          Num=fix(1/2*N*((Max_Gen-t)/Max_Gen))+1; 



%  m = 2 * sin(r + pi / 2) * (Max_Gen - t) / Max_Gen;

%          s = randsample(1:N, 1, true, exp(-obj)); 
%          Fitmax=obj(s);
%          ori_value = randperm(D);
%          cauchy_value = tan((ori_value - 0.5) * pi);
minobj=zeros(1,N); 
newtrees = zeros(N, D);    % 存储新的树木位置
minobj = feval(fhd, newtrees', varargin{:}); % 计算初始适应度值
% m =  2 * sin(r + pi / 2);
m = 2 * sin(r + pi / 2) ; % 根据时间动态调整迁移步长
% 计算每个个体的适应度值和迁移步长
fitness = 1 ./ (1 + obj);  % 将适应度值调整为较大的值，以便适应度较低的个体具有更大的迁移步长
% step_lengths = m * fitness; % 根据适应度值计算迁移步长
% step_lengths = 0.5 + 0.5 * fitness;
step_lengths = m * fitness; % 动态调整迁移步长，适应度较低的个体具有更大的步长
% directions = 2 * randi([0, 1], N, D) - 1;
for i = 1:N
    % 随机选择一个个体作为目标
    target_index = randi([1, N]);
%  % 根据问题特性生成柯西分布参数
%     ori_value = rand(1, D);
%     cauchy_value = tan((ori_value - 0.5) * pi);
%  % 随机生成迁移方向
%         direction = randn(1, D); % 生成服从标准正态分布的随机方向向量
% 
%         % 根据个体适应度调整方向
%         direction = direction .* fitness(i); % 适应度越低，方向的影响越小    
%     
    % 根据个体的适应度值调整迁移步长
    if obj(i) < fitness(target_index)
        step_length = step_lengths(i);
%         direction = direction * fitness(i); 
%          newtrees(i, :) = trees(i, :) + step_length * cauchy_value.* direction.* (best_params - trees(i, :));
    else
        step_length = step_lengths(target_index);
%         newtrees(i, :) = trees(i, :) + step_length * cauchy_value.* direction.* (trees(i, :) - best_params);
%         direction = direction * fitness(target_index);
    end

    % 根据问题特性生成柯西分布参数
    ori_value = rand(1, D);
    cauchy_value = tan((ori_value - 0.5) * pi);
 % 随机生成迁移方向
        direction = randn(1, D); % 生成服从标准正态分布的随机方向向量

        % 根据个体适应度调整方向
        direction = direction .* fitness(i); % 适应度越低，方向的影响越小
     newtrees(i, :) = trees(i, :) + step_length * cauchy_value.* direction;
% directions = 2 * randi([0, 1], N, D) - 1;
    % 更新树木位置
%     newtrees(i, :) = trees(i, :) + step_length * cauchy_value;
%     newtrees(i, :) = trees(i, :) + step_length * cauchy_value.* direction;
    % 边界检查
    newtrees(i, :) = max(newtrees(i, :), dmin);
    newtrees(i, :) = min(newtrees(i, :), dmax);
end

% 根据适应度值比较新位置和原位置，并更新树木位置
for i = 1:N
    if minobj(i) < obj(i)
        trees(i, :) = newtrees(i, :);
%         obj(i)=minobj(i);
         
    end
end


 %%领头者需要修改，Xpos需要修改，Select the optimal fitness value需要思考是否加入，还有num值也要考虑
 
     %%%%%%%%%%%%%%
        

        %%% The determination of best tree location obtained so far.
        [~,indis]=min(obj);
        indis=indis(end);
        temp_best=obj(indis);   
        if(temp_best<best)
            best=temp_best;
            best_params=trees(indis,:);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 记录convergence情况
        convergence(t)=best;
        t=t+1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  end  
    gbest = best_params; %返回最佳值
    gbestval =best; %返回最佳位置
    fitcount=fes; %返回函数评估次数
end