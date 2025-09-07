function [gbest,gbestval,convergence,fitcount]= AMTSA_func(fhd,Dimension,Particle_Number,Max_Gen,VRmin,VRmax,varargin)
%[gbest,gbestval,fitcount]= PSO_func('f8',3500,200000,30,30,-5.12,5.12)

N=Particle_Number; % initial number of trees in an area.
low=ceil(0.1*N);          % is used by determining the number of seeds produced for a tree
high=ceil(0.25*N);        % is used by determining the number of seeds produced for a tree
D= Dimension; % Dimensionality of the problem
ST=0.1;           % Search tendency parameter used for equation selection �������Ʋ���
% ST_initial = 0.9;  % ��ʼֵ����̽��
% ST_final = 0.1;    % ����ֵ���߿���
% ST_initial = 0.9;  % ��ʼֵ����̽��
% ST_final = 0.1;    % ����ֵ���߿���

dmin=VRmin;        % The lower bound of search space
dmax=VRmax;         % The upper bound of search space
fes=N;              % is the fes counter.
% max_fes=500000;    % is the termination condition, ��PSO�ĺ��������������
trees=zeros(N,D);    % tree population
obj=zeros(1,N);      % objective function value for each tree
t=1;
convergence = ones(1,Max_Gen);

% ���� Weibull �ֲ��Ĳ���
%     lambda = (dmax - dmin) / 2; % �߶Ȳ���
%     k = 1 + sqrt(D / 10); % ��״����

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
%        r = mod(2 * cos(2 * pi * j / prime_number_min) * i, 1); % ��Ӧά�ȵ� r
%        trees(i, j) = dmin + r * (dmax - dmin); % ��������λ��
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
%����ע�⵽����ȡobjʱ�����㷨��STASA��һ��
[values,indis]=min(obj);
indis=indis(end);
best=obj(indis);
best_params=trees(indis,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  while  (t <= Max_Gen) %(fes<max_fes)  % Max_Gen ����������δʹ��
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
 
%һ��
%      beta = 0.2+(1.2-0.2)*(1-(t/Max_Gen)^3)^2;                        % Eq.(14.2)
%      alpha = abs(beta.*sin((3*pi/2+sin(3*pi/2*beta))));              % Eq.(14.1)

%����
% ��ʼ���� alpha ��Ϊ 1.0
% alpha = 1.0;
% alpha = 1.0 + (1.0 - 1 / (1 + exp(-Dimension / 10))); 


% �������� beta ��ȡֵ
% beta = 0.2 + (1.2 - 0.2) * (1 - (t / Max_Gen)^3)^2;

% % ������ʼ���� alpha
% initial_alpha = 1.0; % ��ʼ����
% % �������� beta
% initial_beta = 1.0; % ��ʼ beta ֵ
% beta = initial_beta;
% % ��̬������������ alpha
% alpha = abs(beta * sin((3 * pi / 2 + sin(3 * pi / 2 * beta)))); % ���ݹ�ʽ(14.1)���� alpha
% alpha = max(alpha, initial_alpha); % ��ʼ����Ӧ�㹻�������������ֲ����Ž�


% �Ľ���� alpha �� beta ����
% ʹ��ָ��������̬���� alpha �� beta
% z=5 * ( Dimension / 50);
% exp_factor = 5; % ָ��������ϵ�������Ƶ�������
% beta = 0.2 + (1.2 - 0.2) * (1 - (t / Max_Gen)^3)^2;
% alpha = abs(beta .* sin((3 * pi / 2 + sin(3 * pi / 2 * beta)))) * exp(-exp_factor * t / Max_Gen);


% % ��������Ӧ���ӣ�����������չ��̬���������仯����
% adapt_factor = 1 + exp(-exp_factor * t / Max_Gen); % ����Ӧ���ӣ����Ʋ����仯����
% beta = beta * adapt_factor;
% alpha = alpha * adapt_factor;


%      ro = alpha.*(2*rand-1);
%         mean_obj = mean(obj);
%         std_obj = std(obj);
%         lambda = lambda * (1 + std_obj / mean_obj);  % Adjust lambda based on fitness variance
%         k = k * (1 + (mean_obj - best) / mean_obj); % Adjust k based on fitness improvement

        for i=1:N
            ns = ceil(wblrnd(lambda, k)); %weibull�ֲ�
            ns=ns(1);
            seeds=zeros(ns,D);
           
            for j=1:ns %��������
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
                    %��һ���е�seeds(j,d)��Ҫ�޸ģ��ر���Ҫ
%                      step1 = 0.5 * sign(rand(1,D) - 0.5) .* norm(best_params -trees(r,d));
%                      u = rand; 
%                      wbld = B * (-log(1 - u))^(1/A); 
%                      seeds(j,d) = trees(i,d) + step1* wbld ;


%                        step_length = alpha / (t^beta);
% �Ľ���� step_length ����
% ��̬���� step_length������������Ӧ���Ӻ������
% step_exp_factor = 3; % ָ��������ϵ�������Ƶ�������
% adapt_factor = 1 + exp(-step_exp_factor * t / Max_Gen); % ����Ӧ���ӣ����Ʋ����仯����

% % ���㶯̬������� step_length
%                        step_length = alpha / (t^beta)+ randn * 0.1; % �������Ŷ�����������
%                        seeds(j,d) = trees(i,d) + step_length * (best_params(d) - trees(r,d));


%                        step_length = norm(best_params(d)-trees(r,d)) * (rand * 0.5 + 0.5); % ʹ��һ����Χ�ڵ��������
%                         seeds(j,d)=trees(i,d)+(best_params(d)-trees(r,d))*ro; 
%                        exponential_value = exprnd(1); % ����ָ���ֲ��������
%                          seeds(j, d) = trees(i, d) + step_length * (2 * rand - 1) * exponential_value;
%                         seeds(j,d)=trees(i,d)+(best_params(d)-trees(r,d))*(rand-0.5)*2; 


%                        adaptive_step_length = step_length / (1 + 0.1 * abs(trees(i, d) - trees(r, d))); % ���ݵ�ǰλ�ö�̬��������
%                        seeds(j, d) = trees(i, d) + adaptive_step_length * (2 * rand - 1)*(best_params(d)-trees(r,d));
                       
                        

                        step_length = alpha / (t^beta);
                        seeds(j,d) = trees(i,d) + step_length * (best_params(d) - trees(r,d)) ;
                       
                       
                       
%                       seeds(j,d)=trees(i,d)+(best_params(d)-trees(r,d))*step_size;     
%                          competition_strength = rand * k; % ���� Weibull �ֲ�����״������������ǿ��
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
     %%%%%�޸�
        
            fes=fes+ns;            
        end
%      [~, indis] = min(obj);
%     indis = indis(end);
%     temp_best = obj(indis);
%     if (temp_best < best)
%         best = temp_best;
%         best_params = trees(indis, :);
%     end
     
     %   �Ľ�����
%         [fVal,SortIndex]=sort(objs);     %min��sortҲҪ����
       
        %��һ�µĳ������
%          newtrees=zeros(N,D);    % tree population
        
      
%          Num=fix(1/2*N*((Max_Gen-t)/Max_Gen))+1; 



%  m = 2 * sin(r + pi / 2) * (Max_Gen - t) / Max_Gen;

%          s = randsample(1:N, 1, true, exp(-obj)); 
%          Fitmax=obj(s);
%          ori_value = randperm(D);
%          cauchy_value = tan((ori_value - 0.5) * pi);
minobj=zeros(1,N); 
newtrees = zeros(N, D);    % �洢�µ���ľλ��
minobj = feval(fhd, newtrees', varargin{:}); % �����ʼ��Ӧ��ֵ
% m =  2 * sin(r + pi / 2);
m = 2 * sin(r + pi / 2) ; % ����ʱ�䶯̬����Ǩ�Ʋ���
% ����ÿ���������Ӧ��ֵ��Ǩ�Ʋ���
fitness = 1 ./ (1 + obj);  % ����Ӧ��ֵ����Ϊ�ϴ��ֵ���Ա���Ӧ�Ƚϵ͵ĸ�����и����Ǩ�Ʋ���
% step_lengths = m * fitness; % ������Ӧ��ֵ����Ǩ�Ʋ���
% step_lengths = 0.5 + 0.5 * fitness;
step_lengths = m * fitness; % ��̬����Ǩ�Ʋ�������Ӧ�Ƚϵ͵ĸ�����и���Ĳ���
% directions = 2 * randi([0, 1], N, D) - 1;
for i = 1:N
    % ���ѡ��һ��������ΪĿ��
    target_index = randi([1, N]);
%  % ���������������ɿ����ֲ�����
%     ori_value = rand(1, D);
%     cauchy_value = tan((ori_value - 0.5) * pi);
%  % �������Ǩ�Ʒ���
%         direction = randn(1, D); % ���ɷ��ӱ�׼��̬�ֲ��������������
% 
%         % ���ݸ�����Ӧ�ȵ�������
%         direction = direction .* fitness(i); % ��Ӧ��Խ�ͣ������Ӱ��ԽС    
%     
    % ���ݸ������Ӧ��ֵ����Ǩ�Ʋ���
    if obj(i) < fitness(target_index)
        step_length = step_lengths(i);
%         direction = direction * fitness(i); 
%          newtrees(i, :) = trees(i, :) + step_length * cauchy_value.* direction.* (best_params - trees(i, :));
    else
        step_length = step_lengths(target_index);
%         newtrees(i, :) = trees(i, :) + step_length * cauchy_value.* direction.* (trees(i, :) - best_params);
%         direction = direction * fitness(target_index);
    end

    % ���������������ɿ����ֲ�����
    ori_value = rand(1, D);
    cauchy_value = tan((ori_value - 0.5) * pi);
 % �������Ǩ�Ʒ���
        direction = randn(1, D); % ���ɷ��ӱ�׼��̬�ֲ��������������

        % ���ݸ�����Ӧ�ȵ�������
        direction = direction .* fitness(i); % ��Ӧ��Խ�ͣ������Ӱ��ԽС
     newtrees(i, :) = trees(i, :) + step_length * cauchy_value.* direction;
% directions = 2 * randi([0, 1], N, D) - 1;
    % ������ľλ��
%     newtrees(i, :) = trees(i, :) + step_length * cauchy_value;
%     newtrees(i, :) = trees(i, :) + step_length * cauchy_value.* direction;
    % �߽���
    newtrees(i, :) = max(newtrees(i, :), dmin);
    newtrees(i, :) = min(newtrees(i, :), dmax);
end

% ������Ӧ��ֵ�Ƚ���λ�ú�ԭλ�ã���������ľλ��
for i = 1:N
    if minobj(i) < obj(i)
        trees(i, :) = newtrees(i, :);
%         obj(i)=minobj(i);
         
    end
end


 %%��ͷ����Ҫ�޸ģ�Xpos��Ҫ�޸ģ�Select the optimal fitness value��Ҫ˼���Ƿ���룬����numֵҲҪ����
 
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
        % ��¼convergence���
        convergence(t)=best;
        t=t+1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  end  
    gbest = best_params; %�������ֵ
    gbestval =best; %�������λ��
    fitcount=fes; %���غ�����������
end