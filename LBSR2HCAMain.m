clc
clear all
format long
tic;
DNASet = [];
%DNASet = code2DNA(DNASet);
X_min=0; %区间大小
X_max=3;
D=20;
iter_max=50;
popsize=20;
M=5;
lower = repmat(X_min, 1, D);
upper = repmat(X_max, 1, D);
s=size(DNASet,1)+1;
file = fopen('times');
while s <= 7
    t0 = cputime;
    fprintf("strand: %d\n",s)
    X = round(Piecewise(popsize, lower, upper));
    X = repair(popsize,D,X);
    Population.X = X;
    
    objs=zeros(popsize,M);
    for p=1:popsize
        f = fit6([DNASet; X(p,:)]);
        TmVar = var(f(:,5));
        objs(p,1:M) = mean(f(:,1:M),1);
        objs(p,5) = TmVar;
    end
    Population.objs = objs;

    num_vec = 100;
    r = 1 + 1/getH(M,popsize);
    [W, num_vec] = UniformVector(num_vec, M);
    PopObj = Population.objs;
    fmax = max(PopObj, [], 1);
    fmin = min(PopObj,[],1);
    PopObj  = (PopObj-repmat(fmin,size(PopObj,1),1))./repmat(fmax-fmin,size(PopObj,1),1);
    PopObj = [PopObj;zeros(1,M)];
    tensor = InitializeUtilityTensor(popsize,PopObj,W,r,num_vec);
    worst = popsize+1;

    if size(DNASet,1)==7
        DNASet(1,:) = [];
    end
    for it = 1:iter_max
    fprintf("%d :",it);
    for i=1:popsize
        fprintf("%d ", i);
        select = randi(popsize,1,2);
        if i~=round(popsize/2)
            Offspring = OperatorGAhalf(Population.X(select,:),lower, upper);
            Offspring=round(Offspring);
        else
            Offspring = round(Piecewise(1, lower, upper));
        end
        
        Offspring = repair(1,20,Offspring);
        if any(sum(Offspring==Population.X,2)==D)
            continue;
        end
        Offspring = reshape(Offspring',1,D);

        if any(sum(Offspring==Population.X,2)==D)
            continue;
        end
        if worst ==1 
            Population1 = [Offspring; Population.X];
        elseif worst == popsize+1
            Population1 = [Population.X; Offspring];
        else
            Population1 = [Population.X(1:worst-1,:); Offspring; Population.X(worst:end,:)];
        end

        Population.X = Population1;
        
        objs = zeros(popsize, M);
        for i = 1:popsize+1
            f = fit6([DNASet; Population.X(i,:)]);
            TmVar = var(f(:,5));
            objs(i,1:M) = mean(f(:,1:M),1);
            objs(i,5) = TmVar;
        end
        Population.objs = objs;
        [Population,worst,tensor] = Reduce(popsize,M,Population,W,worst,tensor,r,num_vec);
        end
    fprintf("\n");
    end
    % 选择一条DNA加入DNAset
    objs = zeros(popsize, M);
    for i = 1:popsize
        f = fit6([DNASet; Population.X(i,:)]);
        TmVar = var(f(:,5));
        objs(i,1:M) = mean(f(:,1:M),1);
        objs(i,5) = TmVar;
    end
    Zmin = min(objs,[],1);
    objs = objs+1/M;
    Zmax = max(objs,[],1);
    objs = (objs-Zmin) ./ (Zmax-Zmin);
    pro = prod(objs,2);
    best = find(pro==min(pro));
    best = best(1);
    DNASet = [DNASet; Population.X(best,:)];
    [DNAcode2(DNASet) string(fit6(DNASet))]
    iterationTime = cputime - t0;
    iterationTime 
    fprintf(file, "DNA order %d, running time: %f(s)\n", s, iterationTime)
    s = s+1;
end
runningTime = toc;
fprintf(file, "total running time %f(s)", runningTime);
f = fit6(DNASet);
for i = 1:5
    fprintf("%f ", mean(f(:,i)));
end
fprintf("(%f)\n", var(f(:,5)));