function Pop = Piecewise(N, lower, upper)
D = size(lower,2);
Pop = zeros(N, D);
Pop(1,:) = rand(1,D);
P = 0.4;
for i = 1:N
    ind = Pop(i,:)>=0 & Pop(i,:)<P;
    Pop(i+1, ind) = Pop(i,ind)/P;
    
    ind = Pop(i,:)>=P & Pop(i,:)<0.5;
    Pop(i+1, ind) = (Pop(i,ind)-P)/(0.5-P);

    ind = Pop(i,:)>=0.5 & Pop(i,:)<1-P;
    Pop(i+1, ind) = (1-P-Pop(i, ind))/(0.5-P);

    ind = Pop(i,:)>=1-P & Pop(i,:)<1;
    Pop(i+1, ind) = (1-Pop(i,ind))/P;
end
Pop(1,:) = [];
Pop = Pop .* (upper-lower) + lower;
end