function Offspring = OperatorCross(Parent, lower, upper)
    Parent1 = Parent(1,:);
    D = size(Parent1,2);
    Offspring = Parent(2,:);
    dis = abs(Parent1-rand*Offspring);
    c1 = rand(1,D);
    F = 2*rand(1,D) ;
    G = c1 + rand(1,D) - F;
    M = abs(1 + (c1*1));
    A = G./M;
    f=1.2;

    if rand>0.5
        Offspring = Offspring + f.*cosh(dis./(A+1)).*exp(-cosh(dis./(A+1))/f);
        % 
    else
        Offspring = Offspring - f.*cosh(dis./(A+1)).*exp(-cosh(dis./(A+1))/f);
    end

%     s = randi(D,1,2);
%     if rand>0.5
%         Offspring(s) = Offspring([s(2), s(1)]);
%     else 
%         Offspring(s) = mod((5-Offspring(s)),3);
%     end
    if rand>0.5
        s = randi(D,1,round(0.2*D)); %0.2
        Offspring(s) = Offspring(s(end:-1:1));
    else
        Site = rand(1, length(Offspring))<0.2; % 0.3
        Offspring(Site) = mod(5-Offspring(Site),4);
    end

    Offspring = max(lower,min(Offspring, upper));
end