function np=code2DNA(p)
p = char(p);
[m,l]=size(p);
for i=1:m
    for j=1:l %j的值是：初始值是1,以2为步长，到py结束
        switch p(i,j) 
            case 'C'
                np(i,j)= 0;
            case 'T'
                np(i,j)= 1;
            case 'A'
                np(i,j)= 2;
            case 'G'
                np(i,j)= 3;
        end
    end
end

