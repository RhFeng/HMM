function [t,p]=tp(num,x,index)

% index
% 1 for downward -1 for upward

n=length(x)-1;
% m=max(unique(x));
m=num;
p=zeros(m,m);

for i = 1:n
  p(x(i), x(i + 1)) = p(x(i), x(i + 1)) + 1;
end

if index == 1
    p = p;
end

if index == -1
    p = p';
end
    
t=p;

for i = 1:m
  p(i, :) = p(i, :) / sum(p(i, :));
end

for i=1:m
    for j=1:m
        if isnan(p(i,j))
           p(i,j)=0;
        end
    end
end

end