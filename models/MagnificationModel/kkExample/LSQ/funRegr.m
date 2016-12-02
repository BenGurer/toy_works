function [M,b] = funRegr(x,y)

meY = mean(y); N = length(x);
M = (sum(x.*meY)-sum(x)*sum(meY)/N)/(sum(x.^2)-sum(x)^2/N);
b = zeros(size(y,1),1);
for I = 1:size(y,1)
    b(I) = mean(y(I,:))-M*mean(x);
end
