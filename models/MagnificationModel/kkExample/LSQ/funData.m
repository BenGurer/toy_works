function y = funData(soa,df) 
% soa = SOA in ms, and df = frequnecy separation in cents

f = inline('exp(b).* exp(m.*soa)','b','m','soa'); % m and b from log-lin regression of data. 
m = -1.044321/1000;
b = arrayfun(@(x)lcfB(x),df);
y = f(b,m,soa);


function B = lcfB(DF)

df = [0 200 600 1800];
b = [4.378866 4.211515 4.082468 3.848298];
B = b(DF==df);

