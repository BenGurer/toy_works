function erb = funCalERB(ctip,ctail,w)

SF = 25/2; N = 2^15; DF = SF/N; f = DF*(0:N-1); 
FC = 1; ERB1 = funErb(FC);
ptip = ctip*4/ERB1; ptail = ctail*4/ERB1; w = 10.^(w/10);

g = abs((f-FC)/FC);
erb = zeros(1,length(ptip));
for I = 1:length(ptip)
    fw = (1-w(I))*(1+ptip(I)*g).*exp(-ptip(I)*g)+w(I)*(1+ptail(I)*g).*exp(-ptail(I)*g);
    erb(I) = sum(fw*DF);
end



