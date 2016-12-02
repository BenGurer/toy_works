function x = spacing

n = 10;
y = 0.5:0.1:12;
for i = 1:n
    
    x(i) = sqrt(2).*y(i);
    
    plot (x)
end

end