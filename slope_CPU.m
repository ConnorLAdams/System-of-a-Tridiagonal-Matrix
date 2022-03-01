% Calculating slope and intercept for project MAT 5080

function[slope,intercept] = slope_CPU(y,x,sample_size,n)


random = randi([150,n],[sample_size,1]);

slope = sum((log(y(random)) - log(y(random-50))) /...
    (log(x(random)) - log(x(random-50))));


% y = ax^b   a = y/x^b    log(y) = b*log(ax)      log(y) = b*log(x) + b*log(a) 
for i = 1:sample_size
a(i) = y(random(i))/(x(random(i))^slope);
end

intercept = sum(a)/sample_size;


end