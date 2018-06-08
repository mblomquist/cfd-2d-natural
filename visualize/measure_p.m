
dP_x = zeros(100,99);
dP_y = zeros(99,100);

for i = 1:1:99
    for j = 1:1:100
        dP_y(i,j) = p_field(i+1,j)-p_field(i,j);
        
    end 
end

for i = 1:1:100
    for j = 1:1:99
        dP_x(i,j) = p_field(i,j+1)-p_field(i,j);
        
    end 
end

figure
contourf(dP_y)
title 'dP_y'
colorbar
grid minor

figure
contourf(dP_x)
title 'dP_x'
colorbar
grid minor