for j = 1:1:100
    for i = 1:1:100
        
        p_field(j,i) = P(i+(j-1)*100);
        t_field(j,i) = T(i+(j-1)*100);
        
    end
end

for i = 1:1:100
    for j = 1:1:101
        u_field(i,j) = u(j+(i-1)*101);
    end
end

for i = 1:1:101
    for j = 1:1:100
        v_field(i,j) = v(j+(i-1)*100);        
    end
end

% u-vel
for j = 1:1:100
    for i = 1:1:100
        u_vel(i,j) = (u_field(i,j)+u_field(i,j+1))/2;
    end
end

% v-vel
for j = 1:1:100
    for i = 1:1:100
        v_vel(i,j) = (v_field(i,j)+v_field(i+1,j))/2;
    end
end

%p_field = (p_field + 101132100)/101132100;

p_field_fig = figure
contourf(p_field,'ShowText','on')
colorbar
title 'Pressure Field (100x100 Grid)'
grid minor
saveas(p_field_fig, 'p_field.png')

t_field_fig = figure
contourf(t_field,'ShowText','on')
colorbar
title 'Temperature Field (100x100 Grid)'
grid minor
saveas(t_field_fig, 't_field.png')

v_field_fig = figure
quiver(u_vel,v_vel)
axis([0 101 0 101])
title 'Velocity Field (100x100 Grid)'
grid minor

startx = 0:10:100;
starty = ones(size(startx));
streamline(u_vel, v_vel, startx, startx)
saveas(v_field_fig, 'v_field.png')

u_con_fig = figure
contour(u_vel,'ShowText','on')
colorbar
title 'U-Velocity Contour (100x100 Grid)'
grid minor
saveas(u_con_fig, 'u_con.png')

v_con_fig = figure
contour(v_vel,'ShowText','on')
colorbar
title 'V-Velocity Contour (100x100 Grid)'
grid minor
saveas(v_con_fig, 'v_con.png')