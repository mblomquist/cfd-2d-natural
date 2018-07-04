% Calculate Nusselt Number
Nu = zeros(size(t_field));

for i = 1:1:length(t_field)
    for j = 1:1:length(t_field)
        Nu(i,j) = (1-t_field(i,j))/(j/(length(t_field)+1));
    end
end

Nu_field = figure
contourf(Nu)
colorbar
title 'Nussult Number (100x100 Grid)'
grid minor
saveas(Nu_field, 'Nu_field.png')

for i = 1:1:length(t_field)
    Nu_bar(i) = sum(Nu(:,i))*(1/101);
end

Nu_distance = figure 
plot(1:1:length(t_field),Nu_bar)
title 'Nussult Number vs. Distance from Wall (100x100 Grid)'
grid minor
xlabel 'Distance from Hot Wall (dx)'
ylabel 'NuL'
saveas(Nu_distance,'Nu_distance.png')