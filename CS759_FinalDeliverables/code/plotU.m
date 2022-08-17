load uFinal.txt; %loads data into file called uFinal
d = size(uFinal);

dx = 1.0/(d(1)-1);
dy = 1.0/(d(2)-1);

xg = 0:dx:1;
yg = 0:dy:1;

[X,Y] = meshgrid(xg,yg); % create a mesh for plotting

surf(X,Y,uFinal); %make a surface plot
xlabel('x','FontSize',26);
ylabel('y','FontSize',26);
zlabel('u (x,y)','FontSize',26);
set(gca,'FontSize',24);


