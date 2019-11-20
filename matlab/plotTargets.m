function plotTargets()

figure
hold all
circles(0,0,3.5,'facecolor','none','linewidth',1)
theta = [pi/2, -pi/2, pi, 0, 3*pi/4, -3*pi/4, pi/4, -pi/4];
x = .25*7*cos(theta);
y = .25*7*sin(theta);
for i = 1:8
    circles(x(i),y(i),.4,'facecolor','none','linewidth',1)
end

end