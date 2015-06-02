function [value isterminal direction]=stopODE(t,y)

value=y(end)/10^10-.01;
isterminal=1;
direction = -1;