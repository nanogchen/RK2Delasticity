clear;
clc;

% number of nodes
NodeNumbers= importdata('NodeNumbers.txt');
np = NodeNumbers(1);
npx = NodeNumbers(2);
npy = NodeNumbers(3);

X = importdata('xpts.txt'); % xpts(npx) x_start:dx:x_end;  
Y = importdata('ypts.txt'); % ypts(npx) y_start:dy:y_end;

[xx,yy] = meshgrid(X,Y);

xynodes = importdata('xynodes.txt');

dispV = importdata('dispV.txt');
dispV = reshape(dispV,[npx,npy]);

exact_dispV = importdata('exact_dispV.txt');
exact_dispV = reshape(exact_dispV,[npx,npy]);

subplot(1,2,1)
contourf(xx',yy',dispV,'ShowText','on')
title('Approx dispV')
xlabel('x')
ylabel('y')
axis equal

subplot(1,2,2)
contourf(xx',yy',exact_dispV,'ShowText','on')
title('Exact dispV')
xlabel('x')
ylabel('y')
axis equal

hold on 