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

dispU = importdata('dispU.txt');
dispU = reshape(dispU,[npx,npy]);

exact_dispU = importdata('exact_dispU.txt');
exact_dispU = reshape(exact_dispU,[npx,npy]);

subplot(1,2,1)
contourf(xx',yy',dispU,'ShowText','on')
title('Approx dispU')
xlabel('x')
ylabel('y')
axis equal

subplot(1,2,2)
contourf(xx',yy',exact_dispU,'ShowText','on')
title('Exact dispU')
xlabel('x')
ylabel('y')
axis equal