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

U = importdata('solU.txt');

exactU = importdata('ExactsolU.txt');
exactU = reshape(exactU,[npx,npy]);

subplot(1,2,1)
contourf(xx',yy',U,'ShowText','on')
title('Approx U')
xlabel('x')
ylabel('y')
axis equal

subplot(1,2,2)
contourf(xx',yy',exactU,'ShowText','on')
title('Exact U')
xlabel('x')
ylabel('y')
axis equal