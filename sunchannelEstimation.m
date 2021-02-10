clc; clear all; close all;
sigma_sq_dB = -10:1:10; %dB noise variance
sigma_sq = 10.^(sigma_sq_dB/10);  %Noise variance
h = [1;1;1;1];  %True parameter
N = 8; %Number of observations
num_exp = 1000;  %number of experiments
X1 = randn(8,4);
X2 = orth(X1);
norm1 = norm(X1,'fro');
X1 = X1/norm1;
norm2 = norm(X2, 'fro');
X2 = X2/norm2;

for i=1:length(sigma_sq)
    y1 = X1*h + sqrt(sigma_sq(i))*randn(N,num_exp);
    hhat1 = (inv(transpose(X1)*X1))*transpose(X1)*y1;
    total_error1 = sum((hhat1 - h).^2,1);
    error1 = sum(total_error1);
    mse1(i) = error1/num_exp;
    theoretical_mse_gauss(i) = sigma_sq(i)*trace(inv(transpose(X1)*X1));
end 


for i=1:length(sigma_sq)
    y2 = X2*h + sqrt(sigma_sq(i))*randn(N,num_exp);
    hhat2 = (inv(transpose(X2)*X2))*transpose(X2)*y2;
    total_error2 = sum((hhat2 - h).^2,1);
    error2 = sum(total_error2);
    mse2(i) = error2/num_exp;
    theoretical_mse_orth(i) = sigma_sq(i)*trace(inv(transpose(X2)*X2));
end

plot(sigma_sq_dB, mse1, 'b','linewidth',2.5);
hold on;
scatter(sigma_sq_dB, theoretical_mse_gauss, 'r','linewidth',2.5);
hold on;

plot(sigma_sq_dB, mse2, 'r','linewidth',2.5);
hold on;
scatter(sigma_sq_dB, theoretical_mse_orth, 'b','linewidth',2.5);