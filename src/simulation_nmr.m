%% compare the true distribution and sampling distribution.
% Without loss of generality, we take Gaussian as example.
close all;
clear all;
clc;

%% Parameters
SampleSize = [10000,10000];
type = 'Uniform'; 

tic;
switch type
    case 'Uniform'
        data = rand(SampleSize);
        [Recover_Sample,pdf,binSize] = nmr(data);
        N_X = (length(Recover_Sample)-1)/2;
        Recover_Ex = ones(2*N_X+1,1)/(2*N_X+1);
        % for display, add zeros for uniform
        add_zero = round(N_X*0.8);  
        Recover_Ex = [zeros(add_zero,1);Recover_Ex;zeros(add_zero,1)];
        Recover_Sample = [zeros(add_zero,1);Recover_Sample;zeros(add_zero,1)];
        pdf = [zeros(2*add_zero,1);pdf;zeros(2*add_zero,1)];
        N_X = N_X + add_zero;
        pdf_Ex = conv(Recover_Ex,Recover_Ex); 
    case 'Gaussian'
        data = randn(SampleSize);
        [Recover_Sample,pdf,binSize] = nmr(data);
        N_X = (length(Recover_Sample)-1)/2;
        Recover_Ex = zeros(2*N_X+1,1);
        for i = -N_X:N_X
            Recover_Ex(i+N_X+1) = 1/2/pi*exp(-(i*binSize)^2/2);
        end
        Recover_Ex = Recover_Ex/sum(Recover_Ex);
        pdf_Ex = conv(Recover_Ex,Recover_Ex);
    case 'Cauchy'
        data = tan(pi*(rand(SampleSize)-1/2));
        data = data(abs(data)<40);    % in case of too large range
        [Recover_Sample,pdf,binSize] = nmr(data);
        N_X = (length(Recover_Sample)-1)/2;
        Recover_Ex = zeros(2*N_X+1,1);
        for i = -N_X:N_X
            Recover_Ex(i+N_X+1) = 1/pi/(1+(i*binSize)^2);
        end
        Recover_Ex = Recover_Ex/sum(Recover_Ex);
        pdf_Ex = conv(Recover_Ex,Recover_Ex);
    case 'Laplacian'
        data = rand(SampleSize);
        data = sign(0.5-data).*(1/sqrt(2)).*log(2*min(data,1-data));
        [Recover_Sample,pdf,binSize] = nmr(data);
        N_X = (length(Recover_Sample)-1)/2;
        Recover_Ex = zeros(2*N_X+1,1);
        for i = -N_X:N_X
            Recover_Ex(i+N_X+1) = 1/sqrt(2)*exp(-sqrt(2)*(abs(i)*binSize));
        end
        Recover_Ex = Recover_Ex/sum(Recover_Ex);
        pdf_Ex = conv(Recover_Ex,Recover_Ex);
end
toc;
N0 = 2*N_X;

%% Display
f1 = figure;
plot([-N0:N0]*binSize,pdf,'r','Linewidth',1.5);
hold on;
plot([-N0:N0]*binSize,pdf_Ex,'b--','Linewidth',1.5);
legend([type,' From sampling'],[type,' From expression']);
set(f1,'Position',[200,300,560,420]);
grid on;

f2 = figure;
subplot(1,2,1)
plot([-N_X:N_X]*binSize,Recover_Sample,'r','Linewidth',1.5);
hold on;
plot([-N_X:N_X]*binSize,Recover_Ex,'b--','Linewidth',1.5);
legend({'$f_{\hat{X}}(t)$','$f_{X}(t)$'},'Interpreter','latex','FontName', 'Calibri','FontSize',15);
set(f2,'Position',[760,300,1200,420]);
grid on;

subplot(1,2,2)
plot([-2*N_X:2*N_X]*binSize,conv(Recover_Sample,Recover_Sample),'r','Linewidth',1.5);
hold on;
plot([-N0:N0]*binSize,pdf,'b--','Linewidth',1.5);
legend({'$f_{\hat{X}}(t) * f_{\hat{X}}(t)$','$f_Y(t)$'},'Interpreter','latex','FontName', 'Calibri','FontSize',15);
grid on;