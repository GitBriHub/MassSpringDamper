%Final Project
%Author: Gabrielle (Bri) Jackson

%%Clear Cache
clear all %#ok<*CLALL>
close all
clc



%% Initial conditions
x0 = 1;       
v0 = 0;       
a0 = 2.5;  

dt = 1/300;  
t0 = 0;
tf = 10;  
n_steps = round((tf - t0) / dt);
time = linspace(t0, tf, n_steps);  

trials = [3, 200, 2; 4, 50, 45; 5, 125, 50];

%%Main Code
for trial = 1:1:size(trials, 1)
    m = trials(trial, 1); %mass
    k = trials(trial, 2); %spring constant
    c = trials(trial, 3); %damping ratio
   
    xhomo = zeros(1, n_steps); vhomo = zeros(1, n_steps);
    xinhomo = zeros(1, n_steps); vinhomo = zeros(1, n_steps);
    
    xhomo(1) = x0; xinhomo(1) = x0;
    vhomo(1) = v0; vinhomo(1) = v0;

    %%Video Writing
    video_file = VideoWriter(sprintf('Jackson_0540161_video_%d.mp4', trial), 'MPEG-4');
    video_file.FrameRate = 30;
    open(video_file);
    figure;

    framestep = round(1 / (video_file.FrameRate * dt));

     for step = 2:1:n_steps
        %homogeneous
        [xhomo(step), vhomo(step)]= VibrationPosition(xhomo(step-1), vhomo(step-1), m, k, c, 0, dt, 2);
        
        %inhomogeneous
        forcing = a0 * sin(time(step) / (2 * pi));
        [xinhomo(step), vinhomo(step)] = VibrationPosition(xinhomo(step-1), vinhomo(step-1), m, k, c, forcing, dt, 2);
        
        %Plots
        if mod(step,framestep) == 0
            subplot(1, 2, 1);
            plot(time(1:1:step), xhomo(1:1:step), 'b');
            title(['Homogeneous Response - Trial #', num2str(trial)]);
            xlabel('Time (s)'); ylabel('Position (m)');
            grid on;
    
            subplot(1, 2, 2);
            plot(time(1:1:step), xinhomo(1:1:step));
            title(['Inhomogeneous Response - Trial #', num2str(trial)]);
            xlabel('Time (s)'); ylabel('Position (m)');
            grid on;
    
            % Video
            current_frame = getframe(gcf);
            writeVideo(video_file, current_frame); 
        end
     end
    close(video_file);
end






%%Required Function
function [x, v] = VibrationPosition(xk, vk, m, k, c, f, dt, type)

    omega = sqrt(k / m); 
    E = c / (2 * sqrt(1 / m * k));

    % Forward Euler
    if type == 1 
        v = vk + dt * ((-2 * E * omega * vk) - omega^2 * xk + f);
        x = xk + dt * vk;

    % Runge-Kutta
    elseif type == 2 
        cx1=dt* vk;
        cv1=dt*((-2*E*omega*vk)-omega^2 *xk+f);

        cx2=dt *(vk +0.5 * cv1);
        cv2=dt *((-2 *E*omega *(vk +0.5 * cv1))-omega^2*(xk+ 0.5*cx1));

        cx3=dt*(vk+ 0.5*cv2);
        cv3=dt*((-2* E* omega*(vk + 0.5 * cv2))-omega^2*(xk + 0.5 * cx2));

        cx4=dt*(vk+cv3);
        cv4=dt*((-2*E*omega*(vk + cv3))-omega^2*(xk + cx3));

        x=xk+(1/6)*(cx1+2*cx2+2*cx3+cx4);
        v=vk+(1/6)*(cv1+2*cv2+2*cv3+cv4);
    else
        error('Invalid method type! Please enter 1 for F.E or 2 for R.K');
    end
end
