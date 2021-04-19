% ----------------------------------------------------------------------- %
% Example of use of the funcion MDOMOPSO.m, which performs a              %
% Multi-Objective Particle Swarm Optimization (MOPSO),                    %
% based on the following references                                       %
% ----------------------------------------------------------------------- %
%   Author:  Davide Fioriti                                               %
%                                                                         %
%   Code modified from code by Victor Martinez Cagigal Copyright (c) 2017 %
%   Original code available at:                                           %
%   https://it.mathworks.com/matlabcentral/fileexchange/                  %
%        62074-multi-objective-particle-swarm-optimization-mopso          %
%                                                                         %
%   Date:    2/11/2020                                                    %
%   E-mail:  fioritidavidesubs (at) gmail (dot) com                       %
% ----------------------------------------------------------------------- %
%   References:                                                           %
%    [1]Fioriti, D., Lutzemberger, G., Poli, D., Duenas-Martinez, P.,     %
%       Micangeli, A., Coupling economic multi-objective optimization     %
%       and multiple design options: a business-oriented approach to      %
%       optimize an off-grid hybrid microgrid, International Journal of   %
%       Electrical Power and Energy Systems, 2021                         %
%                                                                         %
%    [2]Coello, C. A. C., Pulido, G. T., & Lechuga, M. S. (2004). Handling%
%       multiple objectives with particle swarm optimization. IEEE Tran-  %
%       sactions on evolutionary computation, 8(3), 256-279.              %
%                                                                         %
%    [3]Sierra, M. R., & Coello, C. A. C. (2005, March). Improving PSO-   %
%       based multi-objective optimization using crowding, mutation and ?-%
%       dominance. In International Conference on Evolutionary Multi-Crite%
%       rion Optimization (pp. 505-519). Springer Berlin Heidelberg.      %
% ----------------------------------------------------------------------- %
clear all; clc;

%% initialization
% Poloni function

fun = @(x) opt_test(x);
nVar = 2;
var_min = -pi.*ones(1,nVar);
var_max = pi.*ones(1,nVar);

%% execution of paretosearch for comparison purposes
rng default % for reproducibility
opt = optimoptions('paretosearch', 'Display', 'iter','PlotFcn','psplotparetof', 'ParetoSetSize', 100);
x = paretosearch(fun,nVar,[],[],[],[],var_min,var_max, [], opt);

%% execution of MDO-MOPSO
rng default % for reproducibility
% MOPSO
%option file with only the vectorized parameter set to false, the other default parameters will be used
options = struct('vectorized', false, 'Np', 500);

[REP, MDO] = MDOMOPSO(fun, nVar, var_min, var_max, options);
% Display info
disp('Repository fitness values are stored in REP.fit');
disp('Repository particles positions are store in REP.pos');
disp('History of particles positions are store in REP.history.pos');
disp('History of fitness values are store in REP.history.fit');
disp('History of auxiliary values are store in REP.history.aux');
disp('MDO particles positions are store in MDO.pos');
disp('MDO fitness values are store in MDO.fit');
disp('MDO auxiliary values are store in MDO.aux');


% Plot of pareto frontier and MDO points
[~, sort_ids] = sort(REP.fit(:,1));
figure
plot(REP.history.fit(:,1), REP.history.fit(:,2), '.k');
hold on
plot(MDO.fit(:,1), MDO.fit(:,2), 'xb');
plot(REP.fit(sort_ids,1), REP.fit(sort_ids,2), '-r', 'linewidth', 2);
legend('History points', 'MDO points', 'Pareto frontier');
xlim([0 18]);
ylim([0, 25]);
xlabel('f1');
ylabel('f2');
set(gca, 'Fontsize', 14);
title('MDO-MOPSO')


%% Example function: Poloni function
function [f_val, aux] = opt_test(x_data)
        A1 = 0.5*sin(1)-2*cos(1)+sin(2)-1.5*cos(2);
        A2 = 1.5*sin(1)-cos(1)+2*sin(2)-0.5*cos(2);
        B1 = @(x,y) 0.5.*sin(x)-2.*cos(x)+sin(y)-1.5.*cos(y);
        B2 = @(x,y) 1.5.*sin(x)-cos(x)+2.*sin(y)-0.5.*cos(y);
        f1 = @(x,y) 1+(A1-B1(x,y)).^2+(A2-B2(x,y)).^2;
        f2 = @(x,y) (x+3).^2+(y+1).^2;
        f_val = [f1(x_data(:,1),x_data(:,2)), f2(x_data(:,1),x_data(:,2))];
        aux = [f_val(:, 2), f_val(:, 1)];
end
    