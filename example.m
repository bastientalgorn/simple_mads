%-------------------------------------------------------------------------------------%
%  simple_mads                                                                        %
%                                                                                     %
%  A simple matlab version of the Mesh Adaptive Direct Search algorithm               %
%  for constrained derivative free optimization.                                      %
%  Version 2.0.3                                                                      %
%                                                                                     %
%  Copyright (C) 2012-2017  Bastien Talgorn - McGill University, Montreal             %
%                                                                                     %
%  Author: Bastien Talgorn                                                            %
%  email: bastientalgorn@fastmail.com                                                 %
%                                                                                     %
%  This program is free software: you can redistribute it and/or modify it under the  %
%  terms of the GNU Lesser General Public License as published by the Free Software   %
%  Foundation, either version 3 of the License, or (at your option) any later         %
%  version.                                                                           %
%                                                                                     %
%  This program is distributed in the hope that it will be useful, but WITHOUT ANY    %
%  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A    %
%  PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.   %
%                                                                                     %
%  You should have received a copy of the GNU Lesser General Public License along     %
%  with this program. If not, see <http://www.gnu.org/licenses/>.                     %
%                                                                                     %
%  You can find information on simple_mads at                                         %
%  https://github.com/bastientalgorn/simple_mads                                      %
%-------------------------------------------------------------------------------------%

close all
clear all
clc
disp('========== An exemple for simple_mads.m ==============');

% Dimension of the problem
N = 2;

% blackbox handle
% 1st component is the objective function
% next components are the constraints
bb_handle = @(x) [norm(x-1) x(2)-x(1)+2];

% Bounds:
lb = -5*ones(1,N);
ub = +5*ones(1,N);

% Starting points (one line per starting points)
px0 = 10;
LH = lhsdesign(px0,N);
X0 = ones(px0,1)*lb + LH.*( ones(px0,1)*(ub-lb) );

% MADS options
mads_options.budget      = 100;
mads_options.display     = true;
mads_options.check_cache = true;

% Run the MADS algorithm
[Xmin,fmin,output] = simple_mads(X0,bb_handle,lb,ub,mads_options);

% Display the results
disp('Best feasible point:')
disp(['Xmin = [ ' num2str(Xmin,6) ' ]'])
disp(['fmin = ' num2str(fmin,6)])

disp('Best possible point:')
Xstar = [2 0];
BBOstar = bb_handle(Xstar);
fstar = BBOstar(1);
disp(['X* = [ ' num2str(Xstar,6) ' ]'])
disp(['f* = ' num2str(fstar,6)])
