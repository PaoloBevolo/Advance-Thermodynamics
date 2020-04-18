%This script executes all the functions for the different boundary
%conditions. It is used mainly for testing and to show the capabilities of
%the library.

%Version 2017.1
%2014-2017 Paolo Bardella

clear;
clear classes;
Unattended=true;

if ~exist('shape','class')
    addpath(genpath('..'));
end

Tests=[
    {'Showing basic shapes operations','BasicShapesOperations();'};     %1
    {'Homogeneous Dirichlet B.C.s','dirichletHomo_Main();'};            %2
    {'Homogeneous Dirichlet B.C.s; Diffusion, Transport, Reaction',...  %3
        'dirichletHomo_DiffTransReact_Main();'};
    {'Non - Homogeneous Dirichlet B.C.s','dirichletNonHomo_Main();'};   %4
    {'Non - Homo Dirichlet B.C.s; Diffusion, Transport, Reaction',...   %5
        'dirichletNonHomo_DiffTrans_Main();'};
    {'Non - Homogeneous Neumann B.C.s','neumannNonHomo_Main();'};       %6
    {'Robin B.C.s','robin_Main();'};                                    %7
    {'Periodic B.C.s','periodic_Main();'};                              %8
    {'Coupled Systems','coupledDirichlet_Main();'};                     %9
    {'Coupled Systems with Neumann BCs','coupledNeumann_Main();'};      %10
    {'Darcy','darcy_Main();'};                                          %11
    {'Cylindrical coordinates','cylindrical_Main();'};                  %12
    {'Chimney','chimney_Main();'};                                      %13
    {'Heat Equation (external source, temporal evolution)',...
                'heatEquationExtForce_Main();'};                        %14
    {'Heat Equation (time dependent B.C.s)',...
                'heatEquationVariableDirichlet_Main();'};               %15
    {'Elastic Equation (temporal evolution)','elasticMembrane_Main();'};%16    
    {'Elastic Equation (temporal evolution, lumped mass matrix)',...    %17
        'elasticMembraneLumping_Main();'};    
    ];
TestToExecute=1:length(Tests);

for k=TestToExecute,        
if ~Unattended
    fprintf('Press ENTER to start example %d : %s\n',k,Tests{k,1});
    pause;
end
    close all;
    fprintf('Starting example %d : %s\n',k,Tests{k,1});
    try        
        eval(Tests{k,2});
        drawnow();
    catch exc
         fprintf('Example %d : %s ABNORMALLY terminated\n',k,Tests{k,1});
         rethrow(exc);
    end
end
