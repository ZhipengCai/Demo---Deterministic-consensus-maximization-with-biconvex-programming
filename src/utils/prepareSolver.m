% Prepare standard optimization solver.
% solver.LP = 'linprog'/'sedumi'/'gurobi' : Use linprog/sedumi/gurobi as LP solver.
% solver.SOCP = 'sedumi' : Use sedumi as SOCP solver. 
% If Gurobi is used, need to set gurobiPath to the installation folder of Gurobi

function solverFun = prepareSolver(solver)

addpath('src/utils');

%preparing SOCP solver
addpath('src/utils/sedumi/');

gurobiPath = '/opt/gurobi751/linux64/'; % Please specify another path if you would like to use a different Gurobi installation
%preparing LP solver
if (~strcmp(solver.LP,'gurobi'))

    disp(['Using ', solver.LP,' as the default solver....']);
    disp(['If you have Gurobi installed with a proper license, please set solver to gurobi ']);
    disp(['and specify the path to gurobi installation folder using the variable gurobiPath ']);
    disp(['in prepareSolver.m ']);
    
    if strcmp(solver.LP, 'sedumi')
        solverFun.LP = @sedumiLinProg;
    else
        solverFun.LP = @linprog;
    end
else
    
    solverFun.LP = @gurobiLinProg;
    disp('Using Gurobi as the default solver. Please make sure that you have a proper Gurobi license installed on your system.');
    disp([ 'Gurobi path is set to ' gurobiPath]);
    addpath([gurobiPath 'matlab/']);
    gurobi_setup;    
end
       
end