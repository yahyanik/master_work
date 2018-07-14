function [x,fval,exitFlag,output,population,scores] = ga
(fun,nvars,Aineq,bineq,Aeq,beq,lb,ub,nonlcon,intcon,options)
%GA    Constrained optimization using genetic algorithm.
%   GA attempts to solve problems of the following forms:
%       min F(X)  subject to:  A*X  <= B, Aeq*X  = Beq (linear constraints)
%        X                     C(X) <= 0, Ceq(X) = 0 (nonlinear constraints)
%                              LB <= X <= UB 
%                              X(i) integer, where i is in the index
%                              vector INTCON (integer constraints)
%
%   Note: If INTCON is not empty, then no equality constraints are allowed.
%   That is:-
%   * Aeq and Beq must be empty
%   * Ceq returned from NONLCON must be empty
%
%   X = GA(FITNESSFCN,NVARS) finds a local unconstrained minimum X to the
%   FITNESSFCN using GA. NVARS is the dimension (number of design
%   variables) of the FITNESSFCN. FITNESSFCN accepts a vector X of size
%   1-by-NVARS, and returns a scalar evaluated at X.
%
%   X = GA(FITNESSFCN,NVARS,A,b) finds a local minimum X to the function
%   FITNESSFCN, subject to the linear inequalities A*X <= B. Linear
%   constraints are not satisfied when the PopulationType option is set to
%   'bitString' or 'custom'. See the documentation for details.
%
%   X = GA(FITNESSFCN,NVARS,A,b,Aeq,beq) finds a local minimum X to the
%   function FITNESSFCN, subject to the linear equalities Aeq*X = beq as
%   well as A*X <= B. (Set A=[] and B=[] if no inequalities exist.) Linear
%   constraints are not satisfied when the PopulationType option is set to
%   'bitString' or 'custom'. See the documentation for details.
%
%   X = GA(FITNESSFCN,NVARS,A,b,Aeq,beq,lb,ub) defines a set of lower and
%   upper bounds on the design variables, X, so that a solution is found in
%   the range lb <= X <= ub. Use empty matrices for lb and ub if no bounds
%   exist. Set lb(i) = -Inf if X(i) is unbounded below;  set ub(i) = Inf if
%   X(i) is unbounded above. Linear constraints are not satisfied when the
%   PopulationType option is set to 'bitString' or 'custom'. See the 
%   documentation for details.
%
%   X = GA(FITNESSFCN,NVARS,A,b,Aeq,beq,lb,ub,NONLCON) subjects the
%   minimization to the constraints defined in NONLCON. The function
%   NONLCON accepts X and returns the vectors C and Ceq, representing the
%   nonlinear inequalities and equalities respectively. GA minimizes
%   FITNESSFCN such that C(X)<=0 and Ceq(X)=0. (Set lb=[] and/or ub=[] if
%   no bounds exist.) Nonlinear constraints are not satisfied when the
%   PopulationType option is set to 'bitString' or 'custom'. See the 
%   documentation for details.
%
%   X = GA(FITNESSFCN,NVARS,A,b,Aeq,beq,lb,ub,NONLCON,options) minimizes
%   with the default optimization parameters replaced by values in the
%   structure OPTIONS. OPTIONS can be created with the GAOPTIMSET function.
%   See GAOPTIMSET for details.
%
%   X = GA(FITNESSFCN,NVARS,A,b,[],[],lb,ub,NONLCON,INTCON) requires that
%   the variables listed in INTCON take integer values. Note that GA does
%   not solve problems with integer and equality constraints. Pass empty
%   matrices for the Aeq and beq inputs if INTCON is not empty.
%
%   X = GA(FITNESSFCN,NVARS,A,b,[],[],lb,ub,NONLCON,INTCON,options)
%   minimizes with integer constraints and the default optimization
%   parameters replaced by values in the structure OPTIONS. OPTIONS can be
%   created with the GAOPTIMSET function. See GAOPTIMSET for details.
%
%   X = GA(PROBLEM) finds the minimum for PROBLEM. PROBLEM is a structure
%   that has the following fields:
%       fitnessfcn: <Fitness function>
%            nvars: <Number of design variables>
%            Aineq: <A matrix for inequality constraints>
%            bineq: <b vector for inequality constraints>
%              Aeq: <Aeq matrix for equality constraints>
%              beq: <beq vector for equality constraints>
%               lb: <Lower bound on X>
%               ub: <Upper bound on X>
%          nonlcon: <Nonlinear constraint function>
%           intcon: <Index vector for integer variables>
%          options: <Options structure created with GAOPTIMSET>
%         rngstate: <State of the random number generator>
%
%   [X,FVAL] = GA(FITNESSFCN, ...) returns FVAL, the value of the fitness
%   function FITNESSFCN at the solution X.
%
%   [X,FVAL,EXITFLAG] = GA(FITNESSFCN, ...) returns EXITFLAG which
%   describes the exit condition of GA. Possible values of EXITFLAG and the
%   corresponding exit conditions are
%
%     1 Average change in value of the fitness function over
%        options.StallGenLimit generations less than options.TolFun and
%        constraint violation less than options.TolCon.
%     3 The value of the fitness function did not change in
%        options.StallGenLimit generations and constraint violation less
%        than options.TolCon.
%     4 Magnitude of step smaller than machine precision and constraint
%        violation less than options.TolCon. This exit condition applies
%        only to nonlinear constraints.
%     5 Fitness limit reached and constraint violation less than
%        options.TolCon. 
%     0 Maximum number of generations exceeded.
%    -1 Optimization terminated by the output or plot function.
%    -2 No feasible point found.
%    -4 Stall time limit exceeded.
%    -5 Time limit exceeded.
%
%   [X,FVAL,EXITFLAG,OUTPUT] = GA(FITNESSFCN, ...) returns a
%   structure OUTPUT with the following information:
%             rngstate: <State of the random number generator before GA started>
%          generations: <Total generations, excluding HybridFcn iterations>
%            funccount: <Total function evaluations>
%        maxconstraint: <Maximum constraint violation>, if any
%              message: <GA termination message>
%
%   [X,FVAL,EXITFLAG,OUTPUT,POPULATION] = GA(FITNESSFCN, ...) returns the
%   final POPULATION at termination.
%
%   [X,FVAL,EXITFLAG,OUTPUT,POPULATION,SCORES] = GA(FITNESSFCN, ...) returns
%   the SCORES of the final POPULATION.
%
%
%   Example:
%     Unconstrained minimization of 'rastriginsfcn' fitness function of
%     numberOfVariables = 2
%      x = ga(@rastriginsfcn,2)
%
%     Display plotting functions while GA minimizes
%      options = gaoptimset('PlotFcns',...
%        {@gaplotbestf,@gaplotbestindiv,@gaplotexpectation,@gaplotstopping});
%      [x,fval,exitflag,output] = ga(@rastriginsfcn,2,[],[],[],[],[],[],[],options)
%
%   An example with inequality constraints and lower bounds
%    A = [1 1; -1 2; 2 1];  b = [2; 2; 3];  lb = zeros(2,1);
%    % Use mutation function which can handle constraints
%    options = gaoptimset('MutationFcn',@mutationadaptfeasible);
%    [x,fval,exitflag] = ga(@lincontest6,2,A,b,[],[],lb,[],[],options);
%
%     FITNESSFCN can also be an anonymous function:
%        x = ga(@(x) 3*sin(x(1))+exp(x(2)),2)
%
%   If FITNESSFCN or NONLCON are parameterized, you can use anonymous
%   functions to capture the problem-dependent parameters. Suppose you want
%   to minimize the fitness given in the function myfit, subject to the
%   nonlinear constraint myconstr, where these two functions are
%   parameterized by their second argument a1 and a2, respectively. Here
%   myfit and myconstr are MATLAB file functions such as
%
%        function f = myfit(x,a1)
%        f = exp(x(1))*(4*x(1)^2 + 2*x(2)^2 + 4*x(1)*x(2) + 2*x(2) + a1);
%
%   and
%
%        function [c,ceq] = myconstr(x,a2)
%        c = [1.5 + x(1)*x(2) - x(1) - x(2);
%              -x(1)*x(2) - a2];
%        % No nonlinear equality constraints:
%         ceq = [];
%
%   To optimize for specific values of a1 and a2, first assign the values
%   to these two parameters. Then create two one-argument anonymous
%   functions that capture the values of a1 and a2, and call myfit and
%   myconstr with two arguments. Finally, pass these anonymous functions to
%   GA:
%
%     a1 = 1; a2 = 10; % define parameters first
%     % Mutation function for constrained minimization
%     options = gaoptimset('MutationFcn',@mutationadaptfeasible);
%     x = ga(@(x)myfit(x,a1),2,[],[],[],[],[],[],@(x)myconstr(x,a2),options)
%
%   Example: Solving a mixed-integer optimization problem
%   An example of optimizing a function where a subset of the variables are
%   required to be integers:
%   
%   % Define the objective and call GA. Here variables x(2) and x(3) will
%   % be integer. 
%   fun = @(x) (x(1) - 0.2)^2 + (x(2) - 1.7)^2 + (x(3) -5.1)^2; 
%   x = ga(fun,3,[],[],[],[],[],[],[],[2 3])
%          
%   See also GAOPTIMSET, FITNESSFUNCTION, GAOUTPUTFCNTEMPLATE, PATTERNSEARCH, @.

%   Copyright 2003-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.9 $  $Date: 2011/05/09 00:50:31 $

% If the first arg is not a gaoptimset, then it's a fitness function followed by a genome
% length. Here we make a gaoptimset from the args.
defaultopt = struct('PopulationType', 'doubleVector', ...
    'PopInitRange', [0;1], ...
    'PopulationSize', 20, ...
    'EliteCount', 2, ...
    'CrossoverFraction', 0.8, ...
    'MigrationDirection','forward', ...
    'MigrationInterval',20, ...
    'MigrationFraction',0.2, ...
    'Generations', 100, ...
    'TimeLimit', inf, ...
    'FitnessLimit', -inf, ...
    'StallGenLimit', 50, ...
    'StallTimeLimit', inf, ...
    'TolFun', 1e-6, ...
    'TolCon', 1e-6, ...
    'InitialPopulation',[], ...
    'InitialScores', [], ...
    'InitialPenalty', 10, ...
    'PenaltyFactor', 100, ...
    'PlotInterval',1, ...
    'CreationFcn',@gacreationuniform, ...
    'FitnessScalingFcn', @fitscalingrank, ...
    'SelectionFcn', @selectionstochunif, ...
    'CrossoverFcn',@crossoverscattered, ...
    'MutationFcn',{{@mutationgaussian 1  1}}, ...
    'HybridFcn',[], ...
    'Display', 'final', ...
    'PlotFcns', [], ...
    'OutputFcns', [], ...
    'Vectorized','off', ...
    'UseParallel', 'never');

% Check number of input arguments
errmsg = nargchk(1,11,nargin);
if ~isempty(errmsg)
    error(message('globaloptim:ga:numberOfInputs', errmsg));
end

% If just 'defaults' passed in, return the default options in X
if nargin == 1 && nargout <= 1 && isequal(fun,'defaults')
    x = defaultopt;
    return
end

if nargin < 11, options = [];
    if nargin < 10,  intcon = [];
        if nargin < 9,  nonlcon = [];
            if nargin < 8, ub = [];
                if nargin < 7, lb = [];
                    if nargin <6, beq = [];
                        if nargin <5, Aeq = [];
                            if nargin < 4, bineq = [];
                                if nargin < 3, Aineq = [];
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

% Is third argument a structure
if nargin == 3 && isstruct(Aineq) % Old syntax
    options = Aineq; Aineq = [];
end

% Is tenth argument a structure? If so, integer constraints have not been
% specified
if nargin == 10 && isstruct(intcon)
    options = intcon;
    intcon = [];
end

% One input argument is for problem structure
if nargin == 1
    if isa(fun,'struct')
        [fun,nvars,Aineq,bineq,Aeq,beq,lb,ub,nonlcon,intcon,rngstate,options] = separateOptimStruct(fun);
        % Reset the random number generators
        resetDfltRng(rngstate);
    else % Single input and non-structure. 
        error(message('globaloptim:ga:invalidStructInput'));
    end
end

% If fun is a cell array with additional arguments get the function handle
if iscell(fun)
    FitnessFcn = fun{1};
else
    FitnessFcn = fun;
end

% Only function handles or inlines are allowed for FitnessFcn
if isempty(FitnessFcn) ||  ~(isa(FitnessFcn,'inline') || isa(FitnessFcn,'function_handle'))
    error(message('globaloptim:ga:needFunctionHandle'));
end

% We need to check the nvars here before we call any solver
valid =  isnumeric(nvars) && isscalar(nvars)&& (nvars > 0) ...
    && (nvars == floor(nvars));
if ~valid
    error(message('globaloptim:ga:notValidNvars'));
end

% Specific checks and modification of options for mixed integer GA
if ~isempty(intcon)   
    % Check whether the user has specified options that the mixed integer
    % solver will either ignore or error.
    gaminlpvalidateoptions(options);    
    % If a user doesn't specify PopInitRange, we want to set it to the
    % bounds when we create the initial population. Need to store a flag
    % that indicates whether the user has specified PopInitRange so we can
    % do this in the creation function.
    UserSpecPopInitRange = isa(options, 'struct') && ...
        isfield(options, 'PopInitRange') && ~isempty(options.PopInitRange);
    % Change the default options for PopulationSize and EliteCount here.
    defaultopt.PopulationSize = max(min(10*nvars, 100), 40);
    defaultopt.EliteCount = floor(0.05*defaultopt.PopulationSize);
end

user_options = options;
% Use default options if empty
if ~isempty(options) && ~isa(options,'struct')
        error(message('globaloptim:ga:optionsNotAStruct'));
elseif isempty(options)
    options = defaultopt;
end
% Take defaults for parameters that are not in options structure
options = gaoptimset(defaultopt,options);

% All inputs should be double
try
    dataType = superiorfloat(nvars,Aineq,bineq,Aeq,beq,lb,ub);
    if ~isequal('double', dataType)
        error(message('globaloptim:ga:dataType'))
    end
catch
    error(message('globaloptim:ga:dataType'))
end

[x,fval,exitFlag,output,population,scores,FitnessFcn,nvars,Aineq,bineq,Aeq,beq,lb,ub, ...
    NonconFcn,options,Iterate,type] = gacommon
(nvars,fun,Aineq,bineq,Aeq,beq,lb,ub,nonlcon,intcon,options,user_options);

if exitFlag < 0
    return;
end

% Turn constraints into right size if they are empty.
if isempty(Aineq)
    Aineq = zeros(0,nvars);
end
if isempty(bineq)
    bineq = zeros(0,1);
end
if isempty(Aeq)
    Aeq = zeros(0,nvars); 
end
if isempty(beq)
    beq = zeros(0,1);
end

% Call appropriate single objective optimization solver
if ~isempty(intcon)   
    [x,fval,exitFlag,output,population,scores] = gaminlp(FitnessFcn,nvars, ...
        Aineq,bineq,Aeq,beq,lb,ub,NonconFcn,intcon,options,output,Iterate,...
        UserSpecPopInitRange);
else    
    switch (output.problemtype)
        case 'unconstrained'
            [x,fval,exitFlag,output,population,scores] = gaunc(FitnessFcn,nvars, ...
                options,output,Iterate);
        case {'boundconstraints', 'linearconstraints'}
            [x,fval,exitFlag,output,population,scores] = galincon(FitnessFcn,nvars, ...
                Aineq,bineq,Aeq,beq,lb,ub,options,output,Iterate);
        case 'nonlinearconstr'
            [x,fval,exitFlag,output,population,scores] = gacon(FitnessFcn,nvars, ...
                Aineq,bineq,Aeq,beq,lb,ub,NonconFcn,options,output,Iterate,type);
    end
end