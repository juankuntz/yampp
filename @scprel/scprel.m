classdef scprel < matlab.mixin.SetGet

% A class to store and display the (finite dimensional) relaxations of an 
% scp. It essentially just is a structure with some specific fields and an
% overloaded display method.

% Juan Kuntz, 15/03/2015, last edited 16/04/2015.

properties
    
    % Specification data.
    
    rtyp = [];      % Type of relaxation D,DD,SDD,FKW, or PSD.
    rord = [];      % Order of relaxation.
    objf = [];      % Function specifying objective.
    objs = [];      % String specifying whether we're maximising or 
                    % minimising; either inf or sup.
    mult = [];      % Multiplier.
    FW = [];        % Factor width.
    
    
    % Yalmip relaxation data.
    
    yvar = [];      % Yalmip variables.
    zvar = []       % z-variable = shift(mult,yvar), in case of non-unity multipler.
    ycons = [];     % Yalmip constraints.
    yobj = [];      % Yalmip objective.
    info = [];      % Yalmip diagnostics.
    ops = [] ;      % Yalmip options.
    
    % Conic program data; used to compute dual residues.
    
    F = [];     % Equality constraints: Fy == f.
    f = [];
    
    b = [];     % Objective <b,y>.
    
    A = [];     % Cone constraints: C-A^T(y) >= 0.
    C = [];
    
    % Solution data
    
    pval = [];  % Optimal value of the primal problem.
    dval = [];  % Optimal value of the dual problem.
    
    ppoint = [];    % Optimal point of the primal problem.
    X = [];
    lambda = [];    % (lambda,X) is the optimal point of the dual problem.
    
    pres = [];  % Primal residues.
    dres = [];  % Dual residues.
    
end

methods
    

end    



end

