classdef scp < matlab.mixin.SetGet
    
    % Class definition for scp: meant to model moment problems (scp is 
    % short for "sequence conic program").
    
    % Juan Kuntz, 16/02/2015, last edited 16/02/2015.
    
    properties 
        nvar = [];          % Dimension of underlying space.
        var = [];           % Structure containing symbols and number of components, like for @pol.
        
        obj = [];           % Objective, or vector of objectives, to be specified via polynomials.
        minmax = [];        % Specifies whether the objective should be min(obj) or max(obj).
        
        mass = [];          % The mass of the measure.
        supcon = [];        % Ineqquality constraints on the support of the measure, array of polnomials.
        seqeqcon = [];          % Linear equality constraints on sequence, array of sequences.
        seqineqcon = [];        % Linear inequality constraints on sequence, array of sequences.
        
        reltype = [];   % Type of relaxation {D,DD,SDD,FKW,PSD}
        relorder = [];  % Order of relaxation
        FW = [];        % Specifies the k of the FkW relaxation
        
        status = []; % Solved, unsolved, etc.
        
        sol = []; % A cell array, each entry of which contains the info on some solution.
        ops = []; % Solver options.
        
        yvar = []; % Stores the yalmip variables.
        ycons = []; % Stores the above constraints translated into yalmip constraints.
        yobj = []; % Stores the objective/s re-written in terms of the yalmip variables.
        
        A = [];
        b = [];
        F = [];
        f = [];
        
    end
    
    methods
        
        % Class constructor
        
        function obj = scp(varargin)
            if nargin ~= 0
                obj.nvar = varargin{1};
                obj.relorder = varargin{2};
            end
        end
        
        % Property sets with error messages.
        
        function obj = set.obj(obj,p)
            if ~isa(p,'pol')
                disp('Error: You must specify the objective/s constraints using polynomials.');
                return
            end
            obj.obj = p;
        end
        
        function obj = set.supcon(obj,p)
            if ~isa(p,'pol')
                disp('Error: You must specify the support constraints using polynomials.');
                return
            end
            obj.supcon = p;
        end
        
        function obj = set.seqeqcon(obj,p)
            if ~isa(p,'pol')
                disp('Error: You must specify the sequence constraints using polynomials.');
                return
            end
            obj.seqeqcon = p;
        end
        
        function obj = set.seqineqcon(obj,p)
            if ~isa(p,'pol')
                disp('Error: You must specify the sequence constraints using polynomials.');
                return
            end
            obj.seqineqcon = p;
        end
        
    end
end
