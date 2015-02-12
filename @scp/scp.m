classdef scp < matlab.mixin.SetGet
    
    properties 
        nvar = [];          % Dimension of underlying space.
        supeqcon = [];        % Equality constraints on the support of the measure, array of polnomials.
        supineqcon = [];        % Ineqquality constraints on the support of the measure, array of polnomials.
        seqeqcon = [];          % Linear equality constraints on sequence, array of sequences.
        seqineqcon = [];        % Linear inequality constraints on sequence, array of sequences.
        
        
        reltype = []; % Type of relaxation {LP,SOCP,FKW1,...,FKWD,SDP}
        relorder = []; % Order of relaxation
        
        status = []; % Solved, unsolved, etc.
        
        sol = []; % A cell array, each entry of which contains the info on some solution.
        ops = []; % Solver options.
        
    end
    
    methods
        function obj = scp(varargin)
            if nargin ~= 0
                obj.relorder = varargin{1};
            end
        end
        
        function obj = set.supeqcon(obj,p)
            if ~isa(p,'pol')
                disp('Error: You must specify support constraints using polynomials.');
                return
            end
            obj.supeqcon = p;
        end
        
        function obj = set.supineqcon(obj,p)
            if ~isa(p,'pol')
                disp('Error: You must specify support constraints using polynomials.');
                return
            end
            obj.supineqcon = p;
        end
        
        function obj = set.seqeqcon(obj,s)
            if ~isa(p,'seq')
                disp('Error: You must specify sequence constraints using sequences.');
                return
            end
            obj.seqeqcon = s;
        end
        
        function obj = set.seqineqcon(obj,s)
            if ~isa(p,'seq')
                disp('Error: You must specify sequence constraints using sequences.');
                return
            end
            obj.seqineqcon = s;
        end
        
    end
end
