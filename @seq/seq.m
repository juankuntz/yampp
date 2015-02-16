classdef seq
    
    % Missing comments.
    
    % Juan Kuntz, 16/02/2015, last modified 16/02/2015.

    properties
        coef = []; % 2 x length(coef) array. The first row contains a nonzero coefficient of the polynomial and the corresponding entry of the second row contains the grlex (in sum(var.ncomp) variables) rank.
        dim = []; % Number of variables of underlying space.
        ord = []; % Order of sequence.
        choose = []; % Choose table to speed up lookups.
    end

    methods 
        function obj = seq(varargin) 
            if isempty(varargin)
                disp('Error: you must specify the dimension of the underlying space');
                return
            end
            
            if nargin >= 1
                obj.dim = varargin{1};
            end
            
            if nargin >= 2
                obj.ord = varargin{2};
                obj.choose = ncktab(obj.dim+varargin{2});
            end
            
            if nargin == 3  % Convert vector of doubles into sequence object.
                if numel(varargin{3}) ~= nchoosek(obj.dim+obj.ord,obj.ord)
                    disp('Error the length of the sequence does not match up with the specified dimension and order.')
                    return
                end
                obj.coef = [varargin{3}(:)';1:nchoosek(obj.dim+obj.ord,obj.ord)];
            end
            
        end
    end
end

