classdef seq
    
    %SEQ Summary of this function goes here
    %Detailed explanation goes here

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
            elseif nargin == 2
                obj.ord = varargin{2};
                obj.choose = ncktab(varargin{1}+varargin{2});
            end
            obj.dim = varargin{1};
        end
    end
end

