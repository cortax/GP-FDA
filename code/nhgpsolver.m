classdef nhgpsolver < matlab.mixin.Copyable
    
    properties
    end
    
    methods
        function obj = nhgpsolver()
            
        end
        
        function problem = E_step(~,problem, data)
            problem.logexpectations = log(problem.mixture.proportions) + cell2mat(arrayfun(@(gp) gp.logpdf(data{:,:}),...
                problem.mixture.gp_component,'UniformOutput',false));
            problem.logexpectations = problem.logexpectations - logsumexp(problem.logexpectations,2);

            problem.expectations = exp(problem.logexpectations);
        end
        
        function problem = M_step(~, problem, data)
            problem.mixture.proportions = mean(problem.expectations,1);
            
        end
    end
end



