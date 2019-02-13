function markov_chain = pypgp(hyper, x, y, MCMC_N, burn_in_ratio, thinning_ratio, temperature_ladder, verbose)
% PYPGP ( hyper, x, y, verbose, MCMC_N, burn_in_ratio, thinning_ratio )
% hyper - parameters
% x - equidistant on [0,1] - T values
% y - time series - N fonctions of T values
% verbose - to print and draw each sample of the MCMC
% MCMC_N - number of MCMC samples

	if nargin < 4
        MCMC_N = 1000;
    end
    if nargin < 5
        burn_in_ratio = 0.5;
    end
    if nargin < 6
        thinning_ratio = 1.0;
    end
    if nargin < 7
        temperature_ladder = 1;
    end
    if nargin < 8
        verbose = 0;
    end
    
    N = size(y,2);
    
    assert(thinning_ratio <= 1.0);
    thinning_schedule = round(1:(1/thinning_ratio):MCMC_N);
    
    assert(burn_in_ratio < 1.0);
    burnin_schedule = 1:round(burn_in_ratio*MCMC_N);

    E = length(temperature_ladder);
    chain = cell(1,E);
    
    for e = 1:E
        % Initial assignment of the functions to their component: 
        %display('Drawing initial state from G0');
        chain{e}.Z = CRP(N,hyper.alpha);
        %Z = ones(N,1);
        chain{e}.sub_y = cell(1,size(chain{e}.Z,2));
        for k = 1:size(chain{e}.Z,2)
            chain{e}.sub_y{k} = y(:, chain{e}.Z(:,k)>0);
        end
        C = size(chain{e}.Z,2);
        chain{e}.component = cell(1,C);
        for c = 1:C
            chain{e}.component{c} = sample_G0( hyper, x, chain{e}.sub_y{c});
            chain{e}.component{c}.label = c;
            chain{e}.component{c} = gradient(chain{e}.component{c}, 1000, 0);
        end
    end
    

    markov_chain = cell(MCMC_N,1);
    for sweep = 1:MCMC_N
        for e = 1:length(temperature_ladder)
            chain{e}.component = resample_components( chain{e}.component, chain{e}.Z, x, y, 0);
            [chain{e}.Z, chain{e}.component] = resample_latent_clustering( chain{e}.Z, chain{e}.component, hyper, x, y );
            hyper = update_alpha(chain{e}.Z , chain{e}.hyper);
            %hyper = update_d(Z , hyper);
        end

                    
            
        sample.pi = sum(Z,1);
        sample.Z = Z;
        sample.component = component;
        markov_chain{sweep} = sample;
        drawSweep(x,y,Z,sweep,MCMC_N);
        sample = struct;
    end

end

% for verbose:
function drawSweep(x,y,Z,sweep,MCMC_N)
        fprintf('---------------------------------------------------\n');
        fprintf('Sweep %i/%i \n',sweep,MCMC_N);
        fprintf('found : ');
        fprintf('%i ',sort(sum(Z,1),'descend'));
        fprintf('\n');
        
        global Z_truth;
        fprintf('orig  : ');
        fprintf('%i ',sort(sum(Z_truth,1),'descend'));
        
        fprintf('\n \n');
        
        figure(1);
        clf;
        hold on;
        grid on;
        for k = 1:size(Z,2)
            plot(x, y(:,find(Z(:,k))), 'Color',unifrnd(0.1,1.0, 1,3));
        end
        drawnow;
end
