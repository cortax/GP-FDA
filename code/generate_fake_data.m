function [ y, component, Z ] = generate_fake_data( N, x, hyper )
    y = [];
    component = {};
    Z = CRP(N,hyper.alpha);
    for k = 1:size(Z,2)  
        component{k} = sample_G0( hyper, x, zeros(length(x),0) );
        component{k}.Ky = nsgausskernel(x, x, component{k}.loglambda, component{k}.loglambda, component{k}.loggamma, component{k}.loggamma, component{k}.logeta);
        y = [y, mvnrnd(component{k}.m, component{k}.Ky, sum(Z(:,k)))'];
    end
end

