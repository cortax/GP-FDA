function gpprior_m = make_gpprior_m(x_timegrid)
    T = length(x_timegrid);
    m = zeros(1,T); 
    kernel = gausskernel(x_timegrid, log(1.0), log(0.05));
    gpprior_m = gpmodel(x_timegrid, m, kernel);
end

