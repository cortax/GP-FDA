preprocessing;

phoneme_label = 5;

k_fperiodic = make_fperiodickernel(x_timegrid);            
k_fnoise = make_fnoisekernel(x_timegrid);
kernels = {k_fperiodic, k_fnoise};

