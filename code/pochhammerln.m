function v = lnpochhammer( x , n , a )
    if nargin < 3; a = 1; end;
    v = 0;
    for i = 0:(n-1)
        v = v + log(x + i*a);
    end
end

