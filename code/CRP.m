function Z = CRP( N , alpha )
 
    ClassCounts = 0;
    Z = [];

    for i = 1:N
        %Number of occupied tables
        K = size(ClassCounts,2) - 1;

        %defining probabilities
        proba = ClassCounts(1:K)./(i-1+alpha);
        proba(K+1) = alpha/(i-1+alpha);

        %Custumer sits at a table
        sample = mnrnd(1,proba);
        Z(i,find(sample==1)) = 1;
        ClassCounts = ClassCounts + sample;

        %making additionnal table available
        if ClassCounts(end) > 0; ClassCounts(end+1) = 0; end;
    end

end
