function fractions = computeFraction(labels, G)
    unilb = unique(labels);
    k = length(unilb);
    h = size(G,2);
    fractions = zeros(k,h);
    for i = 1:k
        idx = labels == unilb(i);
        ci = sum(idx);
        fprintf("----cluster #: %d\n----", i);
        for s = 1:h
            fprintf("--group #: %d\n--", s);
            g = G(:,s);
            vfandci = sum(g(idx));
            fprintf("--vfandci #: %d\n--", vfandci);
            fractions(i,s) = vfandci/ci;
        end
    end
end

