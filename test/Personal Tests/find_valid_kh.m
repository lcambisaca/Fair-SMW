function results = find_valid_kh(n_values, max_diff)
    % Finds (k,h) pairs such that:
    % 1. n / (k * h) is an integer
    % 2. k > 1, h > 1
    % 3. n - h + 1 >= k
    % 4. |k - h| <= max_diff (optional balance filter)

    if nargin < 2
        max_diff = 10; % no balance filtering by default
    end

    results = [];

    for i = 1:length(n_values)
        n = n_values(i);

        for k = 2:9
            for h = 2:n
                if mod(n, k * h) ~= 0
                    continue;
                end
                if n - h + 1 < k
                    continue;
                end
                if abs(k - h) > max_diff
                    continue;  % too unbalanced
                end

                results = [results; k, h];
            end
        end
    end

    results = unique(results, 'rows');
end
