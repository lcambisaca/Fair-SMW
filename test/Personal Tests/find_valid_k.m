function results = find_valid_k(n, h)
    % Finds (k,h) pairs such that:
    % 1. n / (k * h) is an integer
    % 2. k > 1, h > 1
    % 3. n - h + 1 >= k
    % 4. |k - h| <= max_diff (optional balance filter)
    results = [];
    
    for k = 2:10
        if mod(n, k * h) ~= 0
            continue;
        end
        if n - h + 1 < k
            continue;
        end
        results = [results; k];
        
    end
end



