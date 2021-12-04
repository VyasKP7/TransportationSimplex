function X = sudoku(num, n)
    
    X = zeros(n, 1);
    
    iterations = 0;
    while num > 0
        q = floor(num/3);
        r = rem(num,3);
        X(n - iterations,1) = r;
        num = q;
        iterations = iterations + 1;
    end
    X = X';
    
end % sudoku