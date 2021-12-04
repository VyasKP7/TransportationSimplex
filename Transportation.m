%Load all variables


clear variables;


prompt = '\nPlease provide name of input file for a transportation problem: ';
x = input(prompt, 's');

start = tic;
fid = fopen(x);



if fid < 0
    fprintf('\nInput file cannot be opened...exiting\n');
    return
end

oname = strcat('Output', x);
fileID = fopen(oname,'w');

tline = fgetl(fid); %Skip num_sources line
tline = fgetl(fid);
%Extract the m, n from TopLine
names = split(tline);
dim = cellfun(@str2num,names);

m = dim(1);  %originals
n = dim(2);
    
tline = fgetl(fid); %Skip
    
%populate a
tline = fgetl(fid);
names = split(tline);
a = cellfun(@str2double,names);

tline = fgetl(fid); %Skip

%populate b
tline = fgetl(fid);
names = split(tline);
b = cellfun(@str2double,names);


tline = fgetl(fid); %Skip

c = zeros(m, n);

%populate C
for i=1:m
   tline = fgetl(fid);
   namei = split(tline);
   constrainti = cellfun(@str2double,namei);
   for j=1:n
       c(i, j) = constrainti(j);
   end
end
fclose(fid); %Close Input


fprintf(fileID,'The transportation problem is:\n');
fprintf(fileID, "a = \n\t");
transposea = a.';
transposeb = b.';
for ii = 1:size(transposea,1)
    fprintf(fileID,'%g\t',transposea(ii,:));
    fprintf(fileID,'\n\t');
end
fprintf(fileID, "\nb = \n\t");
for ii = 1:size(transposeb,1)
    fprintf(fileID,'%g\t',transposeb(ii,:));
    fprintf(fileID,'\n\t');
end

fprintf(fileID, "\nC = \n\t");
for ii = 1:size(c,1)
    fprintf(fileID,'%g\t',c(ii,:));
    fprintf(fileID,'\n\t');
end
clear transposea tranposeb;


%NWCR
remA = a;  %y
remB = b;  %x
bfs = zeros(m,n);
rowsdone = zeros(m);
columnsdone = zeros(n);
basicvariables = zeros(m, n);

for i=1:m
    for j=1:n
        if rowsdone(i) == 1 || columnsdone(j) == 1
            continue
        end
        if remA(i) <= remB(j)
            bfs(i,j) = remA(i);
            basicvariables(i,j) = 1;
            rowsdone(i) = 1;
        else
            bfs(i,j) = remB(j);
            basicvariables(i,j) = 1;
            columnsdone(j) = 1;
        end
        remA(i) = remA(i) - bfs(i,j);
        remB(j) = remB(j) - bfs(i,j);
    end
end
fprintf('After NWCR the Initial bfs is:\n')
bfs
z = sum(sum(bfs.*c));
fprintf('Initial Z = %f\n\n', z);
fprintf(fileID, "\nInitial BFS after NWCR = \n\t");
for ii = 1:size(bfs,1)
    fprintf(fileID,'%g\t',bfs(ii,:));
    fprintf(fileID,'\n\t');
end
fprintf(fileID, '\n');
fprintf(fileID, 'Initial Z = %f\n\n', z);


subC = c;
iteration = 0;

while 1
    iteration = iteration + 1;
    
    fprintf('Simplex iteration %d\n', iteration);
    fprintf(fileID, 'Simplex iteration %d\n', iteration);
    
    %Obtain U and V
    u = zeros(m, 1);
    v = zeros(n, 1);
    
    uFilled = zeros(m,1);
    vFilled = zeros(n,1);
    vFilled(n) = 1;
    while 1
        
        bt = 0;
        for i=1:m
            for j=1:n
                if basicvariables(i,j) == 1
                    if subC(i,j) ~= u(i) +v(j)
                        bt = 1;
                        break;
                    end
                end
            end
            if bt == 1
                break
            end
        end
        
        if bt == 0
            break
        end
        
        for i=1:m
            for j=1:n
                if basicvariables(m-i+1,n-j+1) == 1
                    if uFilled(m-i+1) == 0 && vFilled(n-j+1) == 1
                        u(m-i+1) = subC(m-i+1, n-j+1) - v(n-j+1);
                        uFilled(m-i+1) = 1;
                    end
                    if vFilled(n-j+1) == 0 && uFilled(m-i+1) == 1
                     v(n-j+1) = subC(m-i+1, n-j+1) - u(m-i+1);
                     vFilled(n-j+1) = 1;
                    end
                end
            end
        end
                
    end
    

    r = zeros(m,n);
    %Obtain r for non-basic
    for i=1:m
        for j=1:n
            if basicvariables(i,j) == 0
                r(i,j) = subC(i,j) - (u(i) + v(j));
            end
        end
    end
    
 
    bt = 0;
    smallest = 0;
    for i=1:m
        for j=1:n
            if r(i,j) < 0
                iB = i;
                jB = j;
                smallest = r(i,j);
                bt = 1;
                break
            end

        end
        if bt ==1
            break
        end
        
    end
    

    if smallest>=0   %Optimal
        fprintf('Optimal Solution Reached\n\n');
        fprintf(fileID, 'Optimal Solution Reached\n\n');
        break
    end
    
    fprintf('X%d%d will enter basis\n',iB,jB);
    fprintf(fileID, 'X%d%d will enter basis\n',iB,jB);

    %Determine what leaves basis
    newbasic = zeros(m,n);
    newbasic(iB,jB) = 1; %This is (+)
    
    immutables = zeros(m, n);
    numim = 0;
    
    updatedbasic = basicvariables;
    updatedbasic(iB,jB) = 1;
    for i=1:m
        for j=1:n
            if basicvariables(i,j) == 0
                immutables(i,j) = 1;
                continue
            end
            lol = sum(updatedbasic);     %c
            lolt = sum(transpose(updatedbasic)); %r
            if (lol(j) == 1 || lolt(i) == 1)
                immutables(i,j) = 1;
                continue
            end
            numim = numim + 1; 
        end
    end

    
    permutation = zeros(numim, 1);
    mutablei = zeros(numim, 1);
    mutablej = zeros(numim, 1);
    count = 1;
    %we have 3^numim permutations to try
    for i=1:m
        for j=1:n
            if immutables(i,j) == 0
                mutablei(count) = i;
                mutablej(count) = j;
                count = count + 1;
            end
        end
    end

    % 0 is (0), 1 is (-) and 2 is (+)  
    
    it = 0;
    while 1 %Loop to fill in the basic variables
        continueV = 0;
        columnSums = sum(newbasic);
        rowSums = sum(transpose(newbasic));
        for i=1:m
            if rowSums(i) ~= 0
                continueV = 1;
                break
            end
        end
        
        for j=1:n
            if columnSums(j) ~=0
                continueV = 1;
                break;
            end
        end
        
        if continueV == 0
            break
        end
        
        it = it + 1;
        permutation = perms(it, numim);
        
        for i=1:numim
            if permutation(i) == 0
                newbasic(mutablei(i), mutablej(i)) = 0;
            elseif permutation(i) == 1
                newbasic(mutablei(i), mutablej(i)) = -1;
            else
                newbasic(mutablei(i), mutablej(i)) = 1;
            end
        end
    end
    
    numneg = 0;
    for i=1:m
        for j=1:n
            if newbasic(i,j) == -1
                numneg = numneg + 1;
            end
        end
    end
    thetavector = zeros(numneg, 1);
    counter = 0;
    for i=1:m
        for j=1:n
            if newbasic(i,j) == -1
                thetavector(counter + 1, 1) = bfs(i,j);
                counter = counter + 1;
            end
        end
    end
    theta = min(thetavector);
    
    for i=1:m
        for j=1:n
            if newbasic(i,j) == 1 %add
                bfs(i,j) = bfs(i,j) + theta;
            elseif newbasic(i,j) == -1 %subtract
                bfs(i,j) = bfs(i,j) - theta;
            else
                continue
            end
        end
    end
    
    bt = 0;
    for i=1:m
        for j=1:n
            if basicvariables(i,j) == 1
                if bfs(i,j) == 0
                    basicvariables(i,j) = 0;
                    fprintf('X%d%d will leave basis\n\n',i,j);
                    fprintf(fileID, 'X%d%d will leave basis\n\n',i,j);
                    bt = 1;
                    break
                end
            end
        end
        if bt == 1
            break
        end
    end
    
    basicvariables(iB, jB) = 1;
    
    
    fprintf('Current bfs is:\n')
    bfs
    z = sum(sum(bfs.*c));
    fprintf('Current Z = %f\n\n', z);
    fprintf(fileID, "\nCurrent bfs = \n\t");
    for ii = 1:size(bfs,1)
        fprintf(fileID,'%g\t',bfs(ii,:));
        fprintf(fileID,'\n\t');
    end
    fprintf(fileID, '\n');
    fprintf(fileID, 'Current Z = %f\n\n', z);
    
    
end

fprintf('Optimal Solution:\n');
optimalsolution = bfs
fprintf(fileID, 'Optimal Solution:\n');
for ii = 1:size(bfs,1)
    fprintf(fileID,'\t%g',bfs(ii,:));
    fprintf(fileID,'\n');
end

fprintf('\nConstraint Check:\n');
fprintf(fileID, '\nConstraint Check:\n');
noneless = 0;
for i=1:m
    for j=1:n
        if bfs(i,j) < 0
            fprintf('x%d%d < 0, violated\n', i, j);
            fprintf(fileID, 'x%d%d < 0, violated\n', i, j);
            noneless = 1;
            break
        end
    end
    if noneless == 1
        break
    end
end

if noneless == 0
    fprintf('All Xij >= 0\n');
    fprintf(fileID, 'All Xij >= 0\n');
end

allequal = 0;
for i=1:m   %check all a
    fprintf('source %d = %f, sum of bfs in row %d = %f', i, a(i), i, sum(bfs(i, :)));  
    fprintf(fileID, 'source %d = %f, sum of bfs in row %d = %f', i, a(i), i, sum(bfs(i, :)));  
    if a(i) ~= sum(bfs(i, :))
        fprintf(', This is violated.\n');
        fprintf(fileID, ', This is violated.\n');
        allequal = 1;
    else
        fprintf('\n');
        fprintf(fileID, '\n');
    end
end

for i=1:n   %check all a
    fprintf('destination %d = %f, sum of bfs in column %d = %f', i, b(i), i, sum(bfs(:, i)));  
    fprintf(fileID, 'destination %d = %f, sum of bfs in column %d = %f', i, b(i), i, sum(bfs(:, i)));  
    if b(i) ~= sum(bfs(:, i))
        fprintf(', This is violated.\n');
        fprintf(fileID, ', This is violated.\n');
        allequal = 1;
    else
        fprintf('\n');
        fprintf(fileID, '\n');
    end
    
end

if allequal == 1 || noneless == 1
    fprintf('\nThe solution is infeasible.\n');
    fprintf(fileID, '\nThe solution is infeasible.\n');
else
    fprintf('\nThe solution is feasible.\n');
    fprintf(fileID, '\nThe solution is feasible.\n');
end

    


fprintf('Optimal Cost:\n');
z = sum(sum(bfs.*c))
fprintf(fileID,'\nOptimal Cost= %f\n', z);
Elapsed_time = toc(start);
fprintf(fileID, "\n\nTime Taken: %f seconds\n", Elapsed_time);
fprintf("\n\nTime Taken: %f seconds\n", Elapsed_time);
fclose(fileID);
clear variables;