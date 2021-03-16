function DATA = analyseError(N, m)
% DATA = analyseError(N, m) 
%
% Reads in a table containing the ID of the GRB light curves with the corresponding
% angle values and the T90 values, and calculates the variability of those light curves.
% It uses random light curves generated with given rate errors to calculate
% the error of the variability, and save the errors of the random curves.
% Andor Budai (2020); Eötvös University, Institute of Physics, 1117 Budapest, Hungary; email: arandras@caesar.elte.hu
%
% Input:
%  N - number of iteration (has to be a multiple of m)
%  m - number of memory reduction iterations
%
% Output:
%  DATA - table containing the name of GRBs, the variabilities and the T90 values.
%
% called function:
%  varMeasureMod


%% Variables
folder = './Data' % folder containing the light curves
table = './GRBDataFinal.csv' % path to the data table 
rng('shuffle')
n = N/m; % n is the number of random light curves in one iteration. 
         % m iterations are ran.

%% Reading in the data table
opts = detectImportOptions(table); % properties of the table
opts = setvartype(opts, {'ID'}, 'char'); % set the property of the ID column
                                         % to read the ID numbers as
                                         % characters

DATA = readtable(table, opts); % reading the table with character ID
h = height(DATA); % number of GRB angles
DATA.Var = zeros(h, 1); % adding the Variability column to the table
DATA.Var_Err_Neg = zeros(h, 1); % adding the lower Variability Error column to the table
DATA.Var_Err_Poz = zeros(h, 1); % adding the upper Variability Error column to the table
varErr = zeros(1, N); % every element will be the variability of a random light curve

%% Listing the contents of the folder containing the light curves
files = dir(strcat(folder,'/*.csv')); % listing csv files with information
                                      % (output: struct array)
files = {files.name}; % creating an array of the names (output: cell array)

%% Calculating the variabilities 
i = 0; % number of the calculated variabilities
for grb=files % grb is a cell
    grbname = grb{1}; % converting grbname to string
    grbname = char(grbname); % converting grbname to a char array
    s = strcat(folder,'/', grbname); % path to the light curve
    RATE = readtable(s); % reading the light curve
    grbname = grbname(3:8); % getting the ID from the light curve name
    rate = RATE.rate(:)'; % turning the 'rate' column to an array
    error = RATE.error(:)';% turning the 'error' column to an array 
    t90 = DATA.T90(strcmp(string(DATA.ID),grbname));
    t90 = t90(1); % t90 from literature
    
    l = length(rate); % length of the light curve
    RATE = repmat(rate, n, 1); % creating Rate matrix
                               % (every row is the light curve)
    
    ERROR = repmat(error, n, 1);% creating Error matrix
                                % (every row contains the light curve errors)
    clear rate; % liberate memory
    clear error;% liberate memory
    
    for k = 1:m % iteration is needed to reduce the memory usage
    
        L = ERROR.*randn(n, l) + RATE; % every row is a random light curve
    
        for j =1:n
            varErr((k-1)*n+j) = varMeasureMod(L(j, :))/t90;% calculating the 
                                                       % the variability of
                                                       % the random light
                                                       % curve
        end
        
        clear L; % liberate memory
    end
    
    
    medErr = round(median(varErr), 2);
    nErr = round(medErr - prctile(varErr, 16), 2);
    pErr = round(prctile(varErr, 84) - medErr, 2);
    
    
    DATA.Var(strcmp(string(DATA.ID),grbname)) = medErr; % Variability
    DATA.Var_Err_Neg(strcmp(string(DATA.ID),grbname)) = nErr;% Variability Error lower limit
    DATA.Var_Err_Poz(strcmp(string(DATA.ID),grbname)) = pErr;% Variability Error upper limit 
    
    save(strcat('VarErr/', grbname), 'varErr'); % saving the variability errors
    
    i = i+1 
end

DATA.Var = round(DATA.Var, 2);
DATA.Var_Err_Neg = round(DATA.Var_Err_Neg, 2);
DATA.Var_Err_Poz = round(DATA.Var_Err_Poz, 2);

%% creating result table
writetable(DATA, 'resultTableErr.csv') % WARNING: overwrites the table
end % end of function

% Andor Budai (2020) - arandras@caesar.elte.hu