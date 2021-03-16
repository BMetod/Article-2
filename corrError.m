function corrTable = corrError(N, m)
% corrTable = corrError(N)
%
% Reads in a table containing the ID of the GRB light curves with the corresponding
% angle values and errors, and calculates Spearman's rho and p with lower and
% upper limits between the variabilities and angles. It uses the variabilities 
% of the random curves to calculate the errors. The values are saved for
% the 16 GRBs with fixed angles, and for the 8 possible samples containing
% 16 + 3 GRBs.
% Andor Budai (2020); Eötvös University, Institute of Physics, 1117 Budapest, Hungary; email: arandras@caesar.elte.hu
%
% Input:
%  N - number of iteration (has to be a multiple of m)
%  m - number of memory reduction iterations


%% Variables
VarErr_folder = './VarErr' % folder containing the variabilities
table = './resultTableErr.csv' % path to the data table
rng('shuffle')
n = N/m; % n is the number of random light curves in one iteration. 
         % m iterations are ran.

%% Reading the data table
opts = detectImportOptions(table); % properties of the table
opts = setvartype(opts, {'ID'}, 'char'); % set the property of the ID column
                                         % to read the ID numbers as
                                         % characters

DATA = readtable(table, opts); % reading the table with character ID
h = 19; % max number of GRBs

angle19 = DATA.Theta(:); % angles of the 19 GRBs 
ang_error19 = DATA.Theta_Err(:); % angle errors of the 19 GRBs 
ID19 = DATA.ID(:); % ID of the 19 GRBs
h19 = length(angle19); % height of the table (22) 

GRBrows16 = all(DATA{:,9:16}'); % the index of the rows containing the same angle in every permutation
angle16 = DATA.Theta(GRBrows16); % angles of the 16 GRBs
ang_error16 = DATA.Theta_Err(GRBrows16); % angle errors of the 19 GRBs 
ID16 = DATA.ID(GRBrows16); % ID of the 19 GRBs
h16 = length(angle16); % 16

rv = zeros(N, 1);
pv = zeros(N, 1);



%% Listing the contents of the folder containing the variabilities
files = dir(strcat(VarErr_folder,'/*.mat')); % listing mat files with information
                                      % (output: struct array)
files = {files.name}; % creating an array of the names (output: cell array)

%% 16 GRBs
VAR16 = zeros(h16, N); % variability matrix -- all 16 rows correspond to a GRB 

for j = 1:length(ID16) % looping over the 16 GRBs 
    id = ID16{j}; % ID number
    s = strcat(VarErr_folder,'/',id,'.mat');
    load(s); % opening the corresponding mat file 
    VAR16(j, :) =  varErr; % every column contains 10000 variabilities 
end
    

for i=1:m % iteration is needed to reduce the memory usage
    
    Angle16 = repmat(angle16, 1, n); % angle matrix -- every column contains the 16 angles
    Ang_Error16 = repmat(ang_error16, 1, n); % angle error matrix -- every column contains the 16 error values

    AngleRand16 = Ang_Error16.*randn(h16, n) + Angle16; % random angle matrix -- every row corresponds to the same GRB
    while ~all(all(AngleRand16>=0)) % repeating until every angle is positive
        len = length(AngleRand16(AngleRand16 < 0));
        AngleRand16(AngleRand16 < 0 ) = Ang_Error16(AngleRand16 < 0).*randn(len, 1) + Angle16(AngleRand16 < 0);
    end

    clear Angle16; % liberate memory
    clear Ang_Error16; % liberate memory


    [R, P] = corr(AngleRand16, VAR16(:, (i-1)*n+1 : i*n), 'Type', 'Spearman'); % corr matrix
    rv((i-1)*n+1 : i*n) = diag(R); % Spearman's rho
    pv((i-1)*n+1 : i*n) = diag(P); % p value
    
    clear('R','P') % liberate memory 

end
   
clear VAR16 % liberate memory

mr = round(median(rv), 2); % rho - end result
mp = round(median(pv), 2); % p - end result 
pdr = round(prctile(rv, 84) - mr, 2); % rho - upper limit
ndr = round(mr - prctile(rv, 16), 2); % rho - lower limit
pdp = round(prctile(pv, 84) - mp, 2); % p - upper limit
ndp = round(mp - prctile(pv, 16), 2); % p - lower limit

% creating the table containing the values
RP16 = [mr, ndr, pdr, mp, ndp, pdp];    
RP16 = array2table(RP16, 'VariableNames', {'median_r', 'neg_dr', 'pos_dr', 'median_p', 'neg_dp', 'pos_dp'});
writetable(RP16, 'RP16Sper.csv')

clear('mr', 'mp', 'pdr', 'ndr', 'pdp', 'ndp', 'RP16') % liberate memory

%% 19 GRBs
VAR19 = zeros(h19, N); % variability matrix -- all 22 rows correspond to a GRB

mr = zeros(8, 1);
pdr = zeros(8, 1);
ndr = zeros(8, 1);
mp = zeros(8, 1);
pdp = zeros(8, 1);
ndp = zeros(8, 1);

for j = 1:length(ID19) % looping over the 19 GRBs 
    id = ID19{j}; % ID
    s = strcat(VarErr_folder,'/',id,'.mat');
    load(s); % opening the corresponding mat file
    VAR19(j, :) =  varErr; % every column contains 10000 variabilities
end


for k = 9:16 % looping over every possible sample
    
    index = logical(DATA{:, k}); % the index of the rows

    for i=1:m % iteration is needed to reduce the memory usage
        
        Angle19 = repmat(angle19(index), 1, n); % angle matrix -- every column contains the 19 angles
        Ang_Error19 = repmat(ang_error19(index), 1, n); % angle error matrix -- everz column contains the 19 error values

        AngleRand19 = Ang_Error19.*randn(h, n) + Angle19; % random angle matrix -- evry row corresponds to the same GRB
        while ~all(all(AngleRand19>=0)) % repeating until every angle is positive
            len = length(AngleRand19(AngleRand19 < 0));
            AngleRand19(AngleRand19 < 0 ) = Ang_Error19(AngleRand19 < 0).*randn(len, 1) + Angle19(AngleRand19 < 0);
        end

        clear Angle19; % liberating memory
        clear Ang_Error19; % liberating memory


        [R, P] = corr(AngleRand19, VAR19(index, (i-1)*n+1 : i*n), 'Type', 'Spearman'); % corr matrix
        rv((i-1)*n+1 : i*n) = diag(R); % Spearman's rho
        pv((i-1)*n+1 : i*n) = diag(P); % p value
        clear('R','P') % liberating memory
        
    end

    
    mr(k-8) = round(median(rv), 2); % rho - end result
    mp(k-8) = round(median(pv), 2); % p - end result
    pdr(k-8) = round(prctile(rv, 84) - mr(k-8), 2); % r - upper limit
    ndr(k-8) = round(mr(k-8) - prctile(rv, 16), 2); % r - lower limit
    pdp(k-8) = round(prctile(pv, 84) - mp(k-8), 2); % p - upper limit
    ndp(k-8) = round(mp(k-8) - prctile(pv, 16), 2); % p - lower limit
    
    k-8 % number of the iteration

end

% creating the table containing the values
RP19 = [mr, ndr, pdr, mp, ndp, pdp];
RP19 = array2table(RP19, 'VariableNames', {'median_r', 'neg_dr', 'pos_dr', 'median_p', 'neg_dp', 'pos_dp'});
writetable(RP19, 'RP19Sper.csv')

end % end of function

% Andor Budai (2020) - arandras@caesar.elte.hu
