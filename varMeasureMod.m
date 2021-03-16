function  G = varMeasureMod(Y)
% G = varMeasureMod(Y)
%
% Returns the sum of heights of the peaks of a curve, without the T90 compensation.
% The height is measured by adding the absolute value of the differences between the
% ordinates when it is negative.
% Andor Budai (2020); Eötvös University, Institute of Physics, 1117 Budapest, Hungary; email: arandras@caesar.elte.hu
% 
% Input:
% Y [W] - a vector containing the ordinates of the curve
%
% Output:
% G [1/s] - the variability of the curve

% Pseudocode:
% 1. Checking whether the vector is empty or a null vector.
% 2. Joining zeros to the front end the end of the Y vector, and normalization.
% 3. Calculating the difference between the elements of Y.
% 4. Calculating the sum of the heights of the peaks.
% 5. Returning the variability.


%%
% 1. Checking whether the vector is empty or a null vector.
if(max(Y) == 0 | isempty(Y))
    G = 0;
    return;
end

%%
% 2. Joining zeros to the front end the end of the Y vector, and normalization.
Y = Y/max(Y);
Y = Y;


% 3. Calculating the difference between the elements of Y.
d = Y(2:end) - Y(1:end-1);


% 4. Calculating the sum of the heights of the peaks.
G = abs(sum(d(d<0)));


% 5. Returning the variability.
G = (G - 1); % This makes G = 0, if the number of peaks is one.



%% NUMBER OF PEAKS---------------------------------------------------------
%%% Uncomment, if you want to calculate the number of peaks.
%--------------------------------------------------------------------------
%%% Every time d turns from positive to negative, the number of peaks grows
%%% by one.
%
%%% Finding the peaks of Y.
% p = (d(1:end-1) >= 0 & d(2:end) < 0);
%
%%% Calculating the number of peaks.
% cs = sum(p);
%
%--------------------------------------------------------------------------

end % end of function

% Andor Budai (2020) - arandras@caesar.elte.hu