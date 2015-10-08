%% BBC_MAIN
% This is the entry point script for the Bayesian reanalysis of the
% Reproducibility Project: Psychology
clear all

%% Get the data
% Data from 60 studies whose results could be translated to a t-value with
% associated degrees of freedom

rppdata = [ ...
% ----------------------------------------------------------------------- %
%                  original                         replica               %
%  study    ----------------------------    ----------------------------- %
%  number    df    t-value     eff size       df    t-value     eff size  %
% --------  ----------------------------    ----------------------------- %
     1   ,   13    2.666458    1.479085  ,    28    0.793725    0.300000    
     2   ,   23    3.702702    1.544133  ,    23    1.131371    0.471814    
     3   ,   24    2.300000    0.938971  ,    31    1.227192    0.440820    
     4   ,  190    3.238827    0.469938  ,   268    0.100000    0.012217    
     5   ,   31    2.894823    0.519926  ,    47    0.932738    0.136054    
     6   ,   23    3.550000    0.740226  ,    31    2.400000    0.431053    
     7   ,   99   10.180000    1.023128  ,    14    0.496000    0.132562    
     8   ,   37    4.126742    1.356864  ,    31    0.619677    0.222595    
    10   ,   28    5.166237    1.952654  ,    29    6.728298    2.498827    
    11   ,   21    4.159327    1.815279  ,    29    2.839718    1.054645    
    15   ,   94    1.929000    0.198961  ,   241    3.955000    0.254764    
    19   ,   31    3.768289    1.353609  ,    19    1.913374    0.877916    
    20   ,   94    2.229350    0.459880  ,   106    0.200000    0.038851    
    24   ,  152    4.814146    0.390479  ,    48    2.054264    0.296507    
    26   ,   94    1.581139    0.326164  ,    92    1.396424    0.291175    
    27   ,   31    2.273763    0.816760  ,    70    3.432637    0.820557    
    28   ,   31    2.024846    0.727346  ,    90    0.984886    0.207632    
    33   ,   39    3.770000    0.603683  ,    39    2.080000    0.333067    
    36   ,   20    4.559605    2.039117  ,    20    4.165333    1.862794    
    37   ,   11    2.190890    1.321157  ,    17    1.539480    0.746758    
    44   ,   67    3.080000    0.752564  ,   176    2.016000    0.303923    
    49   ,   34    2.383275    0.817457  ,    86    0.282843    0.060999    
    52   ,  131    2.406242    0.420469  ,   111    0.994987    0.188880    
    53   ,   31    2.267157    0.814387  ,    73    0.657267    0.153855    
    56   ,   99    4.076763    0.819460  ,    38   -0.260000   -0.084355    
    58   ,  182    2.289105    0.339359  ,   278    0.613188    0.073553    
    61   ,  108   -2.340000   -0.450333  ,   220    0.070000    0.009439    
    63   ,   68    2.349468    0.569830  ,   145    0.891067    0.147998    
    65   ,   41    3.065942    0.957639  ,   131    0.134164    0.023444    
    68   ,  116    2.037155    0.378290  ,   222    0.044721    0.006003    
    71   ,  373    4.400000    0.455647  ,   175    0.973000    0.147104    
    72   ,  257    3.402940    0.424539  ,   247    0.700000    0.089080    
    81   ,   90    2.641969    0.556976  ,   137    1.195826    0.204333    
    87   ,   51    3.075711    0.861371  ,    47    0.089443    0.026093    
    94   ,   26    1.870000    0.733474  ,    59    2.325000    0.605378    
    97   ,   73    3.491418    0.817279  ,  1486    1.424781    0.073921    
   106   ,   34    2.408319    0.826047  ,    45    1.534000    0.457350    
   107   ,   84    2.090000    0.456075  ,   156    1.318000    0.211049    
   110   ,  278   11.107655    1.332386  ,   142    1.090871    0.183088    
   111   ,   55    2.622975    0.353682  ,   116    2.495997    0.231747    
   112   ,    9    2.949576    1.966384  ,     9    3.405877    2.270585    
   113   ,  124   10.360000    0.930355  ,   175    15.64000    1.182273    
   114   ,   30    3.806573    1.389964  ,    30    4.719110    1.723175    
   115   ,   31    3.230000    0.580125  ,     7   -1.400000    0.287687    
   118   ,  111    2.304561    0.437478  ,   158    0.615630    0.097954    
   122   ,    7    2.760000    1.043182  ,    16   -9.590000   -2.397500    
   124   ,   34    2.426932    0.832431  ,    68    0.282843    0.068599    
   127   ,   28    4.980000    0.941132  ,    25   -3.103000   -0.620600    
   129   ,   26    2.042058    0.800961  ,    64    0.141421    0.035355    
   133   ,   23    2.387467    0.995643  ,    37    2.842534    0.934619    
   135   ,  562   -0.110000   -0.004640  ,  3511   -6.310000   -0.212979    
   145   ,   76   10.475686    2.403287  ,    36    5.173007    1.724336    
   146   ,   14    3.200000    1.710472  ,    11    1.900000    0.572872    
   148   ,  194    2.675818    0.192113  ,   259    0.485798    0.030186    
   150   ,   13    3.768289    2.090271  ,    18    0.900000    0.212132    
   151   ,   41    2.794638    0.872898  ,   124    0.031623    0.005680    
   153   ,    7    4.450000    1.681942  ,     7    0.320000    0.120949    
   158   ,   38    2.491987    0.808507  ,    93    4.352011    0.902565    
   161   ,   44    3.663332    1.104536  ,    44    1.198749    0.361437    
   167   ,   17    3.054505    0.740826  ,    21    1.204159    0.525538   ...
% ----------------------------------------------------------------------- %
];

n = size(rppdata, 1);  % number of studies

dfo = rppdata(:,2);    % degrees of freedom in original
to  = rppdata(:,3);    % t statistic of original

dfr = rppdata(:,5);    % degrees of freedom in replication
tr  = rppdata(:,6);    % t statistic of replication

%% Compute the Bayes factors
% This is the part that takes some time

B_ori = zeros(size(to));  % pre-allocate memory
B_mit = zeros(size(to));
B_rep = zeros(size(to));

% Loop over studies
parfor c = 1:n
    fprintf('Starting study number %i...\n', c)
    B_ori(c) = bbc_t(to(c), dfo(c), [1 1 0 0 0 0 0 0] / 2);  % face-value
    B_mit(c) = bbc_t(to(c), dfo(c), [1 1 1 1 1 1 1 1] / 8);  % mitigation
    B_rep(c) = bbc_t(tr(c), dfr(c), [1 1 0 0 0 0 0 0] / 2);  % face-value
end

save rpp


%% Output some results tables
% "Critical" Bayes factors
B_crit = log([0 1/10 1/3 3 10 Inf]);

load rpp

original   = histc(log(B_ori), B_crit)';  % frequencies for original
mitigated  = histc(log(B_mit), B_crit)';  % frequencies for mitigated
replicate  = histc(log(B_rep), B_crit)';  % frequencies for replicate
transition = hist3(log([B_mit, B_rep]), 'edges', { B_crit B_crit })';

    % Strip that weird final bin that hist* adds...
    original(end)     = [];
    mitigated(end)    = [];
    replicate(end)    = [];
    transition(:,end) = [];
    transition(end,:) = [];

% Print table
fprintf ' --------------------------------------------\n'
fprintf '  Bayes factor frequencies\n'
fprintf '                <1/10 <1/3    -    >3   >10\n'
fprintf '  original:  ', disp(original)
fprintf '  mitigated: ', disp(mitigated)
fprintf '  replicate: ', disp(replicate)
fprintf ' --------------------------------------------\n'
fprintf '  Transition frequencies\n'
fprintf '  Rows for mitigated, columns for replicate\n'
disp(transition)
fprintf ' --------------------------------------------\n'

%% Make plot

load rpp

mark = rppdata(:,end);

nratio = round(100*(rppdata(:,5) ./ rppdata(:,2)));     % sample size ratio
large = max(abs(log10(B_mit)), abs(log10(B_rep))) > 1;  % find large BFs

% Plot small-BF studies
plot(log10(B_mit(~large)), log10(B_rep(~large)), 'kx', ...
    'MarkerSize', 4, 'LineWidth', 1)
hold on

% Plot large-BF studies
x = log10(B_mit(large));
y = log10(B_rep(large));
s = nratio(large);
for ctr = 1:sum(large)
    scatter(x(ctr), y(ctr), s(ctr), ...
        'MarkerEdgeColor', 'k', 'LineWidth', 2, 'MarkerFaceColor', 'w')
end

% Plot reference circle
rectangle('Position', [3.25 -1.75 3 .5], ...
    'FaceColor', 'w', 'clipping', 'off')
scatter(3.5, -1.5, 100, ...
    'MarkerEdgeColor', 'k', 'LineWidth', 2, 'MarkerFaceColor', 'w')
hold off

% Draw custom grid lines
set(gca, 'ytick', log10([1/10 1/3 3 10]))
set(gca, 'xtick', log10([1/10 1/3 3 10]))
grid on
set(gca, 'xticklabel', {'1/10' '1/3' '3' '10'})    
set(gca, 'yticklabel', {'1/10' '1/3' '3' '10'})
box on

% Zoom in a little
% (this causes a few faraway data points to be clipped!)
axis([-1.5 6.5 -1.5 6.5])

% Labels will be placed in LaTeX because they'll look nicer
% xlabel 'original (mitigated) BF pro-H_A'
% ylabel 'replication BF pro-H_A'

set(gcf, 'Position', [94   100   613   471])