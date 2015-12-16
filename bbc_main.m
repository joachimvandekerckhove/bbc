%% BBC_MAIN
% This is the entry point script for the Bayesian reanalysis of the
% Reproducibility Project: Psychology
clear all

%% Get the data
% Data from 72 studies whose results could be translated to a t-value with
% associated degrees of freedom

rppdata = [ ...
% -------------------------------------------------------------------------- %
%             original           replica             log10(Bayes factor)     %
%  study   ------------------ ------------------ --------------------------- %
%  number    df   t-value       df   t-value        ori      mit      rep    %
% -------- ------------------ ------------------ --------------------------- %
     1   ,   13   2.666458  ,   28   0.793725  ,   0.6441   0.1108  -0.3291
     2   ,   23   3.702702  ,   23   1.131371  ,   1.5810   0.8406  -0.1752
     3   ,   24   2.300000  ,   31   1.227192  ,   0.4905  -0.0513  -0.1799
     4   ,  190   3.238827  ,  268   0.100000  ,   1.3270   0.5640  -0.8735
     5   ,   31   2.894823  ,   47   0.932738  ,   0.9802   0.3126  -0.3767
     6   ,   23   3.550000  ,   31   2.400000  ,   1.4556   0.7298   0.5672
     7   ,   99  10.180000  ,   14   0.496000  ,  13.9830  13.1400  -0.2826
     8   ,   37   4.126742  ,   31   0.619677  ,   2.2564   1.4573  -0.3953
    10   ,   28   5.166237  ,   29   6.728298  ,   3.0554   2.2314   4.5391
    11   ,   21   4.159327  ,   29   2.839718  ,   1.8813   1.1144   0.9241
    15   ,   94   1.929000  ,  241   3.955000  ,   0.0730  -0.1539   2.3326
    19   ,   31   3.768289  ,   19   1.913374  ,   1.7949   1.0271   0.2423
    20   ,   94   2.229350  ,  106   0.200000  ,   0.3236  -0.2641  -0.7108
    24   ,  152   4.814146  ,   48   2.054264  ,   3.7866   2.9503   0.2679
    26   ,   94   1.581139  ,   92   1.396424  ,  -0.1753  -0.3290  -0.2842
    27   ,   31   2.273763  ,   70   3.432637  ,   0.4696  -0.0854   1.6253
    28   ,   31   2.024846  ,   90   0.984886  ,   0.2879   0.0638  -0.4829
    29   ,    7   2.892000  ,   14   3.708000  ,   0.5192   0.0899   1.2588
    32   ,   36   4.783304  ,   37   3.334666  ,   2.9577   2.1352   1.4245
    33   ,   39   3.770000  ,   39   2.080000  ,   1.8938   1.1147   0.3089
    36   ,   20   4.559605  ,   20   4.165333  ,   2.1323   1.3475   1.8434
    37   ,   11   2.190890  ,   17   1.539480  ,   0.3697   0.1730   0.0476
    44   ,   67   3.080000  ,  176   2.016000  ,   1.2134   0.4831   0.0398
    48   ,   92  -2.220000  ,  192  -0.725548  ,   0.3186  -0.2666  -0.7393
    49   ,   34   2.383275  ,   86   0.282843  ,   0.5528  -0.0301  -0.6593
    52   ,  131   2.406242  ,  111   0.994987  ,   0.4373  -0.1962  -0.5215
    53   ,   31   2.267157  ,   73   0.657267  ,   0.4646  -0.0891  -0.5524
    56   ,   99   4.076763  ,   38  -0.260000  ,   2.5232   1.7072  -0.4970
    58   ,  182   2.289105  ,  278   0.613188  ,   0.2790  -0.3413  -0.8382
    61   ,  108  -2.340000  ,  220   0.070000  ,   0.4038  -0.2116  -0.8509
    63   ,   68   2.349468  ,  145   0.891067  ,   0.4744  -0.1317  -0.6207
    65   ,   41   3.065942  ,  131   0.134164  ,   1.1730   0.4637  -0.7584
    68   ,  116   2.037155  ,  222   0.044721  ,   0.1246  -0.4201  -0.8525
    71   ,  373   4.400000  ,  175   0.973000  ,   2.9768   2.1537  -0.6328
    72   ,  257   3.402940  ,  247   0.700000  ,   1.5031   0.7229  -0.8005
    81   ,   90   2.641969  ,  137   1.195826  ,   0.7253   0.0539  -0.4730
    87   ,   51   3.075711  ,   47   0.089443  ,   1.2009   0.4805  -0.5511
    89   ,   26   0.720000  ,   26   0.150000  ,  -0.3374  -0.3756  -0.4331
    93   ,   83   3.050000  ,   68  -1.124000  ,   1.1752   0.4430  -0.3673
    94   ,   26   1.870000  ,   59   2.325000  ,   0.1981   0.0087   0.4679
    97   ,   73   3.491418  , 1486   1.424781  ,   1.7004   0.9244  -0.9051
   106   ,   34   2.408319  ,   45   1.534000  ,   0.5730  -0.0149  -0.0775
   107   ,   84   2.090000  ,  156   1.318000  ,   0.2209  -0.3317  -0.4358
   110   ,  278  11.107655  ,  142   1.090871  ,  20.9560  20.1120  -0.5321
   111   ,   55   2.622975  ,  116   2.495997  ,   0.7462   0.0937   0.5443
   112   ,    9   2.949576  ,    9   3.405877  ,   0.6473   0.1493   0.8169
   113   ,  124  10.360000  ,  175  15.640000  ,  15.4490  14.6050  31.4220
   114   ,   30   3.806573  ,   30   4.719110  ,   1.8159   1.0472   2.7092
   115   ,   31   3.230000  ,    8  -1.426000  ,   1.2825   0.5684   0.0489
   116   ,  172   3.940000  ,  139   4.020000  ,   2.3526   1.5385   2.4698
   118   ,  111   2.304561  ,  158   0.615630  ,   0.3665  -0.2416  -0.7251
   120   ,   29   2.212316  ,   41   1.653280  ,   0.4258  -0.1123   0.0091
   122   ,    7   2.760000  ,   16  -9.590000  ,   0.4803   0.0636   3.8735
   124   ,   34   2.426932  ,   68   0.282843  ,   0.5880  -0.0035  -0.6110
   127   ,   28   4.980000  ,   25  -3.103000  ,   2.8817   2.0618   1.1170
   129   ,   26   2.042058  ,   64   0.141421  ,   0.3105   0.0892  -0.6110
   133   ,   23   2.387467  ,   37   2.842534  ,   0.5513  -0.0038   0.9505
   134   ,  115   2.303000  ,  234   8.836000  ,   0.3596  -0.2489  13.5090
   135   ,  562  -0.110000  , 3511  -6.310000  ,  -0.9042  -0.9044   5.8207
   136   ,   28   3.040000  ,   56  -0.770000  ,   1.0900   0.4085  -0.4664
   145   ,   76  10.475686  ,   36   5.173007  ,  13.1250  12.2820   3.3895
   146   ,   14   3.200000  ,   11   1.900000  ,   0.9709   0.3512   0.2339
   148   ,  194   2.675818  ,  259   0.485798  ,   0.6592  -0.0370  -0.8442
   149   ,  194   2.675818  ,  314   0.324037  ,   0.6592  -0.0370  -0.8799
   150   ,   13   3.768289  ,   18   0.900000  ,   1.2348   0.5677  -0.2222
   151   ,   41   2.794638  ,  124   0.031623  ,   0.9115   0.2428  -0.7509
   153   ,    7   4.450000  ,    7   0.320000  ,   0.8838   0.3464  -0.2037
   154   ,   68   3.927512  ,   14   0.414095  ,   2.2479   1.4437  -0.2958
   155   ,   51   2.328555  ,   70   0.284629  ,   0.4848  -0.1069  -0.6167
   158   ,   38   2.491987  ,   93   4.352011  ,   0.6405   0.0299   2.9206
   161   ,   44   3.663332  ,   44   1.198749  ,   1.8164   1.0407  -0.2521
   167   ,   17   3.054505  ,   21   1.204159  ,   0.9613   0.3306  -0.1315 ...
% -------------------------------------------------------------------------- %
];

n = size(rppdata, 1);  % number of studies

dfo = rppdata(:,2);    % degrees of freedom in original
to  = rppdata(:,3);    % t statistic of original

dfr = rppdata(:,4);    % degrees of freedom in replication
tr  = rppdata(:,5);    % t statistic of replication

%% Compute the Bayes factors
% This is the part that takes some time

B_ori = zeros(size(to));  % pre-allocate memory
B_mit = zeros(size(to));
B_rep = zeros(size(to));

% Loop over studies
parfor c = 1:n
    fprintf('Starting study number %02i ...\n', c)
    B_ori(c) = bbc_t(to(c), dfo(c), [1 1 0 0 0 0 0 0] / 2);  % face-value
    B_mit(c) = bbc_t(to(c), dfo(c), [1 1 1 1 1 1 1 1] / 8);  % mitigation
    B_rep(c) = bbc_t(tr(c), dfr(c), [1 1 0 0 0 0 0 0] / 2);  % face-value
end

save rpp

%% Output some results tables
% "Critical" Bayes factors
B_crit = log10([0 1/10 1/3 3 10 Inf]);

load rpp

B_ori = log10(B_ori);
B_mit = log10(B_mit);
B_rep = log10(B_rep);

original   = histc(B_ori, B_crit)';  % frequencies for original
mitigated  = histc(B_mit, B_crit)';  % frequencies for mitigated
replicate  = histc(B_rep, B_crit)';  % frequencies for replicate
transition = hist3([B_mit, B_rep], 'edges', { B_crit B_crit })';

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
nratio = round(100*(dfr ./ dfo));     % sample size ratio
large = max(abs(B_mit), abs(B_rep)) > 1;  % find large BFs

% Plot small-BF studies
plot(B_mit(~large), B_rep(~large), 'kx', 'MarkerSize', 4, 'LineWidth', 1)
hold on

% Plot large-BF studies
x = B_mit(large);
y = B_rep(large);
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
set(gca, 'ytick', B_crit(2:end-1))
set(gca, 'xtick', B_crit(2:end-1))
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

% savepic bbc_fig eps