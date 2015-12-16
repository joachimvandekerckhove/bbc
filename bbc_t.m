function bayes_factor = bbc_t(t, df, prior)
% BBC_T Compute a t-test Bayes factor, taking into account publication bias
%   This function is a prototype and is not robust against misuse. Do not
%   use this function unless you know what you are doing!

%% Some quick error checking
if ~isfinite(t) || isnan(t) || ~isscalar(t)
    error('bbc:bbc_t:bad_input_t', ...
        'First input (t) must be a finite scalar.')
end
if ~isfinite(df) || isnan(df) || ~isscalar(df)
    error('bbc:bbc_t:bad_input_df', ...
        'Second input (df) must be a finite scalar.')
end
if numel(prior)~=8
    error('bbc:bbc_t:bad_input_prior', ...
        'Third input (prior) must have eight elements.')
end
if sum(prior)~=1
    prior = prior / sum(prior);
end
    

%% Compute classical significance for behavioral model
alpha          = 0.05;
crit_t         = tinv(1 - alpha/2, df);
is_significant = abs(t) >= crit_t;


%% Rescale unit information prior for noncentrality parameter
scale_n = sqrt(df) / 2;       % scale factor to get delta from t
scale_p = 1;                  % std of delta
sigma   = scale_p * scale_n;  % std of the NCP


%% Compute evidence for each model
% integrate each model's likelihood over its parameter space
evidence = zeros(size(prior));

% p(t|df, M1+)
if prior(1) > 0
    evidence(1) = bbc_gquad('normal', ...
        @(delta)bbc_m1_ha_likelihood(delta), 0, sigma);
end

% p(t|df, M1-)
if prior(2) > 0
    evidence(2) = bbc_m1_h0_likelihood();
end

% p(t|df, M2+)
if prior(3) > 0
    evidence(3) = bbc_gquad('normal', ...
        @(delta)bbc_m2_ha_likelihood(delta), 0, sigma);
end

% p(t|df, M2-)
if prior(4) > 0
    evidence(4) = bbc_m2_h0_likelihood();
end

% p(t|df, M3+)
if prior(5) > 0
    evidence(5) = bbc_gquad('normal', ...
        @(delta)bbc_gquad('uniform', ...
        @(pr)bbc_m3_ha_likelihood(pr, delta), 0, 1), 0, sigma);
end

% p(t|df, M3-)
if prior(6) > 0
    evidence(6) = bbc_gquad('uniform', ...
        @(pr)bbc_m3_h0_likelihood(pr), 0, 1);
end

% p(t|df, M4+)
if prior(7) > 0
    evidence(7) = bbc_gquad('normal', ...
        @(delta)bbc_gquad('exponential', ...
        @(lambda)bbc_m4_ha_likelihood(lambda, delta), 5), 0, sigma);
end

% p(t|df, M4-)
if prior(8) > 0
    evidence(8) = bbc_gquad('exponential', ...
        @(lambda)bbc_m4_h0_likelihood(lambda), 5);
end

% p(M|t, df) = p(t|df, M) * p(M) / p(t)
posterior       = (evidence / sum(evidence .* prior)) .* prior ;

% pro-HA prior and posterior ratios
prior_ratio     = sum(prior(1:2:end)) ./ sum(prior(2:2:end));
posterior_ratio = sum(posterior(1:2:end)) ./ sum(posterior(2:2:end));

% pro-HA Bayes factor
bayes_factor    = posterior_ratio / prior_ratio;

% sum(prior(1:2:end).*evidence(1:2:end)) / sum(prior(2:2:end).*evidence(2:2:end))


%% -- Subroutines ------------------------------------------------------ %%
    % ------------------------------------------------------------------- %
    % -- Internal functions for likelihood computation ------------------ %
    % ------------------------------------------------------------------- %

    % -- p(t|df, M1-) --------------------------------------------------- %
    function f = bbc_m1_h0_likelihood()
        f = tpdf(t, df);
    end

    % -- p(t|df, delta, M1+) -------------------------------------------- %
    function f = bbc_m1_ha_likelihood(de)
        f = nctpdf(t, df, de);
    end

    % -- p(t|df, M2-) --------------------------------------------------- %
    function f = bbc_m2_h0_likelihood()
        if is_significant
            f = tpdf(t, df) ./ alpha;
        else
            f = 0;
        end
    end

    % -- p(t|df, delta, M2+) -------------------------------------------- %
    function f = bbc_m2_ha_likelihood(de)
        if is_significant
            A = 1 - nctcdf(crit_t, df, de) + nctcdf(-crit_t, df, de);
            f = nctpdf(t, df, de) ./ A;
        else
            f = 0;
        end
    end

    % -- p(t|df, pr, M3-) ----------------------------------------------- %
    function f = bbc_m3_h0_likelihood(pr)
        A = alpha + pr * (1 - alpha);
        
        if is_significant
            f = tpdf(t, df) ./ A;
        else
            f = pr .* tpdf(t, df) ./ A;
        end
        
    end

    % -- p(t|df, pr, delta, M3+) ---------------------------------------- %
    function f = bbc_m3_ha_likelihood(pr, de)
        A = (pr - 1) .* ...
            (nctcdf(crit_t, df, de) - nctcdf(-crit_t, df, de)) + 1;
        
        if is_significant
            f = nctpdf(t, df, de) ./ A;
        else
            f = pr .* nctpdf(t, df, de) ./ A;
        end
    end

    % -- p(t|df, lambda, M4-) ------------------------------------------- %
    function f = bbc_m4_h0_likelihood(la)
        scale = @(x) ...
            min(max(exp(-(la.*((2 * tcdf(-abs(x), df)) - alpha))), 0), 1);
        
        A = alpha + ...
            integral(@(x) scale(x) .* tpdf(x, df), -crit_t,      0) + ...
            integral(@(x) scale(x) .* tpdf(x, df),       0, crit_t);
        
        if is_significant
            f = tpdf(t, df) ./ A;
        else
            f = scale(t) .* tpdf(t, df) ./ A;
        end
    end

    % -- p(t|df, lambda, delta, M4+) ------------------------------------ %
    function f = bbc_m4_ha_likelihood(la, de)
        scale = @(x) ...
            min(max(exp(-(la.*((2 * tcdf(-abs(x), df)) - alpha))), 0), 1);
        
        fun0 = @(x) scale(x) .* nctpdf(x, df, de);   
        fun1 = @(x) nctpdf(x, df, de);   
        
        A1 = integral(fun1,     -Inf,  -crit_t);
        A2 = integral(fun0,  -crit_t,        0);
        A3 = integral(fun0,        0,   crit_t);
        A4 = integral(fun1,   crit_t,     +Inf);
        
        A = A1 + A2 + A3 + A4;
        
        if is_significant
            f = nctpdf(t, df, de) ./ A;
        else
            f = scale(t) .* nctpdf(t, df, de) ./ A;
        end
    end


    % ------------------------------------------------------------------- %
    % -- Internal function for numerical integration -------------------- %
    % ------------------------------------------------------------------- %

    function [s, f, x, w] = bbc_gquad(type, fcn, varargin)
    % BBC_GQUAD  Gaussian quadrature for norm, unif, or exp factors
    %    S = BBC_GQUAD(TYPE, FCN, ...) where FCN is a function handle and 
    %    TYPE is 'normal', 'exponential', or 'uniform', computes the 
    %    integral int{FCN(x)*TYPE(x)}dx over the range of TYPE.  Additional 
    %    arguments are mean and sd for 'normal', rate for 'exponential', 
    %    and lower and upper boundaries for 'uniform'.
    %
    %    BBC_GQUAD uses Gauss-Hermite (41 nodes), Gauss-Legendre (21 
    %    nodes), or Gauss-Laguerre (11 nodes) integration, as appropriate.
    
    if numel(type) < 3
        error('bbc:bbc_gquad:tooShortFactor', ...
            'Unknown factor distribution "%s".', type)
    end
    
    switch lower(type(1:3))
        case {'exp', 'lag'}
            [s, f, x, w] = glaint(varargin{:});
        case {'uni', 'leg'}
            [s, f, x, w] = gleint(varargin{:});
        case {'nor', 'gau', 'her'}
            [s, f, x, w] = ghint41(varargin{:});
        otherwise
            error('bbc:bbc_gquad:unknownFactor', ...
                'Unknown factor distribution "%s".', type)
    end
    
        % -- Gauss-Hermite integration ---------------------------------- %
        function [s, f, x, w] = ghint41(m, s)
            n = 41;
            x = [...
               -4.4851177450803652;-4.2534584199474237;-4.0230272504840592; 
               -3.7937340917372753;-3.5654937723650511;-3.3382255503269582; 
               -3.1118526303762439;-2.8863017340447339;-2.6615027143383929; 
               -2.4373882085955390;-2.2138933239623531;-1.9909553507579987; 
               -1.7685134996706828;-1.5465086592745176;-1.3248831708071220; 
               -1.1035806175170300;-0.8825456261915653;-0.6617236787209613; 
               -0.4410609317513270;-0.2205040426343026; 0; 
                0.2205040426343017; 0.4410609317513250; 0.6617236787209607;
                0.8825456261915643; 1.1035806175170286; 1.3248831708071209; 
                1.5465086592745159; 1.7685134996706822; 1.9909553507579976; 
                2.2138933239623517; 2.4373882085955385; 2.6615027143383920; 
                2.8863017340447339; 3.1118526303762439; 3.3382255503269587; 
                3.5654937723650506; 3.7937340917372744; 4.0230272504840583; 
                4.2534584199474237; 4.4851177450803661];
            w = [...
                0.0000000004262538; 0.0000000032096075; 0.0000000215027034; 
                0.0000001284750487; 0.0000006860664083; 0.0000032808733242; 
                0.0000140755898509; 0.0000542630994535; 0.0001882540944954; 
                0.0005885255267376; 0.0016599244925103; 0.0042284362004214; 
                0.0097376438372711; 0.0202896840527878; 0.0382790477154379; 
                0.0654308936994394; 0.1013836009016291; 0.1424623145688359; 
                0.1816025879081800; 0.2100559329700009; 0.2204952403727216; 
                0.2100559329699994; 0.1816025879081804; 0.1424623145688367; 
                0.1013836009016292; 0.0654308936994396; 0.0382790477154378; 
                0.0202896840527879; 0.0097376438372711; 0.0042284362004215; 
                0.0016599244925103; 0.0005885255267376; 0.0001882540944954; 
                0.0000542630994535; 0.0000140755898509; 0.0000032808733242; 
                0.0000006860664083; 0.0000001284750487; 0.0000000215027034; 
                0.0000000032096075; 0.0000000004262538];
            
            f = zeros(1, n);
            for g = 1:n
                f(g) = fcn(x(g) * s * sqrt(2) + m);
            end
            
            s = (f * w) / sqrt(pi);
            
        end
    
                
        % -- Gauss-Laguerre integration --------------------------------- %
        function [s, f, x, w] = glaint(l)
            n = 11;
            x = [...
                 0.125796442187968;  0.665418255839227;  1.647150545872169;
                 3.091138143035252;  5.029284401579828;  7.509887863806616;
                10.605950999546968; 14.431613758064190; 19.178857403214682;
                25.217709339677555; 33.497192847175540];
            w = [...
                 0.284933212894201;  0.389720889527850;  0.232781831848991;
                 0.076564453546197;  0.014393282767351;  0.001518880846485;
                 0.000085131224355;  0.000002292403880;  0.000000024863537;
                 0.000000000077126;  0.000000000000029];
            
            f = zeros(1, n);
            for g = 1:n
                f(g) = fcn(x(g) * l);
            end
            
            s = (f * w);
            
        end        
        
        % -- Gauss-Legendre integration --------------------------------- %
        function [s, f, x, w] = gleint(a, b)
            n = 21;
            x = [...
               -0.993752170620389;-0.967226838566306;-0.920099334150401;
               -0.853363364583317;-0.768439963475678;-0.667138804197412;
               -0.551618835887220;-0.424342120207439;-0.288021316802401;
               -0.145561854160896; 0.000000000000000; 0.145561854160895;
                0.288021316802401; 0.424342120207439; 0.551618835887220;
                0.667138804197412; 0.768439963475678; 0.853363364583317;
                0.920099334150401; 0.967226838566306; 0.993752170620390];
            w = [...
                0.016017228257775; 0.036953789770852; 0.057134425426856;
                0.076100113628380; 0.093444423456034; 0.108797299167148;
                0.121831416053729; 0.132268938633337; 0.139887394791073;
                0.144524403989970; 0.146081133649690; 0.144524403989970;
                0.139887394791074; 0.132268938633337; 0.121831416053729;
                0.108797299167149; 0.093444423456034; 0.076100113628379;
                0.057134425426858; 0.036953789770852; 0.016017228257775];
            
            f = zeros(1, n);
            for g = 1:n
                f(g) = fcn(x(g) * (b - a) / 2 + (a + b) / 2);
            end
            
            s = (f * w) / 2 ;
            
        end
        
    end

%% -- End subroutines -------------------------------------------------- %%
end
