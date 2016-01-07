function [values XF DF] = genRandValues(distName, numOfSamples, valMean, valStd, minVal, maxVal)
    if nargin<=4
        minVal = 0;
        maxVal = inf;
    elseif nargin<=5
        maxVal = inf;
    end

    distList = {'beta', 'birnbaumsaunders', 'exponential', ...
        'extreme value', 'gamma', 'generalized extreme value', ...
        'generalized pareto', 'inversegaussian', 'logistic', 'loglogistic', ...
        'lognormal', 'nakagami', 'normal', ...
        'rayleigh', 'rician', 'tlocationscale', 'weibull'};
    if strcmp(distName, 'uniform')
        values = 2*(valMean-minVal)*rand(numOfSamples,1) + minVal;
    elseif strcmp(distName, 'beta')
        values = (valMean-minVal)*betarnd(5,5,numOfSamples,1)+valMean;
    elseif strcmp(distName, 'generalized extreme value') || strcmp(distName,'tlocationscale') || strcmp(distName,'normal')
        pd = makedist(distName); 
        valTemp = random(pd,numOfSamples,1); 
        xCutOff = 3;
        valTemp = valTemp(valTemp>-xCutOff);
        valTemp = valTemp(valTemp<xCutOff); 
        [probabilities,xf] = ksdensity(valTemp);  
        probabilities = probabilities(xf>-xCutOff);
        xf = xf(xf>-xCutOff);       
        [values XF DF] = genValFromProb(xf, probabilities, numOfSamples, valMean, minVal, maxVal, valStd);
    elseif strcmp(distName,'solar') || strcmp(distName,'wind')|| strcmp(distName,'price') || strcmp(distName,'workload')
        % probability density function is based on this data
        hour24 = 13;
        if strcmp(distName,'solar')
            load(['data/probabilities/solar/probability_' int2str(hour24) '.mat']);
        elseif strcmp(distName,'wind')
            load(['data/probabilities/wind/probability_' int2str(hour24) '.mat']);
        elseif strcmp(distName,'price')
            load(['data/probabilities/price/probability_' int2str(hour24) '.mat']);
        else
            load(['data/probabilities/workload/probability_' int2str(hour24) '.mat']);
        end     
%         [values] = genValFromProb(probabilities, numOfSamples, valMean, valStd); 
        [values XF DF] = genValFromProb(xf, probabilities, numOfSamples, valMean, minVal, maxVal, valStd); 
    elseif strcmp(distName,'weird') 
        disp('weird distribution - not done yet');
        normValues =zeros(numOfSamples,1);
        percentages = 0.8;
        minIdx = ceil(0.8*numOfSamples);
        normValues(1:minIdx) = minVal;
        restMean = 1/(1-percentages)*valMean;
        restIdx = minIdx+1:numOfSamples;
        normValues(minIdx+1:numOfSamples) = randn(length(restIdx),1) + restMean; 
        randomIdx = ceil(rand(numOfSamples,1)*numOfSamples);
        curStd = std(normValues(randomIdx)); 
        values = valStd*(normValues./curStd);
    else
        error('unknown distribution');
    end 
    
    % mapping to right values       
%     minVal = max(minVal, min(values));
%     maxVal = min(maxVal, max(values));
%     minList = values<minVal;
%     maxList = values>maxVal;    
%     infeasibleList = minList+maxList;
%     feasibleValues= values(~infeasibleList);
%     feasibleIdx = ceil(rand(sum(infeasibleList),1)*length(feasibleValues));
%     values(infeasibleList==1) = feasibleValues(feasibleIdx);
end

function [values XF DF] = genValFromProb(xf, probabilities, numOfSamples, valMean, minVal, maxVal, valStd)
    pdf = probabilities;
    % normalize the function to have an area of 1.0 under it
    pdf = pdf / sum(pdf);        
    % the integral of PDF is the cumulative distribution function
    cdf = cumsum(pdf);
    % these two variables holds and describes the CDF    %xq:x %cdf: P(x)
    % remove non-unique elements
    [cdf, mask] = unique(cdf);
    xf = xf(mask);
    uniValues = rand(1, numOfSamples);
    % inverse interpolation to achieve P(x) -> x projection of the random values
    normValues = interp1(cdf, xf, uniValues, 'linear','extrap');
    normStd = std(normValues);    
    valScale = valStd/normStd;
    
%     valueRange = (normValues-mean(normValues))*valScale  + valMean;    
    valueRange = xf*valScale  + valMean;    
    
    truncatedList = valueRange>=minVal & valueRange<=maxVal;
    truncatedXf = xf(truncatedList);
    truncatedCdf = cdf(truncatedList);    
    truncatedPdf = pdf(truncatedList);
    truncatedUniValues = (max(truncatedCdf)-min(truncatedCdf))* rand(1, numOfSamples) + min(truncatedCdf);
    truncatedValues = interp1(truncatedCdf, truncatedXf, truncatedUniValues, 'linear','extrap');
    
    XF = truncatedXf; DF = truncatedPdf;
    
    values = truncatedValues*valScale + valMean;    
    std(values)
end

function [values] = genValFromProb2(probabilities, numOfSamples, valMean, valStd)
    x = 1 : length(probabilities);
    % do spline interpolation for smoothing
    xq = 1 : 0.05 : length(probabilities);
    pdf = interp1(x, probabilities, xq, 'spline');
    % remove negative elements due to the spline interpolation
    pdf(pdf < 0) = 0;
    % normalize the function to have an area of 1.0 under it
    pdf = pdf / sum(pdf);        
    % the integral of PDF is the cumulative distribution function
    cdf = cumsum(pdf);
    % these two variables holds and describes the CDF    %xq:x %cdf: P(x)
    % remove non-unique elements
    [cdf, mask] = unique(cdf);
    xq = xq(mask);
    uniValues = rand(1, numOfSamples);
    % inverse interpolation to achieve P(x) -> x projection of the random values
    normValues = interp1(cdf, xq, uniValues);
    normValues(isnan(normValues)) = 0;
    normValues = normValues - mean(normValues);
    curStd = std(normValues); 
    values = valStd*(normValues./curStd);
    values = values - mean(values) + valMean;
end