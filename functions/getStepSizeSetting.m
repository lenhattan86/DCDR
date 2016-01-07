function [ A, a ] = getStepSizeSetting( Start, End, stepStart, stepEnd, alpha )
    A = (stepEnd^(1/alpha)*End - stepStart^(1/alpha)*Start)/(stepStart^(1/alpha) - stepEnd^(1/alpha));
    aRootAlpha = stepStart^(1/alpha)*(Start + A);
    a = aRootAlpha^(alpha);
    if a <= 0 || A <=-1
        error('This step-size setting does not work!');
    end
end

