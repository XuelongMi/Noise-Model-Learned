function [X_optimal] = pdf_recover(pdf,weightType,nonincreasing)

if(~exist('weightType','var'))
    weightType = 'Average';
end
if(~exist('nonincreasing','var'))
    nonincreasing = 0;
end

%% recover pdf
X = deautocorrelation(pdf);
%% refine
X_optimal = refine(pdf,X,weightType,nonincreasing);

end