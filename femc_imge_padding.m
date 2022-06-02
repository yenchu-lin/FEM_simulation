function [I] = femc_imge_padding(letter,imginfo)

% add zero padding based on npad (default npad = 2)
% if (nargin<=2); npad=[]; end
% if isempty(npad) == 1; npad=2; end
newsize = length(letter) * imginfo.npad;

% uint8 to double
doubleI = im2double(letter);
doubleI = doubleI/max(doubleI(:));

% adjust the contrast of the stimulus
adjI = doubleI * imginfo.cstim;

% add padding
n = (newsize-length(adjI))/2; 
I = [zeros(n,newsize);zeros(length(adjI),n),adjI,zeros(length(adjI),n); ...
    zeros(n,newsize)];

end

