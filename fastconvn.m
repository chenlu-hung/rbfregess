% FASTCONVN - performes convolution in fft domain - like convn(x,h,'same')
% Input:    x - data (output vil be the same size as x)
%           h - filter 
% Output:   y - returns convolution (~ within machine precision)
%
% Author:   Anders Ueland
% Version:  September -15: created the file
% Licence:  Use and edit as you like. Provided 'as-is' - use on your own risk.

function y = fastconvn(x,h)

%Ensure double:
% x = double(x);
% h = double(h);

%Check if input is imag
useComplexOut = 0;
if ~ ( isreal(x) && isreal(h))
    useComplexOut = 1;
end

%Measure input
xSize = size(x);
hSize = size(h);

%Check max dim
if  max( [length(hSize), length(xSize)]) > 9
    error('Only arrays with dimension up to 9 is supported. But it is fairly easy to extend the source code if needed.')
end

%Ensure that xSize and hSize have the same size
%diff = abs(length(xSize)-length(hSize));
xSize = [xSize ones(1, max( [length(hSize)-length(xSize),0]))];
hSize = [hSize ones(1, max( [length(xSize)-length(hSize),0]))];

%Length of fft
fftLen =  xSize + hSize - 1;

%Perform convolution in fft - domain
X = fftn(x,fftLen);
H = fftn(h,fftLen);
Y = X.*H;
y = ifftn(Y);

%Convert to real
if ~useComplexOut
    y = real(y);
end

%Chop out central part of output
extra = min([xSize;hSize],[],1)-1;
cutAtStart = ceil(extra./2);
cutAtEnd = floor(extra./2);

%Add some extra zeros so we can suppport higher dimmensions
cutAtStart = [cutAtStart 0,0,0,0,0,0,0,0,0];
cutAtEnd = [cutAtEnd 0,0,0,0,0,0,0,0,0];

%Make the cut
y = y(  1+cutAtStart(1):end-cutAtEnd(1),...
        1+cutAtStart(2):end-cutAtEnd(2),...
        1+cutAtStart(3):end-cutAtEnd(3),...
        1+cutAtStart(4):end-cutAtEnd(4),...
        1+cutAtStart(5):end-cutAtEnd(5),...
        1+cutAtStart(6):end-cutAtEnd(6),...
        1+cutAtStart(7):end-cutAtEnd(7),...
        1+cutAtStart(8):end-cutAtEnd(8),...
        1+cutAtStart(9):end-cutAtEnd(9) );
end