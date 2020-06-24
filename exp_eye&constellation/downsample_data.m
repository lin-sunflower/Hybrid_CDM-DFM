function [output] = downsample_data(input,baudrate,samplerate,point)

%%% point should be no larger thant 
sample_start=floor(point*samplerate);
n=1:length(input); 
m=baudrate;    
s=1/m;         
xi=1:s:length(input);
output= interp1(n,input,xi,'linear');
output=output(sample_start:samplerate:end);
end

