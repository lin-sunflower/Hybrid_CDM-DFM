function [output] = normalized_power(input)

input=input-mean(input);
power=sum(input.^2)/(length(input));
output=input/sqrt(power);

end

