function [approx] = roundd(num,decimal)
%round to desidered decimal.
%Ex: roundd(0.145,2) -> 0.15

    approx = round(num*10^decimal)/10^decimal;

end

