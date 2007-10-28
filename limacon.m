
function f= limacon(x,R_0,diff_len)

f(1) = x(1)*cos(x(3)) + x(2)*cos(2*x(3));
f(2) = sin(x(3))*(x(1) + x(2)*cos(x(3))) - R_0;
f(3) = x(1) + x(2) - diff_len;
