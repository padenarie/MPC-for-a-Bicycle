function stab = stability(sysd)

A = sysd.A;
B = sysd.B;
C = sysd.C;

if sysd.Ts ==0;
    h = 0;
else
    h = sysd.Ts;
end

stab = isstable(sysd);


end