function bvalues = bcfcn(Ti,Tf,Initial_T,Final_T)
bvalues = zeros(2,1);
bvalues(1,1) = Ti(1,1) - Initial_T;
bvalues(2,1) = Tf(1,1) - Final_T;
end