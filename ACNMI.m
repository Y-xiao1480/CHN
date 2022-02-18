function [NMI,AC]=ACNMI(label,gnd)
NMI = MutualInfo(gnd,label);
label = bestMap(gnd,label);
AC=length(find(gnd==label))/length(gnd);