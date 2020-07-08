" to be used on an output of
" toExternalString with complex number entries
%s/,/,/g
%s/},/},/g
%s/}$/},/g
g/^\s*$/d  "remove blank lines
%s/^\s*//g "remove start blanks
%s/[)(]//g
%s/\*ii}/}/g
%s/[+-]\.[0-9]*}/,&/g  "separate imaginary part by ,
