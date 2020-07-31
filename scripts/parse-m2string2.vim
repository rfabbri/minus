"toExternalString in M2,
"then put in bla, let it beign with point { {..
"then join all lines into a single one bf running this 
%s/point { {//ge
%s/} }//g
%s/toCC\s*(\([^(]*\))/{\1}/g
%s/,/,/g
%s/},/},/g
%s/p53//g
g/^\s*$/d  "remove blank lines
%s/^\s*//g "remove start blanks
g/{/,/^[-\.]/join  " join all lines 1) starting with { and 2) starting with + or .
