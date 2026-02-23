v=[ 11 13 25 27 39 41 53 55 67 68 81 82 95 96 109 110 124 125 138 139 152 153 166 167 180 181 194 195];
r = zeros(v);
c = zeros(v);

col_major=%t;
if col_major
  c = int(v/14); 
  r = modulo(int(v),14);
else
  r = int(v/15); 
  c = modulo(int(v),15);
end

// build matrices
m = zeros(14,15);
for i=1:length(v)
  m(r(i)+1,c(i)+1) = 1;
end

muni = zeros(14,15);
muni(v+1) = 1;
cc=[r' c']
