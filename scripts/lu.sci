s = ''
for i=12:-1:0
  s= s + 'd(' + string(i) + ') -= '
  for k=(i+1):13
    s = s + '+ m(' + string(i) + ','  + string(k) + ')*d(' + string(k) + ')';
  end
  s = s + 'xx'
end
