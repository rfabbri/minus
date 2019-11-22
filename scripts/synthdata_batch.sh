while IFS="" read -r p || [ -n "$p" ]
do
  #printf '%s\n' "$p"
  synthdata $p 
done < 100-configurations-synthdata
