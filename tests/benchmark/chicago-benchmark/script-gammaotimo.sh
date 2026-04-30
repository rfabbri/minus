#!/bin/bash
# Author: Juliana Ventura
# See GEMINI.md in benchmarksdir below

#Script para rodar os primeiros 140 problemas aleatórios gerados a partir dos dados sintéticos
#MAIO 2025

set -x
minusdir=/home/juliana/Documents/doutorado/minus   # TODO: infer this automatically eg from cmake
synthdata=$minusdir/scripts/synthdata/synthdata
minus=$minusdir/bin/minus-chicago
benchmarkdir=$minusdir/tests/benchmark
configurations=$benchmarkdir/FRePTcorrigido.txt

for i in $(seq 1 2)
do 
    mkdir frpt-P$i
    fr1=$(cat $configurations | sed -n "$(($i))p" |cut -d' ' -f1)
    fr2=$(cat $configurations | sed -n "$(($i))p" |cut -d' ' -f2)
    fr3=$(cat $configurations | sed -n "$(($i))p" |cut -d' ' -f3)
    pt1=$(cat $configurations | sed -n "$(($i))p" |cut -d' ' -f4)
    pt2=$(cat $configurations | sed -n "$(($i))p" |cut -d' ' -f5)
    pt3=$(cat $configurations | sed -n "$(($i))p" |cut -d' ' -f6)
    for j in $(seq 1 5) # number of runs
    do 
        script -c "$synthdata $fr1 $fr2 $fr3 $pt1 $pt2 $pt3 0 1 |$minus -i -gt  --prefilter_degeneracy=no > out.txt" temp.txt
        sed -i 's/\r$//' temp.txt

#	cat out.txt | cut -d '(' -f 2| cut -d ")" -f 1| sed -n '1,114p'>> frpt-P$i/vec_gamma-P$i.txt
       grep -E "^LOG [0-9]+ " temp.txt | cut -d' ' -f 3 >> frpt-P$i/sol_step-P$i.txt
       grep "isvalid)" temp.txt | cut -d')' -f2 >> frpt-P$i/realsol-P$i.txt
       echo '313,10'>> frpt-P$i/realsol-P$i.txt
       solution=$(grep index temp.txt | cut -d : -f 2)
       if [ -z "$solution" ]
       then
           solution=313
       fi
       echo $solution>> frpt-P$i/gt_sol-P$i.txt
       cat temp.txt | cut -d' ' -f 2-5| sed -n '2,3p' >> frpt-P$i/frpt-$i.txt
    done
done

# Generate Version 1.1 Summaries
echo "Calculating statistics and generating index.html / summary.txt..."
dos2unix frpt-P*/sol_step-P*.txt
cat frpt-P*/sol_step-P*.txt > all_steps.txt

cat << 'EOF' > update_html.py
import os
import subprocess

def get_stats(file_path):
    s = subprocess.check_output(["st-console", "--sum", "--no-header", file_path]).decode().strip()
    m = subprocess.check_output(["st-console", "--median", "--no-header", file_path]).decode().strip()
    return s, m

html = """<!DOCTYPE html>
<html>
<head>
    <title>MINUS Benchmark Results v1.1</title>
    <style>
        table { border-collapse: collapse; width: 50%; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #f2f2f2; }
    </style>
</head>
<body>
    <h1>MINUS Benchmark Results (Version 1.1)</h1>
"""

txt = "MINUS Benchmark Results Summary\n"
txt += "-------------------------------\n\n"
txt += "Per Configuration:\n"

grand_sum = 0
all_files = []

# Find all frpt-P* directories and their sol_step files
dirs = sorted([d for d in os.listdir('.') if d.startswith('frpt-P')], key=lambda x: int(x.split('-P')[1]))

html += "<h2>Per Configuration Statistics</h2>"
html += "<table><tr><th>Configuration</th><th>Total Steps</th><th>Median Steps</th></tr>"

for d in dirs:
    config_name = d.split('-')[1]
    step_file = os.path.join(d, f"sol_step-{config_name}.txt")
    if os.path.exists(step_file):
        s, m = get_stats(step_file)
        html += f"<tr><td>{config_name}</td><td>{s}</td><td>{m}</td></tr>"
        txt += f"Config {config_name}: Sum = {s}, Median = {m}\n"
        grand_sum += float(s)

# Grand Stats
grand_median = subprocess.check_output(["st-console", "--median", "--no-header", "all_steps.txt"]).decode().strip()

html += "</table>"
html += "<h2>Grand Totals</h2>"
html += f"<p>Total iteration steps TOTAL: <strong>{int(grand_sum)}</strong></p>"
html += f"<p>Grand Median of steps: <strong>{grand_median}</strong></p>"
html += "</body></html>"

txt += f"\nGrand Totals:\n"
txt += f"Total Iteration Steps: {int(grand_sum)}\n"
txt += f"Grand Median Steps: {grand_median}\n"

with open("index.html", "w") as f:
    f.write(html)

with open("summary.txt", "w") as f:
    f.write(txt)
EOF

python3 update_html.py
echo "Done! See index.html and summary.txt."

