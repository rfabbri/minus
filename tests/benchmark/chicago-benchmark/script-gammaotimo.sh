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
import json

def get_stats(file_path):
    s = subprocess.check_output(["st-console", "--sum", "--no-header", file_path]).decode().strip()
    m = subprocess.check_output(["st-console", "--median", "--no-header", file_path]).decode().strip()
    return s, m

txt = "MINUS Benchmark Results Summary\n"
txt += "-------------------------------\n\n"
txt += "Per Configuration:\n"

grand_sum = 0
all_files = []
plot_data = []

# Load all lines from FRePTcorrigido.txt
try:
    with open("../../FRePTcorrigido.txt", "r") as frept_f:
        frept_lines = [line.strip() for line in frept_f.readlines()]
except Exception:
    frept_lines = []

# Find all frpt-P* directories and their sol_step files
dirs = sorted([d for d in os.listdir('.') if d.startswith('frpt-P')], key=lambda x: int(x.split('-P')[1]))

for d in dirs:
    config_name = d.split('-')[1]
    config_id = int(config_name.replace('P', ''))
    step_file = os.path.join(d, f"sol_step-{config_name}.txt")
    
    if os.path.exists(step_file):
        s, m = get_stats(step_file)
        txt += f"Config {config_name}: Sum = {s}, Median = {m}\n"
        grand_sum += float(s)
        
        # Read raw data for boxplot
        with open(step_file, 'r') as sf:
            steps = [int(line.strip()) for line in sf if line.strip().isdigit()]
        
        # Read problem config string from FRePTcorrigido.txt
        config_str = "No config found"
        if 0 < config_id <= len(frept_lines):
            config_str = frept_lines[config_id - 1]
                
        plot_data.append({
            "x": config_name,
            "y": steps,
            "config": config_str
        })

# Grand Stats
grand_median = subprocess.check_output(["st-console", "--median", "--no-header", "all_steps.txt"]).decode().strip()

txt += f"\nGrand Totals:\n"
txt += f"Total Iteration Steps: {int(grand_sum)}\n"
txt += f"Grand Median Steps: {grand_median}\n"

with open("summary.txt", "w") as f:
    f.write(txt)

# Write data.js
with open("data.js", "w") as f:
    f.write(f"const benchmarkData = {json.dumps(plot_data)};\n")
    f.write(f"const grandMedian = {grand_median};\n")

html = """<!DOCTYPE html>
<html>
<head>
    <title>MINUS Benchmark Results v1.2</title>
    <script src="https://cdn.plot.ly/plotly-2.32.0.min.js"></script>
    <style>
        body { font-family: sans-serif; margin: 20px; background-color: #fafafa; }
        #myDiv { background-color: white; border: 1px solid #ccc; box-shadow: 0 4px 8px rgba(0,0,0,0.1); }
        #infoPanel { margin-top: 20px; padding: 15px; background: white; border: 1px solid #ccc; display: none; }
    </style>
</head>
<body>
    <h1>MINUS Benchmark Results (Version 1.2)</h1>
    <div id="myDiv" style="width: 100%; height: 800px;"></div>
    
    <div id="infoPanel">
        <h3>Configuration Details</h3>
        <p><strong>ID:</strong> <span id="configId"></span></p>
        <p><strong>Synthdata:</strong> <span id="configStr"></span></p>
    </div>

    <script src="data.js"></script>
    <script>
        const traces = [];
        
        benchmarkData.forEach(item => {
            // Plotly's native boxplots can't perfectly decouple stem/box widths in pixels easily, 
            // but we can make them very thin, set the outlier symbol, and style the median.
            traces.push({
                y: item.y,
                name: item.x,
                type: 'box',
                boxpoints: 'outliers', // Show only outliers points
                jitter: 0,
                pointpos: 0,
                width: 0.2, // Narrow width
                marker: { 
                    color: 'lightgray', 
                    size: 4, 
                    symbol: 'cross' // '+' sign
                }, 
                line: { 
                    color: 'blue', 
                    width: 1 // 1px stem
                }, 
                fillcolor: 'blue',
                hoverinfo: 'y',
                customdata: [item.config]
            });
            
            // Calculate median to overlay the white/black dot
            item.y.sort((a,b) => a-b);
            let med;
            let mid = Math.floor(item.y.length / 2);
            if (item.y.length % 2 === 0) {
                med = (item.y[mid - 1] + item.y[mid]) / 2;
            } else {
                med = item.y[mid];
            }
            
            // Overlay median dot trace
            traces.push({
                x: [item.x],
                y: [med],
                mode: 'markers',
                marker: {
                    color: 'white',
                    size: 10,
                    line: {
                        color: 'black',
                        width: 1.5
                    }
                },
                hoverinfo: 'skip',
                showlegend: false
            });
        });

        const layout = {
            title: '<b>Time spent in paths leading to undesired roots</b>',
            titlefont: { size: 24 },
            showlegend: false,
            yaxis: {
                title: 'Number of iterations per root',
                titlefont: { size: 20 },
                zeroline: false,
                range: [0, 1050],
                dtick: 100,
                gridcolor: '#eee'
            },
            xaxis: {
                title: 'Problem',
                titlefont: { size: 20 },
                tickangle: -90,
                tickfont: { size: 10 }
            },
            shapes: [
                // Green horizontal line for grand median
                {
                    type: 'line',
                    xref: 'paper',
                    x0: 0,
                    x1: 1,
                    yref: 'y',
                    y0: grandMedian,
                    y1: grandMedian,
                    line: {
                        color: 'lime',
                        width: 3
                    }
                }
            ],
            annotations: [
                {
                    x: 1,
                    y: grandMedian,
                    xref: 'paper',
                    yref: 'y',
                    text: 'Median time',
                    showarrow: false,
                    xanchor: 'right',
                    yanchor: 'bottom',
                    font: { size: 16 },
                    bgcolor: 'white',
                    bordercolor: 'black'
                }
            ],
            plot_bgcolor: 'white',
            paper_bgcolor: '#fafafa',
            margin: { l: 80, r: 50, b: 80, t: 80 }
        };

        Plotly.newPlot('myDiv', traces, layout, {responsive: true});

        // Add Interactivity
        document.getElementById('myDiv').on('plotly_click', function(data){
            if(data.points.length > 0) {
                const pt = data.points[0];
                const configId = pt.data.name;
                const configStr = pt.data.customdata[0];
                
                const panel = document.getElementById('infoPanel');
                document.getElementById('configId').innerText = configId;
                document.getElementById('configStr').innerText = configStr;
                panel.style.display = 'block';
            }
        });
    </script>
</body>
</html>
"""

with open("index.html", "w") as f:
    f.write(html)

EOF

python3 update_html.py
echo "Done! See index.html and summary.txt."

