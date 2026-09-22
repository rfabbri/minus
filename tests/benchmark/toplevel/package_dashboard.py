#!/usr/bin/env python3
import os
import sys
import json
import shutil
import subprocess
import re

def main():
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../.."))
    toplevel_dir = os.path.join(repo_root, "tests", "benchmark", "toplevel")
    history_file = os.path.join(toplevel_dir, "history.json")
    public_dir = os.path.join(repo_root, "public")
    
    # 1. Get current Git details
    try:
        commit_6 = subprocess.check_output(["git", "rev-parse", "--short=6", "HEAD"], cwd=repo_root).decode().strip()[:6]
    except:
        commit_6 = "local"
        
    try:
        update_date = subprocess.check_output(["git", "log", "-1", "--format=%cd", "--date=format:%Y-%m-%d %H:%M UTC"], cwd=repo_root).decode().strip()
    except:
        update_date = "Recently"
        
    problems = ["chicago", "linecircle"]
    medians = {}
    
    # 2. Extract Grand Median Steps from currently generated summary.txt inside tmp/
    for p in problems:
        sum_file = os.path.join(repo_root, f"tests/benchmark/individual/{p}-benchmark/tmp/summary.txt")
        m = None
        if os.path.exists(sum_file):
            with open(sum_file, "r") as f:
                content = f.read()
                # e.g., "Grand Median Steps: 179.5"
                match = re.search(r"Grand Median Steps:\s*([\d\.]+)", content)
                if match:
                    m = float(match.group(1))
        medians[p] = m
        
    # 3. Load & Process History (keep up to 5 latest distinct hashes)
    history = []
    if os.path.exists(history_file):
        try:
            with open(history_file, 'r') as f:
                history = json.load(f)
        except Exception as e: 
            print("Warning: Could not parse history.json", e)
            
    # Remove any existing entry for this commit (to update it), then insert at top
    history = [h for h in history if h["commit"] != commit_6]
    history.insert(0, {
        "commit": commit_6,
        "date": update_date,
        "medians": medians
    })
    history = history[:5] # keep max 5 commits
    
    with open(history_file, 'w') as f:
        json.dump(history, f, indent=2)
        
    valid_hashes = [h["commit"] for h in history]
    
    # 4. Handle Archiving: copy tmp/ to <hash>-run/ and cleanup old ones
    for p in problems:
        p_dir = os.path.join(repo_root, f"tests/benchmark/individual/{p}-benchmark")
        tmp_dir = os.path.join(p_dir, "tmp")
        run_name = f"{commit_6}-run"
        archive_dir = os.path.join(p_dir, run_name)
        
        # Copy current tmp into its static <hash>-run snapshot
        if os.path.exists(tmp_dir):
            if os.path.exists(archive_dir):
                shutil.rmtree(archive_dir)
            shutil.copytree(tmp_dir, archive_dir)
            
        # Cleanup any *-run directory that isn't in valid_hashes
        for item in os.listdir(p_dir):
            if item.endswith("-run"):
                item_hash = item.replace("-run", "")
                if item_hash not in valid_hashes:
                    old_path = os.path.join(p_dir, item)
                    shutil.rmtree(old_path)
                    
        # 5. Inject Navigation, Button Bar and Trend Plot into individual problem HTMLs
        # Reversing history for plot (chronological: oldest to newest left to right)
        plot_history = list(reversed(history))
        x_data = [h["commit"] for h in plot_history]
        y_data = [h["medians"].get(p) or 0 for h in plot_history]
        
        def get_injected_html(current_view_hash, top_href):
            buttons = []
            for h in history:
                c = h["commit"]
                badge = " (Latest)" if c == history[0]["commit"] else ""
                is_active = (current_view_hash == c)
                
                style = "padding:6px 14px; margin-right:10px; font-weight:bold; background:#2563eb; color:white; border:none; border-radius:6px; cursor:pointer; font-size:13px; text-decoration:none; display:inline-block;"
                disabled_style = "padding:6px 14px; margin-right:10px; font-weight:bold; background:#e2e8f0; color:#64748b; border:none; border-radius:6px; cursor:default; font-size:13px; display:inline-block;"
                
                if is_active:
                    buttons.append(f'<span style="{disabled_style}">{c}{badge}</span>')
                else:
                    target_dir = f"{c}-run" if c != history[0]["commit"] else "tmp"
                    buttons.append(f'<a href="../{target_dir}/index.html" style="text-decoration:none;"><button style="{style}">{c}{badge}</button></a>')
                    
            btns_str = "".join(buttons)
            
            top_nav = f"""
            <!-- INJECTED_NAV_START -->
            <div style="display:flex; justify-content:space-between; align-items:center; margin:0 0 25px 0; padding:12px 20px; background:white; border:1px solid #e2e8f0; border-radius:8px; box-shadow:0 2px 6px rgba(0,0,0,0.04); font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
              <a href="{top_href}" style="display:inline-flex; align-items:center; gap:8px; text-decoration:none; color:#1d4ed8; font-weight:600; font-size:14px; padding:6px 14px; background:#eff6ff; border:1px solid #bfdbfe; border-radius:6px; transition:background 0.2s;">
                &larr; Back to Benchmark Dashboard
              </a>
              <span style="font-size:14px; color:#64748b;">
                Problem: <strong style="color:#0f172a; text-transform:capitalize;">{p}</strong> &bull; Viewing Commit <code style="background:#f1f5f9; padding:2px 6px; border-radius:4px; font-weight:bold; color:#0f172a;">{current_view_hash}</code>
              </span>
            </div>
            <!-- INJECTED_NAV_END -->
            """
            
            bottom_ui = f"""
            <!-- INJECTED_UI_START -->
            <div style="margin:40px auto; max-width:1000px; padding:24px; background:white; border:1px solid #e2e8f0; border-radius:10px; box-shadow:0 4px 15px rgba(0,0,0,0.05); font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
              <div style="margin-bottom:20px; text-align:center;">
                <span style="font-size:15px; margin-right:15px; color:#334155; font-weight:600;">Compare Commits:</span>
                {btns_str}
              </div>
              <div id="trendPlot" style="width:100%; height:320px;"></div>
            </div>
            <script>
              if(typeof Plotly === 'undefined') {{
                  var s = document.createElement('script');
                  s.src = "https://cdn.plot.ly/plotly-2.32.0.min.js";
                  s.onload = function() {{ renderPlot(); }};
                  document.head.appendChild(s);
              }} else {{ renderPlot(); }}
              
              function renderPlot() {{
                  var trace = {{
                    x: {json.dumps(x_data)},
                    y: {json.dumps(y_data)},
                    type: 'scatter',
                    mode: 'lines+markers',
                    marker: {{size: 10, color: '#e74c3c'}},
                    line: {{width: 3, color: '#c0392b'}},
                    name: 'Grand Median',
                    hovertemplate: '<b>Commit %{{x}}</b><br>Grand Median: %{{y:.1f}}<extra></extra>'
                  }};
                  Plotly.newPlot('trendPlot', [trace], {{
                    title: '<b>Trend: Grand Median Steps per Commit</b>',
                    margin: {{ t: 40, b: 40, l: 40, r: 20 }},
                    paper_bgcolor: 'rgba(0,0,0,0)',
                    plot_bgcolor: 'rgba(0,0,0,0)',
                    xaxis: {{ title: 'Commit Hash', type: 'category' }},
                    yaxis: {{ title: 'Grand Median Steps' }}
                  }});
              }}
            </script>
            <!-- INJECTED_UI_END -->
            """
            return top_nav, bottom_ui

        dirs_to_process = [("tmp", commit_6)] + [(f"{h}-run", h) for h in valid_hashes]
        for d_name, d_hash in set(dirs_to_process):
            idx_path = os.path.join(p_dir, d_name, "index.html")
            if os.path.exists(idx_path):
                with open(idx_path, "r") as f:
                    html = f.read()
                
                # Strip old injections
                html = re.sub(r'<!-- INJECTED_NAV_START -->.*?<!-- INJECTED_NAV_END -->', '', html, flags=re.DOTALL)
                html = re.sub(r'<!-- INJECTED_UI_START -->.*?<!-- INJECTED_UI_END -->', '', html, flags=re.DOTALL)
                
                top_nav, bottom_ui = get_injected_html(d_hash, "../../../toplevel/www/index.html")
                
                # Inject top nav right after <body>
                if "<body>" in html:
                    html = html.replace("<body>", "<body>\n" + top_nav, 1)
                
                # Upgrade green median line to hoverable trace if not already converted
                if "Grand Median: ' + grandMedian" not in html and "Plotly.newPlot('myDiv'" in html:
                    upgrade_script = """
        if (typeof grandMedian !== 'undefined' && !traces.some(t => t.name === 'Grand Median') && typeof benchmarkData !== 'undefined' && benchmarkData.length > 0) {
            traces.push({
                x: benchmarkData.map(d => d.x),
                y: benchmarkData.map(() => grandMedian),
                mode: 'lines',
                line: { color: 'lime', width: 3.5 },
                name: 'Grand Median',
                hoverinfo: 'text',
                hovertext: benchmarkData.map(() => 'Grand Median: ' + grandMedian),
                hoverlabel: { bgcolor: '#1b5e20', font: { color: '#ffffff', size: 14 } },
                showlegend: false
            });
            if (typeof layout !== 'undefined' && layout.shapes) {
                layout.shapes = layout.shapes.filter(s => !(s.line && s.line.color === 'lime'));
            }
        }
        Plotly.newPlot('myDiv'"""
                    html = html.replace("Plotly.newPlot('myDiv'", upgrade_script, 1)
                
                # Inject bottom UI before </body>
                if "</body>" in html:
                    html = html.replace("</body>", bottom_ui + "\n</body>")
                else:
                    html += bottom_ui
                    
                with open(idx_path, "w") as f:
                    f.write(html)
                    
    # 6. Build the public/ site directory and rewrite toplevel index.html
    if os.path.exists(public_dir):
        shutil.rmtree(public_dir)
    os.makedirs(public_dir, exist_ok=True)
    
    template_in = os.path.join(toplevel_dir, "www", "index.html.in")
    with open(template_in, "r") as f:
        template = f.read()
        
    template = template.replace('{{COMMIT_SHA_6}}', commit_6)
    template = template.replace('{{LAST_UPDATED}}', update_date)
    
    # 6.a Local Generation (write back to www/index.html)
    local_html = template.replace('{{CHICAGO_HREF}}', '../../individual/chicago-benchmark/tmp/index.html')
    local_html = local_html.replace('{{LINECIRCLE_HREF}}', '../../individual/linecircle-benchmark/tmp/index.html')
    with open(os.path.join(toplevel_dir, "www", "index.html"), "w") as f:
        f.write(local_html)
        
    # 6.b Public Generation (write to public/index.html)
    pub_html = template.replace('{{CHICAGO_HREF}}', './chicago-benchmark/tmp/index.html')
    pub_html = pub_html.replace('{{LINECIRCLE_HREF}}', './linecircle-benchmark/tmp/index.html')
    with open(os.path.join(public_dir, "index.html"), "w") as f:
        f.write(pub_html)
        
    # 6.c Copy problem directories cleanly preserving tmp and archives identically
    for p in problems:
        pub_p = os.path.join(public_dir, f"{p}-benchmark")
        os.makedirs(pub_p, exist_ok=True)
        p_dir = os.path.join(repo_root, f"tests/benchmark/individual/{p}-benchmark")
        
        # Copy tmp
        if os.path.exists(os.path.join(p_dir, "tmp")):
            shutil.copytree(os.path.join(p_dir, "tmp"), os.path.join(pub_p, "tmp"))
        
        # Copy historic runs
        for h in valid_hashes:
            rname = f"{h}-run"
            if os.path.exists(os.path.join(p_dir, rname)):
                shutil.copytree(os.path.join(p_dir, rname), os.path.join(pub_p, rname))
                
        # Fix public back links (replace local path with public root path)
        for root, dirs, files in os.walk(pub_p):
            for file in files:
                if file == "index.html":
                    fpath = os.path.join(root, file)
                    with open(fpath, "r") as f:
                        f_html = f.read()
                    f_html = f_html.replace("../../../toplevel/www/index.html", "../../index.html")
                    with open(fpath, "w") as f:
                        f.write(f_html)
                        
    print(f"Dashboard packaged successfully for Hash {commit_6}.")
    print(f"Active History (up to 5): {valid_hashes}")

if __name__ == "__main__":
    main()
