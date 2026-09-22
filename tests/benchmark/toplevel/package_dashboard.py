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
        
    # 3. Load & Process History (keep 3 latest distinct hashes)
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
    history = history[:3] # keep max 3 commits
    
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
                    
        # 5. Inject Button Bar and Trend Plot into individual problem HTMLs
        # Reversing history for plot (chronological: oldest to newest left to right)
        plot_history = list(reversed(history))
        x_data = [h["commit"] for h in plot_history]
        y_data = [h["medians"].get(p) or 0 for h in plot_history]
        
        def get_injected_html(current_view_hash):
            buttons = []
            for h in history:
                c = h["commit"]
                badge = " (Latest)" if c == history[0]["commit"] else ""
                is_active = (current_view_hash == c)
                
                style = "padding:5px 12px; margin-right:10px; font-weight:bold; background:#2196F3; color:white; border:none; border-radius:4px; cursor:pointer;"
                disabled_style = "padding:5px 12px; margin-right:10px; font-weight:bold; background:#cbd5e1; color:#64748b; border:none; border-radius:4px; cursor:default;"
                
                if is_active:
                    buttons.append(f'<button style="{disabled_style}" disabled>{c}{badge}</button>')
                else:
                    # Link points to `-run` variant explicitly to avoid circular/broken paths if browsing history locally
                    target_dir = f"{c}-run" if c != history[0]["commit"] else "tmp"
                    buttons.append(f'<a href="../{target_dir}/index.html"><button style="{style}">{c}{badge}</button></a>')
                    
            btns_str = "".join(buttons)
            
            return f"""
            <!-- INJECTED_UI_START -->
            <div style="margin:40px auto; max-width:1000px; padding:20px; background:white; border:1px solid #e2e8f0; border-radius:8px; box-shadow: 0 4px 15px rgba(0,0,0,0.05); font-family:sans-serif;">
              <div style="margin-bottom:20px; text-align:center;">
                <span style="font-size:16px; margin-right:15px; color:#334155;"><strong>Compare Commits: </strong></span>
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
                    name: 'Grand Median'
                  }};
                  Plotly.newPlot('trendPlot', [trace], {{
                    title: '<b>Trend: Grand Median Steps per Commit</b>',
                    margin: {{ t: 40, b: 40, l: 40, r: 20 }},
                    paper_bgcolor: 'rgba(0,0,0,0)',
                    plot_bgcolor: 'rgba(0,0,0,0)'
                  }});
              }}
            </script>
            <!-- INJECTED_UI_END -->
            """

        dirs_to_process = [("tmp", commit_6)] + [(f"{h}-run", h) for h in valid_hashes]
        for d_name, d_hash in set(dirs_to_process):
            idx_path = os.path.join(p_dir, d_name, "index.html")
            if os.path.exists(idx_path):
                with open(idx_path, "r") as f:
                    html = f.read()
                
                # Strip old injections
                html = re.sub(r'<!-- INJECTED_UI_START -->.*?<!-- INJECTED_UI_END -->', '', html, flags=re.DOTALL)
                
                inj_html = get_injected_html(d_hash)
                
                if "</body>" in html:
                    html = html.replace("</body>", inj_html + "\n</body>")
                else:
                    html += inj_html
                    
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
        
        if os.path.exists(os.path.join(p_dir, "tmp")):
            shutil.copytree(os.path.join(p_dir, "tmp"), os.path.join(pub_p, "tmp"))
        
        for h in valid_hashes:
            rname = f"{h}-run"
            if os.path.exists(os.path.join(p_dir, rname)):
                shutil.copytree(os.path.join(p_dir, rname), os.path.join(pub_p, rname))
                
    print(f"Dashboard packaged successfully for Hash {commit_6}.")
    print(f"Active History: {valid_hashes}")

if __name__ == "__main__":
    main()
