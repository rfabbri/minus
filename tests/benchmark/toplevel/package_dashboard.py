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
    robustness = {}
    
    # 2. Extract Grand Median Steps and Robustness from currently generated summary.txt inside tmp/
    for p in problems:
        sum_file = os.path.join(repo_root, f"tests/benchmark/individual/{p}-benchmark/tmp/summary.txt")
        m = None
        r = None
        if os.path.exists(sum_file):
            with open(sum_file, "r") as f:
                content = f.read()
                # e.g., "Grand Median Steps: 179.5"
                match_m = re.search(r"Grand Median Steps:\s*([\d\.]+)", content)
                if match_m:
                    m = float(match_m.group(1))
                # e.g., "Grand robustness / ground-truth found within all runs: 85%"
                # or "Grand Successs Rate / ground-truth found within all runs & all configs: 95.0%"
                match_r = re.search(r"Grand (?:robustness|Successs? Rate).*?:\s*([\d\.]+)%", content)
                if match_r:
                    r = float(match_r.group(1))
        medians[p] = m
        robustness[p] = r
        
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
        "medians": medians,
        "robustness": robustness
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
        y_data_robustness = [h.get("robustness", {}).get(p) or 0 for h in plot_history]
        has_robustness = any(y > 0 for y in y_data_robustness)
        
        def get_injected_html(current_view_hash, top_href, p=p, x_data=x_data, y_data=y_data, 
                              y_data_robustness=y_data_robustness, has_robustness=has_robustness):
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
            <div style="margin: 0 0 25px 0;">
                <div style="display:flex; justify-content:space-between; align-items:center; padding:12px 20px; background:white; border:1px solid #e2e8f0; border-radius:8px 8px 0 0; font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
                  <a href="{top_href}" style="display:inline-flex; align-items:center; gap:8px; text-decoration:none; color:#1d4ed8; font-weight:600; font-size:14px; padding:6px 14px; background:#eff6ff; border:1px solid #bfdbfe; border-radius:6px; transition:background 0.2s;">
                    &larr; Back to Benchmark Dashboard
                  </a>
                  <span style="font-size:14px; color:#64748b;">
                    Problem: <strong style="color:#0f172a; text-transform:capitalize;">{p}</strong> &bull; Viewing Commit <code style="background:#f1f5f9; padding:2px 6px; border-radius:4px; font-weight:bold; color:#0f172a;">{current_view_hash}</code>
                  </span>
                </div>
                <div style="display:flex; justify-content:center; align-items:center; padding:12px 20px; background:#f8fafc; border:1px solid #e2e8f0; border-top:none; border-radius:0 0 8px 8px; box-shadow:0 4px 6px rgba(0,0,0,0.02); font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
                    <span style="font-size:14px; margin-right:15px; color:#475569; font-weight:600; text-transform:uppercase; letter-spacing:0.5px;">Compare Commits:</span>
                    {btns_str}
                </div>
            </div>
            <!-- INJECTED_NAV_END -->
            """
            
            bottom_ui = f"""
            <!-- INJECTED_UI_START -->
            <div style="margin:40px auto; max-width:1000px; padding:24px; background:white; border:1px solid #e2e8f0; border-radius:10px; box-shadow:0 4px 15px rgba(0,0,0,0.05); font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;">
              <div id="trendPlot" style="width:100%; height:320px;"></div>
            </div>
            <script>
              if(typeof Plotly === 'undefined') {{
                  var s = document.createElement('script');
                  s.src = "../../../toplevel/www/js/plotly-2.32.0.min.js";
                  s.onload = function() {{ renderPlot(); }};
                  document.head.appendChild(s);
              }} else {{ renderPlot(); }}
              
              function renderPlot() {{
                  var traces = [];
                  var hasRobustness = {json.dumps(has_robustness)};
                  
                  traces.push({{
                    x: {json.dumps(x_data)},
                    y: {json.dumps(y_data)},
                    type: 'scatter',
                    mode: 'lines+markers',
                    marker: {{size: 8, color: 'rgba(52, 152, 219, 0.7)'}},
                    line: {{width: 4, color: 'rgba(41, 128, 185, 0.7)'}},
                    name: 'Grand Median Steps',
                    hovertemplate: '<b>Commit %{{x}}</b><br>Grand Median: %{{y:.1f}}<extra></extra>'
                  }});
                  
                  if (hasRobustness) {{
                      traces.push({{
                        x: {json.dumps(x_data)},
                        y: {json.dumps(y_data_robustness)},
                        type: 'scatter',
                        mode: 'lines+markers',
                        marker: {{size: 8, color: 'rgba(230, 126, 34, 0.7)'}},
                        line: {{width: 2.5, color: 'rgba(211, 84, 0, 0.7)'}},
                        name: 'Reliability',
                        yaxis: 'y2',
                        hovertemplate: '<b>Commit %{{x}}</b><br>Reliability: %{{y:.1f}}%<extra></extra>'
                      }});
                  }}
                  
                  var layout = {{
                    title: hasRobustness ? '<b>Trend: Grand Median Steps & Reliability per Commit</b>' : '<b>Trend: Grand Median Steps per Commit</b>',
                    margin: {{ t: 40, b: 50, l: 85, r: hasRobustness ? 85 : 30 }},
                    paper_bgcolor: 'rgba(0,0,0,0)',
                    plot_bgcolor: 'rgba(0,0,0,0)',
                    xaxis: {{ title: 'Commit Hash', type: 'category' }},
                    yaxis: {{ title: {{ text: 'Grand Median Steps', standoff: 20 }}, automargin: true }},
                    showlegend: hasRobustness,
                    legend: hasRobustness ? {{ orientation: 'h', y: -0.2 }} : undefined
                  }};
                  
                  if (hasRobustness) {{
                      layout.yaxis2 = {{
                          title: 'Reliability (%)',
                          overlaying: 'y',
                          side: 'right',
                          range: [0, 105],
                          showgrid: false,
                          automargin: true,
                          titlefont: {{ color: 'rgb(211, 84, 0)' }},
                          tickfont: {{ color: 'rgb(211, 84, 0)' }}
                      }};
                  }}
                  
                  Plotly.newPlot('trendPlot', traces, layout);
              }}
            </script>
            <!-- INJECTED_UI_END -->
            """
            return top_nav, bottom_ui

        dirs_to_process = [("tmp", commit_6)] + [(f"{h}-run", h) for h in valid_hashes]
        for d_name, d_hash in set(dirs_to_process):
            # If historical run, synchronize grandMedian in data.js from history if recorded
            h_entry = next((item for item in history if item.get("commit") == d_hash), None)
            if h_entry and h_entry.get("medians") and h_entry["medians"].get(p) is not None:
                rec_med = h_entry["medians"][p]
                djs_path = os.path.join(p_dir, d_name, "data.js")
                if os.path.exists(djs_path):
                    with open(djs_path, "r") as f:
                        djs = f.read()
                    djs = re.sub(r'const grandMedian\s*=\s*[\d\.]+;', f'const grandMedian = {rec_med};', djs)
                    with open(djs_path, "w") as f:
                        f.write(djs)

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
                
                # Upgrade green median line to shape with continuous anywhere-hover and styled badge
                html = re.sub(r'// UPGRADE_MEDIAN_START.*?// UPGRADE_MEDIAN_END\s*', '', html, flags=re.DOTALL)
                html = re.sub(r'// MEDIAN_HOVER_LISTENER_START.*?// MEDIAN_HOVER_LISTENER_END\s*', '', html, flags=re.DOTALL)
                html = re.sub(r'if \(typeof grandMedian !== \'undefined\' && !traces\.some.*?Plotly\.newPlot\(\'myDiv\'', "Plotly.newPlot('myDiv'", html, flags=re.DOTALL)
                
                upgrade_script = """// UPGRADE_MEDIAN_START
        if (typeof grandMedian !== 'undefined') {
            for (let __i = traces.length - 1; __i >= 0; __i--) {
                if (traces[__i].name === 'Grand Median') {
                    traces.splice(__i, 1);
                }
            }
            if (typeof layout !== 'undefined') {
                layout.shapes = (layout.shapes || []).filter(s => !(s.line && (s.line.color === 'lime' || s.line.color === 'rgb(0, 255, 0)' || s.line.color === '0, 255, 0')));
                layout.shapes.push({
                    type: 'line',
                    xref: 'paper',
                    x0: 0,
                    x1: 1,
                    yref: 'y',
                    y0: grandMedian,
                    y1: grandMedian,
                    line: {
                        color: 'lime',
                        width: 3.5
                    }
                });
                if (layout.annotations) {
                    layout.annotations.forEach(ann => {
                        if (ann.text && (ann.text.includes('Median') || ann.text.includes('median'))) {
                            ann.text = '<b>Grand Median: ' + grandMedian + '</b>';
                            ann.font = { size: 14, color: '#1b5e20' };
                            ann.bgcolor = '#dcfce7';
                            ann.bordercolor = '#22c55e';
                            ann.borderwidth = 1.5;
                        }
                    });
                }
            }
        }
        // UPGRADE_MEDIAN_END
        Plotly.newPlot('myDiv'"""
                html = html.replace("Plotly.newPlot('myDiv'", upgrade_script, 1)

                hover_script = """Plotly.newPlot('myDiv', traces, layout, {responsive: true});
        // MEDIAN_HOVER_LISTENER_START
        (function() {
            const myPlot = document.getElementById('myDiv');
            if (!myPlot) return;
            let pointTip = document.getElementById('pointTooltip');
            if (!pointTip) {
                pointTip = document.createElement('div');
                pointTip.id = 'pointTooltip';
                pointTip.style.position = 'fixed';
                pointTip.style.display = 'none';
                pointTip.style.padding = '6px 12px';
                pointTip.style.borderRadius = '6px';
                pointTip.style.fontSize = '13px';
                pointTip.style.fontFamily = 'monospace';
                pointTip.style.pointerEvents = 'none';
                pointTip.style.zIndex = '99999';
                pointTip.style.boxShadow = '0 4px 12px rgba(0,0,0,0.3)';
                pointTip.style.whiteSpace = 'nowrap';
                document.body.appendChild(pointTip);
            }

            function getMedianPath() {
                const paths = myPlot.querySelectorAll('path');
                for (let i = 0; i < paths.length; i++) {
                    const p = paths[i];
                    const s = (p.getAttribute('style') || '') + ' ' + (p.getAttribute('stroke') || '');
                    if (s.includes('lime') || s.includes('0, 255, 0') || s.includes('0,255,0')) {
                        return p;
                    }
                }
                return null;
            }

            function updateMedianHover(e) {
                if (typeof grandMedian === 'undefined' || !myPlot) return;
                const medianPath = getMedianPath();
                let lineY = null;
                let leftX = 0, rightX = 0;
                if (medianPath) {
                    const r = medianPath.getBoundingClientRect();
                    if (r.width > 0) {
                        lineY = (r.top + r.bottom) / 2;
                        leftX = r.left;
                        rightX = r.right;
                    }
                }
                if (lineY === null && myPlot._fullLayout && myPlot._fullLayout.yaxis) {
                    const b = myPlot.getBoundingClientRect();
                    const fl = myPlot._fullLayout;
                    leftX = b.left + fl.margin.l;
                    rightX = b.left + fl.width - fl.margin.r;
                    lineY = b.top + fl.margin.t + fl.yaxis.c2p(grandMedian);
                }
                if (lineY === null) return;

                const hitZone = 14;
                if (e.clientX >= leftX && e.clientX <= rightX &&
                    e.clientY >= (lineY - hitZone) && e.clientY <= (lineY + hitZone)) {
                    pointTip.innerHTML = '<span style="font-weight:bold; font-size:14px; letter-spacing:0.3px;">Grand Median: ' + grandMedian + '</span>';
                    pointTip.style.left = e.clientX + 'px';
                    pointTip.style.top = (e.clientY - 35) + 'px';
                    pointTip.style.transform = 'none';
                    pointTip.style.background = '#1b5e20';
                    pointTip.style.color = '#ffffff';
                    pointTip.style.border = '1px solid #4ade80';
                    pointTip.style.display = 'block';
                    myPlot.style.cursor = 'pointer';
                } else if (pointTip.innerHTML.includes('Grand Median:')) {
                    pointTip.style.display = 'none';
                    pointTip.style.transform = 'translate(12px, 12px)';
                    pointTip.style.background = 'rgba(33, 37, 41, 0.92)';
                    pointTip.style.border = 'none';
                    myPlot.style.cursor = '';
                }
            }

            myPlot.addEventListener('mousemove', updateMedianHover);
            myPlot.addEventListener('mouseleave', function() {
                if (pointTip.innerHTML.includes('Grand Median:')) {
                    pointTip.style.display = 'none';
                    pointTip.style.transform = 'translate(12px, 12px)';
                    pointTip.style.background = 'rgba(33, 37, 41, 0.92)';
                    pointTip.style.border = 'none';
                    myPlot.style.cursor = '';
                }
            });
        })();
        // MEDIAN_HOVER_LISTENER_END"""
                if "Plotly.newPlot('myDiv', traces, layout, {responsive: true});" in html:
                    html = html.replace("Plotly.newPlot('myDiv', traces, layout, {responsive: true});", hover_script, 1)
                
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
    
    # Helper to calculate estimated time: Grand median steps * complex solutions * 10^-6 (seconds)
    # Then express in microseconds (or ms if >= 1000us)
    def compute_time_str(problem_name, num_solutions):
        med = medians.get(problem_name)
        if not med:
            return ""
        # 1 step ~ 1 microsecond (10^-6 s)
        time_us = med * num_solutions
        if time_us >= 1000:
            return f"Time: ~{time_us / 1000.0:.1f}ms"
        else:
            return f"Time: ~{time_us:.0f}&mu;s"

    chicago_time_str = compute_time_str("chicago", 312)
    linecircle_time_str = compute_time_str("linecircle", 4)

    # 6.a Local Generation (write back to www/index.html)
    local_html = template.replace('{{CHICAGO_HREF}}', '../../individual/chicago-benchmark/tmp/index.html')
    local_html = local_html.replace('{{LINECIRCLE_HREF}}', '../../individual/linecircle-benchmark/tmp/index.html')
    local_html = local_html.replace('{{CHICAGO_RELIABILITY}}', f"{robustness.get('chicago', 0):.1f}")
    local_html = local_html.replace('{{LINECIRCLE_RELIABILITY}}', f"{robustness.get('linecircle', 0):.1f}")
    local_html = local_html.replace('{{CHICAGO_TIME}}', chicago_time_str)
    local_html = local_html.replace('{{LINECIRCLE_TIME}}', linecircle_time_str)
    with open(os.path.join(toplevel_dir, "www", "index.html"), "w") as f:
        f.write(local_html)
        
    # 6.b Public Generation (write to public/index.html)
    pub_html = template.replace('{{CHICAGO_HREF}}', './chicago-benchmark/tmp/index.html')
    pub_html = pub_html.replace('{{LINECIRCLE_HREF}}', './linecircle-benchmark/tmp/index.html')
    pub_html = pub_html.replace('{{CHICAGO_RELIABILITY}}', f"{robustness.get('chicago', 0):.1f}")
    pub_html = pub_html.replace('{{LINECIRCLE_RELIABILITY}}', f"{robustness.get('linecircle', 0):.1f}")
    pub_html = pub_html.replace('{{CHICAGO_TIME}}', chicago_time_str)
    pub_html = pub_html.replace('{{LINECIRCLE_TIME}}', linecircle_time_str)
    with open(os.path.join(public_dir, "index.html"), "w") as f:
        f.write(pub_html)
        
    # Copy figs/ directory to public/figs
    figs_src = os.path.join(toplevel_dir, "www", "figs")
    if os.path.exists(figs_src):
        shutil.copytree(figs_src, os.path.join(public_dir, "figs"))
        
    # Copy js/ directory to public/js
    js_src = os.path.join(toplevel_dir, "www", "js")
    if os.path.exists(js_src):
        shutil.copytree(js_src, os.path.join(public_dir, "js"))
        
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
                    f_html = f_html.replace("../../../toplevel/www/js/", "../../js/")
                    with open(fpath, "w") as f:
                        f.write(f_html)
                        
    print(f"Dashboard packaged successfully for Hash {commit_6}.")
    print(f"Active History (up to 5): {valid_hashes}")

if __name__ == "__main__":
    main()
