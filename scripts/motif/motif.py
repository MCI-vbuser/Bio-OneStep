#!/usr/bin/env python3
"""
Generate a statistical summary HTML report and SQLite database from MEME, DREME and Tomtom XML outputs.
"""

import argparse
import sqlite3
import json
import os
import xml.etree.ElementTree as ET
from datetime import datetime

# ---------- 解析辅助函数 ----------
def parse_alphabet(xml_root):
    """从 <alphabet> 提取核心字母列表 (如 ['A','C','G','T'])"""
    alph = xml_root.find(".//alphabet")
    if alph is None:
        return []
    return [l.get('symbol') for l in alph.findall(".//letter")[:4]]  # 取前4个核心字母

def parse_background(xml_root, tag='background_frequencies'):
    """提取背景频率 (A,C,G,T) 列表"""
    bg = xml_root.find(f".//{tag}")
    if bg is None:
        return []
    arr = bg.find(".//alphabet_array")
    if arr is None:
        return []
    freqs = [float(v.text) for v in arr.findall(".//value")]
    return freqs

def parse_motif_pspm(motif_elem):
    """从 <motif> 中提取 PSPM 矩阵 (列表的列表)"""
    probs = motif_elem.find(".//probabilities")
    if probs is None:
        return None
    matrix = []
    for pos in probs.findall(".//alphabet_array"):
        row = [float(v.text) for v in pos.findall(".//value")]
        matrix.append(row)
    return matrix

def parse_meme(xml_path):
    """解析 meme.xml，返回数据结构"""
    tree = ET.parse(xml_path)
    root = tree.getroot()

    model = root.find(".//model")
    cmd = model.findtext("command_line", "").strip()
    version = root.get("version", "")
    release = root.get("release", "")

    training = root.find(".//training_set")
    primary_seq = training.get("primary_sequences", "")
    primary_count = int(training.get("primary_count", 0))
    primary_positions = int(training.get("primary_positions", 0))
    control_source = training.get("control_sequences", "")
    control_count = int(training.get("control_count", 0))

    bg = parse_background(root, "background_frequencies")
    seq_freq_elem = root.find(".//letter_frequencies/alphabet_array")
    seq_freqs = [float(v.text) for v in seq_freq_elem.findall(".//value")] if seq_freq_elem is not None else []

    alphabet = parse_alphabet(root)

    motifs = []
    for motif in root.findall(".//motifs/motif"):
        m = {
            "id": motif.get("id"),
            "alt": motif.get("alt"),
            "name": motif.get("name"),
            "width": int(motif.get("width")),
            "sites": int(motif.get("sites")),
            "ic": float(motif.get("ic")),
            "re": float(motif.get("re")),
            "llr": float(motif.get("llr")),
            "p_value": float(motif.get("p_value")),
            "e_value": motif.get("e_value"),
            "bayes_threshold": float(motif.get("bayes_threshold")),
            "consensus": motif.findtext("regular_expression", "").strip(),
            "pspm": parse_motif_pspm(motif)
        }
        motifs.append(m)

    return {
        "program": "MEME",
        "version": version,
        "release": release,
        "command_line": cmd,
        "primary_source": primary_seq,
        "primary_count": primary_count,
        "primary_positions": primary_positions,
        "control_source": control_source,
        "control_count": control_count,
        "background": bg,
        "seq_freqs": seq_freqs,
        "alphabet": alphabet,
        "motifs": motifs
    }

def parse_dreme(xml_path):
    """解析 dreme.xml"""
    tree = ET.parse(xml_path)
    root = tree.getroot()

    model = root.find(".//model")
    cmd = model.findtext("command_line", "").strip()
    version = root.get("version", "")
    release = root.get("release", "")

    pos = model.find(".//positives")
    neg = model.find(".//negatives")
    primary_source = pos.get("file", "")
    primary_count = int(pos.get("count", 0))
    control_source = neg.get("from", "shuffled")
    control_count = int(neg.get("count", 0))

    bg_attrs = model.find(".//background")
    if bg_attrs is not None:
        bg = [float(bg_attrs.get(k, 0)) for k in ('A','C','G','T')]
    else:
        bg = []

    alphabet = parse_alphabet(root)

    motifs = []
    for motif in root.findall(".//motifs/motif"):
        m = {
            "id": motif.get("id"),
            "alt": motif.get("alt"),
            "seq": motif.get("seq"),
            "length": int(motif.get("length")),
            "nsites": int(motif.get("nsites")),
            "p": int(motif.get("p")),
            "n": int(motif.get("n")),
            "pvalue": motif.get("pvalue"),
            "evalue": motif.get("evalue"),
            "unerased_evalue": motif.get("unerased_evalue"),
            "pspm": []
        }
        for pos in motif.findall(".//pos"):
            row = [float(pos.get('A',0)), float(pos.get('C',0)), float(pos.get('G',0)), float(pos.get('T',0))]
            m["pspm"].append(row)
        motifs.append(m)

    return {
        "program": "DREME",
        "version": version,
        "release": release,
        "command_line": cmd,
        "primary_source": primary_source,
        "primary_count": primary_count,
        "control_source": control_source,
        "control_count": control_count,
        "background": bg,
        "alphabet": alphabet,
        "motifs": motifs
    }

def parse_tomtom(xml_path):
    """解析 tomtom.xml"""
    tree = ET.parse(xml_path)
    root = tree.getroot()

    model = root.find(".//model")
    cmd = model.findtext("command_line", "").strip()
    version = root.get("version", "")
    release = root.get("release", "")

    bg = parse_background(root, "background")
    alphabet = parse_alphabet(root)

    query_dbs = []
    for db in root.findall(".//query_dbs/db"):
        query_dbs.append({
            "source": db.get("source"),
            "name": db.get("name"),
            "loaded": int(db.get("loaded")),
            "excluded": int(db.get("excluded")),
            "last_modified": db.get("last_mod_date")
        })

    target_dbs = []
    for db in root.findall(".//target_dbs/db"):
        target_dbs.append({
            "source": db.get("source"),
            "name": db.get("name"),
            "loaded": int(db.get("loaded")),
            "excluded": int(db.get("excluded")),
            "last_modified": db.get("last_mod_date")
        })

    queries = []
    for q in root.findall(".//queries/motif"):
        qdata = {
            "db": int(q.get("db")),
            "id": q.get("id"),
            "alt": q.get("alt"),
            "length": int(q.get("length")),
            "nsites": int(q.get("nsites")),
            "evalue": q.get("evalue"),
            "pspm": []
        }
        for pos in q.findall(".//pos"):
            row = [float(pos.get('A',0)), float(pos.get('C',0)), float(pos.get('G',0)), float(pos.get('T',0))]
            qdata["pspm"].append(row)
        queries.append(qdata)

    targets = []
    for t in root.findall(".//targets/motif"):
        tdata = {
            "db": int(t.get("db")),
            "id": t.get("id"),
            "alt": t.get("alt"),
            "length": int(t.get("length")),
            "nsites": int(t.get("nsites")),
            "pspm": []
        }
        for pos in t.findall(".//pos"):
            row = [float(pos.get('A',0)), float(pos.get('C',0)), float(pos.get('G',0)), float(pos.get('T',0))]
            tdata["pspm"].append(row)
        targets.append(tdata)

    matches = []
    for qmatch in root.findall(".//matches/query"):
        qidx = int(qmatch.get("idx"))
        for t in qmatch.findall(".//target"):
            matches.append({
                "query_idx": qidx,
                "target_idx": int(t.get("idx")),
                "rc": t.get("rc") == "y",
                "offset": int(t.get("off")),
                "pvalue": float(t.get("pv")),
                "evalue": float(t.get("ev")),
                "qvalue": float(t.get("qv"))
            })

    return {
        "program": "Tomtom",
        "version": version,
        "release": release,
        "command_line": cmd,
        "background": bg,
        "alphabet": alphabet,
        "query_dbs": query_dbs,
        "target_dbs": target_dbs,
        "queries": queries,
        "targets": targets,
        "matches": matches
    }

# ---------- 数据库操作 ----------
def init_db(db_path):
    conn = sqlite3.connect(db_path)
    c = conn.cursor()

    c.execute('''
        CREATE TABLE IF NOT EXISTS program_info (
            id INTEGER PRIMARY KEY,
            program TEXT,
            version TEXT,
            release TEXT,
            command_line TEXT,
            primary_source TEXT,
            primary_count INTEGER,
            control_source TEXT,
            control_count INTEGER,
            background TEXT,  -- JSON 数组
            alphabet TEXT      -- JSON 数组
        )
    ''')

    c.execute('''
        CREATE TABLE IF NOT EXISTS motifs (
            id INTEGER PRIMARY KEY,
            program TEXT,
            motif_id TEXT,
            alt TEXT,
            length INTEGER,
            nsites INTEGER,
            evalue TEXT,
            pvalue TEXT,
            consensus TEXT,
            pspm TEXT         -- JSON 二维数组
        )
    ''')

    c.execute('''
        CREATE TABLE IF NOT EXISTS tomtom_queries (
            id INTEGER PRIMARY KEY,
            program TEXT,
            motif_id TEXT,
            alt TEXT,
            length INTEGER,
            nsites INTEGER,
            evalue TEXT,
            pspm TEXT
        )
    ''')

    c.execute('''
        CREATE TABLE IF NOT EXISTS tomtom_targets (
            id INTEGER PRIMARY KEY,
            program TEXT,
            motif_id TEXT,
            alt TEXT,
            length INTEGER,
            nsites INTEGER,
            pspm TEXT
        )
    ''')

    c.execute('''
        CREATE TABLE IF NOT EXISTS tomtom_matches (
            id INTEGER PRIMARY KEY,
            query_id INTEGER,
            target_id INTEGER,
            rc INTEGER,
            offset INTEGER,
            pvalue REAL,
            evalue REAL,
            qvalue REAL,
            FOREIGN KEY(query_id) REFERENCES tomtom_queries(id),
            FOREIGN KEY(target_id) REFERENCES tomtom_targets(id)
        )
    ''')

    conn.commit()
    return conn

def insert_program_info(conn, prog_data):
    c = conn.cursor()
    c.execute('''
        INSERT INTO program_info
        (program, version, release, command_line, primary_source, primary_count,
         control_source, control_count, background, alphabet)
        VALUES (?,?,?,?,?,?,?,?,?,?)
    ''', (
        prog_data.get("program"),
        prog_data.get("version"),
        prog_data.get("release"),
        prog_data.get("command_line"),
        prog_data.get("primary_source"),
        prog_data.get("primary_count"),
        prog_data.get("control_source"),
        prog_data.get("control_count"),
        json.dumps(prog_data.get("background", [])),
        json.dumps(prog_data.get("alphabet", []))
    ))
    conn.commit()

def insert_motifs(conn, prog_name, motifs):
    c = conn.cursor()
    for m in motifs:
        c.execute('''
            INSERT INTO motifs
            (program, motif_id, alt, length, nsites, evalue, pvalue, consensus, pspm)
            VALUES (?,?,?,?,?,?,?,?,?)
        ''', (
            prog_name,
            m.get("id") or m.get("seq"),
            m.get("alt"),
            m.get("width") or m.get("length"),
            m.get("sites") or m.get("nsites"),
            m.get("e_value") or m.get("evalue"),
            m.get("p_value") or m.get("pvalue"),
            m.get("consensus") or m.get("seq", ""),
            json.dumps(m.get("pspm"))
        ))
    conn.commit()

def insert_tomtom(conn, tomtom_data):
    c = conn.cursor()
    query_id_map = {}
    for i, q in enumerate(tomtom_data["queries"]):
        c.execute('''
            INSERT INTO tomtom_queries
            (program, motif_id, alt, length, nsites, evalue, pspm)
            VALUES (?,?,?,?,?,?,?)
        ''', (
            tomtom_data["program"],
            q["id"],
            q.get("alt"),
            q["length"],
            q["nsites"],
            q["evalue"],
            json.dumps(q["pspm"])
        ))
        query_id_map[i] = c.lastrowid

    target_id_map = {}
    for i, t in enumerate(tomtom_data["targets"]):
        c.execute('''
            INSERT INTO tomtom_targets
            (program, motif_id, alt, length, nsites, pspm)
            VALUES (?,?,?,?,?,?)
        ''', (
            tomtom_data["program"],
            t["id"],
            t.get("alt"),
            t["length"],
            t["nsites"],
            json.dumps(t["pspm"])
        ))
        target_id_map[i] = c.lastrowid

    for m in tomtom_data["matches"]:
        c.execute('''
            INSERT INTO tomtom_matches
            (query_id, target_id, rc, offset, pvalue, evalue, qvalue)
            VALUES (?,?,?,?,?,?,?)
        ''', (
            query_id_map[m["query_idx"]],
            target_id_map[m["target_idx"]],
            1 if m["rc"] else 0,
            m["offset"],
            m["pvalue"],
            m["evalue"],
            m["qvalue"]
        ))
    conn.commit()

# ---------- 统计摘要生成 ----------
def gather_stats(conn):
    """从数据库收集统计摘要"""
    c = conn.cursor()
    stats = {}

    # 程序概览
    prog_rows = c.execute("SELECT program, primary_count, control_count FROM program_info").fetchall()
    stats["programs"] = [{"name": p[0], "primary": p[1], "control": p[2]} for p in prog_rows]

    # MEME 和 DREME 的 motif 统计
    motifs = c.execute("SELECT program, motif_id, length, nsites, evalue, pvalue, consensus FROM motifs").fetchall()
    # 按程序分组
    prog_motifs = {}
    for m in motifs:
        prog = m[0]
        if prog not in prog_motifs:
            prog_motifs[prog] = []
        prog_motifs[prog].append({
            "id": m[1],
            "len": m[2],
            "sites": m[3],
            "evalue": m[4],
            "pvalue": m[5],
            "consensus": m[6]
        })
    stats["motifs_by_prog"] = prog_motifs

    # 每个程序的 motif 数量
    stats["motif_counts"] = {prog: len(lst) for prog, lst in prog_motifs.items()}

    # 最佳 E-value (最小)
    stats["best_evalue"] = {}
    for prog, lst in prog_motifs.items():
        # 将 evalue 字符串转换为浮点数比较（处理科学计数法）
        def parse_e(e):
            try:
                return float(e)
            except:
                return float('inf')
        best = min(lst, key=lambda x: parse_e(x["evalue"]))
        stats["best_evalue"][prog] = best["evalue"]

    # 长度分布
    stats["len_dist"] = {}
    for prog, lst in prog_motifs.items():
        dist = {}
        for m in lst:
            l = m["len"]
            dist[l] = dist.get(l, 0) + 1
        stats["len_dist"][prog] = sorted(dist.items())

    # Tomtom 统计
    # 每个查询 motif 的匹配数
    q_matches = c.execute("""
        SELECT q.motif_id, q.alt, COUNT(m.id)
        FROM tomtom_queries q
        LEFT JOIN tomtom_matches m ON q.id = m.query_id
        GROUP BY q.id
    """).fetchall()
    stats["tomtom_query_counts"] = [{"id": r[0], "alt": r[1], "count": r[2]} for r in q_matches]

    # 最佳匹配（最小 p-value）每个查询
    best_matches = c.execute("""
        SELECT q.motif_id, t.motif_id, m.pvalue, m.evalue
        FROM tomtom_matches m
        JOIN tomtom_queries q ON q.id = m.query_id
        JOIN tomtom_targets t ON t.id = m.target_id
        GROUP BY q.id
        HAVING m.pvalue = MIN(m.pvalue)
    """).fetchall()
    stats["tomtom_best"] = [{"query": b[0], "target": b[1], "pvalue": b[2], "evalue": b[3]} for b in best_matches]

    # 背景频率（取第一个程序的背景，假设一致）
    bg_json = c.execute("SELECT background FROM program_info LIMIT 1").fetchone()
    if bg_json and bg_json[0]:
        stats["background"] = json.loads(bg_json[0])
    else:
        stats["background"] = []

    return stats

# ---------- HTML 生成 ----------
def generate_html(db_path, output_html):
    conn = sqlite3.connect(db_path)
    stats = gather_stats(conn)
    conn.close()

    template = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>MEME Suite - Statistical Summary</title>
    <style>
        :root {
            --bg-primary: #0d1117;
            --bg-secondary: #161b22;
            --text-primary: #c9d1d9;
            --text-secondary: #8b949e;
            --border: #30363d;
            --accent: #58a6ff;
            --card-bg: #21262d;
        }
        body {
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
            background-color: var(--bg-primary);
            color: var(--text-primary);
            line-height: 1.5;
            margin: 0;
            padding: 2rem;
        }
        .container {
            max-width: 1200px;
            margin: 0 auto;
        }
        h1 {
            color: #c9d1d9; /* 改为与正文一致 */
            border-bottom: 1px solid var(--border);
            padding-bottom: 0.3rem;
        }
        h2, h3 {
            color: var(--accent);
            border-bottom: 1px solid var(--border);
            padding-bottom: 0.3rem;
        }
        .cards {
            display: flex;
            flex-wrap: wrap;
            gap: 1rem;
            margin: 1rem 0;
        }
        .card {
            background: var(--card-bg);
            border: 1px solid var(--border);
            border-radius: 6px;
            padding: 1rem;
            flex: 1 1 250px;
        }
        .card h3 {
            margin-top: 0;
            border-bottom: none;
            color: var(--text-primary);
        }
        .card .stat {
            font-size: 2rem;
            font-weight: bold;
            color: var(--accent);
        }
        .stat-label {
            color: var(--text-secondary);
            font-size: 0.9rem;
        }
        table {
            width: 100%;
            border-collapse: collapse;
            margin: 1rem 0;
            background: var(--bg-secondary);
            border: 1px solid var(--border);
            border-radius: 6px;
            overflow: hidden;
        }
        th {
            background: var(--table-header, #161b22);
            color: var(--text-primary);
            font-weight: 600;
            text-align: left;
            padding: 0.75rem 1rem;
            border-bottom: 1px solid var(--border);
        }
        td {
            padding: 0.75rem 1rem;
            border-bottom: 1px solid var(--border);
            color: var(--text-secondary);
        }
        tr:last-child td {
            border-bottom: none;
        }
        tr:nth-child(even) {
            background: #0d1117;
        }
        tr:nth-child(odd) {
            background: #161b22;
        }
        .badge {
            display: inline-block;
            padding: 0.2rem 0.6rem;
            border-radius: 12px;
            font-size: 0.85rem;
            background: var(--accent);
            color: #000;
            font-weight: 600;
        }
        .evalue {
            font-family: monospace;
        }
        .dist-bar {
            display: flex;
            align-items: center;
            margin: 0.25rem 0;
        }
        .dist-label {
            width: 40px;
        }
        .dist-bar-fill {
            height: 20px;
            background: var(--accent);
            border-radius: 4px;
            margin-left: 10px;
            min-width: 2px;
            max-width: 200px; /* 限制条形图最大宽度 */
        }
        .footer {
            margin-top: 3rem;
            text-align: center;
            color: var(--text-secondary);
            font-size: 0.9rem;
        }
        a {
            color: var(--accent);
            text-decoration: none;
        }
        a:hover {
            text-decoration: underline;
        }
    </style>
</head>
<body>
<div class="container">
    <h1>MEME Suite - Statistical Summary</h1>

    <!-- 概览卡片 -->
    <h2>Program Overview</h2>
    <div class="cards">
        {% for prog in stats.programs %}
        <div class="card">
            <h3>{{ prog.name }}</h3>
            <div class="stat">{{ stats.motif_counts[prog.name] or 0 }}</div>
            <div class="stat-label">motifs found</div>
            <p>Primary: {{ prog.primary }} sequences</p>
            <p>Control: {{ prog.control }} sequences</p>
            <p>Best E-value: <span class="evalue">{{ stats.best_evalue[prog.name] }}</span></p>
        </div>
        {% endfor %}
    </div>

    <!-- 背景频率和长度分布并排，长度分布左，背景频率右 -->
    <div style="display: flex; gap: 2rem; align-items: flex-start;">
        <!-- 左侧：Motif Length Distribution -->
        <div style="flex: 1; max-width: 600px;">
            <h2>Motif Length Distribution</h2>
            <div style="display: flex; flex-direction: column; gap: 1rem;">
                {% for prog, dist in stats.len_dist.items() %}
                <div class="card" style="margin: 0; max-height: 100px; overflow-y: auto;">
                    <h4>{{ prog }}</h4>
                    {% for length, count in dist %}
                    <div class="dist-bar">
                        <span class="dist-label">{{ length }}</span>
                        <div class="dist-bar-fill" style="width: {{ count * 20 }}px;"></div>
                        <span style="margin-left: 8px;">({{ count }})</span>
                    </div>
                    {% endfor %}
                </div>
                {% endfor %}
            </div>
        </div>
        <!-- 右侧：Background Frequencies -->
        <div style="flex: 0 0 auto;">
            <h2>Background Frequencies</h2>
            <table style="width: auto;">
                <tr><th>Base</th><th>Frequency</th></tr>
                {% for base, freq in zip(['A','C','G','T'], stats.background) %}
                <tr><td>{{ base }}</td><td>{{ "%.3f"|format(freq) }}</td></tr>
                {% endfor %}
            </table>
        </div>
    </div>

    <!-- 最佳 Motifs (Top 5 by E-value) -->
    <h2>Top Motifs by E-value</h2>
    {% for prog, motifs in stats.motifs_by_prog.items() %}
    <h3>{{ prog }}</h3>
    <table>
        <thead><tr><th>ID</th><th>Consensus</th><th>Length</th><th>Sites</th><th>E-value</th></tr></thead>
        <tbody>
        {% for m in motifs[:5] %}
        <tr>
            <td>{{ m.id }}</td>
            <td><code>{{ m.consensus }}</code></td>
            <td>{{ m.len }}</td>
            <td>{{ m.sites }}</td>
            <td class="evalue">{{ m.evalue }}</td>
        </tr>
        {% endfor %}
        </tbody>
    </table>
    {% endfor %}

    <!-- Tomtom Matches -->
    <h2>Tomtom Matches</h2>
    <table>
        <thead><tr><th>Query Motif</th><th>Matches</th><th>Best Target</th><th>Best p-value</th></tr></thead>
        <tbody>
        {% for item in stats.tomtom_query_counts %}
        <tr>
            <td>{{ item.id }} {% if item.alt %}({{ item.alt }}){% endif %}</td>
            <td>{{ item.count }}</td>
            <td>
                {% set best = stats.tomtom_best | selectattr("query", "equalto", item.id) | first %}
                {% if best %}{{ best.target }} ({{ "%.2e"|format(best.pvalue) }}){% endif %}
            </td>
            <td>{% if best %}{{ "%.2e"|format(best.pvalue) }}{% endif %}</td>
        </tr>
        {% endfor %}
        </tbody>
    </table>
</div>
</body>
</html>
    """

    from datetime import datetime
    try:
        from jinja2 import Template
        tpl = Template(template)
        html = tpl.render(
            date=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            stats=stats,
            zip=zip,
            selectattr=lambda seq, attr, value: [x for x in seq if x.get(attr) == value]
        )
    except ImportError:
        # 简化输出（不依赖 Jinja2）
        html = "<h1>Install jinja2 for full HTML rendering</h1><pre>" + json.dumps(stats, indent=2) + "</pre>"

    with open(output_html, 'w', encoding='utf-8') as f:
        f.write(html)

# ---------- 主程序 ----------
def main():
    parser = argparse.ArgumentParser(description="Generate statistical summary from MEME, DREME and Tomtom XML.")
    parser.add_argument('--meme', required=True, help='Path to meme.xml')
    parser.add_argument('--dreme', required=True, help='Path to dreme.xml')
    parser.add_argument('--tomtom', required=True, help='Path to tomtom.xml')
    parser.add_argument('--output-dir', default='.', help='Output directory for database and HTML')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    db_path = os.path.join(args.output_dir, 'motifs.db')
    html_path = os.path.join(args.output_dir, 'motif_report.html')

    print("Parsing MEME...")
    meme_data = parse_meme(args.meme)
    print("Parsing DREME...")
    dreme_data = parse_dreme(args.dreme)
    print("Parsing Tomtom...")
    tomtom_data = parse_tomtom(args.tomtom)

    conn = init_db(db_path)
    insert_program_info(conn, meme_data)
    insert_program_info(conn, dreme_data)
    insert_motifs(conn, "MEME", meme_data["motifs"])
    insert_motifs(conn, "DREME", dreme_data["motifs"])
    insert_tomtom(conn, tomtom_data)
    conn.close()
    print(f"Database written to {db_path}")

    generate_html(db_path, html_path)
    print(f"HTML report written to {html_path}")

if __name__ == '__main__':
    main()