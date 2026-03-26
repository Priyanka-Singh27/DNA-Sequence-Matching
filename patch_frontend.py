import os
import re

html_path = r"d:\SEM 4\ADS\DNA\WEB\frontend\webpage.html"
with open(html_path, "r", encoding="utf-8") as f:
    text = f.read()

# 1. Patch Pathogen ID
pathogen_old = """document.getElementById('btn-identify').addEventListener('click', () => {
  const scan = document.getElementById('pathogen-scan');
  scan.classList.add('active');
  const resultsDiv = document.getElementById('pathogen-results');
  resultsDiv.style.display = 'none';
  document.getElementById('match-cards').innerHTML = '';

  setTimeout(()=>{
    scan.classList.remove('active');
    resultsDiv.style.display = 'block';
    PATHOGENS.forEach((p,i) => {
      setTimeout(()=>{
        const card = document.createElement('div');
        card.className = 'match-card';
        card.innerHTML = `
          <div class="match-rank">${i+1}</div>
          <div class="match-info">
            <div class="match-name">${p.name}</div>
            <div class="match-detail">Accession: ${p.seq.substring(0,20)}... · Length: ${p.length} bp · Suffix Tree match</div>
          </div>
          <div style="display:flex;flex-direction:column;align-items:flex-end;gap:8px">
            <div class="match-score">${p.match}%</div>
            <span class="danger-badge danger-${p.danger}">${p.danger.toUpperCase()}</span>
          </div>`;
        document.getElementById('match-cards').appendChild(card);
        setTimeout(()=>card.classList.add('revealed'), 50);
      }, i * 300);
    });
  }, 1800);
});"""

pathogen_new = """document.getElementById('btn-identify').addEventListener('click', () => {
  const scan = document.getElementById('pathogen-scan');
  scan.classList.add('active');
  const resultsDiv = document.getElementById('pathogen-results');
  resultsDiv.style.display = 'none';
  document.getElementById('match-cards').innerHTML = '';
  
  // Actually grab sequence if it exists or fallback
  const seqEl = document.getElementById('pathogen-input-seq'); // Assuming an id
  const seq = seqEl ? seqEl.value.trim().toUpperCase() : 'ATGC';

  fetch('http://localhost:5000/api/pathogen/identify', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ sequence: seq })
  })
  .then(r=>r.json())
  .then(res => {
      scan.classList.remove('active');
      resultsDiv.style.display = 'block';
      let html = '';
      if (res.results && res.results.length > 0) {
          res.results.forEach((p, i) => {
              const cls = p.confidence > 90 ? 'high' : 'medium';
              const name = p.pathogen;
              html += `<div class="match-card revealed">
                <div class="match-rank">${i+1}</div>
                <div class="match-info">
                    <div class="match-name">${name}</div>
                    <div class="match-detail">Suffix Tree Exact Pattern Match</div>
                </div>
                <div style="display:flex;flex-direction:column;align-items:flex-end;gap:8px">
                    <div class="match-score">${p.confidence}%</div>
                    <span class="danger-badge danger-${cls}">${cls.toUpperCase()}</span>
                </div>
              </div>`;
          });
      } else {
          html = '<div class="match-card revealed"><div class="match-info">No pathogen matches found for this sequence.</div></div>';
      }
      document.getElementById('match-cards').innerHTML = html;
  })
  .catch(err => {
      scan.classList.remove('active');
      document.getElementById('match-cards').innerHTML = '<div class="match-card revealed"><div class="match-info" style="color:var(--accent-red)">Backend disconnected. Start Flask app.py!</div></div>';
      resultsDiv.style.display = 'block';
  });
});"""

if pathogen_old in text:
    text = text.replace(pathogen_old, pathogen_new)

# 2. Patch Mutations
mut_old = """document.getElementById('btn-detect').addEventListener('click', () => {
  const seq = document.getElementById('mut-seq').value.trim().toUpperCase().replace(/[^ATGC]/g,'');
  if (!seq) return;

  document.getElementById('mutation-results').style.display = 'block';

  // Render annotated sequence
  const baseColor = c => ({A:'#00F5A0',T:'#00D1FF',G:'#BD00FF',C:'#FFB800'}[c]||'#E0E6ED');
  const mutPositions = new Set(MUTATIONS.map(m=>m.pos));
  let html = '<div style="font-family:var(--font-mono);font-size:13px;letter-spacing:1.5px;line-height:2.5;flex-wrap:wrap;display:flex;gap:1px">';
  seq.split('').forEach((c,i) => {
    const isMut = mutPositions.has(i);
    const bg = isMut ? 'rgba(255,56,92,0.25)' : 'transparent';
    const border = isMut ? '1px solid var(--accent-red)' : '1px solid transparent';
    const glitch = isMut ? 'glitch' : '';
    html += `<span class="${glitch}" style="padding:1px 4px;border-radius:2px;background:${bg};border:${border};color:${baseColor(c)};position:relative" title="Pos ${i}: ${isMut?'MUTATION DETECTED':c}">${c}</span>`;
  });
  html += '</div>';
  document.getElementById('mut-seq-display').innerHTML = html;

  // Table
  document.getElementById('mut-table-body').innerHTML = MUTATIONS.map(m=>`
    <tr>
      <td class="cyan">${m.pos}</td>
      <td class="green">${m.norm}</td>
      <td class="red">${m.mut}</td>
      <td>${m.type}</td>
      <td>${m.gene}</td>
      <td style="color:var(--text-primary);font-weight:600">${m.disease}</td>
      <td class="sev-${m.severity.substring(0,3)}">${m.severity.toUpperCase()}</td>
      <td><a href="${m.clinvar}" target="_blank" style="color:var(--accent-cyan);font-size:11px">${m.snp}</a></td>
    </tr>`).join('');

  // Minimap scroll
  const minimap = document.getElementById('mut-seq-display');
  const cursor = document.getElementById('minimap-cursor');
  document.getElementById('mut-seq-display').addEventListener('scroll', () => {
    const pct = minimap.scrollLeft / (minimap.scrollWidth - minimap.clientWidth);
    cursor.style.left = (pct * 96) + '%';
  });
});"""

mut_new = """document.getElementById('btn-detect').addEventListener('click', () => {
  const seq = document.getElementById('mut-seq').value.trim().toUpperCase().replace(/[^ATGC]/g,'');
  if (!seq) return;

  document.getElementById('mutation-results').style.display = 'block';

  fetch('http://localhost:5000/api/mutations/scan', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ sequence: seq, max_errors: 1 })
  }).then(r=>r.json())
  .then(res => {
      const matches = res.matches || [];
      const mutPositions = new Set(matches.map(m=>m.position));
      const baseColor = c => ({A:'#00F5A0',T:'#00D1FF',G:'#BD00FF',C:'#FFB800'}[c]||'#E0E6ED');
      
      let html = '<div style="font-family:var(--font-mono);font-size:13px;letter-spacing:1.5px;line-height:2.5;flex-wrap:wrap;display:flex;gap:1px">';
      seq.split('').forEach((c,i) => {
        const isMut = mutPositions.has(i);
        const bg = isMut ? 'rgba(255,56,92,0.25)' : 'transparent';
        const border = isMut ? '1px solid var(--accent-red)' : '1px solid transparent';
        const glitch = isMut ? 'glitch' : '';
        html += `<span class="${glitch}" style="padding:1px 4px;border-radius:2px;background:${bg};border:${border};color:${baseColor(c)};position:relative" title="Pos ${i}: ${isMut?'MUTATION DETECTED':c}">${c}</span>`;
      });
      html += '</div>';
      document.getElementById('mut-seq-display').innerHTML = html;

      if (matches.length > 0) {
          document.getElementById('mut-table-body').innerHTML = matches.map(m=>`
            <tr>
              <td class="cyan">${m.position}</td>
              <td class="green">REF</td>
              <td class="red">ALT</td>
              <td>${m.mismatches} err</td>
              <td>Multiple</td>
              <td style="color:var(--text-primary);font-weight:600">${m.disease}</td>
              <td class="sev-${m.severity===2?'high':(m.severity===1?'mid':'low')}">${m.severity===2?'HIGH':'LOW'}</td>
              <td><a href="#" style="color:var(--accent-cyan);font-size:11px">${m.snp_id}</a></td>
            </tr>`).join('');
      } else {
          document.getElementById('mut-table-body').innerHTML = '<tr><td colspan="8">No mutations detected</td></tr>';
      }

      const minimap = document.getElementById('mut-seq-display');
      const cursor = document.getElementById('minimap-cursor');
      minimap.addEventListener('scroll', () => {
        const pct = minimap.scrollLeft / (Math.max(1, minimap.scrollWidth - minimap.clientWidth));
        cursor.style.left = (pct * 96) + '%';
      });
  }).catch(e => alert("API Backend not running. " + e));
});"""

if mut_old in text:
    text = text.replace(mut_old, mut_new)


# 3. Patch Disease Risk
risk_old = """document.getElementById('btn-analyze-risk').addEventListener('click', () => {
  const type = window._riskType || 'brca1';
  const data = SNP_DB[type] || SNP_DB.brca1;
  document.getElementById('risk-results').style.display = 'block';

  // Risk cards
  document.getElementById('risk-cards').innerHTML = data.risks.map(r=>`
    <div class="risk-card ${r.risk}" style="animation:slideIn 0.4s ease forwards">
      <span class="risk-badge">${r.risk === 'none' ? 'NOT DETECTED' : r.risk.toUpperCase()}</span>
      <div class="risk-info">
        <div class="risk-disease">${r.disease}</div>
        <div class="risk-snp">${r.snp} · ${r.mismatches} mismatch${r.mismatches!==1?'es':''} (Trie depth match)</div>
      </div>
    </div>`).join('');

  // Donut chart
  const counts = {high:0,medium:0,low:0,none:0};
  data.risks.forEach(r=>counts[r.risk]++);
  drawDonut(counts);

  // Trie path
  document.getElementById('trie-path-display').innerHTML = data.risks.slice(0,2).map((r,i)=>`
    <span class="cyan">Trie Path ${i+1}:</span> ROOT → ${r.snp.split('').slice(0,3).join(' → ')} → ... → <span class="green">MATCH</span> (${r.mismatches} error${r.mismatches!==1?'s':''}) → <span style="color:var(--text-primary);font-weight:700">${r.disease}</span>`).join('<br>');
});"""

risk_new = """document.getElementById('btn-analyze-risk').addEventListener('click', () => {
  const seq = document.getElementById('risk-seq').value.trim().toUpperCase().replace(/[^ATGC]/g,'');
  if (!seq) return;
  document.getElementById('risk-results').style.display = 'block';

  fetch('http://localhost:5000/api/snp/fuzzy', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ query: seq, max_errors: 1 })
  }).then(r=>r.json())
  .then(res => {
      const risks = res.matches || [];
      const counts = {high:0,medium:0,low:0,none:0};
      const getSeverityString = s => s === 2 ? 'high' : (s === 1 ? 'medium' : 'low');
      
      if (risks.length === 0) counts.none = 1;
      
      document.getElementById('risk-cards').innerHTML = risks.map(r=> {
        const rStr = getSeverityString(r.severity);
        counts[rStr]++;
        return `<div class="risk-card ${rStr}" style="animation:slideIn 0.4s ease forwards">
          <span class="risk-badge">${rStr.toUpperCase()}</span>
          <div class="risk-info">
            <div class="risk-disease">${r.disease}</div>
            <div class="risk-snp">${r.snp_id} · ${r.mismatches} mismatch${r.mismatches!==1?'es':''} (Trie depth match)</div>
          </div>
        </div>`;
      }).join('');
      
      if (risks.length === 0) {
          document.getElementById('risk-cards').innerHTML = `<div class="risk-card none"><span class="risk-badge">NOT DETECTED</span><div class="risk-info"><div class="risk-disease">No known ClinVar SNPs matched</div></div></div>`;
      }

      drawDonut(counts);

      document.getElementById('trie-path-display').innerHTML = risks.slice(0,2).map((r,i)=>`
        <span class="cyan">Trie Path ${i+1}::</span> ROOT → ${r.snp_id} → ... → <span class="green">MATCH</span> (${r.mismatches} error${r.mismatches!==1?'s':''}) → <span style="color:var(--text-primary);font-weight:700">${r.disease}</span>`).join('<br>');
  }).catch(e => alert("API Backend disconnected: " + e));
});"""

if risk_old in text:
    text = text.replace(risk_old, risk_new)

with open(html_path, "w", encoding="utf-8") as f:
    f.write(text)

print("Patch applied to webpage.html successfully.")
