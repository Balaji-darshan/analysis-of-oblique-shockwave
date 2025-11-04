<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>Oblique Shock Wave Analyzer</title>
  <script src="https://cdn.plot.ly/plotly-2.25.2.min.js"></script>
  <style>
    body{font-family:system-ui, -apple-system, 'Segoe UI', Roboto, 'Helvetica Neue', Arial; margin:20px; background:#f7f9fc;color:#222}
    .card{background:#fff;padding:18px;border-radius:10px;box-shadow:0 6px 18px rgba(45,55,72,0.08);margin-bottom:16px}
    h1{margin:0 0 8px 0;font-size:20px}
    label{display:block;margin-top:8px;font-weight:600}
    input[type=number]{width:140px;padding:8px;border-radius:6px;border:1px solid #d6dbe6}
    button{margin-top:12px;padding:10px 14px;border-radius:8px;border:0;background:#2563eb;color:#fff;font-weight:700;cursor:pointer}
    .flex{display:flex;gap:16px;align-items:center}
    .results{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:10px}
    pre{background:#0f172a;color:#e6eef8;padding:10px;border-radius:6px;overflow:auto}
    small{color:#556}
  </style>
</head>
<body>
  <div class="card">
    <h1>Oblique Shock Wave Analyzer</h1>
    <p>Enter upstream Mach number (M₁) and wedge/deflection angle θ (degrees). This tool computes weak & strong oblique shock solutions, downstream Mach number (M₂), pressure, temperature and density ratios, and plots the θ–β–M curves for selected Mach numbers.</p>

    <div class="flex">
      <div>
        <label>M₁ (upstream Mach)</label>
        <input id="machInput" type="number" value="2" step="0.01" min="1.01" />
        <label>θ (wedge/deflection) — degrees</label>
        <input id="thetaInput" type="number" value="10" step="0.1" />
        <label>γ (ratio of specific heats)</label>
        <input id="gammaInput" type="number" value="1.4" step="0.01" />
        <div style="margin-top:8px"><button onclick="runCalc()">Calculate</button></div>
      </div>

      <div style="flex:1">
        <label>Results</label>
        <div class="results" id="resultsBox">
          <!-- populated by JS -->
        </div>
      </div>
    </div>
  </div>

  <div class="card">
    <h2>θ–β–M Plot (multiple Mach numbers)</h2>
    <div id="plot" style="height:520px;width:100%"></div>
    <small>Curves show θ (deg) vs β (deg) for different upstream Mach numbers. A point is marked for the computed β (weak/strong) for the chosen M & θ. If θ &gt; θ<sub>max</sub> for that M, the shock detaches (bow shock).</small>
  </div>

  <script>
    // ---------- Utility functions ----------
    const toRad = d => d*Math.PI/180;
    const toDeg = r => r*180/Math.PI;

    // θ–β–M relation: given M and beta returns theta (radians)
    // tan(theta) = 2 cot(beta) * (M^2 sin^2 beta -1) / (M^2(γ + cos(2β)) + 2)
    function thetaFromBeta(M, beta, gamma){
      const sinb = Math.sin(beta);
      const cos2b = Math.cos(2*beta);
      const numerator = 2 * (Math.cos(beta)/Math.sin(beta)) * (M*M*sinb*sinb - 1);
      const denominator = M*M*(gamma + cos2b) + 2;
      const tanTheta = numerator/denominator;
      return Math.atan(tanTheta); // radians (may be negative if no solution)
    }

    // Normal shock relations using Mn1 (normal component)
    function normalShockRelations(Mn1, gamma){
      const rho_ratio = ((gamma+1)*Mn1*Mn1)/((gamma-1)*Mn1*Mn1 + 2);
      const p_ratio = 1 + (2*gamma/(gamma+1))*(Mn1*Mn1 - 1);
      const T_ratio = p_ratio / rho_ratio;
      const Mn2sq = (1 + ((gamma-1)/2)*Mn1*Mn1) / (gamma*Mn1*Mn1 - ((gamma-1)/2));
      const Mn2 = Math.sqrt(Math.max(Mn2sq, 1e-12));
      return {rho_ratio, p_ratio, T_ratio, Mn2};
    }

    // Given M and theta (rad), find beta solutions (weak & strong) using scanning + bisection
    function solveBetaForTheta(M, theta, gamma){
      // Beta limits: from asin(1/M) to pi/2
      const betaMin = Math.asin(1/M) + 1e-8;
      const betaMax = Math.PI/2 - 1e-8;
      function f(beta){
        return thetaFromBeta(M, beta, gamma) - theta;
      }
      const N = 400; // scan
      const betas = [];
      const fvals = [];
      for(let i=0;i<=N;i++){
        const b = betaMin + (betaMax - betaMin)*i/N;
        betas.push(b);
        fvals.push(f(b));
      }
      // find sign changes (roots)
      const roots = [];
      for(let i=0;i<N;i++){
        if(fvals[i] === 0 || (fvals[i]*fvals[i+1] < 0)){
          // bisection between betas[i] and betas[i+1]
          let a = betas[i], c = betas[i+1];
          let fa = f(a), fc = f(c);
          // refine
          for(let k=0;k<60;k++){
            const m = 0.5*(a+c);
            const fm = f(m);
            if(fa*fm <= 0){ c = m; fc = fm; } else { a = m; fa = fm; }
          }
          const root = 0.5*(a+c);
          // avoid duplicates (close roots)
          if(!roots.some(r=>Math.abs(r-root) < 1e-6)) roots.push(root);
        }
      }
      // sort: weak is lower beta, strong is higher beta
      roots.sort((a,b)=>a-b);
      return roots; // array of beta radians (0,1 or 2 roots)
    }

    function computeOblique(M, thetaDeg, gamma){
      const theta = toRad(thetaDeg);
      if(M<=1){ return {error: 'M must be > 1 (supersonic upstream)'} }
      const betas = solveBetaForTheta(M, theta, gamma);
      // compute for each beta: Mn1, normal shock relations, M2, p_ratio, T_ratio, rho_ratio
      const solutions = betas.map(beta=>{
        const Mn1 = M * Math.sin(beta);
        const {rho_ratio, p_ratio, T_ratio, Mn2} = normalShockRelations(Mn1, gamma);
        // downstream Mach number M2 from Mn2 and geometric relation: M2 = Mn2 / sin(beta - theta)
        const denom = Math.sin(beta - theta);
        const M2 = denom>0 ? Mn2/denom : NaN;
        return {beta, Mn1, Mn2, M2, p_ratio, T_ratio, rho_ratio};
      });

      // Also compute theta_max by sampling beta
      let thetaMax = -Infinity, betaAtThetaMax = null;
      const betaStart = Math.asin(1/M) + 1e-6;
      for(let b=betaStart;b<Math.PI/2-1e-6;b+=0.0008){
        const th = thetaFromBeta(M,b,gamma);
        if(th>thetaMax){ thetaMax = th; betaAtThetaMax = b; }
      }
      return {solutions, thetaMax: toDeg(thetaMax), betaAtThetaMax, betasFound: betas.map(toDeg)};
    }

    // ---------- Plotting function ----------
    function plotThetaBeta(selectedM, selectedThetaDeg, gamma){
      const Ms = [1.5,2,3,5,10];
      // ensure selectedM included
      if(!Ms.includes(selectedM)) Ms.unshift(selectedM);
      const traces = [];
      Ms.forEach(M => {
        const betaMin = Math.asin(1/M) + 1e-6;
        const betaMax = Math.PI/2 - 1e-6;
        const B = [];
        const TH = [];
        for(let b=betaMin;b<=betaMax;b+=0.0008){
          const t = thetaFromBeta(M,b,gamma);
          if(t>0){ B.push(toDeg(b)); TH.push(toDeg(t)); }
        }
        traces.push({x:B, y:TH, mode:'lines', name:`M=${M}`, hoverinfo:'name+x+y'});
        // theta max marker
        if(B.length>0){
          const maxTH = Math.max(...TH);
          const idx = TH.indexOf(maxTH);
          traces.push({x:[B[idx]], y:[TH[idx]], mode:'markers', marker:{size:6}, showlegend:false, hoverinfo:'x+y', text:[`θ_max=${TH[idx].toFixed(2)}°`]});
        }
      });

      // compute selected solution markers
      const res = computeOblique(selectedM, selectedThetaDeg, gamma);
      if(res.solutions.length>0){
        res.solutions.forEach((sol,i)=>{
          traces.push({x:[toDeg(sol.beta)], y:[selectedThetaDeg], mode:'markers+text', marker:{size:10}, name:i===0? 'weak solution':'strong solution', text:[i===0? 'weak':'strong'], textposition:'top center'});
        });
      }

      const layout = {
        xaxis:{title:'β (degrees)', range:[0,90]},
        yaxis:{title:'θ (degrees)', range:[0,45]},
        title:`θ–β curves (γ=${gamma}) — selected M=${selectedM}, θ=${selectedThetaDeg}°`,
        hovermode:'closest',
        margin:{l:60,r:30,t:50,b:60}
      };
      Plotly.newPlot('plot', traces, layout, {responsive:true});
    }

    // ---------- UI and wiring ----------
    function displayResults(obj){
      const box = document.getElementById('resultsBox');
      box.innerHTML = '';
      if(obj.error){ box.innerHTML = `<div style="color:#b91c1c;">${obj.error}</div>`; return; }
      if(obj.solutions.length===0){
        box.innerHTML = `<div><strong>No attached oblique shock solution</strong><div>θ_max for M=${obj.M}: ${obj.thetaMax ? obj.thetaMax.toFixed(3) + '°' : '—'}</div><small>If θ &gt; θ_max, the shock detaches (bow shock).</small></div>`;
        return;
      }
      obj.solutions.forEach((s,idx)=>{
        const title = idx===0 ? 'Weak solution (common)' : 'Strong solution (usually not observed)';
        const html = `
          <div class="card" style="padding:10px">
            <strong>${title}</strong>
            <div>β = ${toDeg(s.beta).toFixed(4)}°</div>
            <div>M₂ = ${isFinite(s.M2) ? s.M2.toFixed(4) : '—'}</div>
            <div>p₂/p₁ = ${s.p_ratio.toFixed(4)}</div>
            <div>T₂/T₁ = ${s.T_ratio.toFixed(4)}</div>
            <div>ρ₂/ρ₁ = ${s.rho_ratio.toFixed(4)}</div>
          </div>
        `;
        const el = document.createElement('div'); el.innerHTML = html; box.appendChild(el);
      });
      // theta max info
      const extra = document.createElement('div');
      extra.innerHTML = `<div style="grid-column:1/-1;margin-top:6px"><strong>θ_max for M=${obj.M}:</strong> ${obj.thetaMax.toFixed(4)}° (β at θ_max = ${toDeg(obj.betaAtThetaMax).toFixed(4)}°)</div>`;
      box.appendChild(extra);
    }

    function runCalc(){
      const M = Number(document.getElementById('machInput').value);
      const thetaDeg = Number(document.getElementById('thetaInput').value);
      const gamma = Number(document.getElementById('gammaInput').value) || 1.4;
      const result = computeOblique(M, thetaDeg, gamma);
      result.M = M; result.thetaDeg = thetaDeg; result.gamma = gamma;
      displayResults(result);
      plotThetaBeta(M, thetaDeg, gamma);
    }

    // initial run
    runCalc();
  </script>
</body>
</html>