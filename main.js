/* ===== ChemSpace — shared scripts ===== */

/* ---------- molecule library ---------- */
const MOLECULES = [
  {name:'Aspirin',      smi:'CC(=O)OC1=CC=CC=C1C(=O)O',                         logP:1.2, use:'Analgesic'},
  {name:'Caffeine',     smi:'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',                     logP:-0.1, use:'CNS stimulant'},
  {name:'Ibuprofen',    smi:'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O',                    logP:3.5, use:'NSAID'},
  {name:'Paracetamol',  smi:'CC(=O)NC1=CC=C(C=C1)O',                            logP:0.5, use:'Analgesic'},
  {name:'Nicotine',     smi:'CN1CCCC1C2=CN=CC=C2',                              logP:1.2, use:'Alkaloid stimulant'},
  {name:'Dopamine',     smi:'C1=CC(=C(C=C1CCN)O)O',                             logP:-1.0, use:'Neurotransmitter'},
  {name:'Serotonin',    smi:'C1=CC2=C(C=C1O)C(=CN2)CCN',                        logP:0.2, use:'Neurotransmitter'},
  {name:'Morphine',     smi:'CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O',        logP:0.9, use:'Opioid analgesic'},
  {name:'Penicillin G', smi:'CC1(C(N2C(S1)C(C2=O)NC(=O)CC3=CC=CC=C3)C(=O)O)C',  logP:1.8, use:'β-lactam antibiotic'},
  {name:'Diazepam',     smi:'CN1C(=O)CN=C(C2=C1C=CC(=C2)Cl)C3=CC=CC=C3',        logP:2.8, use:'Benzodiazepine'},
  {name:'Glucose',      smi:'C(C1C(C(C(C(O1)O)O)O)O)O',                         logP:-2.6, use:'Monosaccharide'},
  {name:'Cholesterol',  smi:'CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C',  logP:8.7, use:'Sterol lipid'},
  {name:'Atorvastatin', smi:'CC(C)C1=C(C(=C(N1CCC(CC(CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4', logP:5.0, use:'Statin (lipid-lowering)'}
];

/* ---------- minimal SMILES parser ---------- */
const MASS = {H:1.008,B:10.81,C:12.011,N:14.007,O:15.999,F:18.998,Na:22.99,Si:28.085,P:30.974,S:32.06,Cl:35.45,K:39.098,Se:78.971,Br:79.904,I:126.904};
const VALENCE = {B:[3],C:[4],N:[3,5],O:[2],P:[3,5],S:[2,4,6],F:[1],Cl:[1],Br:[1],I:[1]};

function parseSMILES(s){
  s = s.trim();
  if(!s) throw new Error('Empty SMILES');
  const atoms=[], bonds=[], stack=[], rings={};
  let prev=-1, pend=null, i=0;
  const addAtom=(a)=>{
    atoms.push(a); const idx=atoms.length-1;
    if(prev>=0){ bonds.push({a:prev,b:idx,o:pend??((a.arom&&atoms[prev].arom)?1.5:1)}); }
    prev=idx; pend=null;
  };
  while(i<s.length){
    const ch=s[i];
    if(ch==='('){ if(prev<0) throw new Error('Branch before atom'); stack.push(prev); i++; continue; }
    if(ch===')'){ if(!stack.length) throw new Error('Unmatched ")"'); prev=stack.pop(); i++; continue; }
    if(ch==='-'){pend=1;i++;continue} if(ch==='='){pend=2;i++;continue}
    if(ch==='#'){pend=3;i++;continue} if(ch===':'){pend=1.5;i++;continue}
    if(ch==='/'||ch==='\\'||ch==='.'){ if(ch==='.') prev=-1; i++; continue; }
    if(/[0-9%]/.test(ch)){
      let n; if(ch==='%'){n=s.substr(i+1,2);i+=3}else{n=ch;i++}
      if(prev<0) throw new Error('Ring digit before atom');
      if(rings[n]){ const r=rings[n]; const arom=atoms[r.at].arom&&atoms[prev].arom;
        bonds.push({a:r.at,b:prev,o:pend??r.o??(arom?1.5:1)}); delete rings[n]; }
      else rings[n]={at:prev,o:pend};
      pend=null; continue;
    }
    if(ch==='['){
      const j=s.indexOf(']',i); if(j<0) throw new Error('Unclosed "["');
      const m=s.slice(i+1,j).match(/^(\d*)([A-Z][a-z]?|[cnospb]|se)(@*)(H\d*)?([+-]+\d*)?/);
      if(!m) throw new Error('Bad bracket atom ['+s.slice(i+1,j)+']');
      let el=m[2], arom=/^[a-z]/.test(el); if(arom) el=el[0].toUpperCase()+el.slice(1);
      const h = m[4] ? (m[4].length>1?+m[4].slice(1):1) : 0;
      addAtom({el,arom,hExp:h}); i=j+1; continue;
    }
    const two=s.substr(i,2);
    if(two==='Cl'||two==='Br'){ addAtom({el:two,arom:false,hExp:null}); i+=2; continue; }
    if('BCNOSPFI'.includes(ch)){ addAtom({el:ch,arom:false,hExp:null}); i++; continue; }
    if('cnosp'.includes(ch)){ addAtom({el:ch.toUpperCase(),arom:true,hExp:null}); i++; continue; }
    throw new Error(`Unexpected "${ch}" at position ${i+1}`);
  }
  if(stack.length) throw new Error('Unmatched "("');
  if(Object.keys(rings).length) throw new Error('Unclosed ring '+Object.keys(rings).join(','));
  if(!atoms.length) throw new Error('No atoms');
  atoms.forEach(a=>{if(!MASS[a.el]) throw new Error('Unsupported element '+a.el)});
  // implicit hydrogens
  atoms.forEach((a,idx)=>{
    if(a.hExp!==null){a.h=a.hExp;return}
    const my=bonds.filter(b=>b.a===idx||b.b===idx);
    const sum=my.reduce((t,b)=>t+(b.o===1.5?1:b.o),0);
    const dbl=my.some(b=>b.o===2);
    const vals=VALENCE[a.el]||[0];
    const v=vals.find(x=>x>=sum)??vals[vals.length-1];
    a.h=Math.max(0, v-sum-(a.arom&&!dbl?1:0));
  });
  return {atoms,bonds};
}

function describe(mol){
  const {atoms,bonds}=mol, cnt={};
  let mw=0,hbd=0,hba=0,heavy=atoms.length;
  atoms.forEach(a=>{
    cnt[a.el]=(cnt[a.el]||0)+1; if(a.h) cnt.H=(cnt.H||0)+a.h;
    mw+=MASS[a.el]+a.h*MASS.H;
    if(a.el==='N'||a.el==='O'){hba++; if(a.h>0) hbd++;}
  });
  const order=['C','H',...Object.keys(cnt).filter(k=>k!=='C'&&k!=='H').sort()];
  const formula=order.filter(k=>cnt[k]).map(k=>k+(cnt[k]>1?cnt[k]:'')).join('');
  // adjacency, ring bonds
  const adj=atoms.map(()=>[]); bonds.forEach((b,k)=>{adj[b.a].push([b.b,k]);adj[b.b].push([b.a,k])});
  const connectedWithout=(s,t,skip)=>{const seen=new Set([s]),q=[s];while(q.length){const u=q.pop();if(u===t)return true;
    for(const [v,k] of adj[u]) if(k!==skip&&!seen.has(v)){seen.add(v);q.push(v)}}return false};
  let comps=0; {const seen=new Set();atoms.forEach((_,i)=>{if(seen.has(i))return;comps++;const q=[i];seen.add(i);
    while(q.length){const u=q.pop();for(const [v] of adj[u]) if(!seen.has(v)){seen.add(v);q.push(v)}}})}
  const rings=bonds.length-atoms.length+comps;
  let rot=0; bonds.forEach((b,k)=>{
    if(b.o!==1) return; if(adj[b.a].length<2||adj[b.b].length<2) return;
    if(connectedWithout(b.a,b.b,k)) return; rot++; });
  return {formula,mw,hbd,hba,heavy,rings,rot,adj};
}

function lipinski(d,logP){
  let v=0; if(d.mw>500)v++; if(d.hbd>5)v++; if(d.hba>10)v++; if(logP!=null&&logP>5)v++; return v;
}

/* ---------- path fingerprint (Daylight-style, 128 bits) ---------- */
function fingerprint(mol,nbits=128,maxLen=4){
  const {atoms,bonds}=mol, adj=atoms.map(()=>[]);
  bonds.forEach(b=>{adj[b.a].push([b.b,b.o]);adj[b.b].push([b.a,b.o])});
  const sym=o=>o===2?'=':o===3?'#':o===1.5?':':'-';
  const lab=i=>atoms[i].el+(atoms[i].arom?'a':'');
  const bits=new Set();
  const hash=str=>{let h=2166136261;for(const c of str){h^=c.charCodeAt(0);h=Math.imul(h,16777619)}return (h>>>0)%nbits};
  const walk=(path,str)=>{
    const rev=str.split(/([-=#:])/).reverse().join('');
    bits.add(hash(str<rev?str:rev));
    if(path.length>maxLen) return;
    const last=path[path.length-1];
    for(const [n,o] of adj[last]) if(!path.includes(n)) walk([...path,n],str+sym(o)+lab(n));
  };
  atoms.forEach((_,i)=>walk([i],lab(i)));
  return bits;
}
const tanimoto=(A,B)=>{let c=0;A.forEach(x=>{if(B.has(x))c++});return c/(A.size+B.size-c||1)};

/* ---------- 3D force-directed layout + renderer ---------- */
const COLORS={C:'#cfd8e8',N:'#5b8cff',O:'#ff5c6c',S:'#ffd23f',P:'#ff9f43',F:'#5dff9d',Cl:'#3ee27a',Br:'#c0703a',I:'#a05cff',H:'#ffffff'};

class Mol3D{
  constructor(canvas,mol,opts={}){
    this.c=canvas; this.ctx=canvas.getContext('2d'); this.o=Object.assign({scale:1,auto:0.006,labels:false,drag:true},opts);
    this.rx=0.4; this.ry=0; this.set(mol);
    if(this.o.drag) this.bindDrag();
    this.resize(); new ResizeObserver(()=>this.resize()).observe(canvas);
    const loop=()=>{this.step();this.draw();requestAnimationFrame(loop)}; loop();
  }
  set(mol){
    this.mol=mol; const n=mol.atoms.length;
    this.p=mol.atoms.map(()=>[(Math.random()-.5)*4,(Math.random()-.5)*4,(Math.random()-.5)*4]);
    this.v=mol.atoms.map(()=>[0,0,0]);
    const adj=mol.atoms.map(()=>new Set()); mol.bonds.forEach(b=>{adj[b.a].add(b.b);adj[b.b].add(b.a)});
    this.springs=mol.bonds.map(b=>[b.a,b.b,1]);
    for(let i=0;i<n;i++) for(const j of adj[i]) for(const k of adj[j]) if(k>i&&!adj[i].has(k)) this.springs.push([i,k,1.65]);
    this.heat=1; this.fit=1;
  }
  bindDrag(){
    let down=false,lx=0,ly=0;
    const pt=e=>e.touches?e.touches[0]:e;
    this.c.addEventListener('pointerdown',e=>{down=true;lx=e.clientX;ly=e.clientY;this.c.setPointerCapture(e.pointerId);this.c.style.cursor='grabbing'});
    this.c.addEventListener('pointermove',e=>{if(!down)return;this.ry+=(e.clientX-lx)*0.01;this.rx+=(e.clientY-ly)*0.01;lx=e.clientX;ly=e.clientY});
    this.c.addEventListener('pointerup',()=>{down=false;this.c.style.cursor='grab'});
  }
  resize(){const r=this.c.getBoundingClientRect(),d=devicePixelRatio||1;this.c.width=r.width*d;this.c.height=r.height*d;this.w=r.width;this.h=r.height;this.ctx.setTransform(d,0,0,d,0,0)}
  step(){
    if(this.heat>0.02){
      const p=this.p,v=this.v,n=p.length;
      for(let it=0;it<3;it++){
        const f=p.map(()=>[0,0,0]);
        for(let i=0;i<n;i++)for(let j=i+1;j<n;j++){
          const d=[p[i][0]-p[j][0],p[i][1]-p[j][1],p[i][2]-p[j][2]],r2=d[0]**2+d[1]**2+d[2]**2+0.01,k=0.35/r2;
          for(let a=0;a<3;a++){f[i][a]+=d[a]*k;f[j][a]-=d[a]*k}}
        for(const [i,j,L] of this.springs){
          const d=[p[j][0]-p[i][0],p[j][1]-p[i][1],p[j][2]-p[i][2]],r=Math.hypot(...d)+1e-6,k=(r-L)*0.5/r;
          for(let a=0;a<3;a++){f[i][a]+=d[a]*k;f[j][a]-=d[a]*k}}
        for(let i=0;i<n;i++)f[i][2]-=p[i][2]*0.06;
        for(let i=0;i<n;i++)for(let a=0;a<3;a++){v[i][a]=(v[i][a]+f[i][a]*0.2)*0.82;p[i][a]+=v[i][a]}
      }
      this.heat*=0.992;
      // centre + fit
      const c=[0,0,0];p.forEach(q=>{c[0]+=q[0];c[1]+=q[1];c[2]+=q[2]});c.forEach((_,a)=>c[a]/=n);
      p.forEach(q=>{q[0]-=c[0];q[1]-=c[1];q[2]-=c[2]});
      const R=Math.max(...p.map(q=>Math.hypot(...q)),1); this.fit+=(R-this.fit)*0.1;
    }
    this.ry+=this.o.auto;
  }
  draw(){
    const {ctx,w,h}=this; ctx.clearRect(0,0,w,h); if(!this.p.length)return;
    const s=Math.min(w,h)/(2.4*this.fit)*this.o.scale, cy=Math.cos(this.ry),sy=Math.sin(this.ry),cx=Math.cos(this.rx),sx=Math.sin(this.rx);
    const pr=this.p.map(([x,y,z])=>{const x1=x*cy+z*sy,z1=-x*sy+z*cy,y1=y*cx-z1*sx,z2=y*sx+z1*cx,persp=4/(4+z2*.25);
      return {x:w/2+x1*s*persp,y:h/2+y1*s*persp,z:z2,k:persp}});
    const items=[];
    this.mol.bonds.forEach(b=>items.push({z:(pr[b.a].z+pr[b.b].z)/2-.01,t:'b',b}));
    pr.forEach((q,i)=>items.push({z:q.z,t:'a',i}));
    items.sort((a,b)=>b.z-a.z);
    const R=Math.max(3,s*0.26);
    for(const it of items){
      if(it.t==='b'){const A=pr[it.b.a],B=pr[it.b.b],dx=B.x-A.x,dy=B.y-A.y,L=Math.hypot(dx,dy)||1,nx=-dy/L,ny=dx/L;
        const order=it.b.o===1.5?2:it.b.o, gap=R*0.35;
        ctx.lineWidth=Math.max(1.2,R*0.22); ctx.strokeStyle=`rgba(160,180,215,${0.35+0.4*A.k})`;
        for(let m=0;m<order;m++){const off=(m-(order-1)/2)*gap; ctx.beginPath();
          if(it.b.o===1.5&&m===1) ctx.setLineDash([R*.3,R*.3]);
          ctx.moveTo(A.x+nx*off,A.y+ny*off);ctx.lineTo(B.x+nx*off,B.y+ny*off);ctx.stroke();ctx.setLineDash([])}
      }else{const q=pr[it.i],a=this.mol.atoms[it.i],col=COLORS[a.el]||'#ff6fa8',r=R*q.k*(a.el==='C'?0.85:1);
        const g=ctx.createRadialGradient(q.x-r*.35,q.y-r*.35,r*.1,q.x,q.y,r);
        g.addColorStop(0,'#fff');g.addColorStop(.25,col);g.addColorStop(1,shade(col));
        ctx.fillStyle=g;ctx.beginPath();ctx.arc(q.x,q.y,r,0,Math.PI*2);ctx.fill();
        if(this.o.labels&&a.el!=='C'&&r>8){ctx.fillStyle='#04121a';ctx.font=`700 ${r*.9}px Inter`;ctx.textAlign='center';ctx.textBaseline='middle';ctx.fillText(a.el,q.x,q.y+1)}
      }
    }
  }
}
function shade(hex){const n=parseInt(hex.slice(1),16);return `rgb(${(n>>16)*.3|0},${(n>>8&255)*.3|0},${(n&255)*.3|0})`}

/* ---------- background: drifting atom network ---------- */
function bgNetwork(){
  const c=document.getElementById('bg'); if(!c) return; const ctx=c.getContext('2d'); let W,H,pts;
  const init=()=>{W=c.width=innerWidth;H=c.height=innerHeight;
    pts=Array.from({length:Math.min(90,W*H/16000|0)},()=>({x:Math.random()*W,y:Math.random()*H,vx:(Math.random()-.5)*.3,vy:(Math.random()-.5)*.3,
      c:['#2ee6c5','#8b7cff','#ff6fa8','#5b8cff'][Math.random()*4|0],r:Math.random()*2+1}))};
  let mx=-999,my=-999; addEventListener('mousemove',e=>{mx=e.clientX;my=e.clientY});
  addEventListener('resize',init); init();
  (function loop(){ctx.clearRect(0,0,W,H);
    for(const p of pts){p.x+=p.vx;p.y+=p.vy;if(p.x<0||p.x>W)p.vx*=-1;if(p.y<0||p.y>H)p.vy*=-1;
      const dm=Math.hypot(p.x-mx,p.y-my); if(dm<120){p.x+=(p.x-mx)/dm*.8;p.y+=(p.y-my)/dm*.8}}
    for(let i=0;i<pts.length;i++)for(let j=i+1;j<pts.length;j++){const a=pts[i],b=pts[j],d=Math.hypot(a.x-b.x,a.y-b.y);
      if(d<130){ctx.strokeStyle=`rgba(120,150,210,${(1-d/130)*.35})`;ctx.lineWidth=1;ctx.beginPath();ctx.moveTo(a.x,a.y);ctx.lineTo(b.x,b.y);ctx.stroke()}}
    for(const p of pts){ctx.fillStyle=p.c;ctx.shadowColor=p.c;ctx.shadowBlur=8;ctx.beginPath();ctx.arc(p.x,p.y,p.r,0,7);ctx.fill()}
    ctx.shadowBlur=0; requestAnimationFrame(loop)})();
}

/* ---------- scroll reveal + counters ---------- */
function reveal(){
  const io=new IntersectionObserver(es=>es.forEach(e=>{if(e.isIntersecting){e.target.classList.add('in');
    e.target.querySelectorAll?.('[data-count]').forEach(countUp);io.unobserve(e.target)}}),{threshold:.15});
  document.querySelectorAll('.reveal,.step').forEach(el=>io.observe(el));
}
function countUp(el){
  const end=+el.dataset.count,suf=el.dataset.suffix||'',dec=+(el.dataset.dec||0),t0=performance.now();
  (function f(t){const k=Math.min(1,(t-t0)/1800),e=1-Math.pow(1-k,3);el.textContent=(end*e).toFixed(dec)+suf;if(k<1)requestAnimationFrame(f)})(t0);
}

function navActive(){
  const page=document.body.dataset.page;
  document.querySelectorAll('.links a').forEach(a=>{if(a.dataset.p===page)a.classList.add('active')});
}

document.addEventListener('DOMContentLoaded',()=>{bgNetwork();reveal();navActive();
  const init={home:initHome,molecules:initMolecules,workflow:initWorkflow}[document.body.dataset.page]; init&&init();
  const y=document.getElementById('year'); if(y) y.textContent=new Date().getFullYear();
});

/* ================= page: home ================= */
function initHome(){
  let idx=1; const nameEl=document.getElementById('molname');
  const viewer=new Mol3D(document.getElementById('mol3d'),parseSMILES(MOLECULES[idx].smi),{labels:true,auto:0.005});
  const show=()=>{viewer.set(parseSMILES(MOLECULES[idx].smi));nameEl.textContent=MOLECULES[idx].name};
  show();
  document.getElementById('nextmol').onclick=()=>{idx=(idx+1)%MOLECULES.length;show()};
  // typing SMILES
  const el=document.getElementById('typed'), lab=document.getElementById('typedname'); let m=0;
  (async function type(){
    while(true){const {name,smi}=MOLECULES[m%MOLECULES.length]; lab.textContent=name+':';
      for(let i=1;i<=smi.length;i++){el.textContent=smi.slice(0,i);await sleep(45)}
      await sleep(1600);
      for(let i=smi.length;i>=0;i--){el.textContent=smi.slice(0,i);await sleep(14)}
      m++;}
  })();
}
const sleep=ms=>new Promise(r=>setTimeout(r,ms));

/* ================= page: molecules ================= */
function initMolecules(){
  const inp=document.getElementById('smiIn'), err=document.getElementById('err'), view=new Mol3D(document.getElementById('anaCanvas'),parseSMILES('C'),{labels:true});
  const out=k=>document.getElementById(k);
  const run=()=>{
    try{const mol=parseSMILES(inp.value), d=describe(mol); err.textContent='';
      view.set(mol);
      const vals={rFormula:d.formula,rMW:d.mw.toFixed(2),rHeavy:d.heavy,rRings:d.rings,rHBD:d.hbd,rHBA:d.hba,rRot:d.rot,rLip:lipinski(d)};
      for(const k in vals){const e=out(k);e.textContent=vals[k];e.parentElement.classList.remove('pop');void e.offsetWidth;e.parentElement.classList.add('pop')}
    }catch(e){err.textContent='⚠ '+e.message}
  };
  document.getElementById('anaBtn').onclick=run;
  inp.addEventListener('keydown',e=>{if(e.key==='Enter')run()});
  const chips=document.getElementById('chips');
  MOLECULES.slice(0,8).forEach(m=>{const b=document.createElement('button');b.className='chip';b.textContent=m.name;
    b.onclick=()=>{inp.value=m.smi;run()};chips.appendChild(b)});
  inp.value=MOLECULES[0].smi; run();

  // gallery
  const g=document.getElementById('gallery');
  MOLECULES.forEach((m,i)=>{
    const mol=parseSMILES(m.smi), d=describe(mol), v=lipinski(d,m.logP);
    const card=document.createElement('div'); card.className='card mol-card reveal d'+(i%3+1);
    card.innerHTML=`<div class="flip">
      <div class="face"><h3>${m.name}</h3><div style="color:var(--muted);font-size:.85rem">${m.use}</div>
        <canvas></canvas><div class="smi mono">${m.smi}</div><div class="hint">Click to flip →</div></div>
      <div class="face back"><h3>${m.name}</h3>
        <div class="kv"><span>Formula</span><b class="mono">${d.formula}</b></div>
        <div class="kv"><span>Mol. weight</span><b class="mono">${d.mw.toFixed(2)}</b></div>
        <div class="kv"><span>logP (approx.)</span><b class="mono">${m.logP}</b></div>
        <div class="kv"><span>H-donors / acceptors</span><b class="mono">${d.hbd} / ${d.hba}</b></div>
        <div class="kv"><span>Rotatable bonds</span><b class="mono">${d.rot}</b></div>
        <div><span class="badge ${v<=1?'ok':'no'}">${v===0?'✓ Rule of 5: passes':v===1?'✓ Rule of 5: 1 violation':'✗ Rule of 5: '+v+' violations'}</span></div></div></div>`;
    card.onclick=()=>card.classList.toggle('flipped');
    g.appendChild(card);
    new Mol3D(card.querySelector('canvas'),mol,{drag:false,scale:1.15,auto:0.01+Math.random()*0.006});
  });
  reveal();
}

/* ================= page: workflow ================= */
function initWorkflow(){
  // timeline fill on scroll
  const tl=document.querySelector('.timeline'), fill=tl.querySelector('.fill');
  const upd=()=>{const r=tl.getBoundingClientRect(),k=Math.min(1,Math.max(0,(innerHeight*0.6-r.top)/r.height));fill.style.height=k*100+'%'};
  addEventListener('scroll',upd,{passive:true}); upd();

  // fingerprint demo
  const sa=document.getElementById('molA'), sb=document.getElementById('molB');
  MOLECULES.forEach((m,i)=>{sa.add(new Option(m.name,i));sb.add(new Option(m.name,i))});
  sa.value=5; sb.value=6;
  const mk=id=>{const box=document.getElementById(id);for(let i=0;i<128;i++){const d=document.createElement('div');d.className='bit';box.appendChild(d)}return [...box.children]};
  const bitsA=mk('bitsA'), bitsB=mk('bitsB');
  const circ=2*Math.PI*60, val=document.querySelector('.gauge .val'); val.style.strokeDasharray=circ; val.style.strokeDashoffset=circ;
  const compare=()=>{
    const A=fingerprint(parseSMILES(MOLECULES[sa.value].smi)), B=fingerprint(parseSMILES(MOLECULES[sb.value].smi));
    [bitsA,bitsB].forEach(arr=>arr.forEach(b=>b.className='bit'));
    for(let i=0;i<128;i++){
      setTimeout(()=>{
        if(A.has(i)) bitsA[i].className='bit on'+(B.has(i)?' both':'');
        if(B.has(i)) bitsB[i].className='bit on b'+(A.has(i)?' both':'');
      },i*8);
    }
    const t=tanimoto(A,B); let c=0; A.forEach(x=>{if(B.has(x))c++});
    setTimeout(()=>{val.style.strokeDashoffset=circ*(1-t);const te=document.getElementById('tani');te.dataset.count=t;te.dataset.dec=2;countUp(te)},300);
    document.getElementById('taniInfo').innerHTML=`Bits in A: <b>${A.size}</b> · Bits in B: <b>${B.size}</b> · Shared: <b>${c}</b><br>T = ${c} / (${A.size} + ${B.size} − ${c})`;
  };
  sa.onchange=compare; sb.onchange=compare; document.getElementById('cmpBtn').onclick=compare;
  const io=new IntersectionObserver(e=>{if(e[0].isIntersecting){compare();io.disconnect()}},{threshold:.3});
  io.observe(document.getElementById('fpSection'));

  chemSpace();
}

function chemSpace(){
  const c=document.getElementById('space'), ctx=c.getContext('2d'), tip=document.getElementById('spaceTip');
  const data=MOLECULES.map(m=>{const d=describe(parseSMILES(m.smi));return {...m,mw:d.mw,v:lipinski(d,m.logP)}});
  let W,H,t0=null,started=false,hover=null; const pad={l:56,r:20,t:20,b:44};
  const X=v=>pad.l+(v-100)/(600-100)*(W-pad.l-pad.r), Y=v=>H-pad.b-(v+4)/(10+4)*(H-pad.t-pad.b);
  const size=()=>{const r=c.getBoundingClientRect(),d=devicePixelRatio||1;W=r.width;H=r.height;c.width=W*d;c.height=H*d;ctx.setTransform(d,0,0,d,0,0)};
  size(); addEventListener('resize',size);
  c.addEventListener('mousemove',e=>{const r=c.getBoundingClientRect(),x=e.clientX-r.left,y=e.clientY-r.top;
    hover=data.find(p=>Math.hypot(X(p.mw)-x,Y(p.logP)-y)<10)||null;
    if(hover){tip.style.opacity=1;tip.style.left=X(hover.mw)+14+'px';tip.style.top=Y(hover.logP)-10+'px';
      tip.innerHTML=`<b>${hover.name}</b><br>MW ${hover.mw.toFixed(1)} · logP ${hover.logP}`}else tip.style.opacity=0});
  c.addEventListener('mouseleave',()=>{hover=null;tip.style.opacity=0});
  new IntersectionObserver(e=>{if(e[0].isIntersecting) started=true},{threshold:.3}).observe(c);
  (function loop(t){
    ctx.clearRect(0,0,W,H);
    if(started&&t0===null)t0=t; const k=t0===null?0:Math.min(1,(t-t0)/1500);
    // grid
    ctx.strokeStyle='#1f2c47';ctx.fillStyle='#8a9bbb';ctx.font='12px Inter';ctx.lineWidth=1;
    for(let m=100;m<=600;m+=100){ctx.beginPath();ctx.moveTo(X(m),pad.t);ctx.lineTo(X(m),H-pad.b);ctx.stroke();ctx.textAlign='center';ctx.fillText(m,X(m),H-pad.b+18)}
    for(let l=-4;l<=10;l+=2){ctx.beginPath();ctx.moveTo(pad.l,Y(l));ctx.lineTo(W-pad.r,Y(l));ctx.stroke();ctx.textAlign='right';ctx.fillText(l,pad.l-8,Y(l)+4)}
    ctx.textAlign='center';ctx.fillText('Molecular weight (Da)',(W+pad.l)/2,H-6);
    ctx.save();ctx.translate(14,(H-pad.b+pad.t)/2);ctx.rotate(-Math.PI/2);ctx.fillText('logP',0,0);ctx.restore();
    // Lipinski box
    const e=1-Math.pow(1-k,3);
    ctx.fillStyle=`rgba(46,230,197,${.07*e})`;ctx.strokeStyle=`rgba(46,230,197,${.6*e})`;ctx.setLineDash([6,5]);
    const bx=X(100),by=Y(5),bw=(X(500)-X(100))*e,bh=(Y(-4)-Y(5));
    ctx.fillRect(bx,by,bw,bh);ctx.strokeRect(bx,by,bw,bh);ctx.setLineDash([]);
    if(k>.6){ctx.fillStyle='rgba(46,230,197,.8)';ctx.textAlign='left';ctx.fillText('Lipinski space (MW ≤ 500, logP ≤ 5)',X(110),Y(5)+16)}
    // points
    data.forEach((p,i)=>{const kk=Math.max(0,Math.min(1,(k*1.6-i*0.05)));if(!kk)return;
      const px=X(p.mw),py=Y(-4)+(Y(p.logP)-Y(-4))*easeBack(kk),col=p.v===0?'#2ee6c5':'#ff6fa8',r=(hover===p?9:6)+Math.sin(t/400+i)*1;
      ctx.fillStyle=col;ctx.shadowColor=col;ctx.shadowBlur=14;ctx.beginPath();ctx.arc(px,py,r,0,7);ctx.fill();ctx.shadowBlur=0;
      if(W>600){ctx.fillStyle='#e6ecf7';ctx.textAlign='left';ctx.font='11px Inter';ctx.fillText(p.name,px+10,py-8)}});
    requestAnimationFrame(loop)})(0);
}
function easeBack(x){const c1=1.70158,c3=c1+1;return 1+c3*Math.pow(x-1,3)+c1*Math.pow(x-1,2)}
