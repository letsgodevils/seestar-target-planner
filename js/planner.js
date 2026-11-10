/*!
 * Seestar Target Planner
 * All visibility math runs in the browser. Code: MIT License (see LICENSE).
 * Object catalog (data/catalog.json): OpenNGC by Mattia Verga, CC BY-SA 4.0 (see NOTICE.md).
 */
(function () {
'use strict';

const RAD=Math.PI/180, sinD=x=>Math.sin(x*RAD), cosD=x=>Math.cos(x*RAD), norm=x=>((x%360)+360)%360;
function sunPos(jd){const n=jd-2451545,L=280.460+0.9856474*n,g=357.528+0.9856003*n,lam=L+1.915*sinD(g)+0.020*sinD(2*g),eps=23.439-0.0000004*n;
 return {ra:norm(Math.atan2(cosD(eps)*sinD(lam),cosD(lam))/RAD),dec:Math.asin(sinD(eps)*sinD(lam))/RAD,lam};}
function moonPos(jd){const T=(jd-2451545)/36525;
 const lam=218.32+481267.881*T+6.29*sinD(135.0+477198.87*T)-1.27*sinD(259.3-413335.36*T)+0.66*sinD(235.7+890534.22*T)+0.21*sinD(269.9+954397.74*T)-0.19*sinD(357.5+35999.05*T)-0.11*sinD(186.5+966404.03*T);
 const bet=5.13*sinD(93.3+483202.02*T)+0.28*sinD(228.2+960400.89*T)-0.28*sinD(318.3+6003.15*T)-0.17*sinD(217.6-407332.21*T);
 const par=0.9508+0.0518*cosD(135.0+477198.87*T)+0.0095*cosD(259.3-413335.36*T)+0.0078*cosD(235.7+890534.22*T)+0.0028*cosD(269.9+954397.74*T);
 const eps=23.439-0.0000004*(jd-2451545);
 const x=cosD(bet)*cosD(lam),y=cosD(eps)*cosD(bet)*sinD(lam)-sinD(eps)*sinD(bet),z=sinD(eps)*cosD(bet)*sinD(lam)+cosD(eps)*sinD(bet);
 return {ra:norm(Math.atan2(y,x)/RAD),dec:Math.asin(z)/RAD,lam,bet,par};}
const lstOf=(jd,lon)=>norm(280.46061837+360.98564736629*(jd-2451545)+lon);
const altOf=(ra,dec,lat,ls)=>Math.asin(sinD(lat)*sinD(dec)+cosD(lat)*cosD(dec)*cosD(ls-ra))/RAD;
const raDeg=s=>{const [h,m,x]=s.split(':').map(Number);return 15*(h+m/60+(x||0)/3600);};
const decDeg=s=>{const sg=s.trim()[0]==='-'?-1:1;const [d,m,x]=s.replace(/[+-]/,'').split(':').map(Number);return sg*(d+m/60+(x||0)/3600);};
const MON=["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"];
const NIGHTS=365, STEP=10;
// Precompute dark samples for a site once; each target then needs only altitude checks.
function prepSite(lat,lon,start){
  const nights=[], illum=[], dates=[];
  const y0=start.getFullYear(), m0=start.getMonth(), d0=start.getDate();
  for(let n=0;n<NIGHTS;n++){
    const day=new Date(Date.UTC(y0,m0,d0+n)); dates.push(day);
    const jd0=(day.getTime()+(12-lon/15)*3600e3)/86400000+2440587.5;
    const ls=[], ok=[];
    for(let k=0;k<144;k++){
      const jd=jd0+k*STEP/1440, l=lstOf(jd,lon), s=sunPos(jd);
      if(altOf(s.ra,s.dec,lat,l)>-18) continue;
      const mo=moonPos(jd), ma=altOf(mo.ra,mo.dec,lat,l)-mo.par, kill=(1-cosD(mo.bet)*cosD(mo.lam-s.lam))/2;
      ls.push(l); ok.push(ma<0||kill<0.25);
    }
    nights.push({ls,ok});
    const mm=moonPos(jd0+0.5), ss=sunPos(jd0+0.5); illum.push((1-cosD(mm.bet)*cosD(mm.lam-ss.lam))/2);
  }
  const newMoons=[]; for(let n=1;n<NIGHTS-1;n++) if(illum[n]<=illum[n-1]&&illum[n]<illum[n+1]&&illum[n]<0.05) newMoons.push(fmtD(dates[n]));
  return {lat,lon,nights,dates,newMoons};
}
const fmtD=d=>`${MON[d.getUTCMonth()]} ${d.getUTCDate()}`;
const fmtW=(a,b)=>a.getTime()===b.getTime()?`${fmtD(a)}, ${a.getUTCFullYear()}`:a.getUTCMonth()===b.getUTCMonth()?`${fmtD(a)}–${b.getUTCDate()}, ${b.getUTCFullYear()}`:`${fmtD(a)} – ${fmtD(b)}, ${b.getUTCFullYear()}`;
function evalTarget(S,ra,dec){
  const maxalt=90-Math.abs(S.lat-dec), thr=maxalt>=35?25:Math.max(12,maxalt-8);
  if(maxalt<12) return {visible:false,maxalt:Math.round(maxalt)};
  const tot=new Float32Array(NIGHTS), mf=new Float32Array(NIGHTS), sd=sinD(dec)*sinD(S.lat), cd=cosD(dec)*cosD(S.lat), st=sinD(thr);
  S.nights.forEach((nt,n)=>{for(let i=0;i<nt.ls.length;i++){ if(sd+cd*cosD(nt.ls[i]-ra)>=st){tot[n]+=STEP/60; if(nt.ok[i]) mf[n]+=STEP/60;} }});
  let pi=0; for(let n=1;n<NIGHTS;n++) if(tot[n]>tot[pi]) pi=n;
  const peakH=tot[pi]; if(peakH<0.5) return {visible:false,maxalt:Math.round(maxalt)};
  const sums=Array(12).fill(0), cnt=Array(12).fill(0);
  for(let n=0;n<NIGHTS;n++){const m=S.dates[n].getUTCMonth(); sums[m]+=tot[n]; cnt[m]++;}
  const good=sums.map((s,m)=>cnt[m]&&s/cnt[m]>=Math.max(2,0.5*peakH)); const gc=good.filter(Boolean).length; let season;
  if(gc===12) season="Year-round"; else if(gc===0) season=MON[S.dates[pi].getUTCMonth()];
  else { let s=good.findIndex((g,m)=>g&&!good[(m+11)%12]), e=s; while(good[(e+1)%12]&&(e+1)%12!==s) e=(e+1)%12; season=s===e?MON[s]:`${MON[s]}–${MON[e]}`; }
  const th=Math.max(2,0.6*peakH), wins=[]; let cur=null;
  for(let n=0;n<=NIGHTS;n++){ if(n<NIGHTS&&mf[n]>=th){ if(!cur) cur=[n,n]; else cur[1]=n; } else if(cur){ if(cur[1]-cur[0]>=2) wins.push(cur); cur=null; } }
  const w2=wins.slice(0,2);
  return {visible:true,maxalt:Math.round(maxalt),peak:fmtD(S.dates[pi]),peakKey:S.dates[pi].getUTCMonth()*100+S.dates[pi].getUTCDate(),
    peakH:Math.round(peakH*10)/10,season,windows:w2.map(w=>fmtW(S.dates[w[0]],S.dates[w[1]])),winStart:w2.length?w2[0][0]:9999};
}
let CAT=[], LISTS=null, CITIES=[];
const IDX=new Map(), MIDX=new Map();
const resolve=r=>{ let g=0; while(r&&r[1]==='Dup'&&g++<3) r=IDX.get(r[2].replace(/\s+/g,'').toUpperCase()); return r; };
function buildIndex(){
  CAT.forEach(r=>{ IDX.set(r[0].replace(/\s+/g,'').toUpperCase(), r); if(r[1]!=='Dup'&&r[8]) MIDX.set(r[8], r); });
  CAT.forEach(r=>{ if(r[1]==='Dup'&&r[3]&&!MIDX.has(r[3])){ const t=IDX.get(r[2].replace(/\s+/g,'').toUpperCase()); if(t) MIDX.set(r[3],t);} });
  if(!MIDX.has('M102')){ const t=IDX.get('NGC5866'); if(t){ const c=t.slice(); c[8]='M102'; MIDX.set('M102',c);} }
}
const TYPES={G:'Galaxy',GPair:'Galaxy pair',GTrpl:'Galaxy triplet',GGroup:'Galaxy group',PN:'Planetary nebula',HII:'HII region',EmN:'Emission nebula',Neb:'Nebula',RfN:'Reflection nebula',SNR:'Supernova remnant','Cl+N':'Cluster + nebula',OCl:'Open cluster',GCl:'Globular cluster','*Ass':'Stellar association',DrkN:'Dark nebula','*':'Star','**':'Double star',Nova:'Nova',Other:'Other'};
function estimate(r){
  const t=r[1], mag=r[6], sb=r[7], maj=r[4];
  const tier=['1–2 h','2–3 h','3–4 h','4–5 h','5–6 h','6 h+'];
  let i, f='UV/IR';
  if(/^G/.test(t)){ i=mag==null?4:mag<=7?0:mag<=8.5?1:mag<=9.5?2:mag<=10.5?3:mag<=11.5?4:5; if(sb!=null&&sb>=14&&i<5) i++; }
  else if(t==='PN'){ f='Dual-band'; i=mag==null?3:mag<=8?0:mag<=10?1:mag<=12?2:4; }
  else if(t==='SNR'){ f='Dual-band'; i=3; }
  else if(['HII','EmN','Neb','Cl+N'].includes(t)){ f='Dual-band'; i=mag==null?2:mag<=6?0:mag<=8?1:2; }
  else if(t==='RfN'){ i=mag!=null&&mag<=7?1:2; }
  else if(t==='DrkN'){ i=3; }
  else if(t==='GCl'){ return {integ:'1–2 h',f}; }
  else if(t==='OCl'||t==='*Ass'){ return {integ:'0.5–1 h',f}; }
  else return {integ:'15–30 min',f};
  return {integ:tier[i],f};
}
function lookupRows(q){
  const s=q.trim(); if(!s) return [];
  const k=s.replace(/\s+/g,'').toUpperCase();
  let m;
  if((m=k.match(/^(?:M|MESSIER)(\d{1,3})$/))){ const r=MIDX.get('M'+(+m[1])); return r?[r]:[]; }
  if((m=k.match(/^(NGC|IC)0*(\d+)([A-Z]?)$/))){ const r=resolve(IDX.get(m[1]+m[2]+m[3])); return r?[r]:[]; }
  const ql=s.toLowerCase(); const out=[];
  for(const r of CAT){ if(r[1]!=='Dup'&&r[9]&&r[9].toLowerCase().includes(ql)){ out.push(r); if(out.length>=8) break; } }
  return out;
}
function curatedFor(r){
  const ids=[r[0].replace(/\s+/g,''),(r[8]||'')].filter(Boolean);
  for(const L of [LISTS.galaxies,LISTS.nebulae]) for(const t of L){
    const first=(t.members.match(/\b(M\d+|NGC \d+|IC \d+)/)||[])[0]; if(first&&ids.includes(first.replace(/\s+/g,''))) return t; }
  return null;
}
function lookTarget(r){
  const cur=curatedFor(r), est=estimate(r);
  const label=r[8]?(/^M0*\d+$/.test(r[0])?r[8]:`${r[8]} · ${r[0]}`):r[0];
  const common=(r[9]||'').split(',').filter(Boolean)[0];
  const px=r[4]?Math.round(r[4]*60/3.7):null;
  let notes=cur?cur.notes:(TYPES[r[1]]||r[1])+(r[4]?`, ${r[4]}′ across (~${px} px${r[4]>255?'; larger than the field, plan a mosaic':r[4]<3?'; tiny, crop heavily':''})`:'')+'.';
  const wq=r[8]?`https://en.wikipedia.org/wiki/Messier_${r[8].slice(1)}`:`https://en.wikipedia.org/w/index.php?search=${encodeURIComponent(r[0])}&go=Go`;
  return {name:cur?cur.name:(common||label), members:(cur||common)?label:(TYPES[r[1]]||r[1]), sizeTxt:r[4]?`${r[4]}${r[5]?'x'+r[5]:''}`:'?', mag:r[6], integ:cur?cur.integ:est.integ, filter:cur?cur.filter:est.f,
    notes, wiki:cur?cur.wiki:[[label,wq]], frame:cur&&cur.frame?cur.frame:(r[4]?[r[4],r[5]||r[4]]:null), raD:r[2], decD:r[3], typ:TYPES[r[1]]||r[1], est:!cur};
}
// ---- UI (all element ids are prefixed ss- to avoid clashing with the site template)
// Settings come from the inert JSON block <script type="application/json" id="ss-config"> in index.html.
const OPT=(()=>{ try{ return JSON.parse(document.getElementById('ss-config').textContent)||{}; }catch(e){ return {}; } })();
const DEF=OPT.defaults||{name:'Rockwood, TX',lat:31.4935,lon:-99.3706,list:'galaxies'};
const $=id=>document.getElementById('ss-'+id);
const esc=s=>String(s??'').replace(/[&<>"]/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}[c]));
const store={get(k){try{return localStorage.getItem(k);}catch(e){return null;}},set(k,v){try{localStorage.setItem(k,v);}catch(e){}}};
let page=1, S=null, kind=DEF.list==='nebulae'?'nebulae':'galaxies', R={}, look=[], customName=null, catPromise=null;
const start=new Date(); start.setHours(0,0,0,0);
const getJSON=url=>fetch(url,{credentials:'same-origin'}).then(r=>{ if(!r.ok) throw new Error('HTTP '+r.status); return r.json(); });
function loadCatalog(){ if(!catPromise) catPromise=getJSON(OPT.catalog).then(d=>{ CAT=d; buildIndex(); }); return catPromise; }
function setFields(l){ $('lat').value=(+l.lat).toFixed(4); $('lon').value=(+l.lon).toFixed(4); const i=CITIES.findIndex(c=>c[0]===l.name); $('city').value=i>=0?i:''; customName=i>=0?null:(l.name||null); }
function setKind(k){ kind=k; page=1; $('bG').setAttribute('aria-pressed',String(k==='galaxies')); $('bN').setAttribute('aria-pressed',String(k==='nebulae')); $('hKind').textContent=k==='galaxies'?'Galaxies':'Nebulae'; store.set('seestarKind',k); render(); }
function inSeason(s,m){ if(s==="Year-round")return true; const [a,b]=s.split('–').map(x=>MON.indexOf(x)+1); const e=b||a; return a<=e?(m>=a&&m<=e):(m>=a||m<=e); }
function run(){
  const lat=parseFloat($('lat').value), lon=parseFloat($('lon').value);
  if(!(lat>=-90&&lat<=90)||!(lon>=-180&&lon<=180)){ $('err').textContent='Enter latitude between -90 and 90 and longitude between -180 and 180 (west is negative).'; return; }
  $('err').textContent='';
  const c=CITIES[+$('city').value];
  const loc={name:$('city').value!==''&&c?c[0]:(customName||'Custom location'),lat,lon};
  store.set('seestarLoc',JSON.stringify(loc));
  $('rows').innerHTML='<p class="ss-busy">Calculating…</p>';
  setTimeout(()=>{
    S=prepSite(lat,lon,start); R={}; page=1;
    for(const k of ['galaxies','nebulae']) R[k]=LISTS[k].map(t=>evalTarget(S,raDeg(t.ra),decDeg(t.dec)));
    $('fLoc').textContent=loc.name;
    $('fLL').textContent=`${Math.abs(lat).toFixed(2)}°${lat>=0?'N':'S'}, ${Math.abs(lon).toFixed(2)}°${lon>=0?'E':'W'}`;
    const end=new Date(start); end.setDate(end.getDate()+364);
    $('fRange').textContent=`${MON[start.getMonth()]} ${start.getDate()} ${start.getFullYear()} → ${MON[end.getMonth()]} ${end.getFullYear()}`;
    $('fNM').textContent=S.newMoons.slice(0,7).join(' · ');
    render(); renderLook();
  },20);
}
const SCOPES={
  S30:{name:'Seestar S30',ap:30,fl:150,f:5.0,w:1920,h:1080,batt:6},
  S30P:{name:'Seestar S30 Pro',ap:30,fl:160,f:5.3,w:3840,h:2160,batt:6},
  S50:{name:'Seestar S50',ap:50,fl:250,f:5.0,w:1920,h:1080,batt:6},
  S50P:{name:'Seestar S50 Pro',ap:50,fl:260,f:5.2,w:3840,h:2160,batt:7}};
const PIX=2.9; // µm, all four telephoto sensors
let scope=SCOPES[DEF.scope]?DEF.scope:'S30P';
function scopeGeom(k){ const sc=SCOPES[k], scale=206.265*PIX/sc.fl; return {...sc, scale, L:sc.w*scale/60, S:sc.h*scale/60}; }
let G=scopeGeom(scope);
function framing(fr){
  if(!fr||!fr[0]) return null;
  const a=Math.max(fr[0],fr[1]||fr[0]), b=Math.min(fr[0],fr[1]||fr[0]), px=Math.round(a*60/G.scale);
  if(a<=0.85*G.S) return px<150?{c:'ss-sm',t:`Small · ${px} px`,k:2,f:a/G.S}:{c:'ss-ok',t:`Fits · ${px.toLocaleString()} px`,k:0,f:a/G.S};
  if(a<=0.92*G.L&&b<=0.92*G.S) return {c:'ss-len',t:'Fits lengthwise',k:1,f:a/G.L};
  const cols=Math.ceil(a/(0.85*G.L)), rows=Math.ceil(b/(0.85*G.S));
  return {c:'ss-mos',t:`${cols}×${rows} mosaic`,k:3+cols*rows,f:0};
}
function scaleInteg(str){
  const k=(G.f/5.3)**2; if(Math.abs(k-1)<0.01) return str;
  const m=String(str).match(/^([\d.]+)\s*(?:–\s*([\d.]+))?\s*(h|min)(\+)?$/); if(!m) return str;
  const toH=v=>m[3]==='min'?v/60:v, lo=toH(+m[1])*k, hi=m[2]?toH(+m[2])*k:null;
  const rm=v=>Math.max(5,Math.round(v*60/5)*5), rh=v=>Math.max(0.5,Math.round(v*2)/2);
  if(lo<1&&(hi===null||hi<1)) return `${rm(lo)}${hi!==null?'–'+rm(hi):''} min${m[4]||''}`;
  return `${rh(lo)}${hi!==null?'–'+rh(hi):''} h${m[4]||''}`;
}
function updateScopeUI(){
  document.querySelectorAll('#ss-scopeSeg button').forEach(b=>b.setAttribute('aria-pressed',String(b.dataset.s===scope)));
  $('fScope').textContent=G.name;
  $('fSim').textContent=`Telescopius simulator setup: ${G.fl} mm focal length · ${G.ap} mm aperture · ${(G.w*PIX/1000).toFixed(2)} × ${(G.h*PIX/1000).toFixed(2)} mm sensor`;
  $('fScopeSpec').textContent=`${G.ap} mm · ${G.fl} mm f/${G.f.toFixed(1)} · ${(G.L/60).toFixed(2)}° × ${(G.S/60).toFixed(2)}° · ${G.scale.toFixed(2)}″/px`;
  const k=(G.f/5.3)**2;
  $('tipInteg').textContent = scope==='S30P'
    ? 'Integration times are set for the S30 Pro. Switch telescopes to rescale them.'
    : `Integration times are rescaled from the S30 Pro by focal ratio (×${k.toFixed(2)} for f/${G.f.toFixed(1)}). Per-pixel brightness depends on focal ratio, not aperture, so the times stay similar; ${G.ap===50?'the 50 mm scopes record finer detail (smaller ″/px) rather than needing much less time.':'the S30 records slightly coarser detail than the S30 Pro.'}`;
  $('tipBatt').textContent=`Split long totals across several nights in the same dark window and stack them. Battery is about ${G.batt} h, so plan one or two targets per night.`;
}
function setScope(k){ if(!SCOPES[k]) return; scope=k; G=scopeGeom(k); store.set('seestarScope',k); updateScopeUI(); render(); renderLook(); }
// Telescopius object page for a catalog id like "M31", "NGC 7000", "IC 1396". Returns null if the label isn't a plain catalog id.
function tsSlug(label){
  const m=String(label||'').trim().match(/^(M|NGC|IC)\s*(\d+)/i); if(m) return `${m[1].toLowerCase()}-${m[2]}`;
  const o=String(label||'').trim().match(/^(Sh2|Abell|Ced)[\s-]*(\d+)/i); if(o) return `${o[1].toLowerCase()}-${o[2]}`;
  return null;
}
function tsUrl(slug){ return `https://telescopius.com/deep-sky-objects/${slug}`; }
// Telescope Simulator centred on the target: ra in hours, dec in degrees (the format Telescopius uses in its own links).
function simUrl(o){
  const ra=o.raD!=null?+o.raD:(o.ra?raDeg(o.ra):null), dec=o.decD!=null?+o.decD:(o.dec?decDeg(o.dec):null);
  if(ra==null||dec==null||isNaN(ra)||isNaN(dec)) return null;
  return `https://telescopius.com/telescope-simulator?ra=${(ra/15).toFixed(6)}&dec=${dec.toFixed(6)}`;
}
function tsLinks(o){
  // [label, url] pairs; wiki labels that aren't catalog ids (e.g. "Leo Triplet") are skipped, falling back to ids from the members text.
  let out=o.wiki.map(w=>[w[0],tsSlug(w[0])]).filter(x=>x[1]);
  if(!out.length){ const ids=(String(o.members).match(/\b(M\d+|NGC\s?\d+|IC\s?\d+)/g)||[]); out=ids.map(i=>[i,tsSlug(i)]).filter(x=>x[1]); }
  return out.map(x=>[x[0],tsUrl(x[1])]);
}
function rowHTML(o,rank,r){
  const tl=tsLinks(o), su=simUrl(o), tl0=tl.length?tl[0]:o.wiki[0];
  const links=tl.map(w=>`<a class="ss-wk" href="${esc(w[1])}" target="_blank" rel="noopener">${esc(w[0])}</a>`);
  const nameHTML=`<a class="ss-wk" href="${esc(tl0[1])}" target="_blank" rel="noopener" title="Open the target page on Telescopius">${esc(o.name)}</a>`;
  const simHTML=su?`<a class="ss-sim" href="${esc(su)}" target="_blank" rel="noopener" title="Open the Telescopius telescope simulator centred on this target">Simulator ↗</a>`:'';
  const flt=`<span class="ss-flt ${o.filter==='Dual-band'?'ss-db':''}">${esc(o.filter)}</span>`;
  const fx=framing(o.frame); const fit=fx?`<span class="ss-fit ${fx.c}">${esc(fx.t)}</span>`:'';
  const size=o.sizeTxt??o.size;
  const head=`<span class="ss-rk">${rank}</span>
<div class="ss-c-nm"><div class="ss-nm">${nameHTML}${flt}${fit}${simHTML}</div><div class="ss-mem">${o.typ?`<span class="ss-typ">${esc(o.typ)}</span> · `:''}${esc(o.members)} · <span class="ss-mono">${esc(size)}${/\d$/.test(String(size))?'′':''}${o.mag!=null?` · mag ${o.mag}`:''}</span></div>${links.length>1?`<div class="ss-wl">Telescopius: ${links.join(' · ')}</div>`:''}<div class="ss-note">${esc(o.notes)}</div></div>`;
  if(!r.visible) return `<div class="ss-row ss-off" role="row">${head}<div class="ss-c-se"><span class="ss-lbl">Season</span><span class="ss-mono">Not visible</span></div><div class="ss-c-pk"></div><div class="ss-c-win"></div><div class="ss-c-alt"><span class="ss-lbl">Max alt</span><span class="ss-mono">${r.maxalt}°</span></div><div class="ss-c-int"><span class="ss-int">—</span></div></div>`;
  const tags=(r.maxalt<30?'<span class="ss-tag ss-low">low</span>':'')+(r.maxalt>80?'<span class="ss-tag ss-zen">near zenith</span>':'');
  const w=r.windows.length?r.windows.map((x,i)=>`<span class="ss-pill ${i===0&&r.winStart<10?'ss-now':''}">${esc(x)}</span>`).join(''):'<span class="ss-mono ss-dim">No moon-dark run</span>';
  return `<div class="ss-row" role="row">${head}
<div class="ss-c-se"><span class="ss-lbl">Season</span><span class="ss-mono">${r.season}</span></div>
<div class="ss-c-pk"><span class="ss-lbl">Peak</span><span class="ss-mono">${r.peak}</span><span class="ss-mono ss-dim">${r.peakH} h</span></div>
<div class="ss-c-win"><span class="ss-lbl">Dark</span><div class="ss-win">${w}</div></div>
<div class="ss-c-alt"><span class="ss-lbl">Max alt</span><span class="ss-mono">${r.maxalt}°</span>${tags}</div>
<div class="ss-c-int"><span class="ss-int">${esc(scaleInteg(o.integ))}</span>${o.est?'<span class="ss-tag ss-zen">estimate</span>':''}</div></div>`;
}
function render(){
  if(!S) return;
  const L=LISTS[kind], RR=R[kind];
  const q=$('q').value.trim().toLowerCase(), m=+$('mo').value, s=$('sort').value;
  let rows=L.map((t,i)=>({...t,rank:i+1,r:RR[i]})).filter(o=>(!q||(o.name+' '+o.members).toLowerCase().includes(q))&&(!m||(o.r.visible&&inSeason(o.r.season,m))));
  const key={fit:o=>{const x=framing(o.frame);return x?x.k*10-x.f:99;},rank:o=>o.rank,win:o=>o.r.visible?o.r.winStart:1e5,peak:o=>o.r.visible?o.r.peakKey:1e5,int:o=>parseFloat(o.integ)};
  rows.sort((a,b)=>{const x=key[s](a),y=key[s](b);return x<y?-1:x>y?1:a.rank-b.rank;});
  if(s==='rank') rows.sort((a,b)=>(b.r.visible-a.r.visible)||a.rank-b.rank);
  const total=rows.length, per=+$('per').value||total||1, pages=Math.max(1,Math.ceil(total/per));
  if(page>pages) page=pages;
  const from=(page-1)*per, shown=rows.slice(from,from+per);
  $('count').textContent=`${total} match · ${RR.filter(r=>r.visible).length} of 50 visible here`;
  renderPager(pages,total,from,shown.length);
  $('rows').innerHTML=shown.map(o=>rowHTML(o,o.rank,o.r)).join('');
}
function renderPager(pages,total,from,n){
  const el=$('pager');
  if(pages<=1){ el.innerHTML=''; return; }
  let h=`<button type="button" data-p="${page-1}" ${page===1?'disabled':''} aria-label="Previous page">‹<span class="ss-pt"> Prev</span></button>`;
  for(let i=1;i<=pages;i++) h+=`<button type="button" data-p="${i}" ${i===page?'aria-current="page"':''} aria-label="Page ${i}">${i}</button>`;
  h+=`<button type="button" data-p="${page+1}" ${page===pages?'disabled':''} aria-label="Next page"><span class="ss-pt">Next </span>›</button>`;
  h+=`<span class="ss-range">Showing ${from+1}–${from+n} of ${total}</span>`;
  el.innerHTML=h;
}
function renderLook(){ if(!S||!look.length) return; $('lookRes').innerHTML=look.map(t=>rowHTML(t,'•',evalTarget(S,t.raD,t.decD))).join(''); }
function doLook(){
  const q=$('lookQ').value;
  if(!q.trim()){ $('lookRes').innerHTML=''; return; }
  $('lookRes').innerHTML='<p class="ss-busy">Loading catalog…</p>';
  loadCatalog().then(()=>{
    const rows=lookupRows(q);
    if(!rows.length){ $('lookRes').innerHTML=`<p class="ss-err">No Messier, NGC or IC object matches “${esc(q)}”. Try a form like M27, NGC 7000 or IC 1805.</p>`; return; }
    look=rows.map(lookTarget); renderLook();
  }).catch(()=>{ $('lookRes').innerHTML='<p class="ss-err">The object catalog could not be loaded. Reload the page and try again.</p>'; });
}
function shortName(n){
  const p=n.split(',').map(s=>s.trim()).filter(s=>s&&!/^\d[\d -]*$/.test(s));
  if(p.length<2) return n;
  const last=p[p.length-1];
  if(/^United States/.test(last)&&p.length>=3) return `${p[0]}, ${p[p.length-2]}`;
  return `${p[0]}, ${last}`;
}
// ---- Town search: OpenStreetMap Nominatim, called directly from the browser.
// Usage policy: https://operations.osmfoundation.org/policies/nominatim/
// Explicit button only (no search-as-you-type), at most one request per second, results cached for the session.
const NOMINATIM='https://nominatim.openstreetmap.org/search';
let lastGeo=0, geoBusy=false;
function geoCache(){ try{ return JSON.parse(sessionStorage.getItem('seestarGeo')||'{}'); }catch(e){ return {}; } }
function geoStore(c){ try{ sessionStorage.setItem('seestarGeo',JSON.stringify(c)); }catch(e){} }
function geocode(q){
  const key=q.toLowerCase(), cache=geoCache();
  if(cache[key]) return Promise.resolve(cache[key]);
  const wait=Math.max(0,1100-(Date.now()-lastGeo));
  return new Promise(res=>setTimeout(res,wait)).then(()=>{
    lastGeo=Date.now();
    const url=`${NOMINATIM}?format=jsonv2&limit=6&addressdetails=0&q=${encodeURIComponent(q)}`;
    return fetch(url,{headers:{'Accept':'application/json','Accept-Language':'en'},referrerPolicy:'strict-origin-when-cross-origin',credentials:'omit'});
  }).then(r=>{
    if(r.status===429) throw new Error('Town search is busy. Try again in a minute.');
    if(!r.ok) throw new Error('Town search is unavailable right now.');
    return r.json();
  }).then(d=>{
    const out=(Array.isArray(d)?d:[]).slice(0,6).map(x=>({name:String(x.display_name||'').slice(0,160),lat:parseFloat(x.lat),lon:parseFloat(x.lon)}))
      .filter(p=>p.lat>=-90&&p.lat<=90&&p.lon>=-180&&p.lon<=180);
    cache[key]=out; geoStore(cache); return out;
  });
}
function doFind(){
  const q=$('findQ').value.trim().replace(/\s+/g,' ').slice(0,100), box=$('findRes');
  if(q.length<2){ box.innerHTML='<p class="ss-err">Type at least 2 characters.</p>'; return; }
  if(geoBusy) return;
  geoBusy=true; $('findBtn').disabled=true;
  box.innerHTML='<p class="ss-busy">Searching…</p>';
  geocode(q).then(res=>{
    if(!res.length){ box.innerHTML=`<p class="ss-err">No place found for “${esc(q)}”. Try adding the state or country.</p>`; return; }
    box.innerHTML=res.map((p,i)=>`<button type="button" class="ss-place" data-i="${i}">${esc(p.name)}<span class="ss-mono ss-dim"> ${p.lat.toFixed(2)}, ${p.lon.toFixed(2)}</span></button>`).join('');
    box.querySelectorAll('.ss-place').forEach(b=>b.addEventListener('click',()=>{
      const p=res[+b.dataset.i]; $('city').value=''; customName=shortName(p.name);
      $('lat').value=p.lat.toFixed(4); $('lon').value=p.lon.toFixed(4); box.innerHTML=''; run();
    }));
  }).catch(e=>{
    box.innerHTML=`<p class="ss-err">${esc(e&&e.message&&/^Town search/.test(e.message)?e.message:'Town search failed. Check your connection and try again.')}</p>`;
  }).finally(()=>{ geoBusy=false; $('findBtn').disabled=false; });
}
function init(){
  if(!$('app')) return;
  Promise.all([getJSON(OPT.lists),getJSON(OPT.cities)]).then(([l,c])=>{
    LISTS=l; CITIES=c;
    $('city').innerHTML='<option value="">Custom coordinates</option>'+CITIES.map((c,i)=>`<option value="${i}">${esc(c[0])}</option>`).join('');
    MON.forEach((m,i)=>{const o=document.createElement('option');o.value=i+1;o.textContent=m;$('mo').appendChild(o);});
    $('city').addEventListener('change',()=>{const c=CITIES[+$('city').value]; if($('city').value!==''&&c){ setFields({name:c[0],lat:c[1],lon:c[2]}); run(); }});
    ['lat','lon'].forEach(id=>$(id).addEventListener('input',()=>{ $('city').value=''; customName=null; }));
    $('locForm').addEventListener('submit',e=>{e.preventDefault(); run();});
    $('bG').addEventListener('click',()=>setKind('galaxies')); $('bN').addEventListener('click',()=>setKind('nebulae'));
    $('lookForm').addEventListener('submit',e=>{e.preventDefault(); doLook();});
    $('lookClear').addEventListener('click',()=>{ $('lookQ').value=''; look=[]; $('lookRes').innerHTML=''; });
    if($('findBtn')){
      $('findBtn').addEventListener('click',doFind);
      $('findQ').addEventListener('keydown',e=>{ if(e.key==='Enter'){ e.preventDefault(); doFind(); } });
    }
    ['q','mo','sort'].forEach(id=>$(id).addEventListener('input',()=>{ page=1; render(); }));
    const sp=store.get('seestarPer'); if(sp!==null&&[...$('per').options].some(o=>o.value===sp)) $('per').value=sp;
    $('per').addEventListener('input',()=>{ store.set('seestarPer',$('per').value); page=1; render(); });
    $('pager').addEventListener('click',e=>{
      const b=e.target.closest('button[data-p]'); if(!b||b.disabled) return;
      page=+b.dataset.p; render();
      document.querySelector('#ss-app [role="table"]').scrollIntoView({behavior:matchMedia('(prefers-reduced-motion: reduce)').matches?'auto':'smooth',block:'start'});
    });
    const k=store.get('seestarKind'); if(k==='galaxies'||k==='nebulae') kind=k;
    $('bG').setAttribute('aria-pressed',String(kind==='galaxies')); $('bN').setAttribute('aria-pressed',String(kind==='nebulae')); $('hKind').textContent=kind==='galaxies'?'Galaxies':'Nebulae';
    document.querySelectorAll('#ss-scopeSeg button').forEach(b=>b.addEventListener('click',()=>setScope(b.dataset.s)));
    const sk=store.get('seestarScope'); if(SCOPES[sk]){ scope=sk; G=scopeGeom(sk); }
    updateScopeUI();
    let saved=null; try{ saved=JSON.parse(store.get('seestarLoc')); }catch(e){}
    setFields(saved&&isFinite(saved.lat)&&isFinite(saved.lon)?saved:DEF);
    run();
  }).catch(()=>{ $('rows').innerHTML='<p class="ss-err">The planner data could not be loaded. Reload the page and try again.</p>'; });
}
if(document.readyState==='loading') document.addEventListener('DOMContentLoaded',init); else init();
})();
