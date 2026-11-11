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
    // For each dark sample: sidereal time, whether the Moon is harmless anyway (down or under 25% lit),
    // and the Moon's position so each target can check its distance from it.
    const ls=[], ok=[], mra=[], msd=[], mcd=[];
    for(let k=0;k<144;k++){
      const jd=jd0+k*STEP/1440, l=lstOf(jd,lon), s=sunPos(jd);
      if(altOf(s.ra,s.dec,lat,l)>-18) continue;
      const mo=moonPos(jd), ma=altOf(mo.ra,mo.dec,lat,l)-mo.par, kill=(1-cosD(mo.bet)*cosD(mo.lam-s.lam))/2;
      ls.push(l); ok.push(ma<0||kill<0.25); mra.push(mo.ra); msd.push(sinD(mo.dec)); mcd.push(cosD(mo.dec));
    }
    nights.push({ls,ok,mra,msd,mcd});
    const mm=moonPos(jd0+0.5), ss=sunPos(jd0+0.5); illum.push((1-cosD(mm.bet)*cosD(mm.lam-ss.lam))/2);
  }
  const newMoons=[]; for(let n=1;n<NIGHTS-1;n++) if(illum[n]<=illum[n-1]&&illum[n]<illum[n+1]&&illum[n]<0.05) newMoons.push(fmtD(dates[n]));
  return {lat,lon,nights,dates,newMoons};
}
const fmtD=d=>`${MON[d.getUTCMonth()]} ${d.getUTCDate()}`;
const fmtW=(a,b)=>a.getTime()===b.getTime()?`${fmtD(a)}, ${a.getUTCFullYear()}`:a.getUTCMonth()===b.getUTCMonth()?`${fmtD(a)}–${b.getUTCDate()}, ${b.getUTCFullYear()}`:`${fmtD(a)} – ${fmtD(b)}, ${b.getUTCFullYear()}`;
// Imaging threshold (Min altitude): visitor setting, default 25°. Targets that never get 10° above it
// use a threshold 8° below their highest point (never below 12°, never above the setting).
let THR=25;
const thrFor=maxalt=>maxalt>=THR+10?THR:Math.min(THR,Math.max(12,maxalt-8));
// Moon distance rule: a bright Moon above the horizon is acceptable when it is far enough from the target:
// at least 60° with the dual-band (light-pollution) filter, at least 90° with the UV/IR-cut filter.
const moonSepFor=f=>f==='Dual-band'?60:90;
function evalTarget(S,ra,dec,filter){
  const cosReq=cosD(moonSepFor(filter)), tsd=sinD(dec), tcd=cosD(dec);
  // Max altitude = highest point reached during astronomical darkness (Sun below -18°) over the 12 months,
  // not the geometric culmination, which can fall in daylight or twilight.
  const sd=sinD(dec)*sinD(S.lat), cd=cosD(dec)*cosD(S.lat);
  let maxV=-1;
  S.nights.forEach(nt=>{ for(let i=0;i<nt.ls.length;i++){ const v=sd+cd*cosD(nt.ls[i]-ra); if(v>maxV) maxV=v; } });
  const maxalt=Math.asin(Math.max(-1,Math.min(1,maxV)))/RAD, thr=thrFor(maxalt);
  if(maxalt<12) return {visible:false,maxalt:Math.round(maxalt)};
  const tot=new Float32Array(NIGHTS), mf=new Float32Array(NIGHTS), st=sinD(thr);
  S.nights.forEach((nt,n)=>{for(let i=0;i<nt.ls.length;i++){ if(sd+cd*cosD(nt.ls[i]-ra)>=st){tot[n]+=STEP/60; if(nt.ok[i]||tsd*nt.msd[i]+tcd*nt.mcd[i]*cosD(ra-nt.mra[i])<=cosReq) mf[n]+=STEP/60;} }});
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

// ---- Favorite locations (saved in this browser)
const FAVK='seestarFavs';
function favGet(){ try{ const a=JSON.parse(store.get(FAVK)||'[]'); return Array.isArray(a)?a.filter(f=>f&&f.name&&isFinite(f.lat)&&isFinite(f.lon)):[]; }catch(e){ return []; } }
const favKey=(lat,lon)=>(+lat).toFixed(3)+','+(+lon).toFixed(3);
function curLoc(){
  const lat=parseFloat($('lat').value), lon=parseFloat($('lon').value);
  if(!(lat>=-90&&lat<=90)||!(lon>=-180&&lon<=180)) return null;
  const c=CITIES[+$('city').value];
  return {name:$('city').value!==''&&c?c[0]:(customName||`${lat.toFixed(2)}, ${lon.toFixed(2)}`),lat,lon};
}
function favRender(){
  const a=favGet(), box=$('favs'), cur=curLoc(), saved=!!cur&&a.some(f=>favKey(f.lat,f.lon)===favKey(cur.lat,cur.lon));
  box.hidden=!a.length;
  box.innerHTML=a.length?'<span class="ss-favl">Favorites</span>'+a.map((f,i)=>`<span class="ss-fav${cur&&favKey(f.lat,f.lon)===favKey(cur.lat,cur.lon)?' ss-fav-on':''}"><button type="button" class="ss-favgo" data-i="${i}" title="Use ${esc(f.name)}">★ ${esc(f.name)}</button><button type="button" class="ss-favx" data-i="${i}" aria-label="Remove ${esc(f.name)} from favorites" title="Remove">×</button></span>`).join(''):'';
  const b=$('favBtn'); b.textContent=saved?'★ Saved':'☆ Save location'; b.setAttribute('aria-pressed',String(saved)); b.title=saved?'Remove this location from favorites':'Save this location to favorites';
}
function favInit(){
  $('favBtn').addEventListener('click',()=>{
    const cur=curLoc(); if(!cur){ $('err').textContent='Enter a valid latitude and longitude first.'; return; }
    const a=favGet(), k=favKey(cur.lat,cur.lon), i=a.findIndex(f=>favKey(f.lat,f.lon)===k);
    if(i>=0) a.splice(i,1); else a.push({name:cur.name,lat:+cur.lat.toFixed(4),lon:+cur.lon.toFixed(4)});
    store.set(FAVK,JSON.stringify(a.slice(0,12))); favRender();
  });
  $('favs').addEventListener('click',e=>{
    const go=e.target.closest('.ss-favgo'), x=e.target.closest('.ss-favx'); if(!go&&!x) return;
    const a=favGet(), f=a[+(go||x).dataset.i]; if(!f) return;
    if(x){ a.splice(+x.dataset.i,1); store.set(FAVK,JSON.stringify(a)); favRender(); }
    else { setFields(f); run(); }
  });
  ['lat','lon'].forEach(id=>$(id).addEventListener('input',favRender));
  favRender();
}
function setFields(l){ $('lat').value=(+l.lat).toFixed(4); $('lon').value=(+l.lon).toFixed(4); const i=CITIES.findIndex(c=>c[0]===l.name); $('city').value=i>=0?i:''; customName=i>=0?null:(l.name||null); }
function setKind(k){ kind=k; page=1; $('bG').setAttribute('aria-pressed',String(k==='galaxies')); $('bN').setAttribute('aria-pressed',String(k==='nebulae')); $('hKind').textContent=k==='galaxies'?'Galaxies':'Nebulae'; store.set('seestarKind',k); render(); }
function inSeason(s,m){ if(s==="Year-round")return true; const [a,b]=s.split('–').map(x=>MON.indexOf(x)+1); const e=b||a; return a<=e?(m>=a&&m<=e):(m>=a||m<=e); }
function recompute(){
  if(!S) return;
  for(const k of ['galaxies','nebulae']) R[k]=LISTS[k].map(t=>evalTarget(S,raDeg(t.ra),decDeg(t.dec),t.filter));
}
function setThr(v){
  v=+v; if(!(v>=10&&v<=50)||v%5) v=25; THR=v; $('thr').value=String(v);
  document.querySelectorAll('#ss-app .ss-thrv').forEach(e=>e.textContent=v+'°');
  document.querySelectorAll('#ss-app .ss-thrlow').forEach(e=>e.textContent=(v+5)+'°');
}
function run(){
  const lat=parseFloat($('lat').value), lon=parseFloat($('lon').value);
  if(!(lat>=-90&&lat<=90)||!(lon>=-180&&lon<=180)){ $('err').textContent='Enter latitude between -90 and 90 and longitude between -180 and 180 (west is negative).'; return; }
  $('err').textContent='';
  const c=CITIES[+$('city').value];
  const loc={name:$('city').value!==''&&c?c[0]:(customName||'Custom location'),lat,lon};
  store.set('seestarLoc',JSON.stringify(loc)); favRender();
  $('rows').innerHTML='<p class="ss-busy">Calculating…</p>';
  setTimeout(()=>{
    S=prepSite(lat,lon,start); R={}; page=1;
    recompute();
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
function scopeGeom(k){ const sc=SCOPES[k], scale=206.265*(sc.pix||PIX)/sc.fl; return {...sc, scale, L:sc.w*scale/60, S:sc.h*scale/60}; }
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
  $('fSim').textContent=`Telescopius simulator setup: ${G.fl} mm focal length · ${G.ap} mm aperture · ${(G.w*(G.pix||PIX)/1000).toFixed(2)} × ${(G.h*(G.pix||PIX)/1000).toFixed(2)} mm sensor`;
  $('fScopeSpec').textContent=`${G.ap} mm · ${G.fl} mm f/${G.f.toFixed(1)} · ${(G.L/60).toFixed(2)}° × ${(G.S/60).toFixed(2)}° · ${G.scale.toFixed(2)}″/px`;
  const k=(G.f/5.3)**2;
  $('tipInteg').textContent = scope==='S30P'
    ? 'Integration times are set for the S30 Pro. Switch telescopes to rescale them.'
    : `Integration times are rescaled from the S30 Pro by focal ratio (×${k.toFixed(2)} for f/${G.f.toFixed(1)}). Per-pixel brightness depends on focal ratio, not aperture, so the times stay similar; ${G.ap===50?'the 50 mm scopes record finer detail (smaller ″/px) rather than needing much less time.':'the S30 records slightly coarser detail than the S30 Pro.'}`;
  $('tipBatt').textContent=`Split long totals across several nights in the same dark window and stack them. Battery is about ${G.batt} h, so plan one or two targets per night.`;
}
// Custom telescope/camera: aperture, focal length, sensor pixels and pixel size, saved in this browser.
function custLoad(){ try{ const c=JSON.parse(store.get('seestarCustomScope')||'null'); if(c&&[c.ap,c.fl,c.w,c.h,c.pix].every(v=>isFinite(v)&&v>0)) SCOPES.CUSTOM={name:String(c.name||'Custom').slice(0,40),ap:+c.ap,fl:+c.fl,f:+c.fl/+c.ap,w:+c.w,h:+c.h,pix:+c.pix,batt:6}; }catch(e){} }
function custLabel(){ const b=$('scopeCust'); if(b) b.textContent=SCOPES.CUSTOM?SCOPES.CUSTOM.name:'Custom'; }
function custInit(){
  const f=$('custForm'), b=$('scopeCust'); if(!f||!b) return;
  custLoad(); custLabel();
  const fill=()=>{ const c=SCOPES.CUSTOM||SCOPES[scope]||SCOPES.S30P; f.elements.name.value=SCOPES.CUSTOM?c.name:''; f.elements.ap.value=c.ap; f.elements.fl.value=c.fl; f.elements.w.value=c.w; f.elements.h.value=c.h; f.elements.pix.value=c.pix||PIX; };
  b.addEventListener('click',()=>{ if(SCOPES.CUSTOM&&scope!=='CUSTOM'){ setScope('CUSTOM'); return; } f.hidden=!f.hidden; if(!f.hidden){ fill(); f.elements.ap.focus(); } });
  f.addEventListener('submit',e=>{ e.preventDefault();
    const v=n=>+f.elements[n].value, ok=['ap','fl','w','h','pix'].every(n=>isFinite(v(n))&&v(n)>0);
    if(!ok){ $('custErr').textContent='Enter a number greater than 0 in every field.'; return; }
    $('custErr').textContent='';
    store.set('seestarCustomScope',JSON.stringify({name:f.elements.name.value.trim()||'Custom',ap:v('ap'),fl:v('fl'),w:v('w'),h:v('h'),pix:v('pix')}));
    custLoad(); custLabel(); setScope('CUSTOM'); f.hidden=true; });
  $('custCancel').addEventListener('click',()=>{ f.hidden=true; });
  $('custClear').addEventListener('click',()=>{ store.set('seestarCustomScope',''); delete SCOPES.CUSTOM; custLabel(); f.hidden=true; if(scope==='CUSTOM') setScope('S30P'); });
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
// ---- Tonight's altitude chart (per target). Time axis is local mean solar time at the chosen site.
// ---- Tonight's altitude chart. Times are shown in the device's time zone.
const CW=300, CHT=70, A0=-10, A1=90, CSTEP=10, NIGHT_H=19;
const CMP=['N','NE','E','SE','S','SW','W','NW'];
const fmtClock=ms=>new Date(ms).toLocaleTimeString([], {hour:'numeric',minute:'2-digit'});
const azOf=(ra,dec,lat,l)=>{ const H=l-ra; return norm(Math.atan2(sinD(H),cosD(H)*sinD(lat)-Math.tan(dec*RAD)*cosD(lat))/RAD+180); };
let nightCache=null;
function tonight(){
  // "Tonight" runs from 3 pm to 10 am local solar time at the site, starting from the most recent site noon.
  const now=Date.now(), off=S.lon/15*3600e3, lmt=new Date(now+off);
  let noon=Date.UTC(lmt.getUTCFullYear(),lmt.getUTCMonth(),lmt.getUTCDate(),12)-off;
  if(now<noon) noon-=864e5;
  const key=`${S.lat},${S.lon},${noon}`;
  if(nightCache&&nightCache.key===key) return nightCache;
  const t0=noon+3*3600e3, n=Math.round(NIGHT_H*60/CSTEP)+1, ls=[], sun=[], moon=[], mpos=[];
  let illum=0;
  for(let k=0;k<n;k++){
    const jd=(t0+k*CSTEP*6e4)/864e5+2440587.5, l=lstOf(jd,S.lon), sp=sunPos(jd), mp=moonPos(jd);
    ls.push(l); sun.push(altOf(sp.ra,sp.dec,S.lat,l)); moon.push(altOf(mp.ra,mp.dec,S.lat,l)-mp.par); mpos.push([mp.ra,mp.dec,(1-cosD(mp.bet)*cosD(mp.lam-sp.lam))/2]);
    if(k===Math.floor(n/2)) illum=(1-cosD(mp.bet)*cosD(mp.lam-sp.lam))/2;
  }
  // Sky background: merge runs of the same twilight phase into one band.
  const tone=a=>a>-0.83?0:a>-6?1:a>-12?2:a>-18?3:4, cols=['#1d3150','#2a2443','#1a1b33','#121424','#07080d'], dx=CW/(n-1);
  let bands='', st=0;
  for(let k=1;k<=n;k++){ if(k===n||tone(sun[k])!==tone(sun[st])){ const x0=Math.max(0,(st-0.5)*dx), x1=Math.min(CW,(k-0.5)*dx); bands+=`<rect x="${x0.toFixed(1)}" y="0" width="${(x1-x0).toFixed(1)}" height="${CHT}" fill="${cols[tone(sun[st])]}"/>`; st=k; } }
  let ticks='';
  for(let ms=Math.ceil(t0/18e5)*18e5; ms<=t0+NIGHT_H*36e5; ms+=18e5){
    const d=new Date(ms); if(d.getMinutes()!==0||d.getHours()%3!==0) continue;
    const x=(ms-t0)/(NIGHT_H*36e5)*CW; if(x<14||x>CW-14) continue;
    ticks+=`<line x1="${x.toFixed(1)}" x2="${x.toFixed(1)}" y1="${CHT}" y2="${CHT+2.5}" stroke="#5c6070" stroke-width=".7"/><text x="${x.toFixed(1)}" y="${CHT+10}" text-anchor="middle" font-size="7.5" fill="#8b90a3">${d.toLocaleTimeString([], {hour:'numeric'})}</text>`;
  }
  return nightCache={key,t0,n,ls,sun,moon,mpos,illum,bands,ticks};
}
// Tonight: closest approach of a bright, risen Moon during the target's usable dark time.
function moonTonight(ra,dec,filter){
  if(!S||!isFinite(ra)||!isFinite(dec)) return null;
  const N=tonight(), maxalt=90-Math.abs(S.lat-dec), thr=thrFor(maxalt), req=moonSepFor(filter);
  let usable=false, minSep=null;
  for(let k=0;k<N.n;k++){
    if(N.sun[k]>-18||altOf(ra,dec,S.lat,N.ls[k])<thr) continue;
    usable=true;
    const [mra,mdec,ill]=N.mpos[k]; if(N.moon[k]<0||ill<0.25) continue;
    const c=sinD(dec)*sinD(mdec)+cosD(dec)*cosD(mdec)*cosD(ra-mra), sep=Math.acos(Math.max(-1,Math.min(1,c)))/RAD;
    if(minSep===null||sep<minSep) minSep=sep;
  }
  if(!usable) return {usable:false,req};
  return minSep===null?{usable:true,clear:true,req}:{usable:true,clear:false,sep:Math.round(minSep),ok:minSep>=req,req};
}
function moonLine(m,filter){
  if(!m||!m.usable) return '';
  const f=filter==='Dual-band'?'dual-band':'UV/IR';
  if(m.clear) return `<div class="ss-altmoon ss-ok">Moon: no interference tonight (down or under 25% lit)</div>`;
  return m.ok?`<div class="ss-altmoon ss-ok">Moon: at least ${m.sep}° away tonight · OK for ${f} (needs ${m.req}°)</div>`
             :`<div class="ss-altmoon ss-bad">Moon: comes within ${m.sep}° tonight · too close for ${f} (needs ${m.req}°)</div>`;
}
function tonightMax(ra,dec){
  if(!S||!isFinite(ra)||!isFinite(dec)) return null;
  const N=tonight(); let m=null;
  for(let k=0;k<N.n;k++){ if(N.sun[k]>-18) continue; const a=altOf(ra,dec,S.lat,N.ls[k]); if(m===null||a>m) m=a; }
  return m; // null if there is no full darkness tonight
}
const raOf=o=>o.raD!=null?+o.raD:raDeg(o.ra), decOf=o=>o.decD!=null?+o.decD:decDeg(o.dec);
function altChart(ra,dec,filter){
  if(!S||!isFinite(ra)||!isFinite(dec)) return '';
  const N=tonight(), n=N.n, dx=CW/(n-1), y=a=>(A1-Math.max(A0,Math.min(A1,a)))/(A1-A0)*CHT;
  const maxalt=90-Math.abs(S.lat-dec), thr=thrFor(maxalt);
  const alts=N.ls.map(l=>altOf(ra,dec,S.lat,l));
  // Peak = highest point during astronomical darkness, which is what matters for imaging.
  let best=-1, good=0;
  alts.forEach((a,k)=>{ if(N.sun[k]<=-18){ if(best<0||a>alts[best]) best=k; if(a>=thr) good+=CSTEP; } });
  if(best<0) best=alts.indexOf(Math.max(...alts));
  const path=arr=>arr.map((a,k)=>`${k?'L':'M'}${(k*dx).toFixed(1)},${y(a).toFixed(1)}`).join('');
  const nowF=(Date.now()-N.t0)/(NIGHT_H*36e5);
  const now=nowF>0&&nowF<1?`<line x1="${(nowF*CW).toFixed(1)}" x2="${(nowF*CW).toFixed(1)}" y1="0" y2="${CHT}" stroke="#4aa8ff" stroke-width=".9" stroke-dasharray="3 2"/><text x="${(nowF*CW+2).toFixed(1)}" y="8" font-size="7" font-weight="700" fill="#4aa8ff">now</text>`:'';
  const up=alts[best]>0&&N.sun[best]<=-18, bx=best*dx, by=y(alts[best]);
  const az=alts.map((a,k)=>Math.round(azOf(ra,dec,S.lat,N.ls[k])));
  const moonPct=Math.round(N.illum*100);
  const cap=up?`Tonight: highest in darkness ${Math.round(alts[best])}° ${CMP[Math.round(az[best]/45)%8]} at ${fmtClock(N.t0+best*CSTEP*6e4)} · ${(good/60).toFixed(1)} h above ${Math.round(thr)}° in darkness · Moon ${moonPct}%`
               :`Not above the horizon in darkness tonight · Moon ${moonPct}%`;
  return `<div class="ss-altw"><svg class="ss-altc" viewBox="0 0 ${CW} ${CHT+13}" role="img" aria-label="${esc(cap)}" data-t0="${N.t0}" data-n="${n}" data-best="${best}" data-up="${up?1:0}" data-a="${alts.map(a=>a.toFixed(1)).join(',')}" data-z="${az.join(',')}">
${N.bands}<line x1="0" x2="${CW}" y1="${y(thr).toFixed(1)}" y2="${y(thr).toFixed(1)}" stroke="#7d869f" stroke-width=".7" stroke-dasharray="3 3" opacity=".75"/><line x1="0" x2="${CW}" y1="${y(0).toFixed(1)}" y2="${y(0).toFixed(1)}" stroke="#3a3f55" stroke-width=".6"/>
<path d="${path(N.moon)}" fill="none" stroke="#9aa0b3" stroke-width="1" stroke-dasharray="2 2" opacity=".6"/>${up?`<path d="${path(alts)}" fill="none" stroke="#ff5a4e" stroke-width="1.8" stroke-linejoin="round"/>`:''}${now}
<line class="ss-acur" x1="${bx.toFixed(1)}" x2="${bx.toFixed(1)}" y1="0" y2="${CHT}" stroke="#fff" stroke-width=".6" stroke-dasharray="2 2" opacity="0"/><circle class="ss-adot" cx="${bx.toFixed(1)}" cy="${by.toFixed(1)}" r="3" fill="#fff" stroke="#ff5a4e" stroke-width="1.6" visibility="${up?'visible':'hidden'}"/>
${N.ticks}</svg><div class="ss-altcap" data-def="${esc(cap)}">${esc(cap)}</div>${moonLine(moonTonight(ra,dec,filter),filter)}</div>`;
}
// Hover or drag across a chart to read altitude, direction and time.
(function(){
  let active=null;
  const reset=sv=>{ const w=sv.parentNode, k=+sv.dataset.best, dx=CW/(+sv.dataset.n-1), a=+sv.dataset.a.split(',')[k];
    w.classList.remove('ss-on'); w.querySelector('.ss-altcap').textContent=w.querySelector('.ss-altcap').dataset.def;
    sv.querySelector('.ss-acur').setAttribute('opacity','0');
    const d=sv.querySelector('.ss-adot'); d.setAttribute('cx',(k*dx).toFixed(1)); d.setAttribute('cy',((A1-Math.max(A0,Math.min(A1,a)))/(A1-A0)*CHT).toFixed(1)); d.setAttribute('visibility',sv.dataset.up==='1'?'visible':'hidden'); };
  document.addEventListener('pointermove',e=>{
    const sv=e.target.closest?e.target.closest('.ss-altc'):null;
    if(active&&active!==sv){ reset(active); active=null; }
    if(!sv) return;
    active=sv;
    if(!sv._a){ sv._a=sv.dataset.a.split(',').map(Number); sv._z=sv.dataset.z.split(',').map(Number); }
    const n=+sv.dataset.n, dx=CW/(n-1), r=sv.getBoundingClientRect(), k=Math.max(0,Math.min(n-1,Math.round((e.clientX-r.left)/r.width*CW/dx)));
    const a=sv._a[k], x=(k*dx).toFixed(1), t=fmtClock(+sv.dataset.t0+k*CSTEP*6e4);
    sv.querySelector('.ss-acur').setAttribute('x1',x); sv.querySelector('.ss-acur').setAttribute('x2',x); sv.querySelector('.ss-acur').setAttribute('opacity','.6');
    const d=sv.querySelector('.ss-adot'); d.setAttribute('cx',x); d.setAttribute('cy',((A1-Math.max(A0,Math.min(A1,a)))/(A1-A0)*CHT).toFixed(1)); d.setAttribute('visibility','visible');
    sv.parentNode.classList.add('ss-on');
    sv.parentNode.querySelector('.ss-altcap').textContent= a>0?`${t} · altitude ${Math.round(a)}° · ${CMP[Math.round(sv._z[k]/45)%8]}`:`${t} · below the horizon`;
  });
  const off=()=>{ if(active){ reset(active); active=null; } };
  document.addEventListener('pointercancel',off);
  document.addEventListener('pointerup',e=>{ if(e.pointerType!=='mouse') off(); });
  document.documentElement.addEventListener('mouseleave',off);
})();
let tFavs=null, favOnly=false;
function tfGet(){ if(!tFavs){ try{ const a=JSON.parse(store.get('seestarFavTargets')||'[]'); tFavs=new Set(Array.isArray(a)?a:[]); }catch(e){ tFavs=new Set(); } } return tFavs; }
function tfSave(){ store.set('seestarFavTargets',JSON.stringify([...tfGet()])); }
function starBtn(name){ const on=tfGet().has(name); return `<button type="button" class="ss-star${on?' ss-on':''}" data-t="${esc(name)}" aria-pressed="${on}" aria-label="${on?'Remove':'Add'} ${esc(name)} ${on?'from':'to'} favorites" title="${on?'Remove from favorites':'Favorite this target'}">${on?'★':'☆'}</button>`; }
function tfBtnUpdate(){ const b=$('favOnly'); if(!b) return; const n=tfGet().size; b.textContent=`${favOnly?'★':'☆'} Favorites${n?' ('+n+')':''}`; b.setAttribute('aria-pressed',String(favOnly)); b.classList.toggle('ss-on',favOnly); }
function tfInit(){
  tfBtnUpdate();
  $('favOnly').addEventListener('click',()=>{ favOnly=!favOnly; tfBtnUpdate(); page=1; render(); });
  document.addEventListener('click',e=>{
    const b=e.target.closest&&e.target.closest('.ss-star'); if(!b) return;
    const n=b.dataset.t, f=tfGet(); f.has(n)?f.delete(n):f.add(n); tfSave(); tfBtnUpdate();
    if(favOnly){ render(); return; }
    document.querySelectorAll('.ss-star').forEach(x=>{ if(x.dataset.t===n){ const on=f.has(n); x.classList.toggle('ss-on',on); x.setAttribute('aria-pressed',String(on)); x.textContent=on?'★':'☆'; x.title=on?'Remove from favorites':'Favorite this target'; } });
  });
}
function rowHTML(o,rank,r){
  const tl=tsLinks(o), su=simUrl(o), tl0=tl.length?tl[0]:o.wiki[0];
  const links=tl.map(w=>`<a class="ss-wk" href="${esc(w[1])}" target="_blank" rel="noopener">${esc(w[0])}</a>`);
  const nameHTML=`<a class="ss-wk" href="${esc(tl0[1])}" target="_blank" rel="noopener" title="Open the target page on Telescopius">${esc(o.name)}</a>`;
  const simHTML=su?`<a class="ss-sim" href="${esc(su)}" target="_blank" rel="noopener" title="Open the Telescopius telescope simulator centred on this target">Simulator ↗</a>`:'';
  const flt=`<span class="ss-flt ${o.filter==='Dual-band'?'ss-db':''}">${esc(o.filter)}</span>`;
  const fx=framing(o.frame); const fit=fx?`<span class="ss-fit ${fx.c}">${esc(fx.t)}</span>`:'';
  const size=o.sizeTxt??o.size;
  const head=`<span class="ss-rk">${rank}${starBtn(o.name)}</span>
<div class="ss-c-nm"><div class="ss-nm">${nameHTML}${flt}${fit}${simHTML}</div><div class="ss-mem">${o.typ?`<span class="ss-typ">${esc(o.typ)}</span> · `:''}${esc(o.members)} · <span class="ss-mono">${esc(size)}${/\d$/.test(String(size))?'′':''}${o.mag!=null?` · mag ${o.mag}`:''}</span></div>${links.length>1?`<div class="ss-wl">Telescopius: ${links.join(' · ')}</div>`:''}<div class="ss-note">${esc(o.notes)}</div></div>`;
  if(!r.visible) return `<div class="ss-row ss-off" role="row">${head}<div class="ss-c-se"><span class="ss-lbl">Season</span><span class="ss-mono">Not visible</span></div><div class="ss-c-pk"></div><div class="ss-c-win"></div><div class="ss-c-alt"><span class="ss-lbl">Max alt</span><span class="ss-mono">${r.maxalt}°</span></div><div class="ss-c-int"><span class="ss-int">—</span></div></div>`;
  const tags=(r.maxalt<THR+5?'<span class="ss-tag ss-low">low</span>':'')+(r.maxalt>80?'<span class="ss-tag ss-zen">near zenith</span>':'');
  const w=r.windows.length?r.windows.map((x,i)=>`<span class="ss-pill ${i===0&&r.winStart<10?'ss-now':''}">${esc(x)}</span>`).join(''):'<span class="ss-mono ss-dim">No moon-dark run</span>';
  return `<div class="ss-row" role="row">${head}
<div class="ss-c-se"><span class="ss-lbl">Season</span><span class="ss-mono">${r.season}</span></div>
<div class="ss-c-pk"><span class="ss-lbl">Peak</span><span class="ss-mono">${r.peak}</span><span class="ss-mono ss-dim">${r.peakH} h</span></div>
<div class="ss-c-win"><span class="ss-lbl">Dark</span><div class="ss-win">${w}</div>${altChart(raOf(o),decOf(o),o.filter)}</div>
<div class="ss-c-alt"><span class="ss-lbl">Max alt</span>${(()=>{ const t=tonightMax(raOf(o),decOf(o)); return t===null?'<span class="ss-mono">—</span><span class="ss-tag ss-zen">no full dark tonight</span>':t>0?`<span class="ss-mono">${Math.round(t)}°</span><span class="ss-tag ss-zen">tonight</span>`:'<span class="ss-mono">—</span><span class="ss-tag ss-zen">below horizon tonight</span>'; })()}<span class="ss-tag ss-zen">year ${r.maxalt}°</span>${tags}</div>
<div class="ss-c-int"><span class="ss-int">${esc(scaleInteg(o.integ))}</span>${o.est?'<span class="ss-tag ss-zen">estimate</span>':''}</div></div>`;
}
function render(){
  if(!S) return;
  const L=LISTS[kind], RR=R[kind];
  const q=$('q').value.trim().toLowerCase(), m=+$('mo').value, s=$('sort').value;
  let rows=L.map((t,i)=>({...t,rank:i+1,r:RR[i]})).filter(o=>(!favOnly||tfGet().has(o.name))&&(!q||(o.name+' '+o.members).toLowerCase().includes(q))&&(!m||(o.r.visible&&inSeason(o.r.season,m))));
  const key={moon:o=>{ if(!o.r.visible) return 1e5; const m=moonTonight(raOf(o),decOf(o),o.filter); if(!m||!m.usable) return 5e4; return m.clear?-1000:(m.ok?0:1000)-m.sep; },alt:o=>{ if(!o.r.visible) return 1e5; const t=tonightMax(raOf(o),decOf(o)); return t===null||t<=0?1e4-o.r.maxalt:-t; },fit:o=>{const x=framing(o.frame);return x?x.k*10-x.f:99;},rank:o=>o.rank,win:o=>o.r.visible?o.r.winStart:1e5,peak:o=>o.r.visible?o.r.peakKey:1e5,int:o=>parseFloat(o.integ)};
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
function renderLook(){ if(!S||!look.length) return; $('lookRes').innerHTML=look.map(t=>rowHTML(t,'•',evalTarget(S,t.raD,t.decD,t.filter))).join(''); }
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
// Back-to-top button: appears after scrolling down; returns to the top and moves focus to the heading.
function initTop(){
  const b=document.getElementById('ss-toTop'); if(!b) return;
  let ticking=false;
  const update=()=>{ b.classList.toggle('ss-show',window.scrollY>600); ticking=false; };
  window.addEventListener('scroll',()=>{ if(!ticking){ ticking=true; requestAnimationFrame(update); } },{passive:true});
  b.addEventListener('click',()=>{
    const reduce=matchMedia('(prefers-reduced-motion: reduce)').matches;
    window.scrollTo({top:0,behavior:reduce?'auto':'smooth'});
    const h=document.querySelector('#ss-app .ss-h1'); if(h){ h.setAttribute('tabindex','-1'); h.focus({preventScroll:true}); }
  });
  update();
}
function init(){
  initTop(); tfInit(); chartModalInit();
  if(!$('app')) return;
  Promise.all([getJSON(OPT.lists),getJSON(OPT.cities)]).then(([l,c])=>{
    LISTS=l; CITIES=c;
    $('city').innerHTML='<option value="">Custom coordinates</option>'+CITIES.map((c,i)=>`<option value="${i}">${esc(c[0])}</option>`).join('');
    MON.forEach((m,i)=>{const o=document.createElement('option');o.value=i+1;o.textContent=m;$('mo').appendChild(o);});
    $('mo').value=String(new Date().getMonth()+1); // default to the current month; "Any month" is still one click away
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
    document.querySelectorAll('#ss-scopeSeg button:not(#ss-scopeCust)').forEach(b=>b.addEventListener('click',()=>setScope(b.dataset.s)));
    custInit();
    $('thr').addEventListener('input',()=>{ setThr($('thr').value); store.set('seestarThr',$('thr').value); recompute(); page=1; render(); renderLook(); });
    const st=store.get('seestarThr'); setThr(st!==null?st:(DEF.thr||25));
    const sk=store.get('seestarScope'); if(SCOPES[sk]){ scope=sk; G=scopeGeom(sk); }
    // Red-light (night-vision) toggle; remembered per browser.
    (function(){ const b=$('redBtn'); if(!b) return; const set=on=>{ document.documentElement.classList.toggle('red-light',on); b.setAttribute('aria-pressed',String(on)); store.set('seestarRed',on?'1':'0'); };
      b.addEventListener('click',()=>set(!document.documentElement.classList.contains('red-light')));
      if(store.get('seestarRed')==='1') set(true); })();
    updateScopeUI();
    try{ $('tipChart').textContent+=` Chart times use this device's time zone (${Intl.DateTimeFormat().resolvedOptions().timeZone}).`; }catch(e){}
    let saved=null; try{ saved=JSON.parse(store.get('seestarLoc')); }catch(e){}
    setFields(saved&&isFinite(saved.lat)&&isFinite(saved.lon)?saved:DEF);
    favInit();
    run();
  }).catch(()=>{ $('rows').innerHTML='<p class="ss-err">The planner data could not be loaded. Reload the page and try again.</p>'; });
}

// Click a chart to open a larger copy (same scrubbing). Esc, the backdrop or the close button dismisses it.
function chartModalInit(){
  const close=()=>{ const m=document.getElementById('ss-modal'); if(m) m.remove(); };
  document.addEventListener('click',e=>{
    if(e.target.closest&&e.target.closest('.ss-modal-x')) return close();
    const bd=e.target.closest&&e.target.closest('.ss-modal'); if(bd){ if(e.target===bd) close(); return; }
    const sv=e.target.closest&&e.target.closest('.ss-altc'); if(!sv) return;
    close();
    const row=sv.closest('.ss-row'), nm=row&&row.querySelector('.ss-nm a'), cap=sv.parentNode.querySelector('.ss-altcap');
    const m=document.createElement('div'); m.id='ss-modal'; m.className='ss-modal'; m.setAttribute('role','dialog'); m.setAttribute('aria-label','Altitude chart');
    m.innerHTML=`<div class="ss-modal-box"><button class="ss-modal-x" type="button" aria-label="Close">×</button><div class="ss-modal-t">${esc(nm?nm.textContent:'Altitude tonight')}</div></div>`;
    const w=sv.parentNode.cloneNode(true); w.classList.add('ss-altw-big'); m.firstChild.appendChild(w);
    (document.getElementById('ss-app')||document.body).appendChild(m);
  });
  document.addEventListener('keydown',e=>{ if(e.key==='Escape') close(); });
}

if(document.readyState==='loading') document.addEventListener('DOMContentLoaded',init); else init();
})();
