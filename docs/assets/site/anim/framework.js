/**
 * Native SVG teaching figures for the ANM framework page.
 * All geometry is schematic, not a measured trajectory or experiment replay.
 * The original prose, experimental tables and claims ledger remain authoritative.
 * A single visibility-aware clock drives the visible figures. Controls support
 * pause, replay, scrubbing, phase selection, keyboard and reduced motion.
 */
const NS = 'http://www.w3.org/2000/svg';
const clamp = (v,lo=0,hi=1) => Math.max(lo,Math.min(hi,v));
const esc = s => String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const text = (x,y,s,cls='',anchor='start') => `<text x="${x}" y="${y}" class="${cls}" text-anchor="${anchor}">${esc(s)}</text>`;
const path = (d,cls='edge',extra='') => `<path d="${d}" class="${cls}" ${extra}/>`;
const box = (x,y,w,h,label,sub='',cls='') => `<rect x="${x}" y="${y}" width="${w}" height="${h}" rx="13" class="node ${cls}"/>${text(x+16,y+29,label)}${sub?text(x+16,y+52,sub,'small'):''}`;
const dot = (x,y,r=6,cls='pulse') => `<circle cx="${x}" cy="${y}" r="${r}" class="${cls}"/>`;
const edge = (x,y,x2,y2,on=false,dash=false) => path(`M${x},${y} L${x2},${y2}`,`edge${on?' active':''}${dash?' dashed':''}`);
const arrow = (x,y,x2,y2,on=false,dash=false) => { const a=Math.atan2(y2-y,x2-x),c=Math.cos(a),s=Math.sin(a);return edge(x,y,x2-6*c,y2-6*s,on,dash)+path(`M${x2-14*c+5*s},${y2-14*s-5*c} L${x2-5*c},${y2-5*s} L${x2-14*c-5*s},${y2-14*s+5*c}`,on?'edge active':'edge'); };
const phaseOn = (p,n) => Math.floor(p)===n?'active':'';
const sweep = (x,y,x2,y2,p,n) => { const t=clamp(p-n);return p>=n&&p<n+1?dot(x+(x2-x)*t,y+(y2-y)*t):''; };
const frame = inner => `<svg xmlns="${NS}" class="fw-svg" viewBox="0 0 720 420" role="img">${inner}</svg>`;
const heading = (a,b) => text(30,42,a,'label')+text(30,65,b,'small');
const footer = s => text(30,397,s,'tiny');

function graph(p){
 const events=[['E1','support','A',132],['E2','contradiction','B',214],['E3','uncertainty ↓','A',296]];
 let s=heading('A typed record becomes inspectable evidence','Illustrative values • fixed open-loop updates');
 s+=box(30,127,143,78,'Source record','type / value / source',phaseOn(p,0));
 s+=box(30,260,143,78,'Schema check','validate record',phaseOn(p,1));
 s+=edge(101,205,101,255,p>=1)+sweep(101,205,101,255,p,0);
 for(const [id,type,target,y] of events){
  const yy=target==='A'?160:286;
  s+=path(`M174,300 C220,300 205,${y} 251,${y}`,'edge'+(p>=1?' active':''));
  s+=box(253,y-28,178,64,`${id}  ${type}`,`target ${target} · provenance`,p>=1?'active':'');
  s+=path(`M434,${y} C480,${y} 468,${yy} 520,${yy}`,'edge'+(p>=2?' active':''));
  if(p>=2&&p<3){let t=p-2;s+=dot(434+(520-434)*t,y+(yy-y)*t);}
 }
 s+=box(523,124,164,74,'Answer A','retained value',p>=2?'active':'');
 s+=box(523,250,164,74,'Answer B','retained value',p>=2?'active':'');
 s+=text(525,355,p>=3?'Separate readout → call':'Fixed graph update','label');
 return s+footer('EVENT TYPE + TARGET + SIGN + VALUE + TIME + PROVENANCE');
}
function protocol(p){
 let s=heading('One controlled push. Two separate measurements.','Source → state and source → output are measured separately');
 s+=box(42,170,162,86,'Source u','controlled ±ε',phaseOn(p,0));
 s+=box(286,104,177,80,'Candidate z','measure Δz',phaseOn(p,1));
 s+=box(286,259,177,80,'Readout O','measure ΔO',phaseOn(p,2));
 s+=path('M205,210 C250,210 243,145 280,145','edge'+(p>=1?' active':''));
 s+=path('M205,210 C250,210 243,300 280,300','edge'+(p>=2?' active':''));
 s+=edge(374,190,374,250,false,true)+text(397,221,'mediation?','small');
 s+=box(520,181,163,96,'Declared bounds','works / exact',phaseOn(p,3));
 s+=arrow(467,300,599,300,p>=3)+edge(599,300,599,281,p>=3);
 s+=dot(90+35*Math.sin(p*2),217,7);
 return s+footer('DASHED RELATION = HYPOTHESIS, NOT ESTABLISHED CAUSAL MEDIATION');
}
function failure(p){
 let s=heading('A final score hides different failure locations','Trace the source, retained information and downstream response');
 for(let r=0;r<2;r++){
  let y=125+r*148;let a=r===0?'State loses a distinction':'Readout ignores a distinction';
  s+=text(30,y-18,a,'small');
  s+=box(30,y,150,75,'Input','same evaluation');
  s+=arrow(185,y+38,244,y+38,p>=1);
  s+=box(248,y,190,75,'Carried state',r===0?'distinction lost':'distinction retained',p>=1?(r===0?'bad':'active'):'');
  s+=arrow(441,y+38,490,y+38,p>=2);
  s+=box(495,y,190,75,'Output','same failure score',p>=2?'bad':'');
  s+=sweep(180,y+38,244,y+38,p,0)+sweep(440,y+38,490,y+38,p,1);
 }
 return s+footer('SAME SCORE ≠ SAME MECHANISM');
}
function inverse(p){
 let s=heading('Forward diagnosis / inverse sufficiency test','The question determines which distinctions the state must preserve');
 s+=text(30,113,'FORWARD','tiny')+text(385,113,'INVERSE','tiny');
 s+=box(30,145,126,74,'Δ state','resolved change',p>=1?'active':'');
 s+=arrow(162,182,206,182,false,true)+box(214,145,113,74,'Δ output','flat');
 s+=text(30,261,'Readout does not use it,','small')+text(30,283,'or its effects cancel.','small');
 s+=path('M382,170 C435,170 425,210 475,210','edge active')+path('M382,275 C435,275 425,210 475,210','edge active');
 s+=text(382,151,'history 1','small')+text(382,304,'history 2','small');
 s+=dot(494,211,18)+text(494,253,'same z','small','middle');
 s+=path('M515,211 C570,211 564,164 615,164','edge'+(p>=2?' active':''))+path('M515,211 C570,211 564,287 615,287','edge'+(p>=2?' active':''));
 s+=text(629,170,'O₁')+text(629,292,'O₂');
 s+=text(385,350,p>=3?'Different outputs → reject z':'Keep downstream conditions fixed','label');
 return s+footer('REJECTION IS RELATIVE TO THE DECLARED OBSERVER AND TESTED CONDITIONS');
}
function history(p){
 let s=heading('Same present. Different past. Different future.','Snapshot sufficiency is a testable hypothesis');
 s+=path('M35,244 C90,66 158,126 222,206 C254,240 290,215 320,210 C397,201 420,110 482,127 L640,127','edge blue');
 s+=path('M35,157 C108,326 179,320 238,246 C263,215 291,208 320,210 C399,210 419,312 482,312 L640,312','edge pink');
 s+=edge(320,89,320,339,false,true)+dot(320,210,9);
 s+=text(320,81,'now','tiny','middle')+text(102,350,'retained history','small')+text(493,350,'continuation','small');
 s+=text(510,108,'readout 1','small')+text(510,338,'readout 2','small');
 if(p>=1)s+=`<rect x="40" y="95" width="279" height="230" rx="15" fill="var(--accent)" opacity=".06"/>`;
 if(p>=2)s+=text(351,190,'same candidate z','small');
 s+=text(30,380,p>=3?'A history-aware candidate must be tested again.':'The candidate merges histories that the future separates.','label');
 return s;
}
function linearity(p){
 let s=heading('Check the measurement before interpreting it','Halve the perturbation; compare slopes and declared boundaries');
 s+=text(35,109,'STEP-HALVING','tiny')+text(404,109,'NESTED READOUTS','tiny');
 s+=edge(52,306,304,306)+edge(52,306,52,134);
 s+=path('M62,293 L291,152','edge active');
 s+=edge(118,260,233,189,true)+dot(118,260)+dot(233,189)+text(94,329,'ε / 2','small')+text(238,329,'ε','small');
 s+=text(60,365,'Matching slopes support local linearity.','small');
 s+=`<circle cx="532" cy="226" r="106" fill="var(--accent-soft)" stroke="var(--accent)"/><circle cx="532" cy="226" r="57" fill="var(--bg-elev)" stroke="var(--accent)"/>`;
 s+=text(532,233,'exact','label','middle')+text(532,312,'works','small','middle');
 s+=dot(532+Math.min(p*27,125),226,8,p>3.8?'bad-fill':'pulse');
 return s+footer('ILLUSTRATIVE BOUNDARIES — VALUES COME FROM THE DECLARED TASK');
}
function response(p){
 let s=heading('Each push is propagated to the endpoint','A later response is the sum of propagated contributions');
 for(let i=0;i<4;i++){
  const x=34+i*167;s+=box(x,178,137,83,`State ${i+1}`,i===3?'read at endpoint':`update A${i+1}`,phaseOn(p,i));
  if(i<3)s+=arrow(x+140,218,x+162,218,p>=i+.5)+sweep(x+140,218,x+162,218,p,i);
  if(i<3)s+=arrow(x+65,124,x+65,166,p>=i)+text(x+65,102,`push B${i+1}δu`,'small','middle');
 }
 s+=text(360,320,'δO = sum of contributions carried to the end','label','middle');
 return s+footer('FINITE-WINDOW RESPONSE • THE REAL PROCESSING ORDER IS PRESERVED');
}
function order(p){
 let s=heading('Use the system’s real processing order','Memory, feedback and delays belong in the declared state');
 ['Clock time','Iteration','Network layer'].forEach((a,r)=>{
  const y=132+r*78;s+=text(30,y+5,a,'small');
  for(let i=0;i<6;i++){const x=218+i*79;s+=dot(x,y,14,i<=p+1?'pulse':'wash');if(i<5)s+=edge(x+17,y,x+60,y,i<p);s+=text(x,y+33,String(i),'tiny','middle');}
 });
 s+=text(30,368,'Arbitrary marker order can create artificial weights.','label');
 return s;
}
function ring(p){
 let s=heading('Delayed coupling needs retained history','Schematic ring and response — not a replay of calibration data');
 for(let i=0;i<64;i++){
  let a=i*Math.PI/32,x=183+113*Math.cos(a),y=236+113*Math.sin(a);let b=a+Math.PI/16;
  s+=edge(x,y,183+113*Math.cos(b),236+113*Math.sin(b));
  s+=dot(x,y,3.5+1.4*Math.sin(a+p*2));
 }
 s+=text(183,231,'64 sites','label','middle')+text(183,255,'delay τ','small','middle');
 s+=edge(366,312,682,312)+edge(366,312,366,133);
 for(let k=0;k<2;k++){
  let d='';for(let x=0;x<=310;x+=3){let y=229-55*Math.sin(x/45-p*1.4-k*.65);d+=(x?' L':'M')+(366+x)+','+y;}
  s+=path(d,k?'edge active':'edge dashed');
 }
 s+=text(378,344,'drive','small')+text(540,344,'delayed mean','label');
 return s+footer('ILLUSTRATIVE AMPLITUDE AND PHASE • REPORTED CALIBRATION IS IN THE TEXT');
}
function flock(p){
 let s=heading('Two readouts need not cross their bounds together','Heading dispersion is illustrative, not a physical simulation');
 for(let i=0;i<64;i++){
  let x=45+(i%8)*33,y=126+Math.floor(i/8)*29,angle=Math.sin(i*9.7)*(.12+p*.4);
  s+=`<g transform="translate(${x} ${y}) rotate(${angle*180/Math.PI})">${path('M-7,-4 L7,0 L-7,4','edge active')}</g>`;
 }
 s+=text(38,385,'increasing angular noise →','small');
 s+=box(365,133,310,86,'Within 30°','share in a fixed wedge',p>=1?'active':'');
 s+=box(365,263,310,86,'Polarization','magnitude of mean heading',p>=2?'active':'');
 s+=text(365,377,'These readouts are not nested.','label');return s;
}
function question(p){
 let s=heading('Declare the observer before the run','The answer key grades the result; it must not enter the evidence');
 const names=[['Question','what is being read?'],['Answers','which outputs are allowed?'],['Boundary','what clears the bar?']];
 names.forEach(([a,b],i)=>{s+=box(30,106+i*90,291,73,a,b,phaseOn(p,i));});
 s+=box(414,146,263,87,'Run → readout','answer or no call',p>=2?'active':'');
 s+=box(414,280,263,76,'Held-out answer key','grades only afterwards',p>=3?'active':'');
 s+=arrow(548,238,548,273,p>=3)+text(395,114,'fixed before running','tiny');return s+footer('NO ANSWER-KEY LEAKAGE INTO THE SOURCE OR CARRIED STATE');
}
function contrast(p){
 let s=heading('Measure the state and output apart','Two push sizes test the stability of the measured response');
 s+=box(30,156,180,104,p>=2?'Push ± ε / 2':'Push ± ε','same source and protocol',phaseOn(p,0));
 for(let i=0;i<2;i++){
  const y=135+i*145;s+=arrow(220,208,292,y+32,p>=1)+box(298,y,364,92,i?'Readout response ΔO':'State response Δz','record separately',phaseOn(p,i+1));
  let w=(i?160:220)*(p>=2?.5:1);s+=`<rect x="317" y="${y+66}" width="${w}" height="6" rx="3" class="teal-fill"/>`;
 }
 return s+footer('SCHEMATIC BAR LENGTHS • A CHANGING STATE DOES NOT ESTABLISH MEDIATION');
}
function inventory(p){
 let s=heading('Preserve request–return associations','Exploratory known case — ordinary diagnosis selects the same join');
 s+=box(30,123,186,74,'Request A','warehouse + ID',p>=1?'active':'')+box(30,256,186,74,'Request B','warehouse + ID',p>=1?'active':'');
 s+=box(282,123,171,74,'Return A','matched by ID',p>=2?'active':'')+box(282,256,171,74,'Return B','matched by ID',p>=2?'active':'');
 s+=arrow(222,160,275,160,p>=2)+arrow(222,293,275,293,p>=2);
 s+=path('M460,160 C501,160 480,226 520,226','edge active')+path('M460,293 C501,293 480,226 520,226','edge active');
 s+=box(529,181,158,94,'After swap',p>=3?'14 → 9':'9 → 9',p>=3?'active':'bad');
 return s+footer('ONE EXPLORATORY CASE • NOT A GENERAL REPAIR OR SUPERIORITY CLAIM');
}
function teddy(p){
 let s=heading('Separate the model, head and decision layer','TEDDY and our head stay frozen during these decision-layer tests');
 const stages=[['RNA','input'],['TEDDY','512-D embedding'],['Our head','134 proteins'],['ANM rule','call / no call']];
 stages.forEach(([a,b],i)=>{const x=30+i*173;s+=box(x,140,147,91,a,b,i===1?'orange':phaseOn(p,i));if(i<3)s+=arrow(x+150,186,x+168,186,p>=i)+sweep(x+150,186,x+168,186,p,i);});
 s+=box(376,294,312,68,'Measured protein','evaluation only',p>=3?'active':'');
 s+=edge(622,236,622,286,p>=3,true);
 return s+footer('PREDICTED PROTEIN ENTERS THE DECISION LAYER • MEASURED PROTEIN GRADES');
}
function teddySteps(p){
 let s=heading('The same protocol on TEDDY’s workflow','Mode B uses the head’s outputs; Mode A inspects internal states');
 const ns=[['Question','declared lineage rule'],['Sources','predicted panel proteins'],['Candidate state','gene-mean / token states'],['Readout','call, graded afterwards']];
 ns.forEach(([a,b],i)=>{let x=30+i%2*350,y=112+Math.floor(i/2)*142;s+=box(x,y,310,106,`${i+1}  ${a}`,b,phaseOn(p,i));});
 return s+footer('DECLARED CANDIDATE ≠ ESTABLISHED SUFFICIENT STATE');
}
function sameCalls(p){
 let s=heading('Rescale the score and the threshold together','The ranking and call stay unchanged — not an accuracy gain');
 const vals=[.83,.5,.24];
 for(let side=0;side<2;side++){
  let x0=side?412:50;s+=text(x0,110,side?'ANM = G × rule':'Declared rule','label');
  for(let i=0;i<3;i++){let x=x0+i*72,h=vals[i]*145*(side?.65:1);s+=`<rect x="${x}" y="${304-h}" width="40" height="${h}" rx="5" class="${side?'teal-fill':'teddy-fill'}"/>`+text(x+20,330,['A','B','C'][i],'small','middle');}
  let yy=304-145*.6*(side?.65:1);s+=edge(x0-12,yy,x0+218,yy,false,true)+text(x0+225,yy+5,'θ','small');
 }
 return s+footer('BAR HEIGHTS ARE SCHEMATIC • THE REPORTED ZERO-MISMATCH RESULT IS IN THE TEXT');
}
function compression(p){
 let s=heading('Near-matched embeddings do not settle where loss occurs','v2 illustration — superseded by the registered v3 E3 result');
 s+=dot(130,205,32,'wash')+dot(150,214,13,'pulse')+dot(174,214,13,'teddy-fill');
 s+=text(120,280,'NK / T near in z','small');
 s+=box(322,133,351,80,'Predicted marker gap','smaller in the v2 prototype',p>=1?'active':'');
 s+=box(322,257,351,80,'Measured marker gap','evaluation reference only',p>=2?'active':'');
 s+=text(30,370,'v3 E3 rejected the proposed repair.','bad-text');return s;
}
function depth(p){
 let s=heading('Reduce the input, then follow the whole workflow','Illustration of RNA thinning — the registered E2 result is separate');
 for(let i=0;i<90;i++){
  const shown=i<90*(1-.85*clamp(p/3));s+=`<circle cx="${40+i%10*22}" cy="${130+Math.floor(i/10)*23}" r="4" fill="var(--accent)" opacity="${shown?1:.12}"/>`;
 }
 s+=text(35,365,'RNA counts retained','small');
 s+=box(332,133,341,80,'Re-embed → predict','same frozen workflow',p>=1?'active':'')+box(332,264,341,80,'Evaluate the call','quality and response checks',p>=2?'active':'');
 return s+footer('ILLUSTRATIVE COUNTS • NO SMOOTH OR LINEAR RESPONSE IS ASSUMED');
}
function layers(p){
 let s=heading('A real layer order makes response tracking possible','The derivative validation gate must pass before interpreting a layer');
 for(let i=0;i<12;i++){
  const x=40+i*54,y=151;s+=`<rect x="${x}" y="${y}" width="28" height="132" rx="6" class="node ${i<=p*3?'active':''}"/>`+text(x+14,310,String(i+1),'tiny','middle');
  if(i<11)s+=edge(x+30,217,x+51,217,i<p*3);
 }
 s+=text(40,118,'input push →','small')+text(487,351,'readout probe','label');
 s+=text(30,385,'v3 E5 stopped at its registered derivative gate.','bad-text');return s;
}

const configs={
 'fw-problem-same-score':[failure,['Trace input','Inspect state','Inspect output','Compare failures'],['The two runs receive the same final score.','One run discards a relevant distinction in its state.','Another retains the distinction but the readout ignores or cancels it.','A final score cannot distinguish these mechanisms.']],
 'fw-push-follow-read':[protocol,['Push','Measure state','Measure output','Compare bounds'],['Control one source while holding the protocol fixed.','Measure how the candidate state responds.','Independently measure the readout response.','Compare the readout with boundaries declared before the run.']],
 'fw-forward-inverse':[inverse,['Declare test','Forward','Inverse','Reject candidate'],['Fix the observer, continuation and window.','A changing state with a flat output need not mean information was lost in the state.','Find histories assigned the same state but distinguished by the readout.','A resolved violation rejects sufficiency under the tested conditions.']],
 'fw-response-sum':[response,['Inject','Propagate','Accumulate','Read endpoint'],['A small push enters at a real processing step.','Later updates carry each contribution forward.','Contributions from different steps accumulate.','Read the endpoint response; the formula in the text defines the calculation.']],
 'fw-linearity-nested':[linearity,['Push ε','Push ε/2','Compare slopes','Check bounds'],['Estimate the response using a declared perturbation size.','Repeat at half the perturbation size.','Agreement supports a local linear approximation; it is not assumed.','The stricter exact region sits inside the works region in this declared nested example.']],
 'fw-step-orders':[order,['Clock','Iteration','Layer','Preserve order'],['Clock time can be the real processing order.','Iterations can also define real steps.','Transformer layers give a real feed-forward order.','An arbitrary ordering of simultaneous evidence can silently impose weights.']],
 'fw-delay-ring':[ring,['Drive','Couple','Retain history','Read response'],['A shared drive perturbs the ring.','Neighbouring sites exchange delayed signals.','The current snapshot alone need not specify the delayed evolution.','Read amplitude and phase; this schematic does not replay the reported calibration.']],
 'fw-flock-readouts':[flock,['Align','Add noise','Measure wedge','Measure mean'],['Particles start near a common heading.','Angular noise broadens the headings in this schematic.','One observer counts headings inside a fixed angular wedge.','Another measures mean-heading magnitude; the two readouts are not nested.']],
 'fw-ai-question':[question,['Question','Answers','Boundary','Grade'],['State which question is being asked.','Declare the allowed outputs before evaluating.','Fix the boundary that separates a call from no call.','Use the held-out answer key only after the readout is complete.']],
 'fw-event-graph':[graph,['Record','Schema','Graph update','Readout'],['Record the event type, target, sign, value, time and provenance.','A schema admits, repairs or rejects the record before it enters the graph.','Apply the declared open-loop retention and propagation rules.','A separate model or program reads the state. This drawing does not establish causal mediation.']],
 'fw-sufficiency-delay':[history,['Two histories','Same snapshot','Continuation','Test revision'],['Two histories differ in the past.','A snapshot candidate maps them to the same present value.','Under fixed continuation conditions, the readouts can diverge.','Reject that candidate, then independently test any history-aware revision.']],
 'fw-measure-contrast':[contrast,['Push','Measure apart','Halve push','Compare'],['Perturb one event value up and down.','Record state and readout changes separately.','Repeat at half the perturbation size; shown bar lengths are illustrative.','Check slope consistency before using a linear response interpretation.']],
 'fw-inventory-join':[inventory,['Requests','Provenance','Join returns','Re-evaluate'],['The known inventory example involves two requests and their warehouse assignments.','Keep request identity and source provenance visible.','Join each request to its own return record.','The known case restored the required 14 → 9 response; ordinary diagnosis chose the same join.']],
 'fw-teddy-pipeline':[teddy,['RNA input','Frozen TEDDY','Our head','Decision'],['RNA is the model input.','Frozen TEDDY produces an embedding; it does not predict proteins.','Our separately trained head predicts proteins from that embedding.','ANM’s decision layer reads predictions; measured protein grades only afterwards.']],
 'fw-teddy-steps':[teddySteps,['Question','Sources','Candidate','Readout'],['Declare the lineage question.','Use the predicted panel proteins as evidence.','Test a declared candidate state; declaration is not proof of sufficiency.','Read a call or no call and grade afterwards.']],
 'fw-teddy-same-calls':[sameCalls,['Rule score','Scale score','Scale threshold','Same call'],['Compute the declared rule directly on the predictions.','The linear graph gives a common score multiplier in this setting.','Scale the threshold by the same multiplier.','Ranking and calls stay identical. The experiment verifies equivalence, not better accuracy.']],
 'fw-nkt-compression':[compression,['Near in z','Predict gap','Reference gap','Read v3 result'],['Near-matched embeddings motivated the v2 diagnostic.','The head’s predicted gap was smaller than the measured gap in the prototype.','Measured proteins are an evaluation reference, not the evidence.','The later registered E3 rejected the repair reading. Keep prototype and confirmed results distinct.']],
 'fw-depth-thinning':[depth,['Thin counts','Re-embed','Follow call','Check response'],['Retain fewer RNA counts in this schematic.','Run the same preprocessing and frozen model again.','Follow the change through the head and the declared decision layer.','E2 tested the added value and linearity separately; a smooth response is not assumed.']],
 'fw-mode-a-layers':[layers,['Input push','Track layers','Validate gate','Interpret limits'],['Use an input perturbation with declared conditions.','Follow the response through the real transformer-layer order.','Validate derivatives against finite differences before interpreting them.','The registered E5 gate failed in float32, so it did not localise a layer.']],
};

function setupFigure(fig){
 const conf=configs[fig.dataset.anim];if(!conf)return null;
 const [draw,steps,notes]=conf, stage=fig.querySelector('.anim-stage');if(!stage)return null;
 fig.classList.add('fw-explainer');stage.setAttribute('role','group');stage.setAttribute('aria-label','Interactive explanatory diagram');
 stage.innerHTML=`<div class="fw-viz-header"><span>ANM / visual protocol</span><span>Schematic · not measured data</span></div><div class="fw-steps">${steps.map((s,i)=>`<button type="button" class="fw-step" data-phase="${i}">${i+1}  ${esc(s)}</button>`).join('')}</div><div class="fw-picture"></div><p class="fw-viz-note"></p><div class="fw-controls"><button type="button" class="fw-play" aria-label="Play animation">Play</button><button type="button" class="fw-replay" aria-label="Replay animation">↺</button><input type="range" min="0" max="3999" step="1" value="0" aria-label="Animation progress"/><output>01 / 04</output></div>`;
 const mq=window.matchMedia('(prefers-reduced-motion: reduce)');
 const state={fig,draw,steps,notes,p:0,playing:!mq.matches,visible:false,ended:false};
 const picture=stage.querySelector('.fw-picture'), note=stage.querySelector('.fw-viz-note'),slider=stage.querySelector('input'),play=stage.querySelector('.fw-play'),out=stage.querySelector('output'),buttons=[...stage.querySelectorAll('.fw-step')];
 function render(){
  const i=Math.min(3,Math.floor(state.p));picture.innerHTML=frame(draw(state.p));
  picture.querySelector('svg').setAttribute('aria-label',notes[i]);
  note.textContent=notes[i];slider.value=Math.round(state.p*1000);slider.setAttribute('aria-valuetext',`Step ${i+1}: ${steps[i]}`);out.textContent=`0${i+1} / 04`;
  buttons.forEach((b,j)=>{if(i===j)b.setAttribute('aria-current','step');else b.removeAttribute('aria-current');});
  play.textContent=state.playing?'Pause':state.ended?'Replay':'Play';play.setAttribute('aria-label',state.playing?'Pause animation':state.ended?'Replay animation':'Play animation');fig.dataset.state=state.playing?'playing':'paused';
 }
 state.render=render;
 play.addEventListener('click',()=>{if(state.ended){state.p=0;state.ended=false;}state.playing=!state.playing;render();kick();});
 stage.querySelector('.fw-replay').addEventListener('click',()=>{state.p=0;state.ended=false;state.playing=true;render();kick();});
 slider.addEventListener('input',()=>{state.p=Number(slider.value)/1000;state.playing=false;state.ended=false;render();});
 buttons.forEach((b,i)=>b.addEventListener('click',()=>{state.p=i;state.playing=false;state.ended=false;render();}));
 mq.addEventListener('change',()=>{if(mq.matches){state.playing=false;render();}});
 render();return state;
}
const records=[...document.querySelectorAll('figure.anim[data-anim]')].map(setupFigure).filter(Boolean);
let raf=0,last=0;
function kick(){if(!raf&&!document.hidden&&records.some(r=>r.visible&&r.playing))raf=requestAnimationFrame(tick);}
function tick(now){raf=0;const dt=last?Math.min((now-last)/1000,.08):0;last=now;
 for(const r of records){if(!r.visible||!r.playing)continue;r.p=Math.min(3.999,r.p+dt/4);if(r.p>=3.999){r.playing=false;r.ended=true;}r.render();}kick();}
if('IntersectionObserver' in window){const io=new IntersectionObserver(es=>{for(const e of es){let r=records.find(r=>r.fig===e.target);if(r)r.visible=e.isIntersecting;}last=0;kick();},{threshold:.12});records.forEach(r=>io.observe(r.fig));}
else{records.forEach(r=>r.visible=true);kick();}
document.addEventListener('visibilitychange',()=>{last=0;if(document.hidden){cancelAnimationFrame(raf);raf=0;}else kick();});
