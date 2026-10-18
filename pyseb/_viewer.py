"""HTML presentation for the optional Plotly structure viewer (lazy imports)."""
from functools import lru_cache
import json


_STYLE = r"""
.pyseb-viewer {position:relative; isolation:isolate;}
.pyseb-viewer:hover .modebar--hover .modebar-group {opacity:1!important;}
.pyseb-viewer .pyseb-settings {position:absolute;top:8px;left:8px;z-index:1001;
 font:13px/1.4 system-ui,sans-serif;color:#e9edf3;max-width:calc(100% - 16px);}
.pyseb-settings button,.pyseb-settings select {font:inherit;}
.pyseb-settings .pyseb-toolbar {display:flex;align-items:center;gap:6px;}
.pyseb-settings button {display:inline-flex;align-items:center;justify-content:center;
 box-sizing:border-box;line-height:1.2;vertical-align:middle;}
.pyseb-settings .pyseb-toggle svg {display:block;width:18px;height:18px;flex:none;}
.pyseb-settings .pyseb-toggle {width:32px;height:30px;padding:0;flex:none;border:1px solid #687386;
 border-radius:5px;background:#232b36ed;color:#f3f6fa;cursor:pointer;}
.pyseb-settings .pyseb-fold {height:30px;}
.pyseb-settings .pyseb-toggle,.pyseb-settings .pyseb-fold {opacity:0;
 transition:opacity .15s;pointer-events:none;}
.pyseb-viewer:hover .pyseb-toggle,.pyseb-viewer:hover .pyseb-fold,
.pyseb-viewer:focus-within .pyseb-toggle,.pyseb-viewer:focus-within .pyseb-fold,
.pyseb-settings[data-open] .pyseb-toggle,.pyseb-settings[data-open] .pyseb-fold,
.pyseb-touch .pyseb-toggle,.pyseb-touch .pyseb-fold {opacity:1;pointer-events:auto;}

.pyseb-settings .pyseb-panel {margin-top:5px;width:272px;max-width:100%;
 max-height:min(420px,65vh);overflow:auto;box-sizing:border-box;padding:12px;
 border:1px solid #687386;border-radius:7px;background:#202833f5;
 box-shadow:0 5px 20px #0005;}
.pyseb-settings [hidden] {display:none!important;}
.pyseb-settings .pyseb-heading {display:flex;align-items:center;justify-content:space-between;
 font-size:14px;font-weight:600;margin-bottom:10px;}
.pyseb-settings button {cursor:pointer;color:inherit;background:#354152;
 border:1px solid #718098;border-radius:4px;padding:3px 8px;}
.pyseb-settings label {display:flex;align-items:center;gap:7px;margin:6px 0;}
.pyseb-settings select {width:100%;min-width:0;background:#303b4b;color:inherit;
 border:1px solid #718098;border-radius:4px;padding:4px;}
.pyseb-settings input[type=range] {min-width:0;flex:1;accent-color:#6cb6ff;}
.pyseb-settings input[type=color] {width:32px;height:24px;padding:1px;}
.pyseb-settings input[type=checkbox] {accent-color:#6cb6ff;}
.pyseb-settings .pyseb-name {overflow:hidden;text-overflow:ellipsis;white-space:nowrap;}
.pyseb-settings .pyseb-component {border-top:1px solid #536074;padding-top:4px;margin-top:5px;}
.pyseb-settings .pyseb-component label {margin:3px 0;}
.pyseb-settings .pyseb-tree-heading {display:flex;align-items:center;gap:5px;}
.pyseb-settings .pyseb-tree-heading label {min-width:0;flex:1;}
.pyseb-settings .pyseb-expand {width:22px;flex:none;padding:0;border:0;background:transparent;}
.pyseb-settings .pyseb-tree-children {margin-left:8px;padding-left:8px;border-left:1px solid #536074;}
.pyseb-settings .pyseb-mixed {font-size:11px;color:#c7d4e7;}
.pyseb-settings .pyseb-actions {display:flex;gap:8px;margin-top:12px;}
@media (hover:none) {.pyseb-settings .pyseb-toggle,.pyseb-settings .pyseb-fold {opacity:1;pointer-events:auto;}}
"""

# The post_script hook runs after Plotly.newPlot, including on exported HTML.
# Labels from model metadata are always assigned with textContent, never HTML.
_SCRIPT = r"""
(function() {
 const gd = document.getElementById('{plot_id}');
 if (!gd || gd.querySelector('.pyseb-settings')) return;
 const Plot = window.Plotly;
 const variants = __PYSEB_VARIANTS__;
 if (!document.getElementById('pyseb-viewer-style')) {
   const style = document.createElement('style');
   style.id = 'pyseb-viewer-style'; style.textContent = __PYSEB_STYLE__;
   document.head.appendChild(style);
 }
 gd.classList.add('pyseb-viewer');
 const root = document.createElement('div'); root.className = 'pyseb-settings';
 const toggle = document.createElement('button'); toggle.className = 'pyseb-toggle';
 toggle.type = 'button'; toggle.title = 'Structure settings';
 const svgNS='http://www.w3.org/2000/svg';
 const gear=document.createElementNS(svgNS,'svg');gear.setAttribute('viewBox','0 0 24 24');
 gear.setAttribute('aria-hidden','true');gear.setAttribute('focusable','false');
 const outline=document.createElementNS(svgNS,'path');
 const points=[];
 for(let tooth=0;tooth<8;tooth++)for(let step=0;step<4;step++){
   const angle=-Math.PI/2+(tooth+step/4)*Math.PI/4;
   const radius=step===1||step===2?10:7.5;
   points.push((12+radius*Math.cos(angle)).toFixed(3)+','+(12+radius*Math.sin(angle)).toFixed(3));
 }
 outline.setAttribute('d','M'+points.join(' L')+' Z M12 8.8 a3.2 3.2 0 1 0 0 6.4 a3.2 3.2 0 1 0 0 -6.4 Z');
 outline.setAttribute('fill','currentColor');outline.setAttribute('fill-rule','evenodd');
 gear.append(outline);toggle.append(gear);
 toggle.setAttribute('aria-label', 'Structure settings'); toggle.setAttribute('aria-expanded','false');
 const panel = document.createElement('section'); panel.className = 'pyseb-panel';
 panel.setAttribute('aria-label','Structure settings'); panel.hidden = true;
 const heading = document.createElement('div'); heading.className = 'pyseb-heading';
 const title = document.createElement('span'); title.textContent = 'Structure settings';
 const close = document.createElement('button'); close.type='button';close.textContent='\u00d7';
 close.setAttribute('aria-label','Close structure settings');
 const toolbar=document.createElement('div');toolbar.className='pyseb-toolbar';
 heading.append(title,close);panel.append(heading);toolbar.append(toggle);root.append(toolbar,panel);gd.append(root);
 const sizePanel=()=>{panel.style.maxHeight=Math.min(420,Math.max(80,gd.clientHeight-56))+'px';};
 sizePanel();if(window.ResizeObserver)new ResizeObserver(sizePanel).observe(gd);
 function open(value) {
   panel.hidden=!value; toggle.setAttribute('aria-expanded',String(value));
   root.toggleAttribute('data-open',value);
 }
 toggle.addEventListener('click',()=>open(panel.hidden));
 close.addEventListener('click',()=>{open(false);toggle.focus();});
 root.addEventListener('keydown',e=>{if(e.key==='Escape'){open(false);toggle.focus();}});
 // Events within the overlay must not reach the orbit/zoom handlers.
 for (const event of ['pointerdown','mousedown','touchstart','wheel','dblclick'])
   root.addEventListener(event,e=>e.stopPropagation());
 gd.addEventListener('pointerdown',e=>{
   if(e.pointerType==='touch' && !root.contains(e.target)) gd.classList.add('pyseb-touch');
 });
 function meta(trace) {return trace.meta && trace.meta.pyseb || {};}
 function isReference(trace) {return meta(trace).role==='reference';}
 function groupName(trace,index) {
   const instance=meta(trace).instance;
   return instance ? String(instance).replace(/\//g,':') : trace.name || 'Component '+(index+1);
 }
 function colorOf(trace) {return trace.type==='mesh3d' ? trace.color :
   (trace.mode && trace.mode.includes('markers') ? (trace.marker||{}).color : (trace.line||{}).color);}
 function hex(color) {
   if(typeof color!=='string') return '#8ca6cc';
   if(/^#[0-9a-f]{6}$/i.test(color)) return color;
   const rgb=color.match(/^rgb\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\)$/);
   return rgb ? '#'+rgb.slice(1).map(v=>Number(v).toString(16).padStart(2,'0')).join('') : '#8ca6cc';
 }
 const referenceLabel=document.createElement('label');
 const references=document.createElement('input');references.type='checkbox';
 referenceLabel.append(references,document.createTextNode('Reference points')); panel.append(referenceLabel);
 const scatteringLabel=document.createElement('label');
 const scatteringOpacity=document.createElement('input');scatteringOpacity.type='checkbox';
 scatteringLabel.append(scatteringOpacity,document.createTextNode('Opacity by scattering'));
 scatteringLabel.title='Scale manual opacity by normalized absolute scattering length or contrast density';
 const scatteringMode=document.createElement('select');scatteringMode.setAttribute('aria-label','Opacity weighting');
 for(const [value,text] of [['beta','Total |β|'],['density','Contrast density |β| / V']]){
   const option=document.createElement('option');option.value=value;option.textContent=text;scatteringMode.append(option);
 }
 scatteringMode.hidden=true;
 const scatteringNote=document.createElement('p');scatteringNote.className='pyseb-mixed';
 scatteringNote.setAttribute('role','status');
 panel.append(scatteringLabel,scatteringMode,scatteringNote);
 function scatteringGroups(){
   const groups=new Map();gd.data.forEach((trace,i)=>{
     if(!isReference(trace)&&!groups.has(groupName(trace,i)))groups.set(groupName(trace,i),meta(trace));
   });return groups;
 }
 function scatteringValue(info){
   if(!Number.isFinite(info.beta))return null;
   if(scatteringMode.value==='beta')return Math.abs(info.beta);
   if(info.beta===0)return 0;
   if(!Number.isFinite(info.volume)||info.volume<=0)return null;
   const value=Math.abs(info.beta)/info.volume;return Number.isFinite(value)?value:null;
 }
 function opacityFactors(){
   const groups=scatteringGroups();const factors=new Map();let maximum=0;let missing=0;
   for(const info of groups.values()){
     const value=scatteringValue(info);if(value===null)missing++;else maximum=Math.max(maximum,value);
   }
   for(const [name,info] of groups){
     const value=scatteringValue(info);
     factors.set(name,!scatteringOpacity.checked||value===null?1:(maximum>0?value/maximum:0));
   }
   scatteringOpacity.disabled=![...groups.values()].some(info=>Number.isFinite(info.beta));
   scatteringMode.options[1].disabled=![...groups.values()].some(info=>Number.isFinite(info.beta)&&Number.isFinite(info.volume)&&info.volume>0);
   scatteringMode.options[1].textContent=scatteringMode.options[1].disabled?
     'Contrast density |β| / V (volume unavailable)':'Contrast density |β| / V';
   scatteringMode.hidden=!scatteringOpacity.checked;
   if(scatteringOpacity.disabled){scatteringNote.textContent='No beta values available for this view.';}
   else if(scatteringOpacity.checked){
     scatteringNote.textContent='Normalized within this structure; opacity sliders are multipliers.'+
       (missing?' '+missing+' subunit'+(missing===1?' has':'s have')+' no '+
         (scatteringMode.value==='density'?'beta/material volume':'beta')+'; manual opacity retained.':'');
   }else{scatteringNote.textContent='';}
   return factors;
 }
 scatteringOpacity.addEventListener('change',()=>{opacityFactors();syncControls();applyStyles();});
 scatteringMode.addEventListener('change',()=>{opacityFactors();syncControls();applyStyles();});
 const components=document.createElement('div');panel.append(components);
 const initial = gd.data.map(trace=>({visible:trace.visible===undefined?true:trace.visible,
   opacity:trace.opacity===undefined?1:trace.opacity, color:colorOf(trace)}));
 const styles = new Map();
 const alternatingLabel=document.createElement('label');
 const alternating=document.createElement('input');alternating.type='checkbox';
 alternatingLabel.append(alternating,document.createTextNode('Alternate connected colors'));
 alternatingLabel.title='Use two colors along connected chains and branches';
 panel.insertBefore(alternatingLabel,referenceLabel);
 const palette=document.createElement('div');palette.hidden=true;
 const defaultAlternatingColors=['#269bd2','#f28e2b'];
 const alternatingColors=[...defaultAlternatingColors];
 const partitions=new Map();const adjacency=new Map();
 gd.data.forEach((trace,i)=>{if(!isReference(trace)){
   const name=groupName(trace,i);
   adjacency.set(name,(meta(trace).neighbors||[]).map(value=>String(value).replace(/\//g,':')));
 }});
 for(const start of adjacency.keys()){
   if(partitions.has(start))continue;
   partitions.set(start,0);const queue=[start];
   for(let index=0;index<queue.length;index++){
     const name=queue[index];
     for(const neighbor of adjacency.get(name)||[]){
       if(!adjacency.has(neighbor)||partitions.has(neighbor))continue;
       partitions.set(neighbor,1-partitions.get(name));queue.push(neighbor);
     }
   }
 }
 alternating.disabled=![...adjacency.values()].some(values=>values.length);
 alternating.checked=!alternating.disabled;
 palette.hidden=!alternating.checked;
 for(let index=0;index<2;index++){
   const label=document.createElement('label');label.textContent='Alternating color '+(index+1);
   const input=document.createElement('input');input.type='color';input.value=alternatingColors[index];
   input.setAttribute('aria-label',label.textContent);
   input.addEventListener('input',()=>{alternatingColors[index]=input.value;applyStyles();});
   label.append(input);palette.append(label);
 }
 panel.insertBefore(palette,referenceLabel);
 alternating.addEventListener('change',()=>{palette.hidden=!alternating.checked;applyStyles();});
 const nodes = new Map();
 const expanded = new Set(['']);
 function updateReferences() {
   const indices=[];const visibility=[];
   gd.data.forEach((trace,i)=>{if(isReference(trace)){
     indices.push(i); const style=styles.get(groupName(trace,i));
     visibility.push(references.checked && (!style || style.visible));
   }});
   if(indices.length) return Plot.restyle(gd,{visible:visibility},indices);
 }
 references.addEventListener('change',updateReferences);
 function applyStyles(names=styles.keys()) {
   const selected=new Set(names);const indices=[];const visible=[];const opacity=[];
   const factors=opacityFactors();
   const colors={color:{indices:[],values:[]},'line.color':{indices:[],values:[]},
     'marker.color':{indices:[],values:[]}};
   gd.data.forEach((trace,i)=>{
     const name=groupName(trace,i);if(isReference(trace)||!selected.has(name))return;
     const state=styles.get(name);if(!state)return;
     indices.push(i);visible.push(state.visible);opacity.push(state.opacity*(factors.get(name)??1));
     const attr=trace.type==='mesh3d'?'color':
       (trace.mode&&trace.mode.includes('markers')?'marker.color':'line.color');
     colors[attr].indices.push(i);colors[attr].values.push(alternating.checked?
       alternatingColors[partitions.get(name)??0]:state.color);
   });
   const work=[];
   if(indices.length)work.push(Plot.restyle(gd,{visible,opacity},indices));
   for(const [attr,batch] of Object.entries(colors))
     if(batch.indices.length)work.push(Plot.restyle(gd,{[attr]:batch.values},batch.indices));
   work.push(updateReferences());return Promise.all(work);
 }
 function syncControls() {
   for(const node of nodes.values()){
     const states=[...node.leaves].map(name=>styles.get(name));
     if(!states.length||!node.controls)continue;
     const controls=node.controls;const first=states[0];
     const allVisible=states.every(state=>state.visible);
     const anyVisible=states.some(state=>state.visible);
     const mixedOpacity=states.some(state=>state.opacity!==first.opacity);
     const mixedColor=states.some(state=>state.color!==first.color);
     controls.visible.checked=allVisible;
     controls.visible.indeterminate=anyVisible&&!allVisible;
     controls.visible.setAttribute('aria-checked',anyVisible&&!allVisible?'mixed':String(allVisible));
     controls.color.value=first.color;
     controls.color.title=mixedColor?'Mixed colors; choose one for this whole group':'Color';
     controls.mixed.hidden=!mixedColor;
     controls.opacity.value=String(first.opacity);
     controls.opacity.title=scatteringOpacity.checked?'Manual opacity multiplier':'Opacity';
     controls.opacity.setAttribute('aria-valuetext',mixedOpacity?'Mixed opacity':Math.round(first.opacity*100)+'%');
     controls.value.textContent=mixedOpacity?'Mixed':Math.round(first.opacity*100)+'%';
   }
 }
 function editNode(node,attribute,value) {
   for(const name of node.leaves)styles.get(name)[attribute]=value;
   const manualColor=attribute==='color'&&alternating.checked;
   if(manualColor){alternating.checked=false;palette.hidden=true;}
   syncControls();return manualColor?applyStyles():applyStyles(node.leaves);
 }
 function buildComponents() {
   components.replaceChildren();nodes.clear();
   const tree={path:'',name:'All subunits',children:new Map(),leaves:new Set()};
   nodes.set('',tree);
   gd.data.forEach((trace,i)=>{if(!isReference(trace)){
     const name=groupName(trace,i);
     if(!styles.has(name))styles.set(name,{visible:trace.visible!==false&&trace.visible!=='legendonly',
       opacity:trace.opacity===undefined?1:trace.opacity,color:hex(colorOf(trace))});
     let node=tree;node.leaves.add(name);let path='';
     for(const part of name.split(':')){
       path=path?path+':'+part:part;
       if(!node.children.has(part)){
         const child={path,name:part,children:new Map(),leaves:new Set()};
         node.children.set(part,child);nodes.set(path,child);
       }
       node=node.children.get(part);node.leaves.add(name);
     }
   }});
   references.disabled=!gd.data.some(isReference);
   function renderNode(node,parent){
     const name=node.path||node.name;
     const row=document.createElement('div');row.className='pyseb-component';
     row.dataset.path=node.path;
     const heading=document.createElement('div');heading.className='pyseb-tree-heading';
     const label=document.createElement('label');const visible=document.createElement('input');
     visible.type='checkbox';visible.setAttribute('aria-label','Show '+name);
     const text=document.createElement('span');text.className='pyseb-name';text.textContent=node.name;
     text.title=name+' ('+node.leaves.size+' subunit'+(node.leaves.size===1?'':'s')+')';
     label.append(visible,text);
     const children=document.createElement('div');children.className='pyseb-tree-children';
     if(node.children.size){
       const expand=document.createElement('button');expand.type='button';expand.className='pyseb-expand';
       function updateExpansion(){
         const isOpen=expanded.has(node.path);children.hidden=!isOpen;
         expand.textContent=isOpen?'\u25be':'\u25b8';
         expand.setAttribute('aria-expanded',String(isOpen));
         expand.setAttribute('aria-label',(isOpen?'Collapse ':'Expand ')+name);
       }
       expand.addEventListener('click',()=>{
         if(expanded.has(node.path))expanded.delete(node.path);else expanded.add(node.path);
         updateExpansion();
       });
       updateExpansion();heading.append(expand);
     }
     heading.append(label);row.append(heading);
     const opacityLabel=document.createElement('label');
     const color=document.createElement('input');color.type='color';
     color.setAttribute('aria-label','Color '+name);
     const opacity=document.createElement('input');opacity.type='range';opacity.min='0';opacity.max='1';
     opacity.step='0.05';opacity.setAttribute('aria-label','Opacity '+name);
     const value=document.createElement('output');
     opacityLabel.append(color,opacity,value);row.append(opacityLabel);
     const mixed=document.createElement('span');mixed.className='pyseb-mixed';mixed.textContent='Mixed colors';
     row.append(mixed);parent.append(row);
     node.controls={visible,color,opacity,value,mixed};
     visible.addEventListener('change',()=>editNode(node,'visible',visible.checked));
     color.addEventListener('input',()=>editNode(node,'color',color.value));
     opacity.addEventListener('input',()=>editNode(node,'opacity',Number(opacity.value)));
     if(node.children.size){
       parent.append(children);
       for(const child of node.children.values())renderNode(child,children);
     }
   }
   if(tree.leaves.size)renderNode(tree,components);
   opacityFactors();syncControls();
 }
 let sceneSelect,layoutButton;
 let activeVariant=0;
 function alternateLayout(index) {
   const current=variants[index];if(!current||!current.layout)return -1;
   const other=current.layout==='random'?'readable':'random';
   return variants.findIndex(variant=>variant.seed===current.seed&&variant.layout===other);
 }
 function updateLayoutButton() {
   if(!layoutButton)return;
   layoutButton.hidden=alternateLayout(activeVariant)<0;
   const readable=variants[activeVariant].layout==='readable';
   layoutButton.textContent=readable?'See model':'See unfolded';
   layoutButton.setAttribute('aria-pressed',String(readable));
   layoutButton.title=readable?
     'Restore the same seeded, freely rotating model realization':
     'Reduce overlap for readability; this is not an excluded-volume model';
 }
 async function switchVariant(index) {
   sceneSelect.disabled=true;if(layoutButton)layoutButton.disabled=true;
   try{
     const variant=variants[index];
     const camera=gd.layout.scene&&gd.layout.scene.camera;
     await Plot.react(gd,JSON.parse(JSON.stringify(variant.data)),gd.layout,gd._context);
     await applyStyles();
     await updateReferences();if(camera)await Plot.relayout(gd,{'scene.camera':camera});
     activeVariant=index;sceneSelect.value=String(index);
   }catch(error){sceneSelect.value=String(activeVariant);throw error;}
   finally{
     sceneSelect.disabled=false;if(layoutButton)layoutButton.disabled=false;
     updateLayoutButton();
   }
 }
 if(variants.length>1){
   const label=document.createElement('label');label.textContent='Realization';panel.insertBefore(label,referenceLabel);
   sceneSelect=document.createElement('select');sceneSelect.setAttribute('aria-label','Realization');
   variants.forEach((variant,i)=>{const option=document.createElement('option');option.value=String(i);
     option.textContent=variant.label;sceneSelect.append(option);});panel.insertBefore(sceneSelect,referenceLabel);
   sceneSelect.addEventListener('change',()=>switchVariant(Number(sceneSelect.value)));
   if(variants.some((variant,i)=>alternateLayout(i)>=0)){
     layoutButton=document.createElement('button');layoutButton.type='button';layoutButton.className='pyseb-fold';
     toolbar.append(layoutButton);updateLayoutButton();
     layoutButton.addEventListener('click',()=>{
       const index=alternateLayout(activeVariant);if(index>=0)return switchVariant(index);
     });
   }
 }
 const actions=document.createElement('div');actions.className='pyseb-actions';
 const resetCamera=document.createElement('button');resetCamera.type='button';resetCamera.textContent='Reset camera';
 const camera=JSON.parse(JSON.stringify((gd.layout.scene||{}).camera||{}));
 resetCamera.addEventListener('click',()=>Plot.relayout(gd,{'scene.camera':camera}));
 const resetStyle=document.createElement('button');resetStyle.type='button';resetStyle.textContent='Reset style';
 resetStyle.addEventListener('click',async()=>{
   if(sceneSelect){activeVariant=0;sceneSelect.value='0';updateLayoutButton();await Plot.react(gd,JSON.parse(JSON.stringify(variants[0].data)),gd.layout,gd._context);}
   styles.clear();references.checked=false;scatteringOpacity.checked=false;scatteringMode.value='beta';
   alternatingColors.splice(0,alternatingColors.length,...defaultAlternatingColors);
   palette.querySelectorAll('input').forEach((input,i)=>{input.value=alternatingColors[i];});
   alternating.checked=!alternating.disabled;palette.hidden=!alternating.checked;
   for(let i=0;i<gd.data.length;i++){
     const state=initial[i];if(!state)continue;
     const trace=gd.data[i];const attr=trace.type==='mesh3d'?'color':
       (trace.mode&&trace.mode.includes('markers')?'marker.color':'line.color');
     const update={visible:state.visible,opacity:state.opacity};if(state.color)update[attr]=state.color;
     await Plot.restyle(gd,update,[i]);
   }
   buildComponents();await applyStyles();
 });
 actions.append(resetCamera,resetStyle);panel.append(actions);buildComponents();applyStyles();
})();
"""


def _script(variants):
    from plotly.utils import PlotlyJSONEncoder

    def safe_json(value):
        return json.dumps(value, cls=PlotlyJSONEncoder).replace('<', '\\u003c').replace('>', '\\u003e').replace('&', '\\u0026')

    return _SCRIPT.replace('__PYSEB_STYLE__', safe_json(_STYLE)).replace(
        '__PYSEB_VARIANTS__', safe_json(variants)
    )


@lru_cache(maxsize=1)
def figure_class():
    """Construct the Figure subclass without importing Plotly at package import."""
    import plotly.graph_objects as go

    class StructureFigure(go.Figure):
        """A Plotly figure with an optional HTML settings overlay."""

        def _html_options(self, kwargs):
            options = dict(kwargs)
            script = options.pop('post_script', None)
            scripts = [] if script is None else ([script] if isinstance(script, str) else list(script))
            variants = getattr(self, '_pyseb_variants', [])
            options['post_script'] = [_script(variants), *scripts]
            return options

        def to_html(self, *args, **kwargs):
            return super().to_html(*args, **self._html_options(kwargs))

        def write_html(self, *args, **kwargs):
            return super().write_html(*args, **self._html_options(kwargs))

        def _repr_html_(self):
            return self.to_html(full_html=False)

        def _repr_mimebundle_(self, include=None, exclude=None, **kwargs):
            if (include is not None and 'text/html' not in include) or (exclude and 'text/html' in exclude):
                return super()._repr_mimebundle_(include=include, exclude=exclude, **kwargs)
            return {'text/html': self._notebook_html()}

        def _notebook_html(self, **kwargs):
            try:
                from IPython.display import display
            except ImportError:
                return self.to_html(full_html=False, **kwargs)

            options = dict(kwargs)
            options['full_html'] = False
            if options.get('include_plotlyjs', True) is True:
                from plotly.offline import get_plotlyjs

                # Quarto extracts HTML containing the Plotly bundle as a page
                # dependency. Emit only the library in that output, otherwise
                # the complete viewer is moved above the article title.
                display({'text/html': '<script type="text/javascript">' +
                         get_plotlyjs() + '</script>'}, raw=True)
                options['include_plotlyjs'] = False
            return self.to_html(**options)

        def show(self, renderer=None, **kwargs):
            import plotly.io as pio

            selected = renderer or pio.renderers.default
            html_renderers = {'browser', 'notebook', 'notebook_connected', 'plotly_mimetype',
                              'jupyterlab', 'vscode', 'colab', 'iframe', 'iframe_connected'}
            if selected and any(name not in html_renderers for name in selected.split('+')):
                return super().show(renderer=renderer, **kwargs)
            try:
                from IPython import get_ipython
                notebook = getattr(get_ipython(), 'kernel', None) is not None
            except ImportError:
                notebook = False
            if notebook and selected != 'browser':
                from IPython.display import display
                display({'text/html': self._notebook_html(**kwargs)}, raw=True)
            else:
                from pathlib import Path
                import tempfile

                # Browser launchers such as xdg-open return before the browser
                # reads the page. Keep this file after Python exits, both for
                # the initial load and for subsequent browser refreshes.
                directory = Path(tempfile.mkdtemp(prefix='pyseb-viewer-'))
                self.write_html(directory / 'structure.html', auto_open=True, **kwargs)

    return StructureFigure
