// Nova 2026-11-23. Run: node flat_counterexample.js.
// Double precision exploration, not a certified root count or quadrature.
const A=(a,b)=>[a[0]+b[0],a[1]+b[1]],S=(a,b)=>[a[0]-b[0],a[1]-b[1]],M=(a,b)=>[a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]],C=(a,t)=>[a[0]*t,a[1]*t],abs=a=>Math.hypot(...a),D=(a,b)=>C(M(a,[b[0],-b[1]]),1/(b[0]*b[0]+b[1]*b[1])),E=a=>C([Math.cos(a[1]),Math.sin(a[1])],Math.exp(a[0])),sq=a=>M(a,a),W=z=>E(C(sq(z),-.5)),H=z=>A([1.5,.5],C(M(S([1,0],sq(z)),W(z)),3));
// b=-3,p_b=0,a=-2,p0=1,g=3
const G=(z,w)=>[S(A(A(M([1.5,.5],z),[1.5,.5]),C(M(z,W(z)),3)),w),S(A(A(M([1.5,.5],w),[1.5,1.5]),C(M(w,W(w)),3)),z)];
const P=(z,w)=>S(A(A(A(C(sq(A(z,[2,0])),.5),M([-.25,.25],sq(A(z,[1,0])))),C(sq(S(w,z)),.5)),A(C(sq(A(w,[3,0])),.5),M([-.25,.25],sq(A(w,[3,0]))))),C(A(W(z),W(w)),3));
let roots=[];
for(let x=-3;x<=3;x+=.5)for(let y=-3;y<=3;y+=.5)for(let u=-3;u<=3;u+=.5)for(let v=-3;v<=3;v+=.5){
 let z=[x,y],w=[u,v];
 for(let i=0;i<50;i++){let [f,g]=G(z,w),h=H(z),k=H(w),det=S(M(h,k),[1,0]);if(abs(f)+abs(g)<1e-9)break;let dz=D(A(M(k,f),g),det),dw=D(A(f,M(h,g)),det);z=S(z,dz);w=S(w,dw);if(!Number.isFinite(abs(z)+abs(w))||abs(z)>30||abs(w)>30)break;}
 let [f,g]=G(z,w);if(abs(f)+abs(g)<1e-8&&Math.max(...z.map(Math.abs),...w.map(Math.abs))<=4&&!roots.some(([a,b])=>abs(S(z,a))+abs(S(w,b))<1e-6))roots.push([z,w]);
}
const sqrt=z=>C([Math.cos(Math.atan2(z[1],z[0])/2),Math.sin(Math.atan2(z[1],z[0])/2)],Math.sqrt(abs(z)));
const results=[];
for(const h of [.1,.05]){
 let I=[0,0],N=900,d=20*Math.sqrt(h)/N;
 for(let i=0;i<=N;i++){let x=-1-10*Math.sqrt(h)+i*d;for(let j=0;j<=N;j++){let y=-3-10*Math.sqrt(h)+j*d,ph=P([x,0],[y,0]);I=A(I,C(E([-ph[1]/h,ph[0]/h]),d*d*((i===0||i===N)?.5:1)*((j===0||j===N)?.5:1)));}}
 let rows=roots.map(([z,w])=>{let ph=P(z,w),a=M([0,-1],H(z)),d=M([0,-1],H(w)),disc=sqrt(S(sq(S(a,d)),[4,0])),ev1=C(A(A(a,d),disc),.5),ev2=C(S(A(a,d),disc),.5),L=D(C(E([-ph[1]/h,ph[0]/h]),2*Math.PI*h),M(sqrt(ev1),sqrt(ev2)));return {z,w,im:ph[1],err:Math.min(abs(S(I,L)),abs(A(I,L)))/abs(I)}}).filter(r=>r.im>=0);
 let best=rows.slice().sort((a,b)=>a.err-b.err)[0],base=best.im,loc=[0,0],min=Infinity;
 for(let x=-8;x<=8;x+=.02)for(let y=-8;y<=8;y+=.02){let z=P([x,best.z[1]],[y,best.w[1]])[1];if(z<min){min=z;loc=[x,y]}}
 results.push({h,I,abs:abs(I),best,lowest:rows.sort((a,b)=>a.im-b.im).slice(0,2),flatMin:min,loc,dip:base-min});
}
console.log(JSON.stringify({roots:roots.length,results},null,2));
