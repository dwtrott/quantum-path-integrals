// Nova, 2026-11-23. Run with node contour_review.js.
// IEEE-754 double precision, evidence only; no interval or quadrature certification.
const c=.553223, theta=Math.PI/8;
const add=(a,b)=>[a[0]+b[0],a[1]+b[1]];
const mul=(a,b)=>[a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]];
const sc=(a,t)=>a.map(x=>x*t);
const div=(a,b)=>sc(mul(a,[b[0],-b[1]]),1/(b[0]**2+b[1]**2));
const exp=a=>[Math.exp(a[0])*Math.cos(a[1]),Math.exp(a[0])*Math.sin(a[1])];
const dp=z=>add(add(z,[-c,0]),sc(mul(z,exp(sc(mul(z,z),-.5))),-2));
const dd=z=>add([1,0],sc(mul(add(mul(z,z),[-1,0]),exp(sc(mul(z,z),-.5))),2));
const phiIm=([x,y])=>(x-c)*y-2*Math.exp((y*y-x*x)/2)*Math.sin(x*y);
const saddles=[[1.515,0],[-.634,-.327]].map(z=>{
  for(let i=0;i<8;i++)z=add(z,sc(div(dp(z),dd(z)),-1));
  return z;
});
const segs=saddles.map((z,i)=>{
  const a=dd(z),t=Math.PI/4-Math.atan2(a[1],a[0])/2,u=[Math.cos(t),Math.sin(t)],rho=i?.21:.51;
  return [add(z,sc(u,-rho)),add(z,sc(u,rho))];
});
const paths=[[[-2.5,-1.3],segs[1][0]],[segs[1][1],[.6,-1.3],segs[0][0]],[segs[0][1],[2.7,1.4]]];
const edges=paths.flatMap(p=>p.slice(1).map((b,i)=>[p[i],b]));
const v=(a,b,t)=>phiIm(add(sc(a,1-t),sc(b,t)));
const minline=(a,b)=>Math.min(...Array.from({length:2001},(_,i)=>v(a,b,i/2000)));
const ray=(a,sgn,t)=>phiIm(add(a,sc([Math.cos(theta),Math.sin(theta)],sgn*t)));
function simpson(f,b,n=4000){let s=f(0)+f(b);for(let i=1;i<n;i++)s+=(i%2?4:2)*f(i*b/n);return s*b/(3*n);}
const l=paths[0][0],r=paths[2][1];
const vals=[.05,.025,.0125,.00625].map(h=>({h,T:edges.reduce((s,[a,b])=>s+Math.hypot(b[0]-a[0],b[1]-a[1])*simpson(t=>Math.exp(-v(a,b,t)/h),1),0)+simpson(t=>Math.exp(-ray(l,-1,t)/h),10)+simpson(t=>Math.exp(-ray(r,1,t)/h),10)}));
console.log(JSON.stringify({saddles,segs,mins:edges.map(([a,b])=>minline(a,b)),rayMins:[Math.min(...Array.from({length:1001},(_,i)=>ray(l,-1,i/100))),Math.min(...Array.from({length:1001},(_,i)=>ray(r,1,i/100)))],vals},null,2));
