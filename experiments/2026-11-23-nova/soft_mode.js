// Nova, 2026-11-23. Run: node soft_mode.js
// Floating-point illustration of reduced integral only; not certified quadrature.
const U=1e-7,B=.25,X=100,P=20000;
const results=[];
for(const N of [3,10,40,160]){
 const T=2*N*Math.sin(Math.PI/(2*N))/Math.sqrt(1+U/2), vs=Array.from({length:N},(_,j)=>Math.sin(Math.PI*j/N)**2);
 const val=(x,y)=>{const re=x*x-y*y,im=2*x*y;let sr=0,si=0;for(const v of vs){const z=Math.exp(-v*re/2),t=-v*im/2;sr+=z*Math.cos(t);si+=z*Math.sin(t);}const pr=T*re/8+T*sr/N,pi=T*im/8+T*si/N;return [Math.exp(-pi)*Math.cos(pr),Math.exp(-pi)*Math.sin(pr)];};
 function sim(f,a,b,n){let r=0,i=0;for(let k=0;k<=n;k++){const z=f(a+(b-a)*k/n),w=k===0||k===n?1:k%2?4:2;r+=w*z[0];i+=w*z[1];}return [r*(b-a)/(3*n),i*(b-a)/(3*n)];}
 const hor=sim(x=>val(x,B),0,X,P),vert=sim(y=>val(0,y),-B,B,2000);const J=[2*hor[0]-vert[1],2*hor[1]+vert[0]];
 results.push({N,T,J,abs:Math.hypot(...J)});
}
console.log(JSON.stringify(results,null,2));
