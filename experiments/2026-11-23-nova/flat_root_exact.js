// Nova, 2026-11-23. Exact integer check of the rational degree-6
// Taylor approximation to f(z0) in correspondence 020.
// Run: node flat_root_exact.js. No floating point is used in the assertion.
const D=1000000n,z=[-284121n,-69441n];
const mul=(a,b)=>[a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]];
const add=(a,b)=>[a[0]+b[0],a[1]+b[1]],scale=(a,k)=>a.map(x=>x*k);
const B=2n*D*D,w=scale(mul(z,z),-1n),Q=720n*B**6n;
let p=[0n,0n],pow=[1n,0n],fact=1n;
for(let k=0n;k<=6n;k++){if(k){pow=mul(pow,w);fact*=k;}p=add(p,scale(pow,720n/fact*B**(6n-k)));}
const f=add(add(scale(mul([5n,1n],z),Q),scale([3n,1n],D*Q)),scale(mul(z,p),6n));
const den=2n*D*Q;
const check=(f[0]*f[0]+f[1]*f[1])*10n**14n<49n*den*den;
if(!check)throw new Error("Residual check failed");
console.log({rationalCheck:check,approx:f.map(x=>Number(x)/Number(den))});
