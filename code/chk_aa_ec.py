import sympy as sp
s,t,nu=sp.symbols('s t nu',positive=True)
c=1/(2**((nu-2)/2)*sp.gamma(nu/2))
f=c*s**(nu-1)*sp.exp(-s**2/2)          # chi_nu density
rho={}
for j in range(1,5):
    d=sp.diff(f,s,j-1)
    r=(2*sp.pi)**sp.Rational(-j,2)*(-1)**(j-1)*d
    r=sp.simplify(r.subs(s,sp.sqrt(t)))
    base=(2*sp.pi)**sp.Rational(-j,2)*c*t**((nu-j)/2)*sp.exp(-t/2)
    poly=sp.expand(sp.simplify(r/base))
    rho[j]=poly
    print(j, sp.collect(sp.expand(poly),t))
# check nu=1 against 2*Gaussian EC density at sqrt(t)
u=sp.symbols('u',positive=True)
H={1:sp.Integer(1),2:u,3:u**2-1,4:u**3-3*u}
for j in range(1,5):
    base=(2*sp.pi)**sp.Rational(-j,2)*c*t**((nu-j)/2)*sp.exp(-t/2)
    lhs=(base*rho[j]).subs(nu,1)
    rhs=2*(2*sp.pi)**sp.Rational(-(j+1),2)*H[j].subs(u,sp.sqrt(t))*sp.exp(-t/2)
    print('nu=1 check D=%d:'%j, sp.simplify(lhs-rhs))
# Worsley general formula check
def worsley(D):
    tot=0
    for j in range(0,(D-1)//2+1):
        for k in range(0,D-1-2*j+1):
            tot+= sp.binomial(nu-1,D-1-k-2*j)*(-1)**(D-1+k+j)*sp.factorial(D-1)/(sp.factorial(j)*sp.factorial(k)*2**j)*t**(j+k)
    return sp.expand(tot)
for D in range(1,5):
    for n in [2,3,4,5,6,7,10]:
        a=sp.expand(worsley(D).subs(nu,n)); b=sp.expand(rho[D].subs(nu,n))
        if sp.simplify(a-b)!=0: print('mismatch',D,n,a,b)
print('worsley general formula agrees for nu in 2..10 (ffunc binomial)')
for n in [6]:
    for D in range(1,5): print('nu=6 D=%d'%D, sp.factor(rho[D].subs(nu,n)))
