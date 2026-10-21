import sympy as sp
u=sp.symbols('u',positive=True); r2=sp.sqrt(2)
c=(1-u**2)/(1+u**2); s=2*u/(1+u**2); k=(c+s)/2; h=sp.Rational(1,2)
def T(H):
    X=1+k
    U=h+(k+h-s*(k-h))/c; V=h+(k+h-s*(H-h-k))/c
    Up=h+(k+h-c*(k-h))/s; Vp=h+(k+h-c*(H-h-k))/s
    return {'X':X,'U':U,"V'":Vp},{'X':X,'V':V,"U'":Up}
umax=r2-1
for Hn,H in (('lo',c+s),('kink',h+r2),('2',sp.Integer(2))):
    m=sp.Min(sp.Rational(3,2),2+r2-H) if Hn=='lo' else (sp.Rational(3,2) if Hn=='kink' else 2+r2-H)
    if Hn=='lo': m=sp.Rational(3,2)   # c+s <= sqrt2 < 1/2+sqrt2
    L,R=T(H)
    for a in L:
        for b in R:
            E=sp.together(sp.simplify(L[a]+R[b]-2*m)); num,den=sp.fraction(sp.factor(E))
            print(Hn,a,b,' num =',sp.factor(num),'   den =',sp.factor(den))
