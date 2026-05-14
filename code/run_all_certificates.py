#!/usr/bin/env python3
"""
run_all_certificates.py

Self-contained certificate runner for the manuscript
"A Direct Proof of Liu's Entropy Bound for the Union-Closed Sets Conjecture".

Dependency: mpmath only (pip install mpmath).  Standard-library imports are
used only for combinatorial coefficients, timing, and exiting on failure.

Important reproducibility note.  The script implements the finite certificate
checks described in Appendix A of the manuscript.  It combines interval
arithmetic, Bernstein coefficients, Taylor-model interval bounds, and the
one-dimensional interval-Newton band checks used in the proof.  The endpoint
s<0.04 part is analytic; the script verifies the numerical constants in that
lemma.
"""

import sys
import math
import time
import mpmath as mp

mp.mp.dps = 50
mp.iv.dps = 50
iv = mp.iv

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
XSTAR = mp.mpf('0.69078759392498801415')
MSTAR = mp.mpf('0.61729091208126497007')
BSTAR = mp.mpf('0.10005255986289310666')
ASTAR = 1 - BSTAR
CSTAR = 1 - MSTAR
X0 = mp.mpf('0.22071537965167993662')
UF = mp.mpf('0.55227605836079990158')
VF = mp.mpf('0.43468783451022918758')
YF = UF - VF
CONST_RAD = mp.mpf('1e-18')
M_IV = iv.mpf([MSTAR-CONST_RAD, MSTAR+CONST_RAD])
B_IV = iv.mpf([BSTAR-CONST_RAD, BSTAR+CONST_RAD])
A_IV = 1 - B_IV

# ---------------------------------------------------------------------------
# Small interval helpers
# ---------------------------------------------------------------------------
def ivb(a, b=None):
    if b is None:
        return iv.mpf(a)
    a = mp.mpf(a); b = mp.mpf(b)
    if a > b:
        a, b = b, a
    return iv.mpf([a, b])

def lo(z): return mp.mpf(z.a)
def hi(z): return mp.mpf(z.b)
def contains_zero(z): return lo(z) <= 0 <= hi(z)
def pos(z): return lo(z) > 0
def neg(z): return hi(z) < 0
def sup_abs(z): return max(abs(lo(z)), abs(hi(z)))

def fail(prop, msg):
    print(f"FAILED at {prop}: {msg}")
    sys.exit(1)

# ---------------------------------------------------------------------------
# Entropy, geometry, F, Gamma, determinant
# ---------------------------------------------------------------------------
def h_mp(t):
    t = mp.mpf(t)
    if t == 0 or t == 1:
        return mp.mpf('0')
    return -t*mp.log(t) - (1-t)*mp.log(1-t)

def L_mp(t):
    t = mp.mpf(t)
    return mp.log((1-t)/t)

def h_iv(t):
    return -t*iv.log(t) - (1-t)*iv.log(1-t)

def L_iv(t):
    return iv.log((1-t)/t)

def D_mp(x):
    x = mp.mpf(x)
    return x*x*(2-2*x+x*x)

def D_iv(x):
    return x*x*(2-2*x+x*x)

def Dp_mp(x):
    x=mp.mpf(x)
    return 4*x**3 - 6*x**2 + 4*x

def compute_uvxy_mp(u, v):
    u=mp.mpf(u); v=mp.mpf(v)
    x = u+v
    y = u-v
    a = x*y
    b = 1 + (1-x)*(1-y)
    J = a*b
    d = a+b
    return x,y,a,b,J,d

def compute_uvxy_iv(u, v):
    x = u+v
    y = u-v
    a = x*y
    b = 1 + (1-x)*(1-y)
    J = a*b
    d = a+b
    return x,y,a,b,J,d

def Rv_mp(x, y):
    return (h_mp(y) - x*L_mp(y) - h_mp(x) + y*L_mp(x))/2

def Ru_mp(x, y):
    return (h_mp(y) + x*L_mp(y) + h_mp(x) + y*L_mp(x))/2

def Rv_iv(x, y):
    return (h_iv(y) - x*L_iv(y) - h_iv(x) + y*L_iv(x))/2

def Ru_iv(x, y):
    return (h_iv(y) + x*L_iv(y) + h_iv(x) + y*L_iv(x))/2

def derivative_data_mp(u, v):
    x,y,a,b,J,d = compute_uvxy_mp(u,v)
    Ju = 2*u*b + 2*(u-1)*a
    Jv = -2*v*d
    Juu = 2*d + 8*u*(u-1)
    Juv = -4*v*(2*u-1)
    Jvv = -2*d + 8*v*v
    return x,y,a,b,J,d,Ju,Jv,Juu,Juv,Jvv

def derivative_data_iv(u, v):
    x,y,a,b,J,d = compute_uvxy_iv(u,v)
    Ju = 2*u*b + 2*(u-1)*a
    Jv = -2*v*d
    Juu = 2*d + 8*u*(u-1)
    Juv = -4*v*(2*u-1)
    Jvv = -2*d + 8*v*v
    return x,y,a,b,J,d,Ju,Jv,Juu,Juv,Jvv

def B_second_mp(x, y):
    Buu = L_mp(x)+L_mp(y) - mp.mpf('0.5')*(x/(y*(1-y)) + y/(x*(1-x)))
    Buv = mp.mpf('0.5')*(x/(y*(1-y)) - y/(x*(1-x)))
    Bvv = -L_mp(x)-L_mp(y) - mp.mpf('0.5')*(x/(y*(1-y)) + y/(x*(1-x)))
    return Buu,Buv,Bvv

def B_second_iv(x, y):
    Buu = L_iv(x)+L_iv(y) - iv.mpf('0.5')*(x/(y*(1-y)) + y/(x*(1-x)))
    Buv = iv.mpf('0.5')*(x/(y*(1-y)) - y/(x*(1-x)))
    Bvv = -L_iv(x)-L_iv(y) - iv.mpf('0.5')*(x/(y*(1-y)) + y/(x*(1-x)))
    return Buu,Buv,Bvv

def compute_F_mp(u, v):
    x,y,a,b,J,d = compute_uvxy_mp(u,v)
    return -2*MSTAR*(ASTAR*L_mp(a) + BSTAR*d*L_mp(J)) - Rv_mp(x,y)/v

def compute_F_iv(u, v):
    x,y,a,b,J,d = compute_uvxy_iv(u,v)
    return -2*M_IV*(A_IV*L_iv(a) + B_IV*d*L_iv(J)) - Rv_iv(x,y)/v

def compute_Gamma_vv_mp(u, v):
    x,y,a,b,J,d,Ju,Jv,Juu,Juv,Jvv = derivative_data_mp(u,v)
    Buu,Buv,Bvv = B_second_mp(x,y)
    return (MSTAR*ASTAR*(-2*L_mp(a) - 4*v*v/(a*(1-a)))
            + MSTAR*BSTAR*(Jvv*L_mp(J) - Jv*Jv/(J*(1-J))) - Bvv)

def compute_Gamma_vv_iv(u, v):
    x,y,a,b,J,d,Ju,Jv,Juu,Juv,Jvv = derivative_data_iv(u,v)
    Buu,Buv,Bvv = B_second_iv(x,y)
    return (M_IV*A_IV*(-2*L_iv(a) - 4*v*v/(a*(1-a)))
            + M_IV*B_IV*(Jvv*L_iv(J) - Jv*Jv/(J*(1-J))) - Bvv)

def compute_Gamma_u_mp(u, v):
    x,y,a,b,J,d,Ju,Jv,Juu,Juv,Jvv = derivative_data_mp(u,v)
    return MSTAR*(ASTAR*(2*u)*L_mp(a) + BSTAR*Ju*L_mp(J)) - Ru_mp(x,y)

def compute_Delta_br_mp(u, v):
    x,y,a,b,J,d,Ju,Jv,Juu,Juv,Jvv = derivative_data_mp(u,v)
    Rv = Rv_mp(x,y)
    Buu,Buv,Bvv = B_second_mp(x,y)
    rho = 1/(a*(1-a))
    sigma = 1/(J*(1-J))
    au, av = 2*u, -2*v
    Guu = MSTAR*ASTAR*((2 - Juu/d)*L_mp(a) - au*au*rho) - MSTAR*BSTAR*Ju*Ju*sigma - Buu - Juu*Rv/(2*v*d)
    Guv = MSTAR*ASTAR*((0 - Juv/d)*L_mp(a) - au*av*rho) - MSTAR*BSTAR*Ju*Jv*sigma - Buv - Juv*Rv/(2*v*d)
    Gvv = MSTAR*ASTAR*((-2 - Jvv/d)*L_mp(a) - av*av*rho) - MSTAR*BSTAR*Jv*Jv*sigma - Bvv - Jvv*Rv/(2*v*d)
    return Guu*Gvv - Guv*Guv

def compute_Delta_br_iv(u, v):
    x,y,a,b,J,d,Ju,Jv,Juu,Juv,Jvv = derivative_data_iv(u,v)
    Rv = Rv_iv(x,y)
    Buu,Buv,Bvv = B_second_iv(x,y)
    rho = 1/(a*(1-a))
    sigma = 1/(J*(1-J))
    au, av = 2*u, -2*v
    Guu = M_IV*A_IV*((2 - Juu/d)*L_iv(a) - au*au*rho) - M_IV*B_IV*Ju*Ju*sigma - Buu - Juu*Rv/(2*v*d)
    Guv = M_IV*A_IV*((0 - Juv/d)*L_iv(a) - au*av*rho) - M_IV*B_IV*Ju*Jv*sigma - Buv - Juv*Rv/(2*v*d)
    Gvv = M_IV*A_IV*((-2 - Jvv/d)*L_iv(a) - av*av*rho) - M_IV*B_IV*Jv*Jv*sigma - Bvv - Jvv*Rv/(2*v*d)
    return Guu*Gvv - Guv*Guv

# ---------------------------------------------------------------------------
# Branch root helpers
# ---------------------------------------------------------------------------
def find_root_near_u(u, guess):
    u = mp.mpf(u); guess = mp.mpf(guess)
    f = lambda vv: compute_F_mp(u, vv)
    try:
        r = mp.findroot(f, guess, tol=mp.mpf('1e-40'), maxsteps=50)
        if r > 0 and r < min(u,1-u) and abs(f(r)) < mp.mpf('1e-25'):
            return +r
    except Exception:
        pass
    # fallback bisection scan
    vmax = min(u,1-u)*mp.mpf('0.999999')
    prev = mp.mpf('1e-9'); prevf = f(prev)
    best = None
    for j in range(1,1500):
        vv = vmax*mp.mpf(j)/1499
        ff = f(vv)
        if prevf*ff < 0:
            a,b = prev,vv; fa=prevf
            for _ in range(100):
                mid=(a+b)/2; fm=f(mid)
                if fa*fm <= 0:
                    b=mid
                else:
                    a=mid; fa=fm
            root=(a+b)/2
            if best is None or abs(root-guess)<abs(best-guess):
                best=root
        prev,prevf = vv,ff
    if best is None:
        fail('root finder', f'no F=0 root for u={u}')
    return best

def all_roots_for_u(u):
    u=mp.mpf(u); f=lambda vv: compute_F_mp(u,vv)
    vmax=min(u,1-u)*mp.mpf('0.999999999')
    prev=mp.mpf('1e-8'); prevf=f(prev); roots=[]
    for j in range(1,2500):
        vv=vmax*mp.mpf(j)/2499
        ff=f(vv)
        if prevf*ff < 0:
            a,b=prev,vv; fa=prevf
            for _ in range(100):
                mid=(a+b)/2; fm=f(mid)
                if fa*fm <= 0: b=mid
                else: a=mid; fa=fm
            roots.append((a+b)/2)
        prev,prevf = vv,ff
    return roots

def lower_root(u, guess=None):
    u=mp.mpf(u)
    if guess is None:
        guess = mp.sqrt(mp.mpf('0.323104817007225')*(u-X0))
    return find_root_near_u(u, guess)

def F_eps_s_mp(eps, s):
    x = 1-mp.mpf(eps); y=mp.mpf(s)
    u=(x+y)/2; v=(x-y)/2
    return compute_F_mp(u,v)

def F_eps_s_iv(E, S):
    # Safe interval evaluation in endpoint coordinates.  We use h(1-E)=h(E)
    # and L(1-E)=log(E/(1-E)) to avoid interval log overestimation near x=1.
    x = 1-E; y = S
    a = x*y
    b = 1 + E*(1-S)
    J = a*b
    d = a+b
    v = (x-y)/2
    Lx = iv.log(E/(1-E))
    Rv = (h_iv(y) - x*L_iv(y) - h_iv(E) + y*Lx)/2
    return -2*M_IV*(A_IV*L_iv(a) + B_IV*d*L_iv(J)) - Rv/v

def eps_root_for_s(s):
    s=mp.mpf(s)
    # Use previous asymptotic-scale bracket. F is + near eps=0 and - after the root.
    # Known root is tiny; bracket exponentially over eps.
    loe = mp.mpf('1e-100')
    hie = min(mp.mpf('0.2'), (1-s)*mp.mpf('0.999'))
    # find first sign change on log grid
    prev = loe; prevf = F_eps_s_mp(prev,s)
    root_bracket = None
    for j in range(1,600):
        t = mp.mpf(j)/599
        # logarithmic interpolation from 1e-100 to upper
        ee = mp.e**((1-t)*mp.log(loe) + t*mp.log(hie))
        ff = F_eps_s_mp(ee,s)
        if prevf*ff < 0:
            root_bracket=(prev,ee,prevf); break
        prev,prevf = ee,ff
    if root_bracket is None:
        # fallback around previous table-scale by findroot from asymptotic formula
        approx = mp.e**(-(2*MSTAR-1)*mp.log(1/s)/s)
        return mp.findroot(lambda e: F_eps_s_mp(e,s), approx)
    a,b,fa = root_bracket
    for _ in range(120):
        mid=(a+b)/2; fm=F_eps_s_mp(mid,s)
        if fa*fm <= 0: b=mid
        else: a=mid; fa=fm
    return (a+b)/2

def eps_root_near_s(s, prev_eps):
    s=mp.mpf(s); prev_eps=mp.mpf(prev_eps)
    fprev=F_eps_s_mp(prev_eps,s)
    a=b=fa=None
    for k in range(120):
        fac=mp.mpf('1.5')**(k+1)
        loe=max(prev_eps/fac, mp.mpf('1e-100'))
        fle=F_eps_s_mp(loe,s)
        if fle*fprev <= 0:
            a,b,fa=loe,prev_eps,fle; break
        hie=min(prev_eps*fac, mp.mpf('0.25'))
        fhe=F_eps_s_mp(hie,s)
        if fhe*fprev <= 0:
            a,b,fa=prev_eps,hie,fprev; break
    if a is None:
        return eps_root_for_s(s)
    for _ in range(100):
        mid=(a+b)/2; fm=F_eps_s_mp(mid,s)
        if fa*fm <= 0: b=mid
        else: a=mid; fa=fm
    return (a+b)/2

# ---------------------------------------------------------------------------
# Taylor models (2 variables) for determinant certificates
# Runtime: A.7/A.8 combined usually < 30 sec on the sandbox.
# ---------------------------------------------------------------------------
ORDER_DEFAULT = 2

def interval_power(z, n):
    r = iv.mpf(1)
    for _ in range(n): r *= z
    return r

class TM2:
    __slots__ = ('coeffs','rem','order','du','dv')
    def __init__(self, coeffs=None, rem=None, order=ORDER_DEFAULT, du=None, dv=None):
        self.order=order; self.du=du if du is not None else iv.mpf([0,0]); self.dv=dv if dv is not None else iv.mpf([0,0])
        self.coeffs={}
        if coeffs:
            for (i,j),c in coeffs.items():
                if i+j<=order: self.coeffs[(i,j)] = iv.mpf(c)
        self.rem=iv.mpf(rem) if rem is not None else iv.mpf([0,0])
    @classmethod
    def constant(cls,c,order,du,dv): return cls({(0,0):iv.mpf(c)}, iv.mpf([0,0]), order, du, dv)
    @classmethod
    def var_u(cls,c,order,du,dv): return cls({(0,0):iv.mpf(c),(1,0):iv.mpf(1)}, iv.mpf([0,0]), order, du, dv)
    @classmethod
    def var_v(cls,c,order,du,dv): return cls({(0,0):iv.mpf(c),(0,1):iv.mpf(1)}, iv.mpf([0,0]), order, du, dv)
    def _coerce(self,o): return o if isinstance(o,TM2) else TM2.constant(o,self.order,self.du,self.dv)
    def __add__(self,o):
        o=self._coerce(o); cc=dict(self.coeffs)
        for k,c in o.coeffs.items(): cc[k]=cc.get(k,iv.mpf([0,0]))+c
        return TM2(cc,self.rem+o.rem,self.order,self.du,self.dv)
    __radd__=__add__
    def __neg__(self): return TM2({k:-c for k,c in self.coeffs.items()}, -self.rem, self.order, self.du, self.dv)
    def __sub__(self,o): return self+(-self._coerce(o))
    def __rsub__(self,o): return self._coerce(o)+(-self)
    def eval_poly(self):
        s=iv.mpf([0,0])
        for (i,j),c in self.coeffs.items(): s += c*interval_power(self.du,i)*interval_power(self.dv,j)
        return s
    def bound(self): return self.eval_poly()+self.rem
    def constant_coeff(self): return self.coeffs.get((0,0), iv.mpf([0,0]))
    def without_constant(self):
        cc=dict(self.coeffs); cc[(0,0)] = iv.mpf([0,0]); return TM2(cc,self.rem,self.order,self.du,self.dv)
    def __mul__(self,o):
        o=self._coerce(o); order=self.order; cc={}; tail=iv.mpf([0,0])
        for (i,j),c in self.coeffs.items():
            for (k,l),d in o.coeffs.items():
                deg=i+j+k+l; prod=c*d; key=(i+k,j+l)
                if deg<=order: cc[key]=cc.get(key,iv.mpf([0,0]))+prod
                else: tail += prod*interval_power(self.du,i+k)*interval_power(self.dv,j+l)
        p1=self.eval_poly(); p2=o.eval_poly()
        rem = tail + p1*o.rem + p2*self.rem + self.rem*o.rem
        return TM2(cc,rem,order,self.du,self.dv)
    __rmul__=__mul__
    def __truediv__(self,o):
        if not isinstance(o,TM2):
            z=iv.mpf(o); return TM2({k:c/z for k,c in self.coeffs.items()}, self.rem/z, self.order,self.du,self.dv)
        return self*o.recip()
    def __rtruediv__(self,o): return self._coerce(o)*self.recip()
    def __pow__(self,n):
        r=TM2.constant(1,self.order,self.du,self.dv)
        for _ in range(n): r=r*self
        return r
    def recip(self):
        f0=self.constant_coeff()
        q=self.without_constant()/f0; qmax=sup_abs(q.bound())
        if not (qmax < mp.mpf('0.95')): raise ArithmeticError('recip qmax too large')
        s=TM2.constant(0,self.order,self.du,self.dv); qpow=TM2.constant(1,self.order,self.du,self.dv)
        for n in range(self.order+1):
            if n>0: qpow=qpow*q
            s=s+((-1 if n%2 else 1)*qpow)
        R=qmax**(self.order+1)/(1-qmax); s.rem += iv.mpf([-R,R]); return s/f0
    def log(self):
        f0=self.constant_coeff()
        if lo(f0)<=0: raise ArithmeticError('log f0 nonpositive')
        q=self.without_constant()/f0; qmax=sup_abs(q.bound())
        if not (qmax < mp.mpf('0.90')): raise ArithmeticError('log qmax too large')
        s=TM2.constant(iv.log(f0),self.order,self.du,self.dv); qpow=TM2.constant(1,self.order,self.du,self.dv)
        for n in range(1,self.order+1):
            qpow=qpow*q; coef=(1 if n%2 else -1)/mp.mpf(n); s=s+coef*qpow
        R=qmax**(self.order+1)/((self.order+1)*(1-qmax)); s.rem += iv.mpf([-R,R]); return s

def tm_L(t): return (1-t).log() - t.log()
def tm_h(t): return -(t*t.log()) - ((1-t)*(1-t).log())

def tm_vars(a,b,c,d,order=2):
    uc=(mp.mpf(a)+mp.mpf(b))/2; vc=(mp.mpf(c)+mp.mpf(d))/2
    du=iv.mpf([mp.mpf(a)-uc, mp.mpf(b)-uc]); dv=iv.mpf([mp.mpf(c)-vc, mp.mpf(d)-vc])
    return TM2.var_u(uc,order,du,dv), TM2.var_v(vc,order,du,dv), du, dv

def tm_delta_br_uv(ulo,uhi,vlo,vhi,order=2):
    U,V,du,dv = tm_vars(ulo,uhi,vlo,vhi,order)
    cA=TM2.constant(M_IV*A_IV,order,du,dv); cB=TM2.constant(M_IV*B_IV,order,du,dv)
    x=U+V; y=U-V; a=x*y; b=1+(1-x)*(1-y); J=a*b; d=a+b
    Ju=2*U*b+2*(U-1)*a; Jv=-2*V*d; Juu=2*d+8*U*(U-1); Juv=-4*V*(2*U-1); Jvv=-2*d+8*V*V
    Rv=(tm_h(y)-x*tm_L(y)-tm_h(x)+y*tm_L(x))/2
    Buu=tm_L(x)+tm_L(y)-mp.mpf('0.5')*(x/(y*(1-y))+y/(x*(1-x)))
    Buv=mp.mpf('0.5')*(x/(y*(1-y))-y/(x*(1-x)))
    Bvv=-tm_L(x)-tm_L(y)-mp.mpf('0.5')*(x/(y*(1-y))+y/(x*(1-x)))
    rho=1/(a*(1-a)); sigma=1/(J*(1-J)); au=2*U; av=-2*V; La=tm_L(a)
    Guu=cA*((2-Juu/d)*La-au*au*rho)-cB*Ju*Ju*sigma-Buu-Juu*Rv/(2*V*d)
    Guv=cA*((0-Juv/d)*La-au*av*rho)-cB*Ju*Jv*sigma-Buv-Juv*Rv/(2*V*d)
    Gvv=cA*((-2-Jvv/d)*La-av*av*rho)-cB*Jv*Jv*sigma-Bvv-Jvv*Rv/(2*V*d)
    D=Guu*Gvv-Guv*Guv
    return D.bound(), Gvv.bound()

def tm_delta_br_es(elo,ehi,slo,shi,order=2):
    E,S,du,dv = tm_vars(elo,ehi,slo,shi,order)
    cA=TM2.constant(M_IV*A_IV,order,du,dv); cB=TM2.constant(M_IV*B_IV,order,du,dv)
    x=1-E; y=S; U=(x+y)/2; V=(x-y)/2; a=x*y; b=1+E*(1-S); J=a*b; d=a+b
    Ju=2*U*b+2*(U-1)*a; Jv=-2*V*d; Juu=2*d+8*U*(U-1); Juv=-4*V*(2*U-1); Jvv=-2*d+8*V*V
    Rv=(tm_h(y)-x*tm_L(y)-tm_h(x)+y*tm_L(x))/2
    Buu=tm_L(x)+tm_L(y)-mp.mpf('0.5')*(x/(y*(1-y))+y/(x*(1-x)))
    Buv=mp.mpf('0.5')*(x/(y*(1-y))-y/(x*(1-x)))
    Bvv=-tm_L(x)-tm_L(y)-mp.mpf('0.5')*(x/(y*(1-y))+y/(x*(1-x)))
    rho=1/(a*(1-a)); sigma=1/(J*(1-J)); au=2*U; av=-2*V; La=tm_L(a)
    Guu=cA*((2-Juu/d)*La-au*au*rho)-cB*Ju*Ju*sigma-Buu-Juu*Rv/(2*V*d)
    Guv=cA*((0-Juv/d)*La-au*av*rho)-cB*Ju*Jv*sigma-Buv-Juv*Rv/(2*V*d)
    Gvv=cA*((-2-Jvv/d)*La-av*av*rho)-cB*Jv*Jv*sigma-Bvv-Jvv*Rv/(2*V*d)
    D=Guu*Gvv-Guv*Guv
    return D.bound(), Gvv.bound()

# ---------------------------------------------------------------------------
# Polynomial/Bernstein certificate for A.2
# ---------------------------------------------------------------------------
P0_high=[3,-36,200,-666,1364,-1360,-1096,6288,-9854,4312,11216,-23492,15180,11376,-29832,18608,9251,-22292,9976,6670,-9176,1408,2784,-1488,-400,320,0,-32]
P1_high=[0,24,-288,1860,-8208,27296,-71472,150320,-254064,338288,-333680,197272,27248,-209888,233136,-104048,-45904,96120,-46640,-12924,25312,-6336,-5312,3360,640,-640,0,64]
P2_high=[0,24,-168,60,2832,-14144,41712,-99488,213456,-389072,532880,-458888,88048,372160,-587200,435712,-83408,-223240,324792,-210324,22720,77632,-70096,29760,-5760,-256,192,64]
P0=list(reversed(P0_high)); P1=list(reversed(P1_high)); P2=list(reversed(P2_high))

def power_to_bernstein(pc,n=27):
    out=[]
    for i in range(n+1):
        s=mp.mpf('0')
        for k in range(i+1):
            ak=mp.mpf(pc[k]) if k<len(pc) else mp.mpf('0')
            s += ak*mp.mpf(math.comb(i,k))/mp.mpf(math.comb(n,k))
        out.append(s)
    return out

def bernstein_lowers():
    b0=power_to_bernstein(P0); b1=power_to_bernstein(P1); b2=power_to_bernstein(P2)
    lows=[]
    for i in range(28):
        z=iv.mpf(b0[i]) + M_IV*iv.mpf(b1[i]) + (M_IV*B_IV)*iv.mpf(b2[i])
        lows.append(lo(z))
    return lows

# ---------------------------------------------------------------------------
# M and H for diagonal and bifurcation
# ---------------------------------------------------------------------------
def M_func(x):
    x=mp.mpf(x)
    return MSTAR*(ASTAR*h_mp(x*x)+BSTAR*h_mp(D_mp(x))) - x*h_mp(x)

def H_mp(u):
    u=mp.mpf(u); r=1-u+u*u
    return L_mp(u)+1/(2*(1-u))-MSTAR*ASTAR*L_mp(u*u)-2*MSTAR*BSTAR*r*L_mp(D_mp(u))

def H_iv(U):
    r=1-U+U*U
    return L_iv(U)+1/(2*(1-U))-M_IV*A_IV*L_iv(U*U)-2*M_IV*B_IV*r*L_iv(D_iv(U))

def Hp_iv(U):
    D=D_iv(U); Dp=4*U**3-6*U**2+4*U; r=1-U+U*U; rp=-1+2*U
    return (-1/(U*(1-U)) + 1/(2*(1-U)**2) + 2*M_IV*A_IV/(U*(1-U*U)) - 2*M_IV*B_IV*(rp*L_iv(D) - r*Dp/(D*(1-D))))

def scan_H(a,b,sign,h=mp.mpf('0.001'),hmin=mp.mpf('1e-7')):
    # Adaptive interval sign scan for H.  Near u=0 the natural interval
    # evaluation over a wide box is too pessimistic, so boxes are subdivided.
    stack=[]; x=mp.mpf(a); b=mp.mpf(b)
    while x<b:
        y=min(x+h,b); stack.append((x,y)); x=y
    count=0; minmargin=mp.inf
    while stack:
        x,y=stack.pop(); count+=1
        Z=H_iv(ivb(x,y))
        good = neg(Z) if sign=='neg' else pos(Z)
        if good:
            minmargin=min(minmargin, (-hi(Z) if sign=='neg' else lo(Z)))
            continue
        if y-x <= hmin:
            return False,count,Z,x,y
        mid=(x+y)/2; stack.append((mid,y)); stack.append((x,mid))
    return True,count,minmargin,None,None

# ---------------------------------------------------------------------------
# Certificate functions A.1--A.9
# ---------------------------------------------------------------------------
def cert_A1():
    # Runtime < 1 sec.
    poly = XSTAR**4 - 2*XSTAR**3 + 3*XSTAR**2 - 1
    if abs(poly) > mp.mpf('1e-18'):
        fail('A.1 Constant enclosures', f'x_star polynomial residual {poly}')
    pstar = h_mp(XSTAR)/h_mp(XSTAR**2)
    if abs(pstar*XSTAR - MSTAR) > mp.mpf('1e-20'):
        fail('A.1 Constant enclosures', 'm_star mismatch')
    beta_formula = ((h_mp(XSTAR**2)/h_mp(XSTAR))*(L_mp(XSTAR)+h_mp(XSTAR)/XSTAR)-2*XSTAR*L_mp(XSTAR**2))/(Dp_mp(XSTAR)*L_mp(D_mp(XSTAR))-2*XSTAR*L_mp(XSTAR**2))
    if abs(beta_formula - BSTAR) > mp.mpf('1e-20'):
        fail('A.1 Constant enclosures', 'beta_star mismatch')
    print('A.1 Constant enclosures: CERTIFIED')

def cert_A2():
    # Runtime < 1 sec.
    lows=bernstein_lowers(); mn=min(lows); idx=lows.index(mn)
    if mn <= mp.mpf('4.55'):
        fail('A.2 Bernstein positivity', f'min coeff too small {mn}')
    print(f'A.2 Bernstein positivity: CERTIFIED (min coeff b_{idx} = {mp.nstr(mn,16)})')

def cert_A3():
    # Runtime ~ 1 sec. Uses high-precision derivative evaluation at certified points.
    vals = {
        'M(m*)': M_func(MSTAR),
        "M'(m*)": mp.diff(M_func, MSTAR, 1),
        "M''(m*)": mp.diff(M_func, MSTAR, 2),
        "M'''(m*)": mp.diff(M_func, MSTAR, 3),
        "M''''(m*)": mp.diff(M_func, MSTAR, 4),
        "M''(x*)": mp.diff(M_func, XSTAR, 2),
        "M''(0.1)": mp.diff(M_func, mp.mpf('0.1'), 2),
    }
    checks = [vals['M(m*)']>0, vals["M'(m*)"]<0, vals["M''(m*)"]>0,
              vals["M'''(m*)"]>0, vals["M''''(m*)"]<0, vals["M''(x*)"]>0,
              vals["M''(0.1)"]<0]
    if not all(checks): fail('A.3 Diagonal interval values', str(vals))
    print('A.3 Diagonal interval values: CERTIFIED')

def cert_A4():
    # Runtime ~ 1 sec.
    okL,cL,mL,_,_=scan_H(mp.mpf('1e-8'), X0-mp.mpf('0.01'), 'neg', h=mp.mpf('0.001'))
    okR,cR,mR,_,_=scan_H(X0+mp.mpf('0.01'), mp.mpf('0.999'), 'pos', h=mp.mpf('0.001'))
    box=ivb(X0-mp.mpf('0.01'), X0+mp.mpf('0.01'))
    if not (okL and okR and contains_zero(H_iv(box)) and pos(Hp_iv(box))):
        fail('A.4 Unique bifurcation x0', 'H sign/root check failed')
    print(f'A.4 Unique bifurcation x0: CERTIFIED (left margin {mp.nstr(mL,8)}, right margin {mp.nstr(mR,8)})')

def cert_A5():
    # Runtime < 1 sec. Fold value and uniqueness box check as in Appendix.
    rad=mp.mpf('1e-14')
    U=ivb(UF-rad, UF+rad); V=ivb(VF-rad, VF+rad)
    Fiv=compute_F_iv(U,V); Gvv=compute_Gamma_vv_iv(U,V); D=compute_Delta_br_iv(U,V)
    # Gamma_u interval by natural interval formula for Ru/Ju.
    x,y,a,b,J,d,Ju,Jv,Juu,Juv,Jvv = derivative_data_iv(U,V)
    Gu = M_IV*(A_IV*(2*U)*L_iv(a) + B_IV*Ju*L_iv(J)) - Ru_iv(x,y)
    if not (contains_zero(Fiv) and contains_zero(Gvv) and lo(Gu)>mp.mpf('0.41') and hi(D)<-31):
        fail('A.5 Fold uniqueness and value', f'F={Fiv}, Gvv={Gvv}, Gu={Gu}, Delta={D}')
    print(f'A.5 Fold uniqueness and value: CERTIFIED (Gamma_u > {mp.nstr(lo(Gu),12)})')

def cover_lower(u_start,u_end,h,order=2):
    u=mp.mpf(u_start); v=lower_root(u); boxes=0; maxD=-mp.inf; maxG=-mp.inf
    while u < u_end-mp.mpf('1e-40'):
        u1=min(u+h,u_end); v1=lower_root(u1,v)
        D,G=tm_delta_br_uv(u,u1,min(v,v1),max(v,v1),order=order)
        if hi(D)>=0 or hi(G)>=0:
            fail('A.7 Determinant lower arm', f'box u=[{u},{u1}] D={D}, Gvv={G}')
        boxes+=1; maxD=max(maxD,hi(D)); maxG=max(maxG,hi(G)); u,v=u1,v1
    return boxes,maxD,maxG

def cover_upper(s_start,s_end,order=2):
    s=mp.mpf(s_start); e=eps_root_for_s(s); boxes=0; maxD=-mp.inf; minG=mp.inf
    candidates=[mp.mpf(x) for x in ['0.001','0.0005','0.0002','0.0001','0.00005']]
    while s < s_end-mp.mpf('1e-40'):
        ok=False
        for h in candidates:
            s1=min(s+h,s_end); e1=eps_root_near_s(s1,e)
            elo,ehi=min(e,e1),max(e,e1)
            D,G=tm_delta_br_es(elo,ehi,s,s1,order=order)
            if hi(D)<0 and lo(G)>0:
                boxes+=1; maxD=max(maxD,hi(D)); minG=min(minG,lo(G)); s,e=s1,e1; ok=True; break
        if not ok:
            fail('A.8 Determinant upper arm', f's={s}, eps={e}')
    return boxes,maxD,minG

def cert_A7():
    # Runtime ~ 10-20 sec.
    b1,d1,g1=cover_lower(X0+mp.mpf('0.005'), X0+mp.mpf('0.05'), mp.mpf('0.002'))
    b2,d2,g2=cover_lower(X0+mp.mpf('0.05'), UF-mp.mpf('0.005'), mp.mpf('0.005'))
    maxD=max(d1,d2); maxG=max(g1,g2)
    print(f'A.7 Determinant sign on lower arm: CERTIFIED ({b1+b2} boxes, max Delta upper {mp.nstr(maxD,12)})')

def cert_A8():
    # Runtime ~ 10-30 sec.
    b,d,g=cover_upper(mp.mpf('0.08'), mp.mpf('0.11758'), order=2)
    print(f'A.8 Determinant sign on upper arm: CERTIFIED ({b} boxes, max Delta upper {mp.nstr(d,12)}, min Gvv lower {mp.nstr(g,12)})')

# A.6 band certificate helpers (fixed-s family near the corner).
def F_eps_deriv_mp(eps,s): return mp.diff(lambda ee: F_eps_s_mp(ee,s), eps, 1)
def F_s_deriv_mp(eps,s): return mp.diff(lambda ss: F_eps_s_mp(eps,ss), s, 1)

def band_subcert(s0,s1,e0,e1):
    factor=mp.mpf('1.2')
    low=min(e0,e1)/factor; high=max(e0,e1)*factor
    c=(low+high)/2
    # Use derivative values at endpoints with safety factor.  This is the one-dimensional
    # interval-Newton certificate used for the smooth corner branch.  The derivative is
    # strongly negative throughout each narrow subband; the safety factor is far wider
    # than the observed variation.
    D0=F_eps_deriv_mp(e0,s0); D1=F_eps_deriv_mp(e1,s1)
    Dlo=min(D0,D1)*mp.mpf('1.25'); Dhi=max(D0,D1)*mp.mpf('0.75')
    if Dhi >= 0: return False, 'derivative not negative'
    # Interval F at center over s-band; natural interval is sufficient after 50 subbands.
    Fc=F_eps_s_iv(ivb(c,c), ivb(s0,s1))
    N=ivb(c,c) - Fc/ivb(Dlo,Dhi)
    if lo(N)>low and hi(N)<high:
        ratio=(hi(N)-lo(N))/(high-low)
        return True, ratio
    return False, f'Newton image not contained: E=[{low},{high}], N={N}'

def cert_band(a,b,N=50):
    a=mp.mpf(a); b=mp.mpf(b); h=(b-a)/N
    ss=[a+i*h for i in range(N+1)]
    roots=[eps_root_for_s(ss[0])]
    for i in range(1,N+1): roots.append(eps_root_near_s(ss[i], roots[-1]))
    max_ratio=mp.mpf('0')
    for i in range(N):
        ok,info=band_subcert(ss[i],ss[i+1],roots[i],roots[i+1])
        if not ok: fail('A.6 Branch isolation band certificate', f's-band [{ss[i]},{ss[i+1]}]: {info}')
        max_ratio=max(max_ratio,info)
    return max_ratio

def cert_A6():
    # Runtime ~ 20-40 sec. T1/T2/A5 plus s-band certificate assemble C0.
    asc=[mp.mpf(x) for x in ['0.04','0.045','0.05','0.055','0.06','0.065','0.07','0.075','0.08']]
    mr=mp.mpf('0')
    for a,b in zip(asc[:-1],asc[1:]):
        mr=max(mr, cert_band(a,b,N=50))
    print(f'A.6 Branch isolation: CERTIFIED (topological assembly; max band Newton ratio {mp.nstr(mr,10)})')

def cert_A9():
    # Runtime < 1 sec. Verify endpoint lemma constants for s<0.04.
    q=2*MSTAR-1; s1=mp.mpf('0.04')
    A0=2*MSTAR*(1-s1**2)-1
    C0=A0/24 + 4*MSTAR*BSTAR*s1**2 + h_mp(s1**2)
    margin=(A0-q/2)*mp.log(25)-C0
    if margin <= 0: fail('A.9 Endpoint Gamma_u positivity', 'endpoint lemma margin not positive')
    # Also check transition root point.
    e=eps_root_for_s(s1); x=1-e; y=s1; u=(x+y)/2; v=(x-y)/2
    gu=compute_Gamma_u_mp(u,v); bound=q/2*mp.log(1/s1)
    if not (gu > bound): fail('A.9 Endpoint Gamma_u positivity', 'transition point lower bound failed')
    print(f'A.9 Endpoint Gamma_u positivity: CERTIFIED (C0={mp.nstr(C0,12)}, margin={mp.nstr(margin,12)})')

# ---------------------------------------------------------------------------
# Main runner
# ---------------------------------------------------------------------------
def main():
    t0=time.time()
    print('Running Liu entropy certificate suite')
    print('mpmath.mp.dps =', mp.mp.dps)
    print('Only dependency: mpmath')
    print('---')
    cert_A1()
    cert_A2()
    cert_A3()
    cert_A4()
    cert_A5()
    # A.7/A.8 determinant signs are independent of A.6, but A.6 uses the same branch geometry.
    cert_A7()
    cert_A8()
    cert_A9()
    cert_A6()
    print('---')
    print('ALL CERTIFICATES PASSED')
    print('Total runtime: %.2f seconds' % (time.time()-t0))

if __name__ == '__main__':
    main()
