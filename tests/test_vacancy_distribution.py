"""Python/browser parity and independent high-precision stationarity checks."""
import json
import subprocess
from pathlib import Path
import mpmath as mp
import pytest
from solidgas.vacancy_distribution import solve

ROOT = Path(__file__).resolve().parents[1]

@pytest.mark.parametrize('inventory,radius,temperature,A,xi', [
    (94,450,873.15,1.31,.5), (13.45,450,873.15,1.31,.5),
    (36.85,450,873.15,1.31,.5), (193.1,450,873.15,1.31,.5),
    (781.2,450,873.15,1.31,.5), (94,50,1100,1.31,.25),
    (94,450,873.15,0,.5), (94,450,873.15,1.31,1),
    (94,1,873.15,1.31,.5),
])
def test_parity(inventory,radius,temperature,A,xi):
    p=locals().copy()
    py=solve(**p)
    js=json.loads(subprocess.check_output(['node','-e',
        'console.log(JSON.stringify(require("./web/distribution.js").solve('+json.dumps(p)+')))'],cwd=ROOT))
    for key in py:
        assert js[key] == pytest.approx(py[key],abs=1e-11)
    assert abs(py['closure']/py['mean']) < 1e-8
    assert 0 < py['surface'] < .25
    assert py['interior'] <= py['surface']
    if A == 0:
        assert py['surface'] == pytest.approx(py['mean'],abs=1e-12)
    if radius <= 2:
        assert py['shell'] == pytest.approx(1,abs=1e-10)

def test_paper_values():
    r=solve()
    assert 100*r['surface'] == pytest.approx(22.5,abs=.05)
    assert 100*r['interior'] == pytest.approx(.338,abs=.0005)
    assert 100*r['shell'] == pytest.approx(11.1,abs=.05)
    assert 100*solve(xi=.25)['shell'] == pytest.approx(6.4,abs=.05)
    assert 100*solve(xi=1)['shell'] == pytest.approx(19.3,abs=.05)

def test_independent_80_digit_oracle():
    # Independent quadrature in depth and high-precision local chemical potential.
    # The production midpoint mesh is accepted to 1e-6 relative inventory error.
    mp.mp.dps=80
    r=solve(); mu=mp.mpf(str(r['mu'])); kt=mp.mpf('8.617333262145e-5')*mp.mpf('873.15')
    def occ(z):
        energy=-mp.mpf('1.31')*mp.exp(-z/mp.mpf('.5'))
        def f(x):
            return energy+kt*(mp.log(x/(1-x))+2*mp.log(4*x/(1-4*x)))-mu
        lo=mp.mpf('1e-60'); hi=mp.mpf('.25')-mp.mpf('1e-60')
        for _ in range(100):
            mid=(lo+hi)/2
            if f(mid)<0: lo=mid
            else: hi=mid
        return (lo+hi)/2
    # Gauss-Legendre quadrature independently resolves the surface gradient.
    nodes,weights=mp.gauss_quadrature(48,'legendre')
    def integral(a,b):
        return (b-a)/2*sum(weights[i]*occ((a+b)/2+(b-a)/2*nodes[i])*3*(450-((a+b)/2+(b-a)/2*nodes[i]))**2/450**3 for i in range(48))
    total=sum(integral(mp.mpf(a),mp.mpf(b)) for a,b in [(0,.5),(.5,2),(2,10),(10,450)])
    assert abs(float(total)/r['mean']-1) < 1e-6
    assert float(occ(0)) == pytest.approx(r['surface'],abs=1e-12)
    assert float(integral(mp.mpf(0),mp.mpf(2))/total) == pytest.approx(r['shell'],abs=1e-6)

@pytest.mark.parametrize('p',[{'inventory':0},{'inventory':7000},{'radius':0},{'xi':0},{'temperature':0},{'A':-1}])
def test_invalid(p):
    with pytest.raises(ValueError): solve(**p)
