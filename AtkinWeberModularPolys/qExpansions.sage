# Compute the first terms of the q-expansion of the normalized Hauptmodul of X_0+(ell) when X_0(ell) has genus 0
from sage.modular.etaproducts import qexp_eta #Dedekind eta function divided by q^(1/24)
T = QQ[['q']]
T.inject_variables()

for l in [2,3,5,7,13]:
    s = 12//(l-1)
    # Low precision is sufficient to verify the first few terms
    prc = 3
    etanorm = qexp_eta(T,prc)
    t = (etanorm/etanorm(q^l))^(2*s)/q^((l-1)*s*2/24)
    tplus = t + 2*s + l^s/t
    print('\nell=' + str(l) + ': t^+ =', str(tplus))
