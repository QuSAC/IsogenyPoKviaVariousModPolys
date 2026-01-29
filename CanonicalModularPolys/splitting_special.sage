# Computational verification of splitting  behavior for special j-invariants in small characteristics
R.<X,j> = PolynomialRing(ZZ)

# These are the canonical modular polynomials Phi_ell^c(X.j)
CanonModPoly = {
    2:  X^3 + 48*X^2 + 768*X + 4096 - X*j ,
    3:  X^4 + 36*X^3 + 270*X^2 + 756*X + 729 - X*j ,
    5:  X^6 + 30*X^5 + 315*X^4 + 1300*X^3 + 1575*X^2 + 750*X + 125 - X*j,
    7:  X^8 + 28*X^7 + 322*X^6 + 1904*X^5 + 5915*X^4 + 8624*X^3 + 4018*X^2 + 748*X + 49 - X*j,
    13: X^14 + 26*X^13 + 325*X^12 + 2548*X^11 + 13832*X^10 + 54340*X^9 + 157118*X^8 + 333580*X^7 + 509366*X^6 + 534820*X^5 + 354536*X^4 + 124852*X^3 + 15145*X^2 + 746*X + 13 - X*j
}

# These are the maximal primes p_ell for ell <= 13
max_prime = {
    2: 13,
    3: 53,
    5: 379,
    7: 1217,
    11: 5101,
    13: 8387
}

for l in [2,3,5,7,13]:
    # In this script we confirm that Phi_ell^c(X, j*) splits over F_p^2 for p <= p_ell if j* in {0, 1728} is supersingular
    p_l = max_prime[l]

    for p in prime_range(p_l + 1):
        if p != l:
            # Define the quadratic extension F_p^2, and view Phi_ell^c(X,j) as a polynomial with coefficients in this field
            F.<a> = GF(p^2)
            S.<Xp, jp> = F[]
            Phi = CanonModPoly[l](X=Xp, j=jp)

            # Recall that j = 0 is supersingular if and only if p != 1 mod 3
            if p % 3 != 1:
                # Check that all irreducible factors of Phi_ell^c(X,0) in F[X] have degree 1
                assert max([h.degree() for (h,e) in Phi(jp=0).factor()]) == 1

            # Recall that j = 1728 is supersingular if and only if p != 1 mod 4
            if p % 4 != 1:
                # Check that all irreducible factors of Phi_ell^c(X,1728) in F[X] have degree 1
                assert max([h.degree() for (h,e) in Phi(jp=1728).factor()]) == 1                

    print(f"All claims are correct for l = {l}!")
