########################################
############# General lifting ##########
########################################

def double_prec_squareroot(A, a, m):
    #input: elements A and a so that a^2 = A (mod eps^(2^(m-1))
    #output: lifts A' = A + eps*A1 + ... +eps^{2^m-1)*A_(2^m-1)
    #and a' = a + eps*a1 + ... +eps^{2^m-1)*a_(2^m-1)
    #so that: a'^2 = A' mod eps^(2^m).
    if m>1:
        a = a.lift_to_precision(2^m)
    term = A-a**2
    b = term/(2*a)
    newa = a + b
    return newa
    
    
def double_prec_montgomery_coefficient(A, j_inv, m):
    #given the Montgomery coefficient of an elliptic curve,
    #E: y^2 = x^3 + A*x^2 + x with A in K[eps]/(eps^(2^m)), with 
    #j(E) = j(E) + lam[0]*eps + ... + lam[2^(m-1) -1]*eps^(2^(m-1)) mod eps^(2^(m-1))
    #Compute the newA so that j(Enew) = .... mod(eps^(2^m))
    
    term = 2^8*(A^2-3)^3-j_inv*(A^2-4)
    B = term/(2*A*(j_inv - 3*2^8*(A^2-3)^2))
    newA = A + B

    return newA


########################################
#######lifting isogeny chains ##########
########################################

def double_prec_chain(invariants,isogeny_chain,m):
    '''
    Input: data on a (2^n,2^n)-isogeny chain between products of elliptic curves 
        over a ring K.<eps>/(eps^(2^(m-1))),
        and a lift of the domain to K.<eps>/(eps^(2^(m)))
        (defining again a product isogeny!)
    Output: 
        a lift of the product chain to K.<eps>/(eps^(2^(m)))
    '''

    [A,A1] = invariants
    
    gluing = isogeny_chain[0]
    psi = lift_gluing_kummer_isogeny(invariants,gluing,m)
    lifted_isogeny_chain = [psi]
    
    for phi in isogeny_chain[1:-1]:
        psi = lift_kummer_isogeny(psi[0],phi[1],m)
        lifted_isogeny_chain.append(psi)
    
    splitting = isogeny_chain[-1]
    psi = lift_splitting_kummer_isogeny(psi[0],splitting[1],m)
    lifted_isogeny_chain.append(psi)

    return lifted_isogeny_chain
    

def double_prec_codomain(invariants,isogeny_chain,m):
    '''
    Input: data on a (2^n,2^n)-isogeny chain between products of elliptic curves 
        over a ring K.<eps>/(eps^(2^(m-1))),
        and a lift of the domain to K.<eps>/(eps^(2^(m)))
    Output: 
        a lift of the codomain of the chain (more precisely, the splitting value delta)
    '''

    [A,A1] = invariants
    
    gluing = isogeny_chain[0]
    psi = lift_gluing_kummer_isogeny(invariants,gluing,m)
    lift_isogeny_chain = [psi]
    
    for phi in isogeny_chain[1:-1]:
        psi = lift_kummer_isogeny(psi[0],phi[1],m)
        lift_isogeny_chain.append(psi)
    
    splitting = isogeny_chain[-1]
    delta = lift_determine_codomain(psi[0],splitting[1],m)

    return delta


########################################
########## Lifting an ell-isogeny#######
########################################

def lift_isogeny(E, E1, ker, ell, n):
    '''
    Input: Elliptic curves E, E1 defined over a finite field K and connected by an ell isogeny
    ker is a group so that 2*kernel defines a (2^n,2^n)-isogeny E x E1 -> E x E2 
    Output: codomain of the lifted isogeny over K[eps]/(eps**(ell+2)) with domain j(E) + eps
    '''
    
    m = ceil(log((ell+2),2))
    K = E.base()
    j00 = E.j_invariant()
    S.<eps> = PowerSeriesRing(K)
    A = S(montgomery_coefficient(E),1)
    A1 = S(montgomery_coefficient(E1),1)

    #we compute the isogeny chain over K and auxiliary values
    phi = short_isogeny_chain_kummer(ker, n)

    #we double the precision at each step i=1,...,m
    for i in range(1,m+1):
        precision = min(2**i, (ell+2))
        A = A.lift_to_precision(precision)
        A = double_prec_montgomery_coefficient(A,S(j00+eps,precision),i)
        A1 = A1.lift_to_precision(precision)
        delta = [0,0]
        for lam_i in [0,1]:
            A1_lift = A1 + lam_i*eps^(2^(i-1))
            delta[lam_i] =  double_prec_codomain([A,A1_lift],phi,i)
        #note: delta_i = 0 + eps**(2**(i-1))*(unit), we cannot divide directly
        A1 = A1 - delta[0]/(delta[1]-delta[0]).valuation_zero_part()
        phi = double_prec_chain([A,A1],phi,i)
        if i == m:
            codomains = phi[-1][0]
            if codomains[0][0] == j00:
                j1_lift = codomains[1]
            elif codomains[1][0] == j00:
                j1_lift = codomains[0]
            else:
                print("codomains computed incorrectly")
    return j1_lift