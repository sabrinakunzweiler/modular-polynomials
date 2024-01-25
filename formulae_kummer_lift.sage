###########################################################################
######### Methods for lifting (2,2)-isogenies of Kummer surfaces ##########
###########################################################################


###############################
######## Generic Case #########
###############################


def double_prec_special_trafo(lifted_invariants, lift_aux, m):
    #we assume that kernel = <((x-alpha),0), x*(x-beta, 0)>,
    [A,B,C,E] = lifted_invariants
    [alpha,beta,scalar] = lift_aux
    #find a lift of alpha, so that alpha + 1/alpha = A
    #alpha^2 - A*alpha + 1 = 0 -> 2*alpha - A = =/- sqrt(A^2-4)
    term1 = alpha + alpha - A
    lifted_term1 = double_prec_squareroot(A*A-4, term1, m)
    alpha = (lifted_term1 + A)*half
    #beta^2 - B*beta + C = 0 -> 2*beta - B = +/- sqrt(B^2 - 4*C)
    term2 = beta + beta - B 
    lifted_term2 = double_prec_squareroot(B*B-4*C, term2, m)
    beta = (lifted_term2 + B)*half
    gamma = B - beta
    #assert C == beta*gamma, "lifting trafo didn't work"

    prod = alpha*(alpha-beta)
    scalar = double_prec_squareroot(prod,scalar,m)
    scalar_inv = 1/scalar
    newA = (-alpha-alpha+beta)*scalar_inv #alphanew + 1/alphanew #
    newB = (A+gamma-alpha-alpha-alpha)*scalar_inv #betanew + gammane
    newC = (A-alpha-alpha)*(gamma-alpha)*scalar_inv*scalar_inv #betanew*gammanew
    newE = E*scalar

    return [[newA,newB,newC,newE], [alpha,beta,scalar]]



def lift_kummer_isogeny(lifted_invariants,lift_aux,m):
    """
    ker4 are 4-torsion points lying above the (2,2)-isogeny that is computed
    phi is the previous isogeny 
    """
    type1_invariants = lifted_invariants
    new_type1_invariants, lift_aux = double_prec_special_trafo(type1_invariants, lift_aux, m)
    new_type2_invariants = invariant_type1_type2(new_type1_invariants)
    codomain_type1_invariants = Richelot_type2_invariants(new_type2_invariants)
    return [codomain_type1_invariants, lift_aux]


###############################
######## Gluing Case ##########
###############################


def double_prec_roots(A,roots,m):
    """
    input: roots are the roots of x^3+A*x^2+x up to precision 2^(m-1)
    output: roots in precision 2^m
    """
    new_roots = []
    for el in roots:
        if el == 0:
            new_roots.append(0)
        else: 
            term = el + el + A
            lifted_term = double_prec_squareroot(A**2-4, term, m)
            new_el = (lifted_term - A)*half
            new_roots.append(new_el)
    return new_roots


def double_prec_elliptic_transformation(new_invariants,aux_lift,m):
    """
    we assume Montgomery form.
    invariants are lifted
    roots1 and roots2 are not lifted
    """

    [roots1,roots2] = aux_lift[1]
    [A1,A2] = new_invariants
    
    #compute r1,r2,r3 Weierstrass points of E1, and s1,s2,s3 Weierstrass points of E2.
    [r1,r2,r3] = double_prec_roots(A1,roots1,m)
    [s1,s2,s3] = double_prec_roots(A2,roots2,m)

    #assert r1*r2*r3 == 0 and r1 + r2 + r3 == -A1, 'elliptic lift 1'
    #assert s1*s2*s3 == 0 and s1 + s2 + s3 == -A2, 'elliptic lift 2'
    
    
    #transformation a_r*x+b_r, e_r*y on E1
    term1 = (r2 - r3)*s1 + (r3 - r1)*s2  + (r1 - r2)*s3
    den1 = 1/(2*(r3 - r1)*(r2 - r1)*(s3 - s2))
    a_r = (term1+term1)*den1
    b_r = 1 - r1*a_r    

    alpha1 = a_r*r2 + b_r 
    beta1 = a_r*r3 + b_r 
   
    denom  = 1/(term1 * (alpha1 - 1) * (beta1 - 1))
    A = 2*(alpha1 + 1)*term1*(beta1-1)*denom
    B = 2*(beta1 + 1)*term1*(alpha1-1)*denom
    E = (alpha1-1)*(beta1-1)
  
    return [[A,B,1,E], [[r1,r2,r3],[s1,s2,s3]]]


def lift_gluing_kummer_isogeny(lifted_invariants,aux_lift,m):
    #unnecessary function, only added for consistency with generic method

    lifted_iso = double_prec_elliptic_transformation(lifted_invariants,aux_lift,m)
    
    return lifted_iso


###############################
######## Splitting Case #######
###############################


def lift_determine_codomain(lifted_invariants,aux_lift,m):

    type1_invariants = lifted_invariants
    [new_type1_invariants, aux_lift] = double_prec_special_trafo(type1_invariants, aux_lift,m)
    [A,B,C,E] = new_type1_invariants    

    return C - 1


def lift_splitting_kummer_isogeny(lifted_invariants,aux_lift,m):
    """
    input is a 4-torsion kernel lying above the next 2-isogeny, 
    and data of the previous (2,2)-isogeny
    We assume that the next (2,2)-kernel is of the form <((x-alpha),0), x*(x-beta, 0)>
    """
    type1_invariants = lifted_invariants
    [new_type1_invariants, aux_lift] = double_prec_special_trafo(type1_invariants, aux_lift,m)
    [A,B,C,E] = new_type1_invariants    

    assert C == 1, "the isogeny is not split"

    #splitting maps as in the paper (Prop. 4.2), 
    #alpha1 = (A+2)/(A-2)
    #alpha2 = (B+2)/(B-2)
    #beta1 = 1/alpha1
    #beta2 = 1/alpha2

    denom = ((B*B - 4)*(A*A - 4)*(A-B))^(-2)
    term = (A^2 - A*B + B^2 + 4)
    j1 = (256) * ((B+2)*(A+2))^2 * (term-A-A-B-B)^3 * denom
    j2 = (256) * ((B-2)*(A-2))^2 * (term+A+A+B+B)^3 * denom

    return [[j1,j2],aux_lift]