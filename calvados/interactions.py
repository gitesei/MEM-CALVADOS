import numpy as np
from openmm import openmm, unit

def genParamsDH(temp,ionic):
    """ Debye-Huckel parameters. """

    kT = 8.3145*temp*1e-3
    # Calculate the prefactor for the Yukawa potential
    fepsw = lambda T : 5321/T+233.76-0.9297*T+0.1417*1e-2*T*T-0.8292*1e-6*T**3
    epsw = fepsw(temp)
    lB = 1.6021766**2/(4*np.pi*8.854188*epsw)*6.02214076*1000/kT
    eps_yu = lB*kT
    # Calculate the inverse of the Debye length
    k_yu = np.sqrt(8*np.pi*lB*ionic*6.02214076/10)
    return eps_yu, k_yu

def init_bonded_interactions():
    """ Define bonded interactions. """

    # harmonic bonds
    hb = openmm.HarmonicBondForce()
    hb.setUsesPeriodicBoundaryConditions(True)

    return hb

def init_ah_interactions(eps,rc,fixed_lambda):
    """ Define Ashbaugh-Hatch interactions. """

    # intermolecular interactions
    energy_expression = f'{eps}*select(step(r-2^(1/6)*s),4*l*((s/r)^12-(s/r)^6-shift),4*((s/r)^12-(s/r)^6-l*shift)+(1-l))'
    #ah = openmm.CustomNonbondedForce(energy_expression+f'; s=0.5*(s1+s2); l=0.5*(l1+l2); shift=(0.5*(s1+s2)/{rc})^12-(0.5*(s1+s2)/{rc})^6')
    ah = openmm.CustomNonbondedForce(energy_expression+f'; l=select((id1+id2)*step(id1+id2),0.5*(l1+l2),{fixed_lambda}); shift=(s/{rc})^12-(s/{rc})^6; s=0.5*(s1+s2)')

    ah.addPerParticleParameter('s')
    ah.addPerParticleParameter('l')
    ah.addPerParticleParameter('id')

    ah.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
    ah.setCutoffDistance(rc*unit.nanometer)
    #ah.setForceGroup(0)

    print('Ashbaugh-Hatch potential between particles with lambda=1 and sigma=0.68 at',rc*unit.nanometer,end=': ')
    print(4*eps*((0.68/rc)**12-(0.68/rc)**6)*unit.kilojoules_per_mole)
    return ah

def init_wca_interactions(eps, rc):
    """ Define WCA (Weeks-Chandler-Andersen) interactions for lipid-lipid and lipid-protein pairs.
    
    WCA is a purely repulsive potential used for pairs that are not hydrophobic.
    
    Energy = 4*eps*lambda*((sigma/r)^12 - (sigma/r)^6) - shift  for r < 2^(1/6)*sigma
           = 0  for r >= 2^(1/6)*sigma
    
    Parameters
    ----------
    eps : float
        Energy scale (usually eps_lj in kJ/mol)
    rc : float
        Cutoff distance in nanometers
    
    Returns
    -------
    wca : openmm.CustomNonbondedForce
        WCA force with particles added via addParticle([sigma, lambda])
    
    Notes
    -----
    Mixing rules:
    - sigma_ij = 0.5 * (sigma_i + sigma_j)
    - lambda_ij = sqrt(lambda_i * lambda_j)
    """
    
    wca_expression = f'{eps}*select(step(r-2^(1/6)*s),0,4*l*((s/r)^12-(s/r)^6-shift))'
    
    wca = openmm.CustomNonbondedForce(
        wca_expression + f'; l=sqrt(l1*l2); shift=(1.0)^12-(1.0)^6; s=0.5*(s1+s2)'
    )
    
    wca.addPerParticleParameter('s')
    wca.addPerParticleParameter('l')
    
    wca.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
    wca.setCutoffDistance(rc*unit.nanometer)
    
    return wca

def init_sa_interactions(eps, rc):
    """ Define SA (Surface-Active) interactions for lipid-lipid and lipid-protein pairs.
    
    SA is an attractive potential for hydrophobic residues (both omega > 0).
    Provides a smooth potential well with cubic spline form.
    
    Energy = eps*lambda*select(rmin + omega - r, 3*u^2 - 2*u^3 - 1, 0)
    where u = (r - rmin) / omega, rmin = 2^(1/6)*sigma
    
    Parameters
    ----------
    eps : float
        Energy scale (usually eps_lj in kJ/mol)
    rc : float
        Cutoff distance in nanometers
    
    Returns
    -------
    sa : openmm.CustomNonbondedForce
        SA force with particles added via addParticle([sigma, lambda, omega])
    
    Notes
    -----
    Mixing rules:
    - sigma_ij = 0.5 * (sigma_i + sigma_j)
    - lambda_ij = sqrt(lambda_i * lambda_j)
    - omega_ij = 0.5 * (omega_i + omega_j)  [spatial extent of potential well]
    """
    
    sa_expression = f'{eps}*l*select(step(rmin+o-r),3*u^2-2*u^3-1,0)'
    
    sa = openmm.CustomNonbondedForce(
        sa_expression + f'; l=sqrt(l1*l2); u=(r-rmin)/o; rmin=2^(1/6)*s; s=0.5*(s1+s2); o=0.5*(o1+o2)'
    )
    
    sa.addPerParticleParameter('s')
    sa.addPerParticleParameter('l')
    sa.addPerParticleParameter('o')
    
    sa.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
    sa.setCutoffDistance(rc*unit.nanometer)
    
    return sa

def classify_lipid_interactions(lipid_indices, omegas):
    """ Classify lipid-lipid and lipid-protein interactions based on omega values.
    
    Uses the omega-based classification from iSoLF to partition pairs into three groups:
    - SA (Surface-Active): both omegas positive → attractive well
    - AH (Ashbaugh-Hatch): omega sum <= -3 → stronger attractive interaction
    - WCA (Weeks-Chandler-Andersen): default → repulsive only
    
    Parameters
    ----------
    lipid_indices : set or list
        Global particle indices of all lipid beads
    omegas : array-like
        Omega values for all particles indexed by particle index.
        Particles not in lipid_indices should have omegas[i] <= -3 (treated as AH).
    
    Returns
    -------
    sa_pairs : list of (i, j) tuples
        Pairs with both omega_i > 0 and omega_j > 0
    ah_pairs : list of (i, j) tuples
        Pairs with omega_i + omega_j <= -3 (and not both omegas > 0)
    wca_pairs : list of (i, j) tuples
        All other pairs
    
    Notes
    -----
    Each pair appears exactly once (i < j to avoid duplicates).
    Both intra-lipid and lipid-protein pairs are classified.
    """
    sa_pairs = []
    ah_pairs = []
    wca_pairs = []
    
    lipid_indices_set = set(lipid_indices)
    
    # Only consider pairs where at least one particle is a lipid
    for i in lipid_indices_set:
        for j in range(i + 1, len(omegas)):  # j > i to avoid duplicates
            omega_i = omegas[i]
            omega_j = omegas[j]
            
            # SA: both omegas positive
            if omega_i > 0 and omega_j > 0:
                sa_pairs.append((i, j))
            # AH: sum of omegas <= -3
            elif omega_i + omega_j <= -3:
                ah_pairs.append((i, j))
            # WCA: default
            else:
                wca_pairs.append((i, j))
    
    return sa_pairs, ah_pairs, wca_pairs

def init_yu_interactions(eps, k, rc):
    """ Define Yukawa interactions. """

    shift = np.exp(-k*rc)/rc
    yu = openmm.CustomNonbondedForce(f'q*{eps}*(exp(-{k}*r)/r-{shift}); q=q1*q2')
    yu.addPerParticleParameter('q')

    print('Debye-Hückel potential between unit charges at',rc*unit.nanometer,end=': ')
    print(eps*shift*unit.kilojoules_per_mole)

    yu.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
    yu.setCutoffDistance(rc*unit.nanometer)
    #yu.setForceGroup(1)

    return yu

def init_angles():
    ha = openmm.HarmonicAngleForce()
    ha.setUsesPeriodicBoundaryConditions(True)
    return ha

def init_wcafene(eps_lj):
    wcafene = init_wcafene_interactions(3*eps_lj)
    return wcafene

def init_restraints(restraint_type):
    """ Initialize restraints. """

    if restraint_type == 'harmonic':
        cs = openmm.HarmonicBondForce()
    if restraint_type == 'go':
        go_expr = 'k*(5*(s/r)^12-6*(s/r)^10)'
        cs = openmm.CustomBondForce(go_expr)#; shift=(0.5*(s)/rc)^12-(0.5*(s)/rc)^6')
        cs.addPerBondParameter('s')
        cs.addPerBondParameter('k')
    cs.setUsesPeriodicBoundaryConditions(True)
    return cs

def init_scaled_LJ(eps_lj,cutoff_lj):
    """ Initialize restraints. """

    energy_expression = 'select(step(r-2^(1/6)*s),n*4*eps*l*((s/r)^12-(s/r)^6-shift),n*4*eps*((s/r)^12-(s/r)^6-l*shift)+n*eps*(1-l))'
    scLJ = openmm.CustomBondForce(energy_expression+'; shift=(s/rc)^12-(s/rc)^6')
    scLJ.addGlobalParameter('eps',eps_lj*unit.kilojoules_per_mole)
    scLJ.addGlobalParameter('rc',float(cutoff_lj)*unit.nanometer)
    scLJ.addPerBondParameter('s')
    scLJ.addPerBondParameter('l')
    scLJ.addPerBondParameter('n')
    scLJ.setUsesPeriodicBoundaryConditions(True)
    return scLJ

def init_scaled_YU(eps_yu,k_yu):
    """ Initialize restraints. """

    shift = np.exp(-k_yu*4.0)/4.0
    scYU = openmm.CustomBondForce(f'n*q*{eps_yu}*(exp(-{k_yu}*r)/r-{shift})')
    scYU.addPerBondParameter('q')
    scYU.addPerBondParameter('n')
    scYU.setUsesPeriodicBoundaryConditions(True)
    return scYU

def init_slab_restraints(box,k):
    """ Define restraints towards box center in z direction. """

    mindim = np.amin(box)
    rcent_expr = 'k*abs(periodicdistance(x,y,z,x,y,z0))'
    rcent = openmm.CustomExternalForce(rcent_expr)
    rcent.addGlobalParameter('k',k*unit.kilojoules_per_mole/unit.nanometer)
    rcent.addGlobalParameter('z0',box[2]/2.*unit.nanometer) # center of box in z
    # rcent.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
    # rcent.setCutoffDistance(mindim/2.*unit.nanometer)
    return rcent

def add_single_restraint(
        cs, restraint_type: str,
        dij: float, k: float,
        i: int, j: int):
    """ Add single harmonic or Go restraint. """

    if restraint_type == 'harmonic':
        cs.addBond(
                i,j, dij*unit.nanometer,
                k*unit.kilojoules_per_mole/(unit.nanometer**2))
    elif restraint_type == 'go':
        cs.addBond(
                i,j, [dij*unit.nanometer,
                k*unit.kilojoules_per_mole])
    restr_pair = [i+1, j+1, dij, k] # 1-based
    return cs, restr_pair

def add_scaled_lj(scLJ, i, j, offset, comp):
    """ Add downscaled LJ interaction. """

    s = 0.5 * (comp.sigmas[i] + comp.sigmas[j])
    l = 0.5 * (comp.lambdas[i] + comp.lambdas[j])
    scLJ.addBond(i+offset,j+offset, [s*unit.nanometer, l*unit.dimensionless, comp.bondscale[i,j]*unit.dimensionless])
    scaled_pair = [i+offset+1, j+offset+1, s, l, comp.bondscale[i,j]] # 1-based
    return scLJ, scaled_pair

def add_scaled_yu(scYU, i, j, offset, comp):
    """ Add downsscaled YU interaction. """

    qij = comp.qs[i] * comp.qs[j] * unit.dimensionless
    scYU.addBond(i+offset, j+offset, [qij, comp.bondscale[i,j]*unit.dimensionless])
    scaled_pair = [i+offset+1, j+offset+1, comp.bondscale[i,j]] # 1-based
    return scYU, scaled_pair

def init_wcafene_interactions(eps):
    """ Define FENE interaction. """

    wca_expression = f'4*{eps}*select(step(r-2^(1/6)*s),0,(s/r)^12-(s/r)^6+1/4)'
    fene_expression = '+ -0.5*kfene*(rinf^2)*log(1-(r/rinf)^2); rinf=1.5*s'
    wcafene = openmm.CustomBondForce(wca_expression+fene_expression)
    wcafene.addPerBondParameter('s')
    wcafene.addPerBondParameter('kfene')
    wcafene.setUsesPeriodicBoundaryConditions(True)
    return wcafene

def init_cosine_interactions(eps):
    """ Define cosine interaction (Cooke and Deserno lipid model, DOI: https://doi.org/10.1063/1.2135785). """

    cosine_expression = f'prefactor*select(step(r-rc-1.5*s),0,select(step(r-rc),-{eps}*(cos({np.pi}*(r-rc)/(2*1.5*s)))^2,-{eps}))'
    cosine = openmm.CustomNonbondedForce(cosine_expression+'; prefactor=select(id1*id2,1-delta(l1*l2),(id1+id2)*l1*l2); rc=2^(1/6)*s; s=0.5*(s1+s2)')
    cosine.addPerParticleParameter('s')
    cosine.addPerParticleParameter('l')
    cosine.addPerParticleParameter('id')
    cosine.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
    cosine.setCutoffDistance((2**(1/6)+1.5)*unit.nanometer)
    cosine.setForceGroup(2)
    return cosine

def init_isolf_interactions(eps,rc):
    """ Define interactions between lipids (iSoLF lipid model, DOI: https://doi.org/10.1063/5.0160417).
    
    **DEPRECATED:** Use init_wca_interactions(), init_sa_interactions(), and classify_lipid_interactions()
    instead for better OpenMM performance via interaction group partitioning.
    
    This function is kept for backward compatibility but should not be used in new code.
    """
    isolf_expression = f'{eps}*select(step(r-rmin),l*select(step(o1)*step(o2),step(rmin+o-r)*(3*u^2-2*u^3-1),is_lj*4*uLJ),4*(uLJ+1/4)-select(step(o1)*step(o2),l,is_lj*l))'
    isolf = openmm.CustomNonbondedForce(isolf_expression+'; is_lj=step(-3-o1-o2); l=sqrt(l1*l2); u=(r-rmin)/o; uLJ=(s/r)^12-(s/r)^6; rmin=2^(1/6)*s; s=0.5*(s1+s2); o=0.5*(o1+o2)')
    isolf.addPerParticleParameter('s')
    isolf.addPerParticleParameter('l')
    isolf.addPerParticleParameter('o')
    isolf.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
    isolf.setCutoffDistance(rc*unit.nanometer)
    return isolf

def init_charge_nonpolar_interactions(eps,rc):
    """ Define charge-nonpolar interaction (lipid model, DOI: https://doi.org/10.1063/1.5058234 and DOI: https://doi.org/10.1073/pnas.2311700120). """

    cn = openmm.CustomNonbondedForce(f'-step(id1+id2)*{eps}*alphaq2R3/2*(1/r-1/{rc}); alphaq2R3=alpha1*q2^2*R31+alpha2*q1^2*R32')
    cn.addPerParticleParameter('R3')
    cn.addPerParticleParameter('alpha')
    cn.addPerParticleParameter('q')
    cn.addPerParticleParameter('id')
    cn.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
    cn.setCutoffDistance(rc*unit.nanometer)
    cn.setForceGroup(1)
    return cn
