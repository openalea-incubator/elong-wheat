# -*- coding: latin-1 -*-

from __future__ import division  # use "//" to do integer division

"""
    elongwheat.parameters
    ~~~~~~~~~~~~~~~~~~~~~~

    The module :mod:`elongwheat.parameters` defines the constant parameters.

    :copyright: Copyright 2015 INRA-ECOSYS, see AUTHORS.
    :license: TODO, see LICENSE for details.

    .. seealso:: Barillot et al. 2015.
"""

"""
    Information about this versioned file:
        $LastChangedBy$
        $LastChangedDate$
        $LastChangedRevision$
        $URL$
        $Id$
"""


# -------------------------------------------------------------------------------------------------------------------
# --- SAM
# -------------------------------------------------------------------------------------------------------------------
PLASTOCHRONE = 76.1/12*24*3600    #: Leaf plastochron (s at 12�C) calculated from Ljutovac 2002 with priodia of 5.10-5 m (76 dd) ; Malvoisin 35dd associated with init 3.10-5m
max_nb_leaves = 11   #: Max number of leaves per axis
delta_TT_GA = PLASTOCHRONE * 5  #: Thermal time between floral transition of SAM and Gibberelin production expressed as a function of plastochron (s at 12�C) ; Malvoisin's data give 7 plastochrons

sowing_depth = 0.05  #: Sowing depth (m) used to define plant emergence

# Parameters for temperature responses
Temp_Tref = 12        # Arbitrary reference temperature (�C)
Temp_Ea_R = 8900      # Parameter Ea/R in Eyring equation from Johnson and Lewin (1946) - Parameter value fitted from Kemp and Blacklow (1982) (K)
Temp_DS_R = 68.432    # Parameter deltaS/R in Eyring equation from Johnson and Lewin (1946) - Parameter value fitted from Kemp and Blacklow (1982) (dimensionless)
Temp_DH_R = 20735.5   # Parameter deltaH/R in Eyring equation from Johnson and Lewin (1946) - Parameter value fitted from Kemp and Blacklow (1982) (K)
Temp_Ttransition = 9  # Below this temperature f = linear function of temperature instead of Arrhenius-like(�C)

# -------------------------------------------------------------------------------------------------------------------
# --- Leaves
# -------------------------------------------------------------------------------------------------------------------
# Exponential elongation
RERmax_dict =  { 5 : 0.000003 , 6 : 0.00000175 , 7 : 0.00000164 , 8 : 0.00000154 , 9 : 0.00000152 , 10 : 0.00000135 , 11 : 0.0000013 } # s-1 at 12�C fit dec 18
# { 5 : 0.00000279 , 6 : 0.00000176 , 7 : 0.00000162 , 8 : 0.00000144 , 9 : 0.00000144 , 10 : 0.00000144 , 11 : 0.00000142 } # Ljutovac 2002 RER (s-1 at 12�C)
# {3: 4.1E-06, 4: 4.1E-06, 5: 6E-06, 6: 4.1E-06, 7: 3.6E-06, 8: 3.2E-06, 9: 3.08E-06, 10: 2.9E-06, 11: 2.75E-06}  # MG2
# {3 : 0.00000279 ,4 : 0.00000279 ,5 : 0.00000279 , 6 : 0.00000176 , 7 : 0.00000162 , 8 : 0.00000154 , 9 : 0.0000015 , 10 : 0.00000146 , 11 : 0.00000146 }#: s-1 at 12�C fit nov 18
# {3: 4.1E-06, 4: 4.1E-06, 5: 4.1E-06, 6: 4.1E-06, 7: 3.6E-06, 8: 3.3E-06, 9: 3.2E-06, 10: 2.9E-06, 11: 2.75E-06}
# {5 : 0.009/3600, 6 : 0.009/3600, 7: 0.0088/3600, 8: 0.00875/3600, 9: 0.00875/3600, 10: 0.0086/3600, 11: 0.008/3600} # RB 2013

RERmax_dict2 =  { 5 : 0.000004 , 6 : 0.00000170 , 7 : 0.00000160 , 8 : 0.00000145 , 9 : 0.0000013 , 10 : 0.0000012 , 11 : 0.0000013 } # s-1 at 12�C fit pour pilotage par metabolisme
RERmax = 2.8E-06 #: s-1 at 12�C # 5.56e-06 Ljutovac 2002 # 4e-06 Anne # 2.43E-06 RB v1
Kc = 300 # 145.6 # 350 #: affinity coefficient of RER to C (�mol g-1)
Kn = 300 #200 #16.64  # 40  #: affinity coefficient of RER to C N (�mol g-1)

# Automate elongation
te = 300 * 3600 * 24 / 12 #: end of leaf elongation in automate growth (s at 12�c); fitted from adapted data from Fournier 2005
tm = 190 * 3600 * 24 / 12 #: time at which leaf elongation rate is maximal in automate growth (s at 12�c); fitted from adapted data from Fournier 2005
tb = -25 * 3600 * 24 / 12 #: beginning of leaf elongation in automate growth (s at 12�c); fitted from adapted data from Fournier 2005
# NB : Previous fit on adapted data from Fournier 2005 in phyllochronic time te = 271, tm=176, tb=-25
L0 = abs((1 + (te / (te - tm))) * (min(1.0, float(-tb) / float(te - tb))**((te - tb) / (te - tm))))  #: Leaf length at t=0 in automate growth (beta function) (m)
FITTED_L0 = 0.01557936             #: Fitted value of leaf length at t=0 after rescaling the beta function with L0 (m); Fournier 2005 sur courbe corrigee
OFFSET_LEAF = FITTED_L0 - L0       #: Offset used for the final fitting of the beta function (m)
SCALING_FACTOR_LEAF = 1/FITTED_L0  #: Scaling factor of the leaf in automate growth (dimensionless)
leaf_Lmax_MAX = 0.7                #: Maximum leaf_Lmax (m)

# Leaf maximal width
leaf_Wmax_dict = { 3 : 0.0040 ,4 : 0.0045 ,5 : 0.0056 , 6 : 0.0075 , 7 : 0.010 , 8 : 0.012 , 9 : 0.013 , 10 : 0.014 , 11 : 0.018 }#: m (Ljutovac 2002)
EC_wmax = 0.3  #: variation de + ou - 15% de maximal leaf width (SU)
leaf_Wmax_MAX = 0.03 #: Maximum leaf_Wmax (m)
leaf_Wmax_MIN = 0.002 #: Manimum leaf_Wmax (m)
Ksslw = 4160   #: Affinite SSLW aux fructanes (�mol C g-1)
min_SSLW = 22  #: g m-2
max_SSLW = 50  #: g m-2
ratio_LSSW_SSLW = 0.003  #: ratio lineic structural mass sheath / specific strucutal mass lamina of the specific structural dry masses (from data of J. Bertheloot, 2004) (m)
#TODO : adapter ce param�tre pour les feuilles adultes car change drastiquement
leaf_SSLW_dict = {1: 15,2: 23,3: 25,4: 24,5: 21,6: 18,7: 16,8: 18,9: 21,10: 26,11: 33} # Manip NEMA 05/06 traitments N+ (from data of J. Bertheloot, 2004) sauf pour F7/F8
leaf_SSLW_dict2 = {1: 15,2: 23,3: 25,4: 18,5: 22,6: 25,7: 20,8: 23,9: 26,10: 28,11: 31} # Test correction
leaf_SSLW_nominal = {1: 15,2: 23,3: 25,4: 23,5: 22,6: 21,7: 19,8: 19,9: 21,10: 22,11: 24}
leaf_SSLW_MIN = 10
leaf_SSLW_MAX = 50
leaf_LSSW_dict = {1: 0.08,2: 0.09,3: 0.11,4: 0.18,5: 0.17,6: 0.21,7: 0.24,8: 0.4,9: 0.5,10: 0.55,11: 0.65}  # Manip NEMA 05/06 Soissons N+ (from data of J. Bertheloot, 2004)
leaf_LSSW_nominal = {1: 0.09, 2: 0.088, 3: 0.11, 4: 0.19, 5: 0.17, 6: 0.21, 7: 0.23, 8: 0.36, 9: 0.45, 10: 0.5, 11: 0.58}

# Share of division zone in the growing leaf.  The model was fitted on litterature data on wheat(Fournier 2005, Beemer and Masle 1996, Schuppler 1998)
ratio_DZ_l_init = 0.065     #: normalized length of the whole leaf
ratio_DZ_l_mid = 0.075      #: normalized length of the whole leaf
ratio_DZ_l_end = 0.9        #: normalized length of the whole leaf

# -------------------------------------------------------------------------------------------------------------------
# --- Internodes
# -------------------------------------------------------------------------------------------------------------------
# Exponential elongation
RERmax_dict_IN = { 3 : 2.48E-06 ,4 : 2.48E-06 ,5 : 2.48E-06 , 6 : 2.48E-06 , 7 : 3.1E-06 , 8 : 1.9E-06 , 9 : 2E-06 , 10 :2.05E-06 , 11 :1.92E-06  }#: s-1 at 12�C FIT dec 18
# { 3 : 2.48E-06 ,4 : 2.48E-06 ,5 : 2.48E-06 , 6 : 2.48E-06 , 7 : 2.9E-06 , 8 : 2.3E-06 , 9 : 2.45E-06 , 10 :2.35E-06 , 11 : 2.2E-06 }#: s-1 at 12�C FIT nov 18
#{ 3 : 2.48E-06 ,4 : 2.48E-06 ,5 : 2.48E-06 , 6 : 2.48E-06 , 7 : 2.48E-06 , 8 : 2.48E-06 , 9 : 2.48E-06 , 10 : 1.9E-06 , 11 : 1.6E-06 }#: s-1 at 12�C
# estimate from Ljutovac 2002 over the period until leaf ligulation i.e. wider than in the model.
# Because i) not enought data if we consider only up to previous leaf ligulation, ii) same exponential like period

SCALING_FACTOR_INT = 59  # 5.2 #: Scaling factor of the internode in automate growth (dimensionless), Malvoisin 1984 II

# Initiation of internode elongation
nb_PLASTO_internode_init = 5  #: Delay between leaf initiation and internode initiation expressed as a number of plastochron. From Malvoisin 1984b, associated with priodia of 5.10-4 m
internode_L_init = 5E-5       #: Initial internode length (m)

# Automate elongation
te_IN = 210 * 3600 * 24 /12  #: end of internode elongation in automate growth; Malvoisin 1984 II
tm_IN = 156 * 3600 * 24 /12  #: time at which internode elongation rate is maximal in automate growth (s); Malvoisin 1984 II
tb_IN = -70 * 3600 * 24 /12  #: beginning of internode elongation in automate growth (s); Malvoisin 1984 II
L0_INT = (1 + (te_IN / (te_IN - tm_IN))) * (min(1.0, float(-tb_IN) / float(te_IN - tb_IN))**((te_IN - tb_IN) / (te_IN - tm_IN)))  #: Internode length at t=0 in automate growth (beta function) (m)
OFFSET_INT = 1 / SCALING_FACTOR_INT - L0_INT

ratio_LSIW_LSSW = 2.5 #: ratio lineic structural internode mass / lineic structural sheath mass  of the specific structural dry masses (from data of J. Bertheloot, 2004)
internode_LSIW_dict = {1:2.8, 2:2.8, 3:2.8, 4:2.8, 5:2.8, 6:2.8, 7:2.8, 8: 2.8, 9: 2.3,10: 1.7,11: 1.6,12: 1.4,13: 0.7} #: Manip Marion GAuthier 2017/18 en accord avec Manip Romain Barillot 2014


class HiddenZoneInit(object):
    """
    Initial values for hidden zones
    """
    def __init__(self):
        self.leaf_is_growing = True
        self.internode_is_growing = False
        self.leaf_pseudostem_length = 4E-5           #: m TODO: ok comme valeur?
        self.delta_leaf_pseudostem_length = 0        #: m
        self.internode_distance_to_emerge = 0        #: m
        self.delta_internode_distance_to_emerge = 0  #: m, needded for growthwheat
        self.leaf_L = 5E-5                       #: m, should be consistent with PLASTOCHRONE
        self.delta_leaf_L = 0                    #: m, needded for growthwheat
        self.internode_L = 0                     #: m
        self.delta_internode_L = 0               #: m, needded for growthwheat
        self.leaf_Lmax = None                    #: m, no calculation before emergence Ln-1
        self.lamina_Lmax = None                  #: m, no calculation before emergence Ln-1
        self.sheath_Lmax = None                  #: m, no calculation before emergence Ln-1
        self.leaf_Wmax_int = 2e-7                #: m, intermediate maximum leaf width
        self.leaf_Wmax = None                    #: m, no calculation before emergence Ln-1
        self.SSLW = None                         #: g m-2, no calculation before emergence Ln-1
        self.LSSW = None                         #: g m-1, no calculation before emergence Ln-1 (about 2)
        self.leaf_is_emerged = False
        self.internode_Lmax = None               #: m, no calculation before ligulation Ln
        self.LSIW = None                         #: g m-1, no calculation before ligulation Ln
        self.internode_is_visible = False
        self.leaf_pseudo_age = 0
        self.internode_pseudo_age = 0
        self.delta_leaf_pseudo_age = 0
        self.delta_internode_pseudo_age = 0
        self.hiddenzone_age = 0
        self.is_over = False
        self.ratio_DZ = 1.0
        self.ratio_EOZ = 0

        # Default values used for RER calculation in elong wheat
        self.sucrose = 5E-6                      #: �mol C
        self.amino_acids = 4E-6                  #: �mol N
        self.fructan = 0 #1.7e-3                    #: �mol C - about 10% DM
        self.leaf_enclosed_mstruct = 1.26E-07    #: g
        self.internode_enclosed_mstruct = 0      #: g
        self.mstruct = self.leaf_enclosed_mstruct + self.internode_enclosed_mstruct  #: g
        self.leaf_enclosed_Nstruct = self.leaf_enclosed_mstruct * 0.005  #: g, parameter value in growth wheat #: g
        self.internode_enclosed_Nstruct = self.internode_enclosed_mstruct * 0.0322   #: g, parameter value in growth wheat
        self.Nstruct = self.leaf_enclosed_Nstruct + self.internode_enclosed_Nstruct  #: g
        self.proteins = 2.6E-03                   #: �mol N - about 9% N
        self.conc_cytokinins = 150                #: AU / g mstruct
        self.integral_conc_sucrose = self.sucrose / self.mstruct #: �mol C / g mstruct


class ElementInit(object):
    """
    Initial values for emerged and growing elements
    """
    def __init__(self):
        self.is_growing = True
        self.length = 0               #: m
        self.senesced_length = 0      #: m
        self.age = 0                  #: Thermal Time
        self.max_proteins = 0         #: �mol N
        self.Nresidual = 0            #: g
        self.sucrose = 0              #: �mol C
        self.amino_acids = 0          #: �mol N
        self.fructan = 0              #: �mol C
        self.proteins = 0             #: �mol N
        self.mstruct = 0              #: g
        self.Nstruct = 0              #: g
        self.cytokinins = 0           #: g
