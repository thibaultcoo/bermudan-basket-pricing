import numpy as np
import instance  # is a general instance of our algorithm: stores the wanted parameters as well as the generated sequences and schemes.
import result  # is an object that will store all our prices and other results for a given structure and its variance reduction enhancements.
import generators  # stores all random number sequence generators
import schemes  # stores all the resolution schemes available
import europeanSingle  # is an object for the specific single asset european call option (contains initialization and three pricing methods)


"""
Vanilla and exotic derivatives pricing tool: M203 Numerical finance project
"""

"""
Inner logic behind the algorithm:
1/ We instantiate an object 'instance' that will encompasses all the algorithm features and sequences. It already has params, and all the generated sequences are filled into this object.
2/ Once we have everything stored within this instance, we fetch what we need from it to price and evaluate all the different structures that we come across.
3/ Those pricing results will be stored in a second large object 'result' which we use to centralize all our pricing results that we then display for the user to see.
"""

print("--- Beginning: instantiation and sequence generation ---")

# Creating an instance of the instance class to store the algorithm object
o_instance = instance.instance()

# Storing the generated sequences and schemes within the general instance object (independent of the option and structure so we run it beforehand, getting them from 'generators.py')
o_instance.genEcuyerCombined = generators.EcuyerCombined
o_instance.genNormalBoxMuller = None
o_instance.genNormalCLT = None
o_instance.genVanDerCorput = None
o_instance.genNormalInvCdf = None

# Storing the generated schemes for regular MC (getting them from 'schemes.py', will be function of the 'o_instance' object)
o_instance.bsEuler1D = None
o_instance.bsMilstein1D = None
o_instance.bsEulerND = None
o_instance.bsMilsteinND = None

# Storing the generated schemes for quasi MC (getting them from 'schemes.py', will be function of the 'o_instance' object)
o_instance.bsEuler1DquasiMC = None
o_instance.bsMilstein1DquasiMC = None
o_instance.bsEulerNDquasiMC = None
o_instance.bsMilsteinNDquasiMC = None

print("--- Ending: instantiation and sequence generation ---")

"""
Now that we instantiated a general object containing all the params and required sequences used to price, we successively run the pricing of our desired options.
"""

print("--- Beginning: pricing single European call option ---")

o_result_europeanSingle = result.result()

"""
(1) We price here all the Euler variations.
"""

o_option_europeanSingleEuler = europeanSingle.EuropeanSingleCall(scheme=o_instance.bsEuler1D, strike=o_instance.strike, rate=o_instance.rate, matu=o_instance.matu)
o_option_europeanSingleEulerQuasiMC = europeanSingle.EuropeanSingleCall(scheme=o_instance.bsEuler1DquasiMC, strike=o_instance.strike, rate=o_instance.rate, matu=o_instance.matu)

# Pricing and directly storing the result into the result object (dealing with the variance and confidence interval later)
o_result_europeanSingle.rawEulerPrice = o_option_europeanSingleEuler.price_raw(nb_simu=o_instance.nb_simu)
o_result_europeanSingle.antitheticEulerPrice = o_option_europeanSingleEuler.price_antithetic(nb_simu=o_instance.nb_simu)
o_result_europeanSingle.controlEulerPrice = o_option_europeanSingleEuler.price_control_variate(nb_simu=o_instance.nb_simu)
o_result_europeanSingle.quasiMCEulerPrice = o_option_europeanSingleEulerQuasiMC.price_quasiMC(nb_simu=o_instance.nb_simu)

"""
(2) We price here all the Milstein variations.
"""

o_option_europeanSingleMilstein = europeanSingle.EuropeanSingleCall(scheme=o_instance.bsMilstein1D, strike=o_instance.strike, rate=o_instance.rate, matu=o_instance.matu)
o_option_europeanSingleMilsteinQuasiMC = europeanSingle.EuropeanSingleCall(scheme=o_instance.bsMilstein1DquasiMC, strike=o_instance.strike, rate=o_instance.rate, matu=o_instance.matu)

# Pricing and directly storing the result into the result object (dealing with the variance and confidence interval later)
o_result_europeanSingle.rawMilsteinPrice = o_option_europeanSingleMilstein.price_raw(nb_simu=o_instance.nb_simu)
o_result_europeanSingle.antitheticMilsteinPrice = o_option_europeanSingleMilstein.price_antithetic(nb_simu=o_instance.nb_simu)
o_result_europeanSingle.controlMilsteinPrice = o_option_europeanSingleMilstein.price_control_variate(nb_simu=o_instance.nb_simu)
o_result_europeanSingle.quasiMCMilsteinPrice = o_option_europeanSingleMilstein.price_quasiMC(nb_simu=o_instance.nb_simu)

print("--- Ending: pricing single European call option ---")

"""
Then: we will simply need to build an efficient display function to exhibit all results containted within each solved option object.
"""