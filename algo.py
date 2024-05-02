import numpy as np

"""
Overall chain of thought for pricing pipeline:
1/ define set of parameters for both the option and the algorithm
2/ choose methodology (what option + what random sequence + what process Euler/Millstein + what enhancements)
3/ prepare the successively introduced steps of the methodology
4/ run the final pipeline and present the results
"""

"""
Bermudan basket option pricer: M203 Numerical finance project
"""

class algorithm():
    def __init__(self) -> None:
        # general option parameters
        self.matu = 2/12
        self.strike = 50
        self.rate = 0.05
        
        # simulation parameters
        self.nb_simu = 10e6

        # single underlying options parameters
        self.spot_single = 50
        self.vol_single = 0.25

        # basket underlying options parameters
        self.spot_basket = [40, 50, 60]
        self.vol_basket = [[0.2, 0.05, 0.08], [0.05, 0.23, 0.12], [0.08, 0.12, 0.31]]
        self.weights = [0.2, 0.5, 0.3]
        self.size_basket = len(self.weights)

        # bermudan feature schedule
        self.bermudan = [1/24, 1/12, 3/24, 2/12]

        # initializing random number sequences
        self.seqEcuyerCombined = None
        self.seqNormalBoxMuller = None
        self.seqNormalCLT = None

        # initializing both schemes
        self.euler = None
        self.milstein = None


    def buildEuropeanSinglePriceables(self):
        """
        Encompasses the overall pricing methodology for European Single call options.
        """

        pass


    def displayPricingParams(self):
        """
        Outputs the option and parameters currently ran for pricing.
        """

        pass


    def displayPricingResults(self):
        """
        Ouputs the option and algorithm pricing results.    
        """

        pass



# runs the global algorithm structure
algorithm().buildBermudanBasketPriceables()