import numpy as np

"""
Overall chain of thought for pricing pipeline:
1/ define set of parameters for both the option and the algorithm
2/ choose methodology (what option + what random sequence + what process Euler/Millstein + what enhancements)
3/ prepare the successively introduced steps of the methodology
4/ run the final pipeline and present the results
"""


class algorithm():
    def __init__(self) -> None:
        # general option parameters
        self.matu = 2/12
        self.strike = 50
        self.rate = 0.05
        
        # simulation parameters
        self.nb_simu = 10e6
        self.method = None
        self.random_sequence = None
        self.resolution_scheme = None

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

        # pricing results
        self.option = None
        self.price = None
        self.variance = None
        self.interval = None

        # enhancements (variance reduction methods, quasi Monte-Carlo)
        self.control_variate = None
        self.antithetic_variable = None
        self.quasi_MC = None


    def displayPricingParams(self):
        """
        Outputs the option and parameters currently ran for pricing.
        """

        print("----- Pricing a " + self.method + " option -----")


    def displayPricingResults(self):
        """
        Ouputs the option and algorithm pricing results.    
        """

        print("----- Pricing results -----")
        print("Option price : " + str(self.price))
        print("Price variance : " + str(self.variance))
        print("Price confidence interval : " + str(self.interval))


    def generatingSequences(self):
        """
        Prepares the different random number generation algorithms.
        """

        # ecuyerCombined = ecuyerCombined()
        # normalBoxMuller = normalBoxMuller()
        # normalCLT = normalCLT()

        pass


    def generatingTrajectories(self):
        """
        Builds the required diffusion algorithms. Generalized to N-dimensions.
        """

        # blackScholesEuler = blackScholesEuler()
        # blackScholesMilstein = blackScholesMilstein()

        pass

    
    def computePrice(self, option=None):
        pass


    def computeVariance(self, option=None):
        pass


    def computeIC(self, option=None):
        pass


    def priceEuropeanSingle(self):
        """
        Encompasses the overall pricing methodology for European Single call options.
        """

        self.method = "European Single Call"
        self.displayPricingParams()
        self.generatingSequences()
        self.generatingTrajectories()

        self.price = self.computePrice(self.option)
        self.variance = self.computeVariance(self.option)
        self.interval = self.computeIC(self.option)

        self.displayPricingResults()

        return None


    def priceEuropeanBasket(self):
        """
        Encompasses the overall pricing methodology for European Basket call options.
        """

        self.method = "European Basket Call"
        self.displayPricingParams()


    def priceBermudanSingle(self):
        """
        Encompasses the overall pricing methodology for Bermudan Single call options.
        """

        self.method = "Bermudan Single Call"
        self.displayPricingParams()


    def priceBermudanBasket(self):
        """
        Encompasses the overall pricing methodology for Bermudan Baket call options.
        """

        self.method = "Bermudan Basket Call"
        self.displayPricingParams()


# runs the global algorithm structure
algorithm().priceEuropeanSingle()