import numpy as np
import schemes

"""
Need to build a framework for a trajectory simulation, considering or not antithetic variable feature.

"""

class EuropeanSingleCall():
    def __init__(self, scheme=None, strike=None, rate=None, matu=None) -> None:
        self.scheme = scheme
        self.strike = strike
        self.rate = rate
        self.matu = matu
        self.v = None


    def price_raw(self, nb_simu=None):
        """
        Computes the price of the European single call option with no further variance enhancement.
        """

        sum = 0.0
        self.v = [0] * nb_simu

        for n in range(nb_simu):
            self.scheme.simulate(0, self.matu, int(self.matu * 365))
            final = max(self.process.get_value(self.matu) - self.strike, 0)

            sum += final
            self.v[n] = final

        price = np.exp(-self.rate * self.matu) * (sum / nb_simu)
        return price


    def price_antithetic(self, nb_simu=None):
        """
        Computes the price of the European single call option considering an antithetic variable.
        """

        sum = 0.0
        self.v = [0] * nb_simu

        for n in range(nb_simu):
            self.process.simulate_antithetic(0, self.matu, int(self.matu * 365))
            final1 = max(self.process.get_value(self.matu) - self.strike, 0)
            final2 = max(self.process.get_value_antithetic(self.matu) - self.strike, 0)

            sum += (final1 + final2) 
            self.v[n] = (final1 + final2) / 2

        price = np.exp(-self.rate * self.matu) * (sum / (nb_simu * 2))
        return price


    def price_control_variate(self, nb_simu=None):
        """
        Computes the price of the European single call option considering a control variate.
        """
        
        sum = 0.0
        self.v = [0] * nb_simu

        for n in range(nb_simu):
            self.process.simulate(0, self.matu, int(self.matu * 365))
            final = max(self.strike - self.process.get_value(self.matu), 0) + np.exp(self.rate * self.matu) * self.process.get_value(0) - self.strike

            sum += final
            self.v[n] = final

        price = np.exp(-self.rate * self.matu) * (sum / nb_simu)
        return price
    

    def price_quasiMC(self, nb_simu=None):
        pass