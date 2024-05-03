class instance():
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

        # storing the generated random number sequences
        self.genEcuyerCombined = None
        self.genNormalBoxMuller = None
        self.genNormalCLT = None
        self.genVanDerCorput = None
        self.genNormalInvCdf = None

        # storing the generated schemes (no quasi-MC)
        self.bsEuler1D = None
        self.bsMilstein1D = None
        self.bsEulerND = None
        self.bsMilsteinND = None

        # storing the generated schemes (quasi-MC)
        self.bsEuler1DquasiMC = None
        self.bsMilstein1DquasiMC = None
        self.bsEulerNDquasiMC = None
        self.bsMilsteinNDquasiMC = None