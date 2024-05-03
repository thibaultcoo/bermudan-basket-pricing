class result():
    def __init__(self) -> None:
        # Euler without no variance reduction pricing result
        self.rawEulerPrice = None
        self.rawEulerVar = None
        self.rawEulerIntervalLeft = None
        self.rawEulerIntervalRight = None

        # Milstein without no variance reduction pricing result
        self.rawMilsteinPrice = None
        self.rawMilsteinVar = None
        self.rawMilsteinIntervalLeft = None
        self.rawMilsteinIntervalRight = None

        # Euler with antithetic variable
        self.antitheticEulerPrice = None
        self.antitheticEulerVar = None
        self.antitheticEulerIntervalLeft = None
        self.antitheticEulerIntervalRight = None

        # Milstein wiht antithetic variable
        self.antitheticMilsteinPrice = None
        self.antitheticMilsteinVar = None
        self.antitheticMilsteinIntervalLeft = None
        self.antitheticMilsteinIntervalRight = None

        # Euler with control variate
        self.controlEulerPrice = None
        self.controlEulerVar = None
        self.controlEulerIntervalLeft = None
        self.controlEulerIntervalRight = None

        # Milstein with control variate
        self.controlMilsteinPrice = None
        self.controlMilsteinVar = None
        self.controlMilsteinIntervalLeft = None
        self.controlMilsteinIntervalRight = None

        # Euler with quasi-MC
        self.quasiMCEulerPrice = None
        self.quasiMCEulerVar = None
        self.quasiMCEulerIntervalLeft = None
        self.quasiMCEulerIntervalRight = None

        # Milstein with quasi-MC
        self.quasiMCMilsteinPrice = None
        self.quasiMCMilsteinVar = None
        self.quasiMCMilsteinIntervalLeft = None
        self.quasiMCMilsteinIntervalRight = None