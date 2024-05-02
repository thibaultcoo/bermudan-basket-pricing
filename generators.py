from math import log, exp, sqrt, pi, cos, sin


class RandomGenerator:
    def __init__(self):
        self.current = 0

    def get_current(self):
        return self.current

    def generate(self):
        raise NotImplementedError("Generate method should be implemented in subclasses")

    def mean(self, nb_sim: int) -> float:
        simulations = [self.generate() for _ in range(nb_sim)]
        mean_value = sum(simulations) / nb_sim
        return mean_value

    def variance(self, nb_sim: int) -> float:
        simulations = [self.generate() for _ in range(nb_sim)]
        mean_value = sum(simulations) / nb_sim
        variance_value = sum([x - mean_value ** 2 for x in simulations]) / nb_sim
        return variance_value


class DiscreteGenerator(RandomGenerator):
    def __init__(self):
        super().__init__()


class HeadTail(DiscreteGenerator):
    def __init__(self):
        super().__init__()
        self.generator = EcuyerCombined()

    def generate(self) -> int:
        return 1 if self.generator.generate() <= 0.5 else 0


class Bernoulli(DiscreteGenerator):
    def __init__(self, p: float):
        super().__init__()
        self.p = p
        self.generator = EcuyerCombined()

    def generate(self) -> int:
        return 1 if self.generator.generate() <= self.p else 0


class Binomial(DiscreteGenerator):
    def __init__(self, p: float, n: int):
        super().__init__()
        self.p = p
        self.n = n
        self.generator = Bernoulli(self.p)

    def generate(self) -> int:
        return sum([self.generator.generate() for _ in range(self.n)])


class UniformGenerator(RandomGenerator):
    def __init__(self):
        super().__init__()


class PseudoGenerator(UniformGenerator):
    def __init__(self, seed=203):
        super().__init__()
        self.seed = seed
        self.current = seed


class LinearCongruential(PseudoGenerator):
    def __init__(self, multiplier=17, increment=43, modulus=100, seed=27):
        super().__init__(seed)
        self.multiplier = multiplier
        self.increment = increment
        self.modulus = modulus

    def generate(self):
        # Generate the next number in the sequence
        self.current = (self.multiplier * self.current + self.increment) % self.modulus
        return self.current / self.modulus  # Return a normalized number between 0 and 1


class EcuyerCombined(PseudoGenerator):
    def __init__(self, seed1=202, seed2=203):
        super().__init__()
        # Initialize two linear congruential generators
        self.first_linear = LinearCongruential(40014, 0, 2147483563, seed1)
        self.second_linear = LinearCongruential(40692, 0, 2147483399, seed2)

    def generate(self):
        # Generate a number from each Linear Congruential generator
        self.first_linear.generate()
        self.second_linear.generate()

        x1 = self.first_linear.get_current()
        x2 = self.second_linear.get_current()

        # Combining the two generators
        m1 = 2147483563
        current = (x1 - x2) % (m1 - 1)

        if current > 0:
            result = current / m1
        elif current < 0:
            result = current / m1 + 1
        else:
            result = (m1 - 1) / m1

        self.current = result
        return result


class ContinuousGenerator(RandomGenerator):
    def __init__(self):
        super().__init__()


class Exponential(ContinuousGenerator):
    def __init__(self, algorithm: str, lambda_: float):
        super().__init__()
        self.algorithm = algorithm
        if lambda_ <= 0:
            raise ValueError("Lambda must be strictly positive")
        self.lambda_ = lambda_
        self.uniform1 = EcuyerCombined()
        self.uniform2 = EcuyerCombined(27, 42)

    def generate(self):
        if self.algorithm == "inverse_distribution":
            return -log(self.uniform1.generate()) / self.lambda_

        elif self.algorithm == "rejection_sampling":
            M = self.lambda_
            while True:
                X = self.uniform1.generate()
                Y = M * self.uniform2.generate()
                if Y <= self.lambda_ * exp(-self.lambda_ * X):
                    return X

        else:
            raise ValueError(f"{self.algorithm} has not yet been implemented")


class Normal(ContinuousGenerator):
    def __init__(self, algorithm: str, mu: float, sigma: float, require_new_simulation: bool):
        super().__init__()
        if sigma <= 0:
            raise ValueError("The variance must be strictly positive")
        self.algorithm = algorithm
        self.mu = mu
        self.sigma = sigma
        self.uniform1 = EcuyerCombined()
        self.uniform2 = EcuyerCombined(27, 42)
        self.require_new_simulation = True
        self.X = None
        self.Y = None

    def generate(self):
        if self.algorithm == "box_muller":
            if self.require_new_simulation:
                R = sqrt(-2 * log(self.uniform1.generate()))
                Theta = 2 * pi * self.uniform2.generate()

                self.X = R * cos(Theta)
                self.Y = R * sin(Theta)

                self.require_new_simulation = False
                return self.mu + self.sigma * self.X
            else:
                self.require_new_simulation = True
                return self.mu + self.sigma * self.Y

        elif self.algorithm == "clt":
            standard_normal = sum([self.uniform1.generate() for _ in range(12)]) - 6
            return self.mu + self.sigma * standard_normal

        elif self.algorithm == "rejection_sampling":
            raise NotImplementedError("Rejection Sampling Algorithm has not yet been implemented")

        else:
            raise ValueError(f"{self.algorithm} has not yet been implemented")

