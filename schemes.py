import numpy as np


class BS():
    pass


class BSEuler():
    def __init__(self, gen, s, r, VarCovar, dim):
        self.gen = gen
        self.s = s
        self.r = r
        self.VarCovar = VarCovar
        self.dim = dim

        # Compute the Cholesky decomposition of the covariance matrix
        if np.linalg.det(VarCovar) == 0:
            eig_vals, eig_vecs = np.linalg.eig(VarCovar)
            D = np.diag(np.maximum(eig_vals, 0))
            P = eig_vecs
            self.B = P @ np.sqrt(D)
        else:
            self.B = np.linalg.cholesky(VarCovar)

    def simulate(self, start_time, end_time, nb_steps):
        dt = (end_time - start_time) / nb_steps
        last = np.array(self.s)

        paths = []

        for _ in range(self.dim):
            paths.append([last[_]])

        for _ in range(nb_steps):
            dW = self.gen.generate_vector(self.dim) * np.sqrt(dt)
            Z = self.r * dt + np.dot(self.B, dW)
            next_ = last + np.multiply(last, Z)
            last = next_

            for j in range(self.dim):
                paths[j].append(last[j])

        return paths


class BSMilstein():
    pass