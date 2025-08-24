import cvxpy as cp
import numpy as np
from scipy.special import expit

class PairwiseComparison:
    def __init__(self, W):
        assert W.shape[0] == W.shape[1]
        self.W = W
        self.n = W.shape[0]
        self.theta_opt = None
        self.obj_value = None
        self.topk_vector = None  # Stores binary vector for top/bottom-k

    def fit(self, gamma, noise):
        n = self.n
        W = self.W
        theta = cp.Variable(n)

        theta_i = cp.reshape(theta, (n, 1))
        theta_j = cp.reshape(theta, (1, n))
        M = theta_i - theta_j

        log_sigmoid = -cp.logistic(-M)
        log_likelihood = cp.sum(cp.multiply(W, log_sigmoid))

        l2_penalty  = 0.5 * gamma * cp.sum_squares(theta)
        linear_term = noise @ theta

        objective   = cp.Minimize(-log_likelihood + l2_penalty + linear_term)
        constraints = [cp.sum(theta) == 0]

        problem = cp.Problem(objective, constraints)
        problem.solve()

        self.theta_opt  = theta.value
        self.obj_value  = problem.value

    def fit_by_count(self, noise, k, select="top"):
        """
        Selects k entries based on the noisy row-wise sum.
        By default, selects the bottom-k. Set select='top' to select top-k.
        """
        row_sum = np.sum(self.W, axis=1)
        noisy_row_sum = row_sum + noise

        if select == "top":
            indices = np.argsort(noisy_row_sum)[-k:]
        elif select == "bottom":
            indices = np.argsort(noisy_row_sum)[:k]
        else:
            raise ValueError("select must be 'top' or 'bottom'")

        # Create binary vector
        topk_vector = np.zeros(self.n, dtype=int)
        topk_vector[indices] = 1

        # Store for access
        self.topk_vector = topk_vector

        return topk_vector
    
    def get_noisy_counts(self, noise):
        """
        Computes the row‐wise sum of W, adds the provided noise, 
        and returns those noisy counts.
        """
        # true row counts
        row_sum = np.sum(self.W, axis=1)
        # add noise
        noisy_row_sum = row_sum + noise
        # store for later inspection if you like
        self.noisy_row_sum = noisy_row_sum
        return noisy_row_sum

    def get_theta(self):
        return self.theta_opt

    def get_objective_value(self):
        return self.obj_value

    def get_topk_vector(self):
        return self.topk_vector
    

class individualDP(PairwiseComparison):
    def __init__(self, n, k, L, theta_true=None, W=None):
        self.n = n
        self.k = k
        self.L = L
        if theta_true:
            self.theta_true = np.asarray(theta_true).flatten()

        if W is None:
            self.samples = self._generate_samples()
            W = self._aggregate_samples()
        else:
            self.samples = None

        super().__init__(W)

    def _generate_samples(self):
        sample_list = []
        pairs = [(i, j) for i in range(self.n) for j in range(i + 1, self.n)]

        for _ in range(self.k):
            mat = np.zeros((self.n, self.n), dtype=int)
            sampled_indices = np.random.choice(len(pairs), size=self.L, replace=True)
            for idx in sampled_indices:
                i, j = pairs[idx]
                delta = self.theta_true[i] - self.theta_true[j]
                prob = expit(delta)
                if np.random.rand() < prob:
                    mat[i, j] += 1
                else:
                    mat[j, i] += 1
            sample_list.append(mat)

        return sample_list

    def _aggregate_samples(self):
        return sum(self.samples)
    
relative_l2 = lambda pred, true: np.log(np.linalg.norm(pred - true) / np.linalg.norm(true))
relative_l_inf = lambda pred, true: np.log(np.linalg.norm(pred - true, ord=np.inf) / np.linalg.norm(true, ord=np.inf))
hit_at_k = lambda pred, true, k: 1 - (np.intersect1d(np.argpartition(pred, k)[-k:], np.argpartition(true, k)[-k:]).size / k)
hit_at_k_by_count = lambda pred_vec, true_vec, k: 1 - (np.intersect1d(np.where(pred_vec == 1)[0], np.where(true_vec == 1)[0]).size / k)

rank = lambda x: np.argsort(np.argsort(x))

avg_rank_diff = lambda pred, true: np.mean(
    np.abs(rank(pred) - rank(true))
)