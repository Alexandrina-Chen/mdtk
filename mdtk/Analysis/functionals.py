def double_exp_1(x, alpha, beta, tau1, tau2):
    return alpha *np.exp(-x / tau1) + beta * np.exp(-x / tau2)

def double_exp_2(x, alpha, beta, tau1, tau2):
    return alpha * tau1 * (1 - np.exp(-x / tau1)) + beta * tau2 * (1 - np.exp(-x / tau2))

def single_exp_1(x, tau):
    return np.exp(-x / tau)

def single_exp_2(x, tau):
    return tau * (1 - np.exp(-x / tau))