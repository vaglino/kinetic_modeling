using PyCall
BayesOpt = pyimport("bayes_opt")
BayesOpt.BayesianOptimization
so = pyimport("scipy.optimize")

function black_box_function(x)
    # Function with unknown internals we wish to maximize.
    return x[1] ^ 2 + x[2] ^2
end

so.minimize(black_box_function,[1.,1.],method="Nelder-Mead")

bds = Dict('x' => (2, 4), 'y' => (-3, 3))
bds = Dict("x"=> (2,2))
bds = Dict(-2,2)
rates = []
Dict(bds)
Dict(i => bd for (i,bd) in enumerate(bds))


optimizer = BayesOpt.BayesianOptimization(
    f=black_box_function,
    pbounds=bds,
    verbose=2, # verbose = 1 prints only when a maximum is observed, verbose = 0 is silent
    random_state=1)

py"""
    from bayes_opt import BayesianOptimization


    def black_box_function(x, y):
        # Function with unknown internals we wish to maximize.
        return -x ** 2 - (y - 1) ** 2 + 1

    # Bounded region of parameter space
    pbounds = {'x': (2, 4), 'y': (-3, 3)}

    optimizer = BayesianOptimization(
        f=black_box_function,
        pbounds=pbounds,
        verbose=2, # verbose = 1 prints only when a maximum is observed, verbose = 0 is silent
        random_state=1,
    )

    optimizer.maximize(
        init_points=2,
        n_iter=3,
    )
"""
optimizer.maximize(
    init_points=2,
    n_iter=3)

py"""
import numpy as np

def sinpi(x):
    return np.sin(np.pi * x)
"""


#----------------------------------------------------------------------------

skopt = pyimport("skopt")
noise_level = 0.1

function f(x,ct)
    x[1]^2 + x[2]^2 + ct
end

g = (x) -> f(x,ct)
ct = 1.0
# Plot f(x) + contours
x = range(-2, stop=2, length=400) #.reshape(-1, 1)

fx = [f(x_i) for x_i in x]
plot(x, fx, label="True (unknown)")
plt.fill(np.concatenate([x, x[::-1]]),
        np.concatenate(([fx_i - 1.9600 * noise_level for fx_i in fx],
                        [fx_i + 1.9600 * noise_level for fx_i in fx[::-1]])),
        alpha=.2, fc="r", ec="None")
plt.legend()
plt.grid()
plt.show()


res = skopt.gp_minimize(g,                  # the function to minimize
                  [(-2.0, 2.0),(-2.0, 2.0)],      # the bounds on each dimension of x
                  acq_func="EI",      # the acquisition function
                  n_calls=50,         # the number of evaluations of f
                  n_random_starts=5,  # the number of random initialization points
                  noise=0.,       # the noise level (optional)
                  random_state=1234)   # the random seed

res["x"]
