import math
from random import gauss

class CancerGrowthModel:
    def __init__(self, model, ys_data, ts_data):
        self.model = model
        self.ys_data = ys_data
        self.ts_data = ts_data

    def __repr__(self):
        return (
            f"CancerGrowthModel(model={self.model}, "
            f"len ys_data={len(self.ys_data)}, "
            f"len ts_data={len(self.ts_data)})"
        )

    def __str__(self):
        return (
            "Cancer Growth Model\n"
            f"Chosen Model = {self.model}\n"

        )

    # make sure numbers are never 0, to prevent errors for log() or when dividing
    def positive(self, x, eps=1e-10):
        return max(x, eps)

    def linear_model(self, t, y, c):
        c = self.positive(c)
        return c

    def exponential_model(self, t, y, c):
        y = self.positive(y)
        c = self.positive(c)
        return c * y

    def mendelsohn_model(self, t, y, c, a):
        y = self.positive(y)
        c = self.positive(c)
        a = self.positive(a)
        return c * y ** a

    def exponential_flattening_model(self, t, y, c, ymax):
        y = self.positive(y)
        c = self.positive(c)
        ymax = self.positive(ymax)
        return c * (ymax - y)

    def logistic_model(self, t, y, c, ymax):
        y = self.positive(y)
        c = self.positive(c)
        ymax = self.positive(ymax)
        return c * y * (1 - y / ymax)

    def montroll_model(self, t, y, c, ymax, a):
        y = self.positive(y)
        c = self.positive(c)
        ymax = self.positive(ymax)
        a = self.positive(a)

        return c * y * (1- y**a/ymax**a)

    def allee_model(self, t, y, c, ymin, ymax):
        y = self.positive(y)
        ymax = self.positive(ymax)
        ymin = self.positive(ymin)
        c = self.positive(c)
        return c * (1 - ymin/y) * (1-y/ymax)

    def linear_limited_model(self, t, y, c, a):
        y = self.positive(y)
        c = self.positive(c)
        a = self.positive(a)
        return c * y/(y+a)

    def surface_limited_model(self, t, y, c, a):
        y = self.positive(y)
        c = self.positive(c)
        a = self.positive(a)
        return c * y/(y+a)**(1/3)

    def bertalanffy_model(self, t, y, c, a):
        y = self.positive(y)
        c = self.positive(c)
        a = self.positive(a)
        return c * y**(2/3) - a*y

    def gompertz_model(self, t, y, c, ymax):
        y = self.positive(y)
        c = self.positive(c)
        ymax = self.positive(ymax)
        return c * y * math.log(ymax/y)

    def positive(self, x, eps=1e-10):
        return max(x, eps)

    def compute_curve(self, solver, **params):
        # Kies ODE
        model_dict = {
            "linear": self.linear_model,
            "exponential_incr": self.exponential_model,
            "mendelsohn": self.mendelsohn_model,
            "exponential_flat": self.exponential_flattening_model,
            "logistic": self.logistic_model,
            "montroll": self.montroll_model,
            "allee": self.allee_model,
            "linear_limited": self.linear_limited_model,
            "surface_limited": self.surface_limited_model,
            "bertalanffy": self.bertalanffy_model,
            "gompertz": self.gompertz_model,
        }

        ODE = model_dict[self.model]
        model_params = self.find_needed_params(params)

        t = 0.0
        dt = 1.0
        y = params["y0"]
        ts, ys = [], []

        for t_target in self.ts_data:
            while t < t_target:
                if solver == "euler":
                    dydt = ODE(t, y, **model_params)
                    t = t + dt
                    y = y + dydt * dt
                elif solver == "heun":
                    dydt1 = ODE(t, y, **model_params)
                    t1 = t + dt
                    y1 = y + dydt1 * dt

                    dydt2 = ODE(t1, y1, **model_params)
                    t = t + dt
                    y = y + (dydt1 + dydt2) / 2.0 * dt

                elif solver == "runge-kutta":
                    k1 = ODE(t, y, **model_params)
                    k2 = ODE(t + 0.5 * dt, y + 0.5 * dt * k1, **model_params)
                    k3 = ODE(t + 0.5 * dt, y + 0.5 * dt * k2, **model_params)
                    k4 = ODE(t + dt, y + dt * k3, **model_params)
                    y = y + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
                    t = t + dt

            ts.append(t_target)
            ys.append(y)

        return ts, ys

    def find_needed_params(self, params):
        model_params_dict = {
            "linear": ["c"],
            "exponential_incr": ["c"],
            "logistic": ["c", "ymax"],
            "gompertz": ["c", "ymax"],
            "mendelsohn": ["c", "a"],
            "exponential_flat": ["c", "ymax"],
            "montroll": ["c", "ymax", "a"],
            "allee": ["c", "ymax", "ymin"],
            "linear_limited": ["c", "a"],
            "surface_limited": ["c", "a"],
            "bertalanffy": ["c", "a"],
        }

        needed_params = model_params_dict[self.model]
        model_params = {k: params[k] for k in needed_params if k in params}

        return model_params

    def MSE(self, solver, **params):
        _, ys_model = self.compute_curve(solver=solver, **params)
        mse = sum((y_data - y_model) ** 2 for y_data, y_model in zip(self.ys_data, ys_model)) / len(ys_model)
        return mse

    def search(self, params, solver, tries=1000):
        mse = self.MSE(solver=solver, **params)

        deltas = {k: 1.0 for k in params}
        max_iter, it = tries, 0
        while max(abs(d) for d in deltas.values()) > 1e-6 and it < max_iter:
            it += 1
            for k in params:
                new_params = params.copy()
                new_params[k] = params[k] + deltas[k]
                new_mse = self.MSE(solver=solver, **new_params)
                if new_mse < mse:
                    params, mse = new_params, new_mse
                    deltas[k] *= 1.2
                    continue

                new_params[k] = params[k] - deltas[k]
                new_mse = self.MSE(solver=solver, **new_params)
                if new_mse < mse:
                    params, mse = new_params, new_mse
                    deltas[k] *= -1.2
                    continue
                deltas[k] *= 0.2
        return params

    def score(self, mse, params, solver):
        _, ys_model = self.compute_curve(solver = solver, **params)
        n, k = len(ys_model), len(params)
        BIC = n * math.log(mse) + k * math.log(n)
        AIC = n * math.log(mse) + k * 2.0
        AICc = n * math.log(mse) + k * 2.0 * n / n - k - 1

        print(f"{BIC  = :.1f}")
        print(f"{AIC  = :.1f}")
        print(f"{AICc = :.1f}")

