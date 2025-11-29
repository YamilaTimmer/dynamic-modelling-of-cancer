import math

class CancerGrowthModel:
    def __init__(self, start_size, t_end, step_size = 1):
        self.start_size = start_size
        self.t_end = t_end
        self.step_size = step_size


    def __repr__(self):
        return (f"CancerGrowthModel("
                f"start_size={self.start_size}, "
                f"t_end={self.t_end}, ")

    def __str__(self):
        return (
            "Hidden Markov Model\n"
            f"Start size tumor:\n{self.start_size}\n\n"
            f"Ending time (in days):\n{self.t_end}"
            f"Growth rate:\n{self.growth_rate}\n"
        )

    def linear_model(self, t, y, c):
        return c

    def exponential_model(self, t, y, c):
        return c * y

    def mendelsohn_model(self, t, y, c, a):
        return c * y ** a

    def exponential_flattening_model(self, t, y, c, ymax):
        return c * (ymax - y)

    def logistic_model(self, t, y, c, ymax):
        return c * y * (1 - y / ymax)

    def montroll_model(self, t, y, c, ymax, a):
        return c * y * (1- y**a/ymax**a)

    def allee_model(self, t, y, c, ymin, ymax):
        return c * (y - ymin) * (ymax-y)

    def linear_limited_model(self, t, y, c, a):
        return c * y/(y+a)

    def surface_limited_model(self, t, y, c, a):
        return c * y/(y+a)**(1/3)

    def bertalanffy_model(self, t, y, c, a):
        return c * y**(2/3) - a*y

    def gompertz_model(self, t, y, c, ymax):
        return c * y * math.log(ymax/y)

    def compute_curve(self, model = "linear", solver = "runge-kutta", **params):
        dt = 1.0 / self.step_size
        t, y = 0.0, params["y0"]
        ts, ys = [t], [y]

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

        ODE = model_dict[model]
        needed_params = model_params_dict[model]
        model_params = {k: params[k] for k in needed_params if k in params}

        if solver == "euler":
            # ODE-solver euler
            for _ in range(self.t_end):
                # Bereken de slope
                slope = ODE(t, y, **model_params)
                # bereken de verandering van de slope
                change_y = dt * slope

                # update the results.
                t = t + dt
                y = y + change_y
                ts.append(t)
                ys.append(y)

        elif solver == "heun":
            # ODE-solver
            for _ in range(self.t_end):
                # Eerste update
                dydt1 = ODE(t, y, **model_params)
                t1 = t + dt
                y1 = y + dydt1 * dt
                # Tweede update volgens Heun
                dydt2 = ODE(t1, y1, **model_params)
                t = t + dt
                y = y + (dydt1 + dydt2) / 2.0 * dt
                # Bewaar tussenresultaten
                ts.append(t)
                ys.append(y)

        elif solver == "runge-kutta":

            for _ in range(self.t_end):
                # Eerste update
                dydt1 = ODE(t, y, **model_params)
                t1 = t + 0.5 * dt
                y1 = y + 0.5 * dydt1 * dt
                # Tweede update
                dydt2 = ODE(t1, y1, **model_params)
                t2 = t + 0.5 * dt
                y2 = y + 0.5 * dydt2 * dt
                # Derde update
                dydt3 = ODE(t2, y2, **model_params)
                t3 = t + dt
                y3 = y + dydt3 * dt
                # Vierde update volgens Runge-Kutta
                dydt4 = ODE(t3, y3, **model_params)
                t = t + dt
                y = y + (dydt1 + 2.0 * dydt2 + 2.0 * dydt3 + dydt4) / 6.0 * dt
                # Bewaar tussenresultaten
                ts.append(t)
                ys.append(y)

        else:
            return None

        return ts, ys

    def MSE(self, model = "linear", solver = "runge-kutta", ys_data=None, **params):

            _, ys_model = self.compute_curve(model=model, solver=solver, **params)
            sum_squared_error = 0.0
            for y_data, y_model in zip(ys_data, ys_model):
                error = y_data - y_model
                sum_squared_error += error * error
            mean_squared_error = sum_squared_error / len(ys_model)
            return mean_squared_error

