import math

from matplotlib import pyplot as plt

class CancerGrowthModel:
    def __init__(self, model, start_size, t_0, t_end, growth_rate, step_size = 1, max_size=0.0, a=0.0, min_size=0.0):
        self.model = model
        self.start_size = start_size
        self.t_0 = t_0
        self.t_end = t_end
        self.growth_rate = growth_rate
        self.step_size = step_size
        self.max_size = max_size
        self.min_size = min_size
        self.a = a

    def __repr__(self):
        return (f"CancerGrowthModel("
                f"model={self.model}, "
                f"start_size={self.start_size}, "
                f"t_0={self.t_0}, "
                f"t_end={self.t_end}, "
                f"growth_rate={self.growth_rate})")

    def __str__(self):
        return (
            "Hidden Markov Model\n"
            f"Chosen model(s):\n{self.model}\n\n"
            f"Start size tumor:\n{self.start_size}\n\n"
            f"Starting time (in days):\n{self.t_0}\n"
            f"Ending time (in days):\n{self.t_end}"
            f"Growth rate:\n{self.growth_rate}\n"
        )

    def linear_model(self):
        return self.start_size + self.growth_rate

    # def linear_model(self):
    #     return self.start_size + self.growth_rate

    # def exponential_model(self):
    #     return self.start_size + self.growth_rate * self.start_size
    #
    # def mendelsohn_model(self):
    #     return self.start_size + self.growth_rate * self.start_size ** self.a
    #
    # def exponential_flattening_model(self):
    #     return self.start_size + self.growth_rate * (1.0 - self.start_size / self.max_size)
    #
    # def logistic_model(self):
    #     return self.start_size + self.growth_rate * self.start_size * (1.0 - self.start_size / self.max_size)
    #
    # def montroll_model(self):
    #     return self.start_size + self.growth_rate * self.start_size * (1.0 - self.start_size ** self.a / self.max_size ** self.a)
    #
    # def allee_model(self):
    #     return self.start_size + self.growth_rate * self.start_size * (1.0 - self.start_size / self.max_size) * (
    #                 1.0 - self.min_size / self.start_size)
    #
    # def linear_limited_model(self):
    #     return self.start_size + self.growth_rate * (self.start_size / (self.start_size + self.a))
    #
    # def surface_limited_model(self):
    #     return self.start_size + self.growth_rate * (self.start_size / ((self.start_size + self.a) ** 1/3))
    #
    # def bertalanffy_model(self):
    #     return self.start_size + self.growth_rate * (1.0 - self.a / (self.start_size ** 2 / 3)) * self.start_size
    #
    # def gompertz_model(self):
    #     return self.start_size + self.growth_rate * self.start_size * math.log(self.max_size/self.start_size)

    def predict(self):
        prediction = [self.start_size]  # add first value

        model_dict = {
            "linear": self.linear_model,
            # "exponential_incr": self.exponential_model,
            # "mendelsohn": self.mendelsohn_model,
            # "exponential_flat": self.exponential_flattening_model,
            # "logistic": self.logistic_model,
            # "montroll": self.montroll_model,
            # "allee": self.allee_model,
            # "linear_limited": self.linear_limited_model,
            # "surface_limited": self.surface_limited_model,
            # "bertalanffy": self.bertalanffy_model,
            # "gompertz": self.gompertz_model,
        }

        model_to_use = model_dict[self.model]

        for day in range(self.t_0, self.t_end, self.step_size):
            end_size = model_to_use()
            prediction.append(round(end_size, 3))
            self.start_size = end_size

        return prediction

    def plot(self, prediction, t_s = [], y_s = []):
        plt.axhline(0.0, color='k')
        plt.axvline(0.0, color='k')
        plt.plot(list(range(self.t_0, self.t_end + 1)), prediction, 'o-', linewidth=2)
        plt.plot(t_s, y_s, 'o-', linewidth=2)
        plt.suptitle('Tumor volume growth')
        plt.title(f'Model = {self.model}')
        plt.xlabel('$t$ (days)')
        plt.ylabel('$V(t)$ (mm³)')
        plt.grid(True)
        plt.box(False)
        plt.show()

    def compute_curve(self, c):

        def ODE(t, y):
            return  c  # <- het model heeft nu 2 parameters

        dt = 1.0 / 1.0
        t, y = 0.0, self.start_size  # <- de beginwaarde y(0) is ook een parameter
        ts, ys = [t], [y]
        for _ in range(self.t_end):
            # Eerste update
            dydt1 = ODE(t, y)
            t1 = t + 0.5 * dt
            y1 = y + 0.5 * dydt1 * dt
            # Tweede update
            dydt2 = ODE(t1, y1)
            t2 = t + 0.5 * dt
            y2 = y + 0.5 * dydt2 * dt
            # Derde update
            dydt3 = ODE(t2, y2)
            t3 = t + dt
            y3 = y + dydt2 * dt
            # Vierde update volgens Runge-Kutta
            dydt4 = ODE(t3, y3)
            t = t + dt
            y = y + (dydt1 + 2.0 * dydt2 + 2.0 * dydt3 + dydt4) / 6.0 * dt
            # Bewaar tussenresultaten
            ts.append(t)
            ys.append(y)
        return ts, ys



    # def ODE(self, t, y):
    #     return self.growth_rate
    #
    #
    #
    # def ode_solver(self, solver):
    #
    #     if solver == "heun":
    #         t, y = 0.0, self.start_size
    #         ts_heun, ys_heun = [t], [y]
    #         dt = 1.0 / self.t_end
    #
    #         for _ in range(self.t_end):
    #             # Eerste update
    #             dydt1 = self.ODE(t, y)
    #             print(dydt1)
    #             t1 = t + dt
    #             y1 = y + dydt1 * dt
    #             # Tweede update volgens Heun
    #             dydt2 = self.ODE(t1, y1)
    #             t = t + dt
    #             y = y + (dydt1 + dydt2) / 2.0 * dt
    #             # Bewaar tussenresultaten
    #             ts_heun.append(t)
    #             ys_heun.append(y)
    #         print(ts_heun)
    #         print(ys_heun)
    #         return ts_heun, ys_heun
    #     return None