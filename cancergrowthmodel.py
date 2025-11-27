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

    def exponential_model(self):
        return self.start_size + self.growth_rate * self.start_size

    def mendelsohn_model(self):
        return self.start_size + self.growth_rate * self.start_size ** self.a

    def exponential_flattening_model(self):
        return self.start_size + self.growth_rate * (1.0 - self.start_size / self.max_size)

    def logistic_model(self):
        return self.start_size + self.growth_rate * self.start_size * (1.0 - self.start_size / self.max_size)

    def montroll_model(self):
        return self.start_size + self.growth_rate * self.start_size * (1.0 - self.start_size ** self.a / self.max_size ** self.a)

    def allee_model(self):
        return self.start_size + self.growth_rate * self.start_size * (1.0 - self.start_size / self.max_size) * (
                    1.0 - self.min_size / self.start_size)

    def linear_limited_model(self):
        return self.start_size + self.growth_rate * (self.start_size / (self.start_size + self.a))

    def surface_limited_model(self):
        return self.start_size + self.growth_rate * (self.start_size / ((self.start_size + self.a) ** 1/3))

    def bertalanffy_model(self):
        return self.start_size + self.growth_rate * (1.0 - self.a / (self.start_size ** 2 / 3)) * self.start_size

    def gompertz_model(self):
        return self.start_size + self.growth_rate * self.start_size * math.log(self.max_size/self.start_size)

    def predict(self):
        prediction = [self.start_size]  # add first value

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

        model_to_use = model_dict[self.model]

        for day in range(self.t_0, self.t_end, self.step_size):
            end_size = model_to_use()
            prediction.append(round(end_size, 3))
            self.start_size = end_size

        return prediction

    def plot(self, prediction):
        plt.axhline(0.0, color='k')
        plt.axvline(0.0, color='k')
        plt.plot(list(range(self.t_0, self.t_end + 1)), prediction, 'o-', linewidth=2)
        plt.suptitle('Tumor volume growth')
        plt.title(f'Model = {self.model}')
        plt.xlabel('$t$ (days)')
        plt.ylabel('$V(t)$ (mm³)')
        plt.grid(True)
        plt.box(False)
        plt.show()
