from matplotlib import pyplot as plt

class CancerGrowthModel:
    def __init__(self, model, start_size, t_0, t_end, growth_rate, max_size=0.0, a=0.0, min_size=0.0):
        self.model = model
        self.start_size = start_size
        self.t_0 = t_0
        self.t_end = t_end
        self.growth_rate = growth_rate
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
        return self.start_size + self.growth_rate * self.start_size ** self.a  # werkt nog niet

    def exponential_flattening_model(self):
        return self.start_size + self.growth_rate * self.start_size * (1.0 - self.start_size / self.max_size)

    def logistic_model(self):
        return self.start_size + self.growth_rate * self.start_size * (1.0 - self.start_size / self.max_size)

    def montroll_model(self):
        return

    def allee_model(self):
        return self.start_size + self.growth_rate * self.start_size * (1.0 - self.start_size / self.max_size) * (
                    1.0 - self.min_size / self.start_size)

    def linear_limited_model(self):
        return self.start_size + self.growth_rate * (self.start_size / self.start_size + self.a)  # werkt nog niet

    def surface_limited_model(self):
        return

    def bertalanffy_model(self):
        return

    def gompertz_model(self):
        return

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

        for day in range(self.t_0, self.t_end):
            end_size = model_to_use()
            prediction.append(round(end_size, 3))
            self.start_size = end_size

        print(self.model, prediction)
        return prediction

    def plot(self, prediction):
        plt.axhline(0.0, color='k')
        plt.axvline(0.0, color='k')
        plt.plot(list(range(self.t_0, self.t_end + 1)), prediction, 'or')
        plt.suptitle('Tumor volume growth')
        plt.title(f'Model = {self.model}')
        plt.xlabel('$t$ (days)')
        plt.ylabel('$V(t)$ (mm³)')
        plt.grid(True)
        plt.box(False)
        plt.show()


def main():
    # Dit allemaal in notebook zetten ipv hier
    linear_model = CancerGrowthModel("linear", 1, 0, 20, 1.2)
    linear_prediction = linear_model.predict()

    exponential_model = CancerGrowthModel("exponential_incr", 3, 0, 20, 1.4)
    exponential_incr_prediction = exponential_model.predict()

    exponential_flattening_model = CancerGrowthModel("exponential_flat", 1, 0, 20, 1.2, 100)
    exponential_flattening_prediction = exponential_flattening_model.predict()

    logistic_model = CancerGrowthModel("logistic", 1, 0, 20, 1.2, 40)
    logistic_model_prediction = logistic_model.predict()

    allee_model = CancerGrowthModel("allee", 10, 0, 20, 1.2, 40, 0, 5)
    allee_prediction = allee_model.predict()

    mendelson_model = CancerGrowthModel("mendelsohn", 1, 0, 20, 1.2, 40, 1.2)
    mendelsohn_prediction = mendelson_model.predict()

    linear_model.plot(linear_prediction)
    exponential_model.plot(exponential_incr_prediction)
    exponential_flattening_model.plot(exponential_flattening_prediction)
    logistic_model.plot(logistic_model_prediction)
    allee_model.plot(allee_prediction)
    mendelson_model.plot(mendelsohn_prediction)

if __name__ == "__main__":
    main()