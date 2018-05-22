from __future__ import division

import random

import igraph as ig
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint


class SirNaive(object):
    """Naive implementation SIR model"""

    def __init__(self, graph: ig.Graph, infection_rate=0.002, recover_rate=0.01,
                 verbose=True):
        """"Initialize a sir model with designated infection and recove rate"""
        self.graph = graph
        self.infection_rate = infection_rate
        self.recover_rate = recover_rate
        self.verbose = verbose

        self.infected = set(random.sample(list(self.graph.vs.indices), k=1))
        self.recovered = set()
        self.susceptible = set(self.graph.vs.indices).difference(self.infected)

        self.links = set()
        self.links.update((i, s) for i in self.infected for s in self.graph.neighbors(i) if s in self.susceptible)

        self.time_history = [0.]
        self.susceptible_history = [len(self.susceptible)]
        self.infected_history = [len(self.infected)]
        self.recovered_history = [len(self.recovered)]

        self.time = 0.
        self.count = 0

        if self.verbose:
            print("Proportion of three types: {:.3f} {:.3f} {:.3f}".format(self.get_susceptible() / self.get_all(),
                                                                           self.get_infected() / self.get_all(),
                                                                           self.get_recovered() / self.get_all()))

    def run(self):
        """Run the simulation of the sir model"""
        while len(self.infected) > 0:
            rate1, rate2 = len(self.links) * self.infection_rate, len(self.infected) * self.recover_rate
            rate = [float(rate1) / (rate1 + rate2), float(rate2) / (rate1 + rate2)]
            if np.random.choice([0, 1], p=rate) == 0:
                self.infect_event()
            else:
                self.recover_event()
            self.susceptible_history.append(self.get_susceptible())
            self.infected_history.append(self.get_infected())
            self.recovered_history.append(self.get_recovered())
            self.count += 1
            self.time += random.expovariate(1 / (rate1 + rate2))
            self.time_history.append(self.time)

            if self.verbose:
                print("Proportion of three types: {:.3f} {:.3f} {:.3f}".format(self.get_susceptible() / self.get_all(),
                                                                               self.get_infected() / self.get_all(),
                                                                               self.get_recovered() / self.get_all()))
        if self.verbose:
            print("Outbreak size:", self.get_all() - self.get_susceptible(), "out of", self.get_all())

    def recover_event(self):
        """A node is recovered"""
        node = random.sample(self.infected, 1)[0]
        self.infected.remove(node)
        self.recovered.add(node)
        self.links.difference_update((node, s) for s in self.graph.neighbors(node))
        if self.verbose:
            print("Node", node, "has recovered. ")

    def infect_event(self):
        """A node is infected"""
        link = random.sample(self.links, 1)[0]
        source, node = link
        self.susceptible.remove(node)
        self.infected.add(node)
        self.links.update((node, s) for s in self.graph.neighbors(node) if s in self.susceptible)
        self.links.difference_update((i, node) for i in self.graph.neighbors(node))
        if self.verbose:
            print("Node", node, "is infect by node", source)

    def get_all(self):
        """Return the number of all nodes"""
        return len(self.graph.vs)

    def get_susceptible(self):
        """Return the number of susceptible nodes at a time"""
        return len(self.susceptible)

    def get_infected(self):
        """Return the number of infected nodes at a time"""
        return len(self.infected)

    def get_recovered(self):
        """Return the number of recovered nodes at a time"""
        return len(self.recovered)

    def plot_time(self):
        """Plot the sir process w.r.t. time"""
        plt.plot(self.time_history, self.susceptible_history, label="Susceptible")
        plt.plot(self.time_history, self.infected_history, label="Infected")
        plt.plot(self.time_history, self.recovered_history, label="Recovered")
        plt.legend()
        plt.show()

    def plot_count(self):
        """Plot the sir process w.r.t count"""
        plt.plot(range(self.count + 1), self.susceptible_history, label="Susceptible")
        plt.plot(range(self.count + 1), self.infected_history, label="Infected")
        plt.plot(range(self.count + 1), self.recovered_history, label="Recovered")
        plt.legend()
        plt.show()


class SirReedFrost(object):
    """Using Reed Frost to solve sir model"""

    def __init__(self, graph: ig.Graph, transmission_rate=0.2, verbose=True):
        self.graph = graph
        self.transmission_rate = transmission_rate
        self.verbose = verbose
        self.infected = set(random.sample(list(self.graph.vs.indices), k=1))
        self.recovered = set()
        self.susceptible = set(self.graph.vs.indices).difference(self.infected)
        self.time = 0
        self.susceptible_history = [self.get_susceptible()]
        self.infected_history = [self.get_infected()]
        self.recovered_history = [self.get_recovered()]

    def run(self):
        """Run the simulation"""
        while len(self.infected) > 0:
            self.time += 1
            infected_set = set()
            for vertex in self.susceptible:
                self.infection_event(vertex, infected_set)
            self.infected.update(infected_set)
            self.susceptible.difference_update(infected_set)
            if self.verbose:
                print("Node", infected_set, "are infected.")

            self.recovered.update(self.infected)
            if self.verbose:
                print("Node", self.infected, "have recovered.")
            self.infected = infected_set

            self.susceptible_history.append(self.get_susceptible())
            self.infected_history.append(self.get_infected())
            self.recovered_history.append(self.get_recovered())
            if self.verbose:
                print("Proportion of three types: {:.3f} {:.3f} {:.3f}".format(
                    self.get_susceptible() / self.get_all(),
                    self.get_infected() / self.get_all(),
                    self.get_recovered() / self.get_all()))
        if self.verbose:
            print("Outbreak size:", self.get_all() - self.get_susceptible(), "out of", self.get_all())

    def infection_event(self, vertex, infected_ls):
        """If vertex is infected, push it into the list"""
        infected_neighbors = len(self.infected.intersection(self.graph.neighbors(vertex)))
        if random.random() > (1 - self.transmission_rate) ** infected_neighbors:
            infected_ls.add(vertex)

    def recover_event(self, vertex, recover_ls, recover_time):
        """Recover a vertex after some time"""
        if recover_time[vertex] == self.time:
            recover_ls.append(vertex)

    def get_all(self):
        """Return the number of all nodes"""
        return len(self.graph.vs)

    def get_susceptible(self):
        """Return the number of susceptible nodes at a time"""
        return len(self.susceptible)

    def get_infected(self):
        """Return the number of infected nodes at a time"""
        return len(self.infected)

    def get_recovered(self):
        """Return the number of recovered nodes at a time"""
        return len(self.recovered)

    def plot(self):
        """Plot the sir process w.r.t. time"""
        plt.plot(range(self.time + 1), self.susceptible_history, label="Susceptible")
        plt.plot(range(self.time + 1), self.infected_history, label="Infected")
        plt.plot(range(self.time + 1), self.recovered_history, label="Recovered")
        plt.legend()
        plt.show()


class SirODE(object):
    """Using ODE to solve sir model"""

    def __init__(self, size=10000, infection_rate=0.02, recover_rate=0.01, stop_time=None,
                 verbose=True):
        """Initialize the graph"""
        self.size = size
        self.infection_rate = infection_rate
        self.recover_rate = recover_rate
        self.verbose = verbose
        if stop_time is None:
            self.stop_time = 10000
        else:
            self.stop_time = stop_time

    def run(self):
        """Run the model"""

        def deriv(y, t, n=self.size, beta=self.infection_rate, gamma=self.recover_rate):
            """Calculate the derivatives"""
            s, i, r = y
            dsdt = -beta * s * i / n
            didt = beta * s * i / n - gamma * i
            drdt = gamma * i
            return dsdt, didt, drdt

        y0 = (self.size - 1, 1, 0)

        self.time_history = np.linspace(0, self.stop_time, 100000)
        self.history = odeint(deriv, y0, self.time_history)

    def plot(self):
        """Plot the graph w.r.t. the time"""
        cut_off = np.argwhere(self.history[:, 1] < 1)[0][0]
        plt.plot(self.time_history[:cut_off], self.history[:cut_off, 0], label="Susceptible")
        plt.plot(self.time_history[:cut_off], self.history[:cut_off, 1], label="Infected")
        plt.plot(self.time_history[:cut_off], self.history[:cut_off, 2], label="Recovered")
        plt.legend()
        plt.show()
