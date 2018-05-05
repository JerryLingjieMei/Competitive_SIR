from __future__ import division

import random

import igraph as ig
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

from model.sir import get_threhsold


class SirODE(object):
    """Using ODE to solve sir model"""

    def __init__(self, size=10000, infection_rate=(0.02, 0.03), recover_rate=(0.01, 0.01), stop_time=None,
                 verbose=True, plot=True):
        """Initialize the graph"""
        self.size = size
        self.infection_rate1, self.infection_rate2 = infection_rate
        self.recover_rate1, self.recover_rate2 = recover_rate
        self.verbose = verbose
        self.plot = plot
        if stop_time is None:
            self.stop_time = 10000
        else:
            self.stop_time = stop_time

    def run(self):
        """Run the model"""

        def deriv(y, t, n=self.size, beta1=self.infection_rate1, beta2=self.infection_rate2, gamma1=self.recover_rate1,
                  gamma2=self.recover_rate2):
            """Calculate the derivatives"""
            s, i1, i2, r1, r2 = y
            dsdt = -beta1 * s * i1 / n - beta2 * s * i2 / n
            di1dt = beta1 * s * i1 / n - gamma1 * i1
            di2dt = beta2 * s * i2 / n - gamma2 * i2
            dr1dt = gamma1 * i1
            dr2dt = gamma2 * i2
            return dsdt, di1dt, di2dt, dr1dt, dr2dt

        y0 = (self.size - 2, 1, 1, 0, 0)

        self.time_history = np.linspace(0, self.stop_time, 100000)
        self.history = odeint(deriv, y0, self.time_history)
        if self.plot:
            self.plot_time()

    def plot_time(self):
        """Plot the graph w.r.t. the time"""
        cut_off = max(np.argwhere(self.history[:, 1] < 1)[0][0], np.argwhere(self.history[:, 2] < 1)[0][0])
        plt.plot(self.time_history[:cut_off], self.history[:cut_off, 0], label="Susceptible")
        plt.plot(self.time_history[:cut_off], self.history[:cut_off, 1], label="Infected by 1")
        plt.plot(self.time_history[:cut_off], self.history[:cut_off, 2], label="Infected by 2")
        plt.plot(self.time_history[:cut_off], self.history[:cut_off, 3], label="Recovered from 1")
        plt.plot(self.time_history[:cut_off], self.history[:cut_off, 4], label="Recovered from 2")
        plt.legend()
        plt.show()


class SirReedFrost(object):
    """Using Reed Frost to solve sir model"""

    def __init__(self, graph: ig.Graph, time_step=(15, 20), transmission_rate=(0.02, 0.02), verbose=True, plot=True):
        self.graph = graph
        self.time_step = time_step
        self.transmission_rate = transmission_rate
        self.infection_rate = 1 - (1 - transmission_rate[0]) ** (1 / time_step[0]), 1 - (1 - transmission_rate[1]) ** (
                1 / time_step[1])
        self.verbose = verbose
        self.plot = plot
        infected1, infected2 = random.sample(list(self.graph.vs.indices), k=2)
        self.infected = {infected1}, {infected2}
        self.recovered = set(), set()
        self.susceptible = set(self.graph.vs.indices).difference(self.infected[0]).difference(
            self.infected[1])
        self.time = 0
        self.susceptible_history = [self.get_susceptible()]
        self.infected_history = [self.get_infected(0)], [self.get_infected(1)]
        self.recovered_history = [self.get_recovered(0)], [self.get_recovered(1)]

    def run(self):
        """Run the simulation"""
        if self.verbose:
            print("Threshold for epidemics:", get_threhsold(self.graph))
            print("Rate of transmission:", self.transmission_rate)
        recover_time = [0] * len(self.graph.vs.indices), [0] * len(self.graph.vs.indices)
        for i in (0, 1):
            for vertex in self.infected[i]:
                recover_time[i][vertex] = self.time + self.time_step[i]
        while len(self.infected[0]) > 0 or len(self.infected[1]) > 0:
            self.time += 1
            infected_ls = [], []
            for vertex in self.susceptible:
                self.infection_event(vertex, infected_ls, recover_time)
            for i in (0, 1):
                self.infected[i].update(infected_ls[i])
            self.susceptible.difference_update(infected_ls[0])
            self.susceptible.difference_update(infected_ls[1])
            if self.verbose:
                print("Node", infected_ls, "are infected.")

            recover_ls = [], []
            for i in (0, 1):
                for vertex in self.infected[i]:
                    self.recover_event(vertex, recover_ls, recover_time)
            for i in (0, 1):
                self.infected[i].difference_update(recover_ls[i])
            for i in (0, 1):
                self.recovered[i].update(recover_ls[i])
            if self.verbose:
                print("Node", recover_ls, "have recovered.")

            self.susceptible_history.append(self.get_susceptible())
            for i in (0, 1):
                self.infected_history[i].append(self.get_infected(i))
                self.recovered_history[i].append(self.get_recovered(i))
            if self.verbose:
                print("Proportion of five types: {:.3f} {:.3f} {:.3f} {:.3f} {:.3f} ".format(
                    self.get_susceptible() / self.get_all(),
                    self.get_infected(0) / self.get_all(),
                    self.get_infected(1) / self.get_all(),
                    self.get_recovered(0) / self.get_all(),
                    self.get_recovered(1) / self.get_all()))
        if self.verbose:
            print("Outbreak size:", self.get_recovered(0), self.get_recovered(1), "out of", self.get_all())

        if self.plot:
            self.plot_time()

    def infection_event(self, vertex, infected_ls, recover_time):
        """If vertex is infected, push it into the list"""
        infected_neighbors = len(self.infected[0].intersection(self.graph.neighbors(vertex))), len(
            self.infected[1].intersection(self.graph.neighbors(vertex)))
        rand = random.random()
        if rand > (1 - self.infection_rate[0]) ** infected_neighbors[0]:
            infected_ls[0].append(vertex)
            recover_time[0][vertex] = self.time + self.time_step[0]
            return
        if rand < 1 - (1 - self.infection_rate[1]) ** infected_neighbors[1]:
            infected_ls[1].append(vertex)
            recover_time[1][vertex] = self.time + self.time_step[1]

    def recover_event(self, vertex, recover_ls, recover_time):
        """Recover a vertex after some time"""
        for i in (0, 1):
            if recover_time[i][vertex] == self.time:
                recover_ls[i].append(vertex)

    def get_all(self):
        """Return the number of all nodes"""
        return len(self.graph.vs)

    def get_susceptible(self):
        """Return the number of susceptible nodes at a time"""
        return len(self.susceptible)

    def get_infected(self, i):
        """Return the number of infected nodes at a time"""
        return len(self.infected[i])

    def get_recovered(self, i):
        """Return the number of recovered nodes at a time"""
        return len(self.recovered[i])

    def plot_time(self):
        """Plot the sir process w.r.t. time"""
        plt.plot(range(self.time + 1), self.susceptible_history, label="Susceptible")
        plt.plot(range(self.time + 1), self.infected_history[0], label="Infected by 0")
        plt.plot(range(self.time + 1), self.infected_history[1], label="Infected by 1")
        plt.plot(range(self.time + 1), self.recovered_history[0], label="Recovered from 0")
        plt.plot(range(self.time + 1), self.recovered_history[1], label="Recovered from 1")
        plt.legend()
        plt.show()
