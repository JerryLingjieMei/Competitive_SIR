from __future__ import division

import random

import igraph as ig
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from copy import copy


class SirODE(object):
    """Using ODE to solve sir model"""

    def __init__(self, size=10000, infection_rate=(0.02, 0.03), recover_rate=(0.01, 0.01), stop_time=None,
                 verbose=True):
        """Initialize the graph"""
        self.size = size
        self.infection_rate1, self.infection_rate2 = infection_rate
        self.recover_rate1, self.recover_rate2 = recover_rate
        self.verbose = verbose
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

    def plot(self):
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

    def __init__(self, graph: ig.Graph, time_step=(1., 2.), transmission_rate=(0.1, 0.2), verbose=True):
        self.graph = graph
        self.time_step = time_step
        self.transmission_rate = transmission_rate
        self.verbose = verbose
        infected1, infected2 = random.sample(list(self.graph.vs.indices), k=2)
        self.infected = [{infected1}, {infected2}]
        self.recovered = [set(), set()]
        self.susceptible = set(self.graph.vs.indices).difference(self.infected[0]).difference(
            self.infected[1])
        self.time = 0.
        self.next_time = [time_step[0], time_step[1]]
        self.time_history = [0.]
        self.susceptible_history = [self.get_susceptible()]
        self.infected_history = [self.get_infected(0)], [self.get_infected(1)]
        self.recovered_history = [self.get_recovered(0)], [self.get_recovered(1)]

    def run(self):
        """Run the simulation"""
        while len(self.infected[0]) > 0 or len(self.infected[1]) > 0:
            self.time = min(self.next_time)
            i = 0 if self.time == self.next_time[0] else 1
            self.next_time[i] += self.time_step[i]
            infected_set = set()
            for vertex in self.susceptible:
                self.infection_event(vertex, infected_set, i)

            if self.verbose:
                print("Node", infected_set, "are infected.")
            self.susceptible.difference_update(infected_set)
            self.recovered[i].update(self.infected[i])

            if self.verbose:
                print("Node", self.infected[i], "have recovered.")
            self.infected[i] = infected_set

            self.time_history.append(copy(self.time))
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


    def infection_event(self, vertex, infected_ls, i):
        """If vertex is infected, push it into the list"""
        infected_neighbors = len(self.infected[i].intersection(self.graph.neighbors(vertex)))
        rand = random.random()
        if rand > (1 - self.transmission_rate[i]) ** infected_neighbors:
            infected_ls.add(vertex)

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

    def plot(self):
        """Plot the sir process w.r.t. time"""
        plt.plot(self.time_history, self.susceptible_history, label="Susceptible")
        plt.plot(self.time_history, self.infected_history[0], label="Infected by 0")
        plt.plot(self.time_history, self.infected_history[1], label="Infected by 1")
        plt.plot(self.time_history, self.recovered_history[0], label="Recovered from 0")
        plt.plot(self.time_history, self.recovered_history[1], label="Recovered from 1")
        plt.legend()
        plt.show()


class ImmunizationReedFrost(SirReedFrost):
    """Immunized Reed Frost model"""
    def __init__(self, graph: ig.Graph, time_step=(1., 2.01), transmission_rate=(0.1, 0.2), immunity_rate=0.0,
                 verbose=True):
        super(ImmunizationReedFrost, self).__init__(graph, time_step, transmission_rate, verbose)
        self.removed = set()
        self.removed.update(random.choices(self.graph.vs.indices, k=int(immunity_rate * len(self.graph.vs.indices))))
        self.removed.difference_update(self.infected[0])
        self.removed.difference_update(self.infected[1])
        self.susceptible.difference_update(self.removed)
        self.susceptible_history = [self.get_susceptible()]
        self.infected_history = [self.get_infected(0)], [self.get_infected(1)]
        self.recovered_history = [self.get_recovered(0)], [self.get_recovered(1)]
