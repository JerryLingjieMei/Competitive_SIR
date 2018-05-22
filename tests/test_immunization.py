from __future__ import division
from competitive_sir import competitive
from tests.data_loader import map_loader
import numpy as np
import json


def test_immunization(network, time_step=(10, 20), transmission_rate=(0.1, 0.2), immunity_rate=.5, iterations=50):
    """Run simulations of competitive sir with immunization and return the number of infected nodes"""
    result = [], []
    for iter in range(iterations):
        epidemic = competitive.ImmunizationReedFrost(network, time_step=time_step, transmission_rate=transmission_rate,
                                                     immunity_rate=immunity_rate, verbose=False)
        epidemic.run()
        result[0].append(epidemic.get_recovered(0))
        result[1].append(epidemic.get_recovered(1))
    return result


if __name__ == "__main__":
    network = map_loader("government")
    x = np.linspace(0, 1., 100)
    infections = list(test_immunization(network, time_step=(10, 20), immunity_rate=xx) for xx in x)
    with open("immunization_prob.json", "w") as out:
        json.dump(infections, out)
