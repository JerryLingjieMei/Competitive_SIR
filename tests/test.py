from competitive_sir import sir, competitive
from tests.data_loader import map_loader


def test_sir_naive(network, infection_rate=.002, recover_rate=.01, verbose=False, plot=False):
    """Run naively on network"""
    epidemic = sir.SirNaive(network, infection_rate=infection_rate, recover_rate=recover_rate, verbose=verbose)
    epidemic.run()
    if plot:
        epidemic.plot_count()


def test_sir_ode(size, infection_rate=.02, recover_rate=.01, verbose=False, plot=False):
    """Run ode on network"""
    epidemic = sir.SirODE(size, infection_rate=infection_rate, recover_rate=recover_rate, verbose=verbose)
    epidemic.run()
    if plot:
        epidemic.plot()


def test_sir_redd_frost(network, transmission_rate=.2, verbose=False, plot=False):
    """Run reedfrost on network"""
    epidemic = sir.SirReedFrost(network, transmission_rate=transmission_rate, verbose=verbose)
    epidemic.run()
    if plot:
        epidemic.plot()


def test_competitive_reed_frost(network, time_step=(1., 2.), transmission_rate=(.1, .2), verbose=False,
                                plot=False):
    """Run competitve reedfrost on network"""
    epidemic = competitive.SirReedFrost(network, time_step=time_step, transmission_rate=transmission_rate,
                                        verbose=verbose)
    epidemic.run()
    if plot:
        epidemic.plot()


def test_competitive_immunization(network, time_step=(1., 2.), transmission_rate=(.1, .2), immunity_rate=.5,
                                  verbose=False,
                                  plot=False):
    """Run competitive reedfrost with immunization on network"""
    epidemic = competitive.ImmunizationReedFrost(network, time_step=time_step, transmission_rate=transmission_rate,
                                                 immunity_rate=immunity_rate, verbose=verbose)
    epidemic.run()
    if plot:
        epidemic.plot()


def test_competitive_ode(size, infection_rate=(.02, .03), recover_rate=(.01, .01), verbose=False, plot=False):
    """Run competitive ode on network"""
    epidemic = competitive.SirODE(size, infection_rate=infection_rate, recover_rate=recover_rate, verbose=verbose)
    epidemic.run()
    if plot:
        epidemic.plot()


if __name__ == "__main__":
    network = map_loader("government")
    for test in [test_sir_naive, test_sir_redd_frost, test_competitive_reed_frost,
                 test_competitive_immunization]:
        test(network, verbose=True, plot=True)
    for test in [test_sir_ode, test_competitive_ode]:
        test(len(network.vs), verbose=True, plot=True)
