import os
import random
import igraph as ig
from model import sir, competitive


def map_loader(name):
    """Load a network"""
    facebook_pages = ["artist", "athletes", "company", "government", "new_sites", "politician", "public_figure",
                      "tv_show"]
    name2file = {name: "test/test_data/facebook_clean_data/" + name + "_edges.txt" for name in facebook_pages}
    other_names = ["small_facebook", "karate", "wikipedia", "cornell"]
    other_files = ["facebook.txt", "karate.graphml", "wikipedia.graphml", "socfb-Cornell5.txt"]
    name2file.update(zip(other_names, ("test/test_data/" + file for file in other_files)))
    if name in name2file:
        graph = ig.Graph.Load(name2file[name])
        print("Number of nodes:", len(graph.vs))
        print("Number of edges:", len(graph.es))
        return graph
    elif name == "school":
        filenames = os.listdir("test/test_data/facebook_schools")
        filename = random.choice(filenames)
        graph = ig.Graph.Load("test/test_data/facebook_schools/" + filename)
        print("Number of nodes:", len(graph.vs))
        print("Number of edges:", len(graph.es))
        return graph


def test_network_naive(network_name, plot_mode="count"):
    """Run naively on network"""
    network = map_loader(network_name)
    epidemic = sir.SirNaive(network, plot_mode=plot_mode)
    epidemic.run()


def test_network_ode(size, plot=True):
    """Run ode on network"""
    epidemic = sir.SirODE(size, plot=plot)
    epidemic.run()


def test_network_redd_frost(network_name, plot=True):
    """Run reddfrost on network"""
    network = map_loader(network_name)
    epidemic = sir.SirReedFrost(network, plot=plot)
    epidemic.run()


def test_network_competitive_redd_frost(network_name, plot=True):
    """Run reddfrost on network"""
    network = map_loader(network_name)
    epidemic = competitive.SirReedFrost(network, plot=plot)
    epidemic.run()


if __name__ == "__main__":
    test_network_competitive_redd_frost("school")
    # test_network_ode(10000)
    # test_network_naive("facebook","time")
