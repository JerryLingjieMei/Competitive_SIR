import os
import random
import igraph as ig
from competitive_sir.property import get_threhsold

def map_loader(name, verbose=True):
    """Load a network"""
    facebook_pages = ["artist", "athletes", "company", "government", "new_sites", "politician", "public_figure",
                      "tv_show"]
    name2file = {name: "tests/test_data/facebook_clean_data/" + name + "_edges.txt" for name in facebook_pages}
    other_names = ["small_facebook", "karate", "wikipedia", "cornell"]
    other_files = ["facebook.txt", "karate.graphml", "wikipedia.graphml", "socfb-Cornell5.txt"]
    name2file.update(zip(other_names, ("tests/test_data/" + file for file in other_files)))

    if name in name2file:
        graph = ig.Graph.Load(name2file[name])
        if verbose:
            print("Number of nodes:", len(graph.vs))
            print("Number of edges:", len(graph.es))
            print("Transmission threshold:", get_threhsold(graph))
        return graph
    elif name == "school":
        filenames = list(f for f in os.listdir("tests/test_data/facebook_schools") if f.endswith(".txt"))
        filename = random.choice(filenames)
        graph = ig.Graph.Load("tests/test_data/facebook_schools/" + filename)
        if verbose:
            print("Number of nodes:", len(graph.vs))
            print("Number of edges:", len(graph.es))
            print("Transmission threshold:", get_threhsold(graph))
        return graph


def power_law(size=7000, edge=90000, alpha=2.6, verbose=True):
    """Generate power law network model"""
    graph = ig.Graph.Static_Power_Law(size, edge, alpha)
    if verbose:
        print("Number of nodes:", len(graph.vs))
        print("Number of edges:", len(graph.es))
        print("Transmission threshold:", get_threhsold(graph))
    return graph


def small_world(size=7000, edges=91000, verbose=True):
    """Generate small world network model"""
    graph = ig.Graph.Watts_Strogatz(1, size, edges // size, .25)
    if verbose:
        print("Number of nodes:", len(graph.vs))
        print("Number of edges:", len(graph.es))
        print("Transmission threshold:", get_threhsold(graph))
    return graph


def erdos_renyi(size=7000, edges=90000, verbose=True):
    """Generate Erdos Renyi network model"""
    graph = ig.Graph.Erdos_Renyi(size, m=edges)
    if verbose:
        print("Number of nodes:", len(graph.vs))
        print("Number of edges:", len(graph.es))
        print("Transmission threshold:", get_threhsold(graph))
    return graph
