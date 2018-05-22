import igraph as ig
import numpy as np

def excess_degree(network: ig.Graph):
    """Return the generating function of excess degree"""
    degree_ls = list(network.degree(vertex) for vertex in network.vs.indices)
    mean = np.mean(degree_ls)
    return lambda x: sum(k * x ** (k - 1) for k in degree_ls) / mean / len(degree_ls)


def excess_degree_deriv(network: ig.Graph):
    """Return the derivative of the generarting function of excess degree"""
    degree_ls = list(network.degree(vertex) for vertex in network.vs.indices)
    mean = np.mean(degree_ls)
    return lambda x: sum(k * (k - 1) * x ** (k - 2) for k in degree_ls) / mean / len(degree_ls)


def get_threhsold(graph: ig.Graph):
    """Return the threshold for the epidemic to spread"""
    expected_degree = sum(graph.degree(vertex) for vertex in graph.vs)
    squared_degree = sum(graph.degree(vertex) ** 2 for vertex in graph.vs)
    return expected_degree / (squared_degree - expected_degree)
