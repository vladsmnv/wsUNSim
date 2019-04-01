#!/usr/bin/env/python
# -*- coding: utf-8 -*-


import timeit
from gui import MapDrawer
from reports import write_results_to_file
from ltm import ltm
from examples import net_1_inhomogeneous, net_2_Yperman, net_3_merge, net_4_Visum_Yperman, net_5_diverge, net_6_mimo_node, net_7_mimo_node

# TODO: Add unit tests and profilers based on examples (what to compare for correctness?)
# net, T, t = net_1_inhomogeneous(15.)  # inhomogeneous node example
net, T, t = net_2_Yperman(15.)  # Yperman example (without disaggregated links)
# net, T, t = net_3_merge(15.)  # merge node example
# net, T, t = net_4_Visum_Yperman(15.)  # visum import (Yperman) (with disggregated links by 500m)
# net, T, t = net_5_diverge(15.)  # diverge node example
# net, T, t = net_6_mimo_node(15.)  # MIMO example
# net, T, t = net_7_mimo_node(15.)  # Simple network MIMO example

# net, totT, dt = net_1_inhomogeneous(15.)  # inhomogeneous node example
# net, totT, dt = net_2_Yperman(15.)  # Yperman example (without disaggregated links)
# net, totT, dt = net_3_merge(15.)  # merge node example
# net, totT, dt = net_4_Visum_Yperman(15.)  # visum import (Yperman) (with disggregated links by 500m)
# net, totT, dt = net_5_diverge(15.)  # diverge node example
# net, T, t = net_6_mimo_node(15.)  # MIMO example
# net, totT, dt = net_7_mimo_node(15.)  # Simple network MIMO example


ltm(net, T, t)
write_results_to_file(net, T, t)
mapdrawer = MapDrawer(net, T, t)
# pyltm(net, totT, dt)
# write_results_to_file(net, totT, dt)
# mapdrawer = MapDrawer(net, totT, dt)
mapdrawer.plot_net()
