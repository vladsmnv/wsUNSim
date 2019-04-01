# -*- coding: utf-8 -*-


from netmodel import Net, Node, Turn, Link, Route


# ----------------------------------------------------------------------------------------------
# ----------------------------------------- EXAMPLE 1 ------------------------------------------
# ----------------------------------------------------------------------------------------------
# +------------------------------+
# |                              |
# |                 ^3           |
# |                 |            |
# |                 |            |
# |                 |            |
# |                 |            |
# |                 |            |
# |                 |            |
# |                 |            |
# |                 |            |
# |                 |            |
# |                 |            |
# |                 |            |
# |    ------------->            |
# |    1             2           |
# +------------------------------+
def net_1_inhomogeneous(t):
    # Main params
    # T = 3600.0 * 2.  # time interval in sec (with extension interval)
    # t = t  # time step in sec - should be smaller or the same as the smallest link travel time!!!

    T = 60.0 * 4  # time interval in sec (with extension interval)
    t = t  # time step in sec - should be smaller or the same as the smallest link travel time!!!

    # network
    net = Net()
    node1 = Node(1, x=0.0, y=10.0)
    node2 = Node(2, x=20.0, y=10.0)
    node3 = Node(3, x=20.0, y=40.0)
    net.nodes.add_node(node1)
    net.nodes.add_node(node2)
    net.nodes.add_node(node3)
    link1_2 = Link(1, node1, node2, length=2.)
    link1_2.update_main_links_params_by_time_step(t)
    # link1_2.cvn_init(t, T)
    net.links.add_link(link1_2)
    link2_3 = Link(2, node2, node3, length=2.)
    link2_3.update_main_links_params_by_time_step(t)
    # link2_3.cvn_init(t, T)
    link2_3.outcapacity *= 0.75
    net.links.add_link(link2_3)
    net.generate_transition_flows(True)
    # routes
    N = [400.0, 800.0, 500.0]  # number of vehicles during each 20min
    route1 = Route(1, node1, node3, N, t, 1200.0)
    route1.path = [node1, link1_2, node2, link2_3, node3]
    net.routes.add_route(route1)
    net.routes.calc_origin_flows(t)

    return net, T, t


# ----------------------------------------------------------------------------------------------
# ----------------------------------------- EXAMPLE 2 ------------------------------------------
# ----------------------------------------------------------------------------------------------
# +------------------------------+
# | 3**                     **4  |
# |    **                 **     |
# |      ** 6           **       |
# |       ***         **         |
# |      **  **     **           |
# |     **     ** **             |
# |    **        * 5             |
# |   **         *               |
# |  2           *               |
# |              *               |
# |              *               |
# |              *               |
# |              *               |
# |              *               |
# |               1              |
# +------------------------------+
def net_2_Yperman(t):
    # Main params
    T = 90.0 * 60. + 3600.  # time interval in sec with extension interval
    t = t  # time step in sec - should be smaller or the same as the smallest link travel time!!!
    # network
    net = Net()
    node_orig1 = Node(1, x=10., y=0.)
    node_orig2 = Node(2, x=-30., y=70.)
    node_dest3 = Node(3, x=-30., y=110.)
    node_dest4 = Node(4, x=50., y=110.)
    node_diverge5 = Node(5, x=10., y=70.)
    node_merge6 = Node(6, x=-10., y=90.)
    net.nodes.add_node(node_orig1)
    net.nodes.add_node(node_orig2)
    net.nodes.add_node(node_diverge5)
    net.nodes.add_node(node_merge6)
    net.nodes.add_node(node_dest3)
    net.nodes.add_node(node_dest4)
    link1_5 = Link(15, node_orig1, node_diverge5, length=7., numlanes=4, capacity=8000, v0=120., w=-15.)
    link5_4 = Link(54, node_diverge5, node_dest4, length=5., numlanes=2, capacity=4000, v0=120., w=-15.)
    link5_6 = Link(56, node_diverge5, node_merge6, length=2.5, numlanes=2, capacity=4000, v0=120., w=-15.)
    link6_3 = Link(63, node_merge6, node_dest3, length=2.5, numlanes=2, capacity=4000, v0=120., w=-15.)
    link2_6 = Link(26, node_orig2, node_merge6, length=2.5, numlanes=2, capacity=4000, v0=120., w=-15.)
    net.links.add_link(link1_5)
    net.links.add_link(link5_4)
    net.links.add_link(link5_6)
    net.links.add_link(link6_3)
    net.links.add_link(link2_6)
    for link in net.links:
        link.update_main_links_params_by_time_step(t)
    net.generate_transition_flows(True)
    # routes
    routeflow13 = [1500.0, 1500.0, 1500.0]  # number of vehicles during each 30min
    routeflow14 = [1500.0, 1500.0, 1500.0]  # number of vehicles during each 30min
    routeflow23 = [1500.0, 0.0, 0.]  # number of vehicles during each 30min
    route13 = Route(1, node_orig1, node_dest3, routeflow13, t, 30.0 * 60.)
    route14 = Route(2, node_orig1, node_dest4, routeflow14, t, 30.0 * 60.)
    route23 = Route(3, node_orig2, node_dest3, routeflow23, t, 30.0 * 60.)
    route13.path = [node_orig1, link1_5, node_diverge5, link5_6, node_merge6, link6_3, node_dest3]
    route14.path = [node_orig1, link1_5, node_diverge5, link5_4, node_dest4]
    route23.path = [node_orig2, link2_6, node_merge6, link6_3, node_dest3]
    net.routes.add_route(route13)
    net.routes.add_route(route14)
    net.routes.add_route(route23)
    net.routes.calc_origin_flows(t)

    return net, T, t


# ----------------------------------------------------------------------------------------------
# ----------------------------------------- EXAMPLE 3 ------------------------------------------
# ----------------------------------------------------------------------------------------------
# +------------------------------+
# |                              |
# |                 ^3           |
# |                 |            |
# |                 |            |
# |                 |            |
# |                 |            |
# |                 |            |
# |    ------------>^ 2          |
# |    1            |            |
# |                 |            |
# |                 |            |
# |                 |            |
# |                 |            |
# |                 | 4          |
# |                              |
# +------------------------------+
# В этом примере задана часовая интенсивность на маршрутах (списочные переменные N1, N2) только для первого часа
# моделирования, тогда как общая длительность моделирования составляет 2 часа (T = 3600.0 * 2).
def net_3_merge(t):
    # Main params
    T = 3600.0 * 2  # time interval in sec (with extension interval)
    # network
    net = Net()
    node1 = Node(1, x=0., y=0.)
    node2 = Node(2, x=20., y=10., signalized=True)
    node3 = Node(3, x=20., y=40.)
    node4 = Node(4, x=20., y=-30.)
    net.nodes.add_node(node1)
    net.nodes.add_node(node2)
    net.nodes.add_node(node3)
    net.nodes.add_node(node4)

    # link1_2 = Link(12, node1, node2, length=2., numlanes=1, capacity=2000, v0=60., w=-20.) #, flow_det={2400:100, 2460:1, 2520:1})
    # link2_3 = Link(23, node2, node3, length=2., numlanes=1, capacity=2000, v0=60., w=-20.) #, flow_det={0:0})
    # link4_2 = Link(42, node4, node2, length=2., numlanes=1, capacity=2000, v0=60., w=-20.) #, flow_det={0:0})

    # ++++++++++++++++++++++++++++Node Model+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    link1_2 = Link(12, node1, node2, length=2., numlanes=1, capacity=2000, v0=60., w=-20., p=0.3) # p is the input links nonzero priority
    link2_3 = Link(23, node2, node3, length=2., numlanes=1, capacity=2000, v0=60., w=-20., p=0)
    link4_2 = Link(42, node4, node2, length=2., numlanes=1, capacity=2000, v0=60., w=-20., p=0.7)
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    net.links.add_link(link1_2)
    net.links.add_link(link2_3)
    net.links.add_link(link4_2)
    # TODO: Add this to every example
    t = min(t, net.links.calc_min_time_step())  # time step in sec - should be smaller or the same as the smallest link travel time!!!
    for link in net.links:
        # link.update_main_links_params_by_time_step(t)
        link.update_main_links_params_by_time_step(t)
    net.generate_transition_flows(True)

    # routes
    N1 = [1500.]  # number of vehicles during each hour
    N2 = [2000.]  # number of vehicles during each hour
    route13 = Route(13, node1, node3, N1, t, 3600.0)
    route13.path = [node1, link1_2, node2, link2_3, node3]
    net.routes.add_route(route13)
    route43 = Route(43, node4, node3, N2, t, 3600.0)
    route43.path = [node4, link4_2, node2, link2_3, node3]
    net.routes.add_route(route43)
    net.routes.calc_origin_flows(t)

    return net, T, t


# ----------------------------------------------------------------------------------------------
# ------------------------- EXAMPLE 4 as 2 but links splitted by 500m --------------------------
# ----------------------------------------------------------------------------------------------
# +------------------------------+
# | 3^                      ^ 4  |
# |    **                 **     |
# |      ** 6           **       |
# |       *^*         **         |
# |      **  **     **           |
# |     **     ** **             |
# |    **        ^ 5             |
# |   **         *               |
# |  2           *               |
# |              *               |
# |              *               |
# |              *               |
# |              *               |
# |              *               |
# |               1              |
# +------------------------------+
def net_4_Visum_Yperman(t):
    # Main params
    T = 90.0 * 60. + 3600.  # time interval in sec with extension interval
    t = t  # time step in sec - should be smaller or the same as the smallest link travel time!!!
    # network
    net = Net()
    net.import_visum_net('test_example_split500m.net')
    for link in net.links:
        link.update_main_links_params_by_time_step(t)
    # routes
    net.import_visum_routes('test_example_split500m_routes.att', (
        {'flow': [1500.0, 1500.0, 1500.0], 'time_step': t, 'time_int': 30. * 60.},
        {'flow': [1500.0, 1500.0, 1500.0], 'time_step': t, 'time_int': 30. * 60.},
        {'flow': [1500.0, 0.0, 0.0], 'time_step': t, 'time_int': 30. * 60.}))
    net.routes.calc_origin_flows(t)
    return net, T, t


# ----------------------------------------------------------------------------------------------
# ----------------------------------------- EXAMPLE 5 ------------------------------------------
# ----------------------------------------------------------------------------------------------
# +------------------------------+
# |                              |
# |                 ^3           |
# |                 |            |
# |                 |            |
# |                 |            |
# |                 |            |
# |                 |            |
# |    <------------^ 2          |
# |    4            |            |
# |                 |            |
# |                 | 5          |
# |                 |            |
# |                 |            |
# |                 | 1          |
# |                              |
# +------------------------------+
def net_5_diverge(t):
    # Main params
    T = 3600.0 * 2  # time interval in sec (with extension interval)
    t = t  # time step in sec - should be smaller or the same as the smallest link travel time!!!
    # network
    net = Net()
    node1 = Node(1, x=0., y=0.)
    node2 = Node(2, x=0., y=40.)
    node3 = Node(3, x=0., y=80.)
    node4 = Node(4, x=-40., y=50.)
    node5 = Node(5, x=0., y=20.)
    net.nodes.add_node(node1)
    net.nodes.add_node(node2)
    net.nodes.add_node(node3)
    net.nodes.add_node(node4)
    net.nodes.add_node(node5)
    link1_5 = Link(15, node1, node5, length=1., numlanes=1, capacity=2000, v0=60., w=-20.)
    link5_2 = Link(52, node5, node2, length=1., numlanes=1, capacity=2000, v0=60., w=-20.)
    link2_3 = Link(23, node2, node3, length=2., numlanes=1, capacity=2000, v0=60., w=-20.)
    link2_4 = Link(24, node2, node4, length=2., numlanes=1, capacity=1000, v0=60., w=-20.)
    link2_4.outcapacity *= 0.5
    net.links.add_link(link1_5)
    net.links.add_link(link5_2)
    net.links.add_link(link2_3)
    net.links.add_link(link2_4)
    for link in net.links:
        link.update_main_links_params_by_time_step(t)
    net.generate_transition_flows(True)
    # routes
    N1 = [1000.]  # number of vehicles during each hour
    N2 = [600.]  # number of vehicles during each hour
    route13 = Route(13, node1, node3, N1, t, 3600.0)
    route13.path = [node1, link1_5, node5, link5_2, node2, link2_3, node3]
    net.routes.add_route(route13)
    route14 = Route(14, node1, node4, N2, t, 3600.0)
    route14.path = [node1, link1_5, node5, link5_2, node2, link2_4, node4]
    net.routes.add_route(route14)
    net.routes.calc_origin_flows(t)

    return net, T, t


# ----------------------------------------------------------------------------------------------
# ----------------------------------------- EXAMPLE 6 ------------------------------------------
# ----------------------------------------------------------------------------------------------
# +-------------------------------+
# |                               |
# |         8<-     ^3            |
# |            -    |             |
# |             -   |             |
# |              -  |             |
# |               - |             |
# |                -|             |
# |    <------------^ 5<--------7 |
# |    4            | <-          |
# |                 |   -         |
# |                 | 2  -        |
# |                 |     -       |
# |                 |      -      |
# |                 | 1     -6    |
# |                               |
# +------------------------------+
def net_6_mimo_node(t):
    # Main params
    T = 3600.0 * 2  # time interval in sec (with extension interval)
    t = t  # time step in sec - should be smaller or the same as the smallest link travel time!!!
    # T = 60.0 * 2  # time interval in sec (with extension interval)
    # t = t  # time step in sec - should be smaller or the same as the smallest link travel time!!!
    # totT = 10 # number of simulation steps
    # dt = 15 # time step
    # network
    net = Net()

    node1 = Node(1, x=0., y=0.)
    node2 = Node(2, x=0., y=20.)
    node3 = Node(3, x=0., y=80.)
    node4 = Node(4, x=-40., y=50.)
    node5 = Node(5, x=0., y=40.)
    node6 = Node(6, x=40., y=25.)
    node7 = Node(7, x=40., y=50.)
    node8 = Node(8, x=-40., y=25.)

    net.nodes.add_node(node1)
    net.nodes.add_node(node2)
    net.nodes.add_node(node3)
    net.nodes.add_node(node4)
    net.nodes.add_node(node5)
    net.nodes.add_node(node6)
    net.nodes.add_node(node7)
    net.nodes.add_node(node8)

    link1_2 = Link(12, node1, node2, length=1., numlanes=1, capacity=2000, v0=60., w=-20.)
    link2_5 = Link(25, node2, node5, length=1., numlanes=1, capacity=2000, v0=60., w=-20.)
    link6_5 = Link(65, node6, node5, length=2., numlanes=1, capacity=2000, v0=60., w=-20.)
    link7_5 = Link(75, node7, node5, length=2., numlanes=1, capacity=2000, v0=60., w=-20.)
    link5_3 = Link(53, node5, node3, length=2., numlanes=1, capacity=2000, v0=60., w=-20.)
    link5_4 = Link(54, node5, node4, length=2., numlanes=1, capacity=2000, v0=60., w=-20.)
    link5_8 = Link(58, node5, node8, length=2., numlanes=1, capacity=2000, v0=60., w=-20.)

    link5_4.outcapacity *= 0.5

    net.links.add_link(link1_2)
    net.links.add_link(link2_5)
    net.links.add_link(link6_5)
    net.links.add_link(link7_5)

    net.links.add_link(link5_3)
    net.links.add_link(link5_4)
    net.links.add_link(link5_8)

    for link in net.links:
        # link.update_main_links_params_by_time_step(t)
        link.update_main_links_params_by_time_step(t)
    net.generate_transition_flows(True)
    # routes
    N1 = [1.]  # number of vehicles during each hour
    N2 = [1.]  # number of vehicles during each hour

    N3 = [1.]
    N4 = [1.]

    # route13 = Route(13, node1, node3, N1, t, 3600.0)
    route13 = Route(13, node1, node3, N1, t, 60.0)
    route13.path = [node1, link1_2, node2, link2_5, node5, link5_3, node3]
    net.routes.add_route(route13)

    # route14 = Route(14, node1, node4, N2, t, 3600.0)
    route14 = Route(14, node1, node4, N2, t, 60.0)
    route14.path = [node1, link1_2, node2, link2_5, node5, link5_4, node4]
    net.routes.add_route(route14)

    # route18 = Route(18, node1, node8, N3, t, 3600.0)
    route18 = Route(18, node1, node8, N3, t, 60.0)
    route18.path = [node1, link1_2, node2, link2_5, node5, link5_8, node8]
    net.routes.add_route(route18)

    # route64 = Route(64, node6, node4, N3, t, 3600.0)
    route64 = Route(64, node6, node4, N3, t, 60.0)
    route64.path = [node6, link6_5, node5, link5_4, node4]
    net.routes.add_route(route64)

    # route73 = Route(73, node7, node3, N4, t, 3600.0)
    route73 = Route(73, node7, node3, N4, t, 60.0)
    route73.path = [node7, link7_5, node5, link5_3, node3]
    net.routes.add_route(route73)

    # route78 = Route(78, node7, node8, N4, t, 3600.0)
    route78 = Route(78, node7, node8, N4, t, 60.0)
    route78.path = [node7, link7_5, node5, link5_8, node8]
    net.routes.add_route(route78)

    # net.routes.calc_origin_flows(t)
    net.routes.calc_origin_flows(t)

    return net, T, t


    # ----------------------------------------------------------------------------------------------
    # ----------------------------------------- EXAMPLE 7 ------------------------------------------
    # ----------------------------------------------------------------------------------------------
def net_7_mimo_node(t):
    # Main params
    T = 3600.0 * 2  # time interval in sec (with extension interval)
    t = t  # time step in sec - should be smaller or the same as the smallest link travel time!!!
    # network
    net = Net()

    node1 = Node(1, x=-40., y=20.)
    node2 = Node(2, x=-20., y=20.)
    node3 = Node(3, x=-20., y=40.)
    node4 = Node(4, x=20., y=40.)
    node5 = Node(5, x=20., y=20.)
    node6 = Node(6, x=40., y=20.)
    node7 = Node(7, x=40., y=-20.)
    node8 = Node(8, x=20., y=-20.)
    node9 = Node(9, x=20., y=-40.)
    node10 = Node(10, x=-20., y=-40.)
    node11 = Node(11, x=-20., y=-20.)
    node12 = Node(12, x=-40., y=-20.)

    net.nodes.add_node(node1)
    net.nodes.add_node(node2)
    net.nodes.add_node(node3)
    net.nodes.add_node(node4)
    net.nodes.add_node(node5)
    net.nodes.add_node(node6)
    net.nodes.add_node(node7)
    net.nodes.add_node(node8)
    net.nodes.add_node(node9)
    net.nodes.add_node(node10)
    net.nodes.add_node(node11)
    net.nodes.add_node(node12)



    link1_2 = Link(12, node1, node2, length=1., numlanes=1, capacity=2000, v0=60., w=-20.)
    link2_3 = Link(23, node2, node3, length=1., numlanes=1, capacity=2000, v0=60., w=-20.)
    link2_5 = Link(25, node2, node5, length=2., numlanes=1, capacity=2000, v0=60., w=-20.)
    link2_11 = Link(211, node2, node11, length=2., numlanes=1, capacity=2000, v0=60., w=-20.)
    link5_4 = Link(54, node5, node4, length=2., numlanes=1, capacity=2000, v0=60., w=-20.)
    link5_6 = Link(56, node5, node6, length=2., numlanes=1, capacity=2000, v0=60., w=-20.)
    link5_8 = Link(58, node5, node8, length=2., numlanes=1, capacity=2000, v0=60., w=-20.)
    link12_11 = Link(1211, node12, node11, length=2., numlanes=1, capacity=2000, v0=60., w=-20.)
    link11_2 = Link(112, node11, node2, length=2., numlanes=1, capacity=2000, v0=60., w=-20.)
    link11_8 = Link(118, node11, node8, length=2., numlanes=1, capacity=2000, v0=60., w=-20.)
    link11_10 = Link(1110, node11, node10, length=2., numlanes=1, capacity=2000, v0=60., w=-20.)

    link8_5 = Link(85, node8, node5, length=2., numlanes=1, capacity=2000, v0=60., w=-20.)
    link8_7 = Link(87, node8, node7, length=2., numlanes=1, capacity=2000, v0=60., w=-20.)
    link8_9 = Link(89, node8, node9, length=2., numlanes=1, capacity=2000, v0=60., w=-20.)

    # link5_6.outcapacity *= 0.5

    net.links.add_link(link1_2)
    net.links.add_link(link2_3)
    net.links.add_link(link2_5)
    net.links.add_link(link5_4)

    net.links.add_link(link5_6)
    net.links.add_link(link5_6)
    net.links.add_link(link5_8)

    net.links.add_link(link12_11)
    net.links.add_link(link11_2)
    net.links.add_link(link11_8)
    net.links.add_link(link11_10)
    net.links.add_link(link8_5)
    net.links.add_link(link8_7)
    net.links.add_link(link8_9)

    for link in net.links:
        link.update_main_links_params_by_time_step(t)
    net.generate_transition_flows(True)
    # routes
    N1 = [300.]  # number of vehicles during each hour
    N2 = [200.]  # number of vehicles during each hour
    N3 = [400.]

    N4 = [100.]
    N5 = [200.]
    N6 = [300.]

    N7 = [400.]
    N8 = [200.]

    route13 = Route(13, node1, node3, N1, t, 3600.0)
    route13.path = [node1, link1_2, node2, link2_3, node3]
    net.routes.add_route(route13)

    route14 = Route(14, node1, node4, N2, t, 3600.0)
    route14.path = [node1, link1_2, node2, link2_5, node5, link5_4, node4]
    net.routes.add_route(route14)

    route16 = Route(16, node1, node6, N3, t, 3600.0)
    route16.path = [node1, link1_2, node2, link2_5, node5, link5_6, node6]
    net.routes.add_route(route16)

    route1210 = Route(1210, node12, node10, N4, t, 3600.0)
    route1210.path = [node12, link12_11, node11, link11_10, node10]
    net.routes.add_route(route1210)

    route129 = Route(129, node12, node9, N5, t, 3600.0)
    route129.path = [node12, link12_11, node11, link11_8, node8, link8_9, node9]
    net.routes.add_route(route129)

    route127 = Route(127, node12, node7, N6, t, 3600.0)
    route127.path = [node12, link12_11, node11, link11_8, node8, link8_7, node7]
    net.routes.add_route(route127)

    route110 = Route(110, node1, node10, N7, t, 3600.0)
    route110.path = [node1, link1_2, node2, link2_11, node11, link11_10, node10]
    net.routes.add_route(route110)

    route124 = Route(124, node12, node4, N8, t, 3600.0)
    route124.path = [node12, link12_11, node11, link11_8, node8, link8_5, node5, link5_4, node4]
    net.routes.add_route(route124)

    net.routes.calc_origin_flows(t)

    return net, T, t
