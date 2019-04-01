# -*- coding: utf-8 -*-


from cmd_helpers import printProgress
from netmodel import Link, Node, Net
import numpy as np
import copy


# TODO: Use assert checks and def docs for type: param: return: etc.
# TODO: Make docs and comments throughout the code
# TODO: Check for some very small values!!! Maybe better to round them or use Decimal!!!
# TODO: Time dependent attrs (so-called events)
# TODO: Change to event-based method (as alternative!)
# TODO: Make check and auto-correction for: 1. outcapa <= incapa 2. t >= min(L/v0 for all links)
# TODO: Make unit tests!!!
# TODO: REALTIME DATA FROM sensors. Change cuminflows with sensors? What to do with speeds? Do in separate def

def ltm(net, time_horizon, time_step):
# def pyltm(net, totT, dt):
    assert isinstance(net, Net)


    def get_inversecumflow(flowval, linkflow_mint, routeflow_mint, linkflow_dt, routeflow_dt):
        if linkflow_dt > 0.:
            coef = (flowval - linkflow_mint) / linkflow_dt
            return max(routeflow_mint + coef * routeflow_dt, 0.)
        else:
            return routeflow_mint


    # def get_timeinterpolatedcumflow(flow, attime, tcur, tstep=time_step, inlinkno=0):
    def get_timeinterpolatedcumflow(flow, attime, tcur, tstep=time_step, inlinkno=0):
        # Addition to Yperman
        coef = divmod((attime - (tcur - tstep)), tstep)
        coef_1 = coef[0]
        coef_2 = coef[1] / tstep
        if tcur > 0. and coef_1 != 1.:
            return max(flow[tcur + coef_1 * tstep - tstep] + coef_2 * (
                flow[tcur + coef_1 * tstep] - flow[tcur + coef_1 * tstep - tstep]), 0.)
        else:
            return flow[tcur]

    # ===========================================================================================
    # Nested function that assigns the origin flow
    # def loadNormalNodes(net_nodes, origins, destinations):
    #     # update origin nodes
    #     normalNodes = list()
    #     for node in net_nodes:
    #         if origins.count(node) == 0 and destinations.count(node) == 0:
    #             normalNodes.append(node)
    #     return normalNodes
    #
    # # def loadOriginNodes(t):
    # #     pass
    #
    # def calculateSendingFlow(inlink, t):
    #     SF = inlink.capacity
    #     time = t - inlink.length / inlink.v0
    #     val = findCVN(inlink, time)
    #     return min(SF, val - inlink.cumflow_at_start[t])
    #
    # def calculateReceivingFlowFQ(outlink, t):
    #     RF = outlink.capacity
    #     time = timeSlices.get(t) - outlink.length / outlink.w
    #     find_cvn = findCVN(outlink, time)
    #     val = find_cvn + outlink.jamdensity * outlink.length
    #     if t == 0:
    #         return min(RF, val - outlink.cumflow_at_start.get(t))
    #     else:
    #         return min(RF, val - outlink.cumflow_at_start.get(t - time_step))
    #
    # def findCVN(inlink, time):
    #     if time <= timeSlices.get(0):
    #         return 0
    #     elif time >= timeSlices.get(int(time_horizon /time_step) - 1):
    #         return inlink.cumflow_at_end.get(int(time_horizon / time_step))
    #     else:
    #         t1 = round(time / time_step)
    #         t2 = t1 + 1
    #         return inlink.cumflow_at_start.get(t1) + (time / time_step - t1 + 1) * \
    #                                                  (inlink.cumflow_at_end.get(t2) - inlink.cumflow_at_end.get(t1))

    #============================================================================================
    # # for MISO model
    # def get_update(inlink_no, outlink_no):
    #     # update the set of unprocessed input links
    #     node.inlinks_set.remove(inlink_no)
    #     # update supply for next turn of output link
    #     for i in node.inlinks_set():
    #         net.turns[(i, outlink_no)].to_link.receiving_flow[t] -= \
    #             net.turns[(inlink_no, outlink_no)].transition_flows[t]


    # for MIMO model
    def get_update_mimo(inlink_no, outlink_no):
        # update the set of unprocessed input links
        # node.inlinks_set.remove(inlink_no)
        # update supply for next turn of output link
        for i in node.inlinks_set():
            net.turns[(i, outlink_no)].to_link.receiving_flow[t] -= \
                net.turns[(inlink_no, outlink_no)].transition_flows[t]


    def min_flow_reduction_factor(inlink_no, outlinks_set, t):
        alpha_min = {outlinks_set: net.turns[(inlink_no, outlinks_set)].to_link.receiving_flow[t] /
                                   net.turns[(inlink_no, outlinks_set)].from_link.sending_flow[t]
                     for outlinks_set in outlinks_set
                     if net.turns[(inlink_no, outlinks_set)].from_link.sending_flow[t] != 0}
        if len(alpha_min.values()) == 0:
            return 'zero_sending_flow'
        else:
            return {k: alpha_min[k] for k in alpha_min if alpha_min[k] == min(alpha_min.values())}


    def oriented_prio_set(node, t):
        sum_sending_flow = {[node.inlinks_set(), node.outlinks_set()]: sum(net.turns[(inlink, outlink)].sending_flow[t] for inlink in node.inlinks_set() for outlink in node.outlinks_set())}
        # sum_sending_flow = sum([node.turns[turn].sending_flow[t] for turn in node.turns])
        tilde_oriented_prio = {[inlink, outlink]: inlink_prio * sum_sending_flow / node.turns[turn].sending_flow[t] for turn in node.turns}


    def prio_fact(inlink_prio, t):
        sum_sending_flow = {(inlink, outlink): sum(net.turns[(inlink, outlink)].from_link.sending_flow[t]
                                                   for inlink in node.inlinks_set()
                                                   for outlink in node.outlinks_set())
                            for inlink in node.inlinks_set() for outlink in node.outlinks_set()}
        oriented_prio = {(inlink, outlink): inlink_prio.get(inlink) * net.turns[inlink, outlink].sending_flow[t]
                                            / sum_sending_flow.get((inlink, outlink))
                         for inlink in node.inlinks_set()
                         for outlink in node.outlinks_set()
                         if sum_sending_flow.get((inlink, outlink)) != 0}
        alpha_min = {(inlink, outlink): net.turns[(inlink.no, outlink.no)].to_link.receiving_flow[t] /
                                            oriented_prio.get((inlink.no, outlink.no))
                         for inlink in node.inlinks_set()
                         for outlink in node.outlinks_set()
                         if len(oriented_prio.values()) != 0}
        return oriented_prio,\
                   {(inlink, outlink): alpha_min[(inlink, outlink)] for inlink in alpha_min for outlink in alpha_min
                    if alpha_min[(inlink, outlink)] == min(alpha_min.values())}


    def flow_reduction_factors(inlink_no, turns_to_link, t):
        alpha_min = {turns_to_link: net.turns[(inlink_no, turns_to_link)].to_link.receiving_flow[t] /
                                    net.turns[(inlink_no,
     turns_to_link)].from_link.sending_flow[t]
                     for turns_to_link in turns_to_link
                     if net.turns[(inlink_no, turns_to_link)].from_link.sending_flow[t] != 0}
        if len(alpha_min) != 0:
            alpha_min = {k: alpha_min[k] for k in alpha_min if alpha_min[k] == min(alpha_min.values())}
            if len(alpha_min.values()) > 1:
                while len(alpha_min.values()) != 1:
                    alpha_min.popitem()
        else:
            alpha_min = {turns_to_link[0]: 1}
        return alpha_min

    # for MIMO model
    def inlink_contrib_tolink(node):
        U_j = dict()
        for j in node.outlinks:
            fromlink = [i.from_link for i in node.turns_tolink(j)] # if i.sending_flow[t] > 0]
            U_j.update({j: fromlink})
        return U_j

    # for MIMO model
    def get_min_j(min_a, min_val):
        for k, v in min_a.items():
            if v == min_val:
                return k

    # timeSlices = {i: i * time_step for i in range(0, int(time_horizon / time_step), 1)}
    total = int(time_horizon) - int(time_step)
    for t in range(0, int(time_horizon) + int(time_step), int(time_step)):
        printProgress(t, total, 'PYLTM time steps')
        # Node sequence is not important because of algorithm logic, the main reason is FIFO rule, which is also
        # ensured by the condition that time step should be smaller or the same as the smallest link travel time!
        # However, FIFO is not "strictly" ensured, see p. 25 and Node Models
        t = float(t)
        for node in net.nodes:
            assert isinstance(node, Node)
            # ----------------------------------------------------------------------------------------------
            # ------------------------------------- STEP 1. Link Model -------------------------------------
            # ----------------------------------------------------------------------------------------------
            # Step 1: For each incoming link i determine the Sending flow Si(t)
            # at the downstream link end (end node) and for each outgoing
            # link j determine the Receiving flow Rj(t) at the upstream link
            # end (start node)
            # TODO: Add piecewise linear fundamental diagram modification as alternative
            for inlink in node.inlinks:
                assert isinstance(inlink, Link)

                # find sending flow
                # inlink.sending_flow[t] = calculateSendingFlow(inlink, t)

                cumveh_i_atstart = get_timeinterpolatedcumflow(inlink.cumflow_at_start,
                                                               max(t + time_step - inlink.length / inlink.v0, 0.), t)
                cumveh_i_atend = inlink.cumflow_at_end[t]
                inlink.sending_flow[t] = min(cumveh_i_atstart - cumveh_i_atend, inlink.outcapacity)  # Changed to OUT!
            for outlink in node.outlinks:
                assert isinstance(outlink, Link)

                # find receiving flow
                # outlink.receiving_flow[t] = calculateReceivingFlowFQ(outlink, t)

                cumveh_j_atend = get_timeinterpolatedcumflow(outlink.cumflow_at_end,
                                                             max(t + time_step + outlink.length / outlink.w, 0.), t)
                jammedveh = outlink.jamdensity * outlink.length
                cumveh_j_atstart = outlink.cumflow_at_start[t]
                outlink.receiving_flow[t] = min(cumveh_j_atend + jammedveh - cumveh_j_atstart,
                                                outlink.capacity)


            # -------------------------- Sending Flows Disaggregated by Turn -------------------------
            # TODO: Can be done only for SOME nodes!
            for turn in node.turns:
                turn.sending_flow[t] = 0.
                for route in net.routes.get_routes_through_link(turn.from_link):
                    # Addition to Yperman
                    # ------------------------- START correct delta calc -------------------------
                    if 'link' in turn.from_link.maxflowdeltas:
                        linkcumflowdt = turn.from_link.maxflowdeltas['link']
                    else:
                        linkcumflowdt = turn.from_link.cumflow_at_start[t] - turn.from_link.cumflow_at_start[
                            max(t - time_step, 0.)]
                    if 'route' in turn.from_link.maxflowdeltas and route.no in turn.from_link.maxflowdeltas['route']:
                        routecumflowdt = turn.from_link.maxflowdeltas['route'][route.no]
                    else:
                        routecumflowdt = turn.from_link.cumrouteflow_at_start[(route.no, t)] - \
                                         turn.from_link.cumrouteflow_at_start[
                                             (route.no, max(t - time_step, 0.))]
                    # -------------------------- END correct delta calc --------------------------
                    turn.sending_flow[t] += (turn.to_link in route.path) * (
                        get_inversecumflow(turn.from_link.cumflow_at_end[t] + turn.from_link.sending_flow[t],
                                           turn.from_link.cumflow_at_start[max(t - time_step, 0.)],
                                           turn.from_link.cumrouteflow_at_start[(route.no, max(t - time_step, 0.))],
                                           linkcumflowdt, routecumflowdt) -
                        turn.from_link.cumrouteflow_at_end[(route.no, t)])


            # ----------------------------------------------------------------------------------------------------------
            # ------------------------- STEP 2. KURZHANSKIY NODE-MODEL -------------------------------------------------
            # ----------------------------------------------------------------------------------------------------------
            # ===================== Terminal nodes processing ==========================================================
            for turn in node.turns:
                # TODO: Remake origin and destination nodes - not only end or begin nodes but any!
                if turn.from_link is None:
                    # origin node
                    # *-
                    if t + time_step in node.Nr:
                        originflow = node.Nr[t + time_step]
                    else:
                        # Addition to Yperman
                        nodenr = max([key for key in node.Nr.keys() if not isinstance(key, tuple)])  # t == 5760
                        originflow = node.Nr[nodenr]
                    turn.transition_flows[t] = min(originflow - turn.to_link.cumflow_at_start[t],
                                                   turn.to_link.receiving_flow[t])

                elif turn.to_link is None:
                    # destination node
                    # -*
                    turn.transition_flows[t] = turn.from_link.sending_flow[t]

# ======================================================================================================================
#                 # elif len(node.turns_fromlink(turn.from_link)) == len(node.turns_tolink(turn.to_link)) == 1:
#                 elif len(node.inlinks) == len(node.outlinks) == 1:
#                     # TODO: Check for no crossing of turns also
#                     # inhomogeneous node
#                     # -*-
#                     turn.transition_flows[t] = min(turn.from_link.sending_flow[t],
#                                                    turn.to_link.receiving_flow[t])
#
#                 elif len(node.turns_fromlink(turn.from_link)) == 2 and len(node.turns_tolink(turn.to_link)) == 1:
#                     # elif len(node.inlinks) == 1 and len(node.outlinks) == 2:
#                     # TODO: Check for no crossing of turns also
#                     # diverge node - without perfect FIFO
#                     # -*<
#                     tsf = turn.sending_flow[t]
#                     ttf_list = [tsf]
#                     for turnj in node.turns_fromlink(turn.from_link):
#                         if turnj.sending_flow[t] != 0.:
#                             ttf_list.append(
#                                 tsf * turnj.to_link.receiving_flow[t] / turnj.sending_flow[t])
#                     turn.transition_flows[t] = min(ttf_list)
# ======================================================================================================================

            # ++++++++++++++++++++++++++++ MIMO Node Model +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # ---------------------------- M-N Node Model --------------------------------------------------------------
            # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            if len(node.inlinks) >= 1 and len(node.outlinks) >= 1:
            # if len(node.inlinks) == 2 and len(node.outlinks) == 1:
                # +++++++++++++++++++++++++++++++++++++++++++++
                for inlink in node.inlinks:                     #
                    sumflows_at_end = 0.0                       #
                    for turn in node.turns_fromlink(inlink):    #
                        turn.transition_flows[t] = 0            #
                # +++++++++++++++++++++++++++++++++++++++++++++
                # ======================================================================================================
                # STEP 1. Initialize
                # ======================================================================================================
                k = 0 # Iteration index
                tilde_R_j = { # Adjusted supply of link j at iteration k:
                    node.outlinks[i].no:
                        node.outlinks[i].receiving_flow[t] for i in range(len(node.outlinks_set()))
                    }
                # tilde_U_j = inlink_contrib_tolink(node) # TODO - delite it from def section
                tilde_U_j = dict()
                for j in node.outlinks: # output link
                    fromlink = list()
                    for i in [link_no.path for link_no in net.routes.get_routes_through_link(j)]: # nodes and paths of route are goes through link j
                        for pathlink in i: # network element of the route
                            for inlinks in [inl.from_link for inl in node.turns_tolink(j)]:
                                if inlinks.no == pathlink.no:
                                    fromlink.append(inlinks)
                    tilde_U_j.update({j: fromlink})

                # tilde_S_ij = dict() # Oriented demand

                # ======================================================================================================
                # STEP 3.
                # ======================================================================================================
                # 3.1. Check that at least one of the unprocessed input links has nonzero priority, otherwise,
                # assign equal positive priorities to all the unprocessed input links
                # +++++++ def inlinks_priority() +++++++++++++++++++++++++++++++
                tilde_p_i = dict()
                # set of prioroties for input links arranged into groups of output links:
                if len(node.inlinks) != 1:
                    for outlink in tilde_U_j.keys():
                        for inlink in tilde_U_j.get(outlink):
                            if inlink.priority != 0:
                                tilde_p_i.update({inlink.no: inlink.priority})
                            else: # otherwise set up: 1 / |Union j in V(k) for Uj(k)|
                                tilde_p_i.update({inlink.no: 1.0 / len(tilde_U_j.get(outlink))})
                else: # if only one input link:
                    if inlink.priority != 0:
                        tilde_p_i.update({inlink.no: inlink.priority})
                    else:
                        tilde_p_i.update({inlink.no: 1.0})

                # 3.2. For each output link j ∈ V (k) and input link i ∈ Uj(k) compute oriented priority:
                # ~p_ij(k) = ˜p_i(k) * betta_ij
                #++++++++ def oriented_priority() ++++++++++++++++++++++++++++++
                tilde_p_ij = dict()
                for outlink in tilde_U_j.keys():
                    for inlink in tilde_U_j.get(outlink):
                        if inlink.sending_flow[t] != 0:
                            betta_ij = net.turns[(inlink.no, outlink.no)].sending_flow[t] / inlink.sending_flow[t]
                            tilde_p_ij.update({(inlink.no, outlink.no): tilde_p_i.get(inlink.no) * betta_ij})
                        else: # if turns.sending_flow == 0:
                            tilde_p_ij.update({(inlink.no, outlink.no): 0})

                # ======================================================================================================
                # STEP 4. For each j ∈ V(k), compute factors and find the smallest of these factors:
                # ======================================================================================================
                # 4.1. Compute the factors:
                #+++++++++ def min_factors() +++++++++++++++++++++++++++++++++++
                a_j = dict() # dict of a_ij factors: {j: sum_[i in U_j](a_ij)}
                for outlink in tilde_U_j.keys():
                    tilde_p_sum_i_Uj = sum(tilde_p_ij[inlink_no.no, outlink.no] for inlink_no in tilde_U_j.get(outlink))
                    if tilde_p_sum_i_Uj != 0:
                        a_j.update({outlink: tilde_R_j.get(outlink.no) / tilde_p_sum_i_Uj})
                    else:
                        a_j.update({outlink: 0})

                # 4.2. Find the smallest of the factors:
                min_a_j = min(a_j.values())
                j_min_a_j = get_min_j(a_j, min_a_j) # find the number of output link for minimum a_j

                # ======================================================================================================
                # STEP 5. Define the set of input links, whose demand does not exceed the allocated supply:
                # ======================================================================================================
                for outlink in tilde_U_j.keys():
                    for inlink in tilde_U_j.get(outlink):
                        # ------------if ~U(k) is NOT empty (free flow)-------------------------------------------------
                        if net.turns[(inlink.no, outlink.no)].from_link.sending_flow[t] <= tilde_p_i.get(inlink.no) * min_a_j:
                            # demand is satisfied
                            net.turns[(inlink.no, outlink.no)].transition_flows[t] = net.turns[(inlink.no, outlink.no)].from_link.sending_flow[t]
                            # +++++++ def get_update() +++++++++++++++++++++++++++++++++++++
                            # get_update_mimo(inlink.no, outlink.no):
                            tilde_R_j.update({outlink.no: tilde_R_j.get(outlink.no) - net.turns[(inlink.no, outlink.no)].transition_flows[t]})
                        # -----------(congested flow)-------------------------------------------------------------------
                        else:
                            try:
                                net.turns[(inlink.no, outlink.no)].transition_flows[t] = net.turns[(inlink.no, outlink.no)].from_link.sending_flow[t] - net.turns[(inlink.no, outlink.no)].from_link.sending_flow[t] * (1 - tilde_p_ij.get((inlink.no, j_min_a_j.no)) * min_a_j / net.turns[(inlink.no, j_min_a_j.no)].from_link.sending_flow[t])
                            except TypeError:
                                print t
                            # +++++++ get_update_mimo(inlink.no, outlink.n+++++++++++++++++++
                            new_tilde_R_j = tilde_R_j.get(outlink.no) - net.turns[(inlink.no, outlink.no)].transition_flows[t]
                            tilde_R_j.update({outlink.no: new_tilde_R_j})
                            # ========= TEMP ===========================================================================
                            if new_tilde_R_j < 0:
                                print (inlink.no, outlink.no, t, new_tilde_R_j)
                # STEP 6. Set k:= k + 1, and return to step 2.
                k += 1
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            # # +++++++++++ SIMO Node Model (for deverge and inhomogeneous nodes) ++++++++++++++++++++++++++++++++++++++
            # # -*< (diverge node) or -*- (inhomogeneous node)
            # # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # if len(node.inlinks) == 1 and len(node.outlinks) == 2 or len(node.inlinks) == len(node.outlinks) == 1:
            #     outlinks_set = node.outlinks_set()  # set of unprocessed outnput links
            #     for inlink in node.inlinks: # the only input link for SIMO
            #         alpha_min_j = min_flow_reduction_factor(inlink.no, outlinks_set, t)
            #         for outlink in node.outlinks: # output links' cycle (the base calculation cycle)
            #             # ------------(free flow)---------------------------------------------------------------------
            #             if alpha_min_j == 'zero_sending_flow':
            #                 net.turns[(inlink.no, outlink.no)].transition_flows[t] = 0
            #             elif alpha_min_j >= 1:
            #                 # demand is satisfied
            #                 net.turns[(inlink.no, outlink.no)].transition_flows[t] = \
            #                    net.turns[(inlink.no, outlink.no)].from_link.sending_flow[t]
            #             # -----------(congested flow)-----------------------------------------------------------------
            #             else:
            #                 # tilde_S1_j_ast = net.turns[(inlink.no, outlink.no)].from_link.sending_flow[t] * \
            #                 #                  sum(alpha_min_j.values())
            #                 net.turns[(inlink.no, outlink.no)].transition_flows[t] = \
            #                     net.turns[(inlink.no, outlink.no)].from_link.sending_flow[t] - \
            #                     net.turns[(inlink.no, outlink.no)].from_link.sending_flow[t] * \
            #                     (1 - sum(alpha_min_j.values()))
            #
            #
            #
            # # +++++++++++ MISO Node Model (for merge and inhomogeneous nodes) ++++++++++++++++++++++++++++++++++++++++
            # # # >*- or -*- (inhomogeneous node)
            # # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            # # TODO: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            # # TODO: Add check of links occupancy and if one of them is empty do something (change priorities or???)
            # # TODO: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            # if len(node.inlinks) == 2 and len(node.outlinks) == 1 or len(node.inlinks) == len(node.outlinks) == 1:
            #     inlinks_set = node.inlinks_set()  # set of unprocessed input links
            #     for inlinkprio, inlink in node.inlinks_sorted_prio():  # sorted by priority input links' cycle
            #         for outlink in node.outlinks:  # output links' cycle
            #             # ------------if ~U(k) is NOT empty (free flow)-----------------------------------------------
            #             if net.turns[(inlink.no, outlink.no)].from_link.sending_flow[t] <= \
            #                                     inlinkprio / node.inlink_sum_prio() * \
            #                             net.turns[(inlink.no, outlink.no)].to_link.receiving_flow[t]:
            #                 # demand is satisfied
            #                 net.turns[(inlink.no, outlink.no)].transition_flows[t] = \
            #                     net.turns[(inlink.no, outlink.no)].from_link.sending_flow[t]
            #                 get_update(inlink.no, outlink.no)
            #             # -----------if ~U(k) IS empty (congested flow)-----------------------------------------------
            #             else:
            #                 # only some of demand is satisfied
            #                 if inlinkprio > 0:  # we have user-defined priorities
            #                     net.turns[(inlink.no, outlink.no)].transition_flows[t] = \
            #                         inlinkprio / node.inlink_sum_prio() * \
            #                         net.turns[inlink.no, outlink.no].to_link.receiving_flow[t]
            #                     get_update(inlink.no, outlink.no)
            #                 else:  # if (inlink.priority == 0) we need to get normalized priorities
            #                     net.turns[(inlink.no, outlink.no)].transition_flows[t] = \
            #                         node.norm_prio() * \
            #                         net.turns[inlink.no, outlink.no].to_link.receiving_flow[t]
            #                     get_update(inlink.no, outlink.no)

            # ----------------------------------------------------------------------------------------------
            # ------------------------------------- STEP 2. Node Model -------------------------------------
            # ----------------------------------------------------------------------------------------------
            # Step 2: Determine the transition flows Gij(t) from incoming
            # links i to outgoing links j, i.e. determine which parts of the
            # Sending and Receiving flows can actually be sent and received.
            # TODO: Capacity + delay for turns -> remember FIFO rule, see Yperman, remember sneaking in Optima

            # for turn in node.turns:
            #     # TODO: Remake origin and destination nodes - not only end or begin nodes but any!
            #     if turn.from_link is None:
            #         # origin node
            #         # *-
            #         if t + time_step in node.Nr:
            #             originflow = node.Nr[t + time_step]
            #         else:
            #             # Addition to Yperman
            #             nodenr = max([key for key in node.Nr.keys() if not isinstance(key, tuple)])
            #             originflow = node.Nr[nodenr]
            #         turn.transition_flows[t] = min(originflow - turn.to_link.cumflow_at_start[t],
            #                                        turn.to_link.receiving_flow[t])
            #
            #     elif turn.to_link is None:
            #         # destination node
            #         # -*
            #         turn.transition_flows[t] = turn.from_link.sending_flow[t]
            #
            #     elif len(node.turns_fromlink(turn.from_link)) == len(node.turns_tolink(turn.to_link)) == 1:
            #         # elif len(node.inlinks) == len(node.outlinks) == 1:
            #         # TODO: Check for no crossing of turns also
            #         # inhomogeneous node
            #         # -*-
            #         turn.transition_flows[t] = min(turn.from_link.sending_flow[t],
            #                                        turn.to_link.receiving_flow[t])
            #
            #     elif len(node.turns_fromlink(turn.from_link)) == 2 and len(node.turns_tolink(turn.to_link)) == 1:
            #         # elif len(node.inlinks) == 1 and len(node.outlinks) == 2:
            #         # TODO: Check for no crossing of turns also
            #         # diverge node - without perfect FIFO
            #         # -*<
            #         tsf = turn.sending_flow[t]
            #         ttf_list = [tsf]
            #         for turnj in node.turns_fromlink(turn.from_link):
            #             if turnj.sending_flow[t] != 0.:
            #                 ttf_list.append(
            #                     tsf * turnj.to_link.receiving_flow[t] / turnj.sending_flow[t])
            #         turn.transition_flows[t] = min(ttf_list)
            #
            #         # if turn.from_link.sending_flow[t] != 0.:
            #         #     dmd_proportional_flow = turn.from_link.sending_flow[t] * turn.to_link.receiving_flow[t] / sum(
            #         #         [inlink.sending_flow[t] for inlink in node.inlinks])
            #         #
            #         #     turn.transition_flows[t] = min(turn.from_link.sending_flow[t],
            #         #                                    dmd_proportional_flow)
            #         # else:
            #         #     turn.transition_flows[t] = 0.
            #
            #     else:
            #         # cross (urban) node
            #         # >*<
            #         if node.signalized:
            #             # TODO: Signalized Cross Node calculation
            #             turn.receiving_flow[t] = turn.capa / sum(
            #                 [turni.capa for turni in node.turns_tolink(turn.to_link)])
            #             ttf_list = [1.]
            #             for turnj in node.turns_fromlink(turn.from_link):
            #                 if turnj.sending_flow[t] != 0.:
            #                     ttf_list.append(turnj.capa / turnj.sending_flow[t])
            #                     ttf_list.append(turnj.receiving_flow[t] / turnj.sending_flow[t])
            #             turn.transition_flows[t] = min(ttf_list) * turn.sending_flow[t]
            #         else:
            #             # TODO: UNSignalized Cross Node calculation
            #             pass

            # ----------------------------------------------------------------------------------------------
            # ---------------------------------- STEP 3. Cumulative Flows ----------------------------------
            # ----------------------------------------------------------------------------------------------
            # Step 3: For the downstream link boundary (at the end node) of each incoming link i and for the upstream
            #  link boundary (at the start node) of each outgoing link j update the cumulative vehicle numbers N(x,t)
            for inlink in node.inlinks:
                sumflows_at_end = 0.0
                for turn in node.turns_fromlink(inlink):
                    sumflows_at_end += turn.transition_flows[t]
                inlink.cumflow_at_end[t + time_step] = inlink.cumflow_at_end[t] + sumflows_at_end
            for outlink in node.outlinks:
                sumflows_at_start = 0.0
                for turn in node.turns_tolink(outlink):
                    sumflows_at_start += turn.transition_flows[t]

                # # ======================= DETECTOR DATA ===================================
                # if t in outlink.flow_det:
                #     sumflows_at_start = outlink.flow_det[t]
                # # ======================= DETECTOR DATA ===================================

                outlink.cumflow_at_start[t + time_step] = outlink.cumflow_at_start[t] + sumflows_at_start

            # --------------------------- Cumulative Flows Disaggregated by Route --------------------------
            # TODO: I think that we can calc cumulative disaggregated flows by route after each cycle of nodes - check
            # TODO: Cumulative route flows - maybe transition flows by routes?
            # Assumption. Flows by route at start of link are in the same proportion as at start of previous link
            for outlink in node.outlinks:
                sumrouteflows = 0.0
                for route in net.routes.get_routes_through_link(outlink):
                    for turn in node.turns_tolink(outlink):
                        if turn.from_link in route.path:
                            sumrouteflows += turn.from_link.cumrouteflow_at_start[(route.no, t)]
                for route in net.routes.get_routes_through_link(outlink):
                    sumrouteflows_at_start = 0.0
                    for turn in node.turns_tolink(outlink):
                        # Addition to Yperman. Yperman didn't take into account (maybe for some reason...) that
                        # we should take transition flows disaggregated by routes??? very strange
                        if turn.via_node == route.path[0]:
                            if t + time_step in node.Nr:
                                nodenr = t + time_step
                            else:
                                nodenr = max([key for key in node.Nr.keys() if not isinstance(key, tuple)])  # t == 5760
                            originflow = node.Nr[nodenr]
                            sumrouteflows_at_start += turn.transition_flows[t] * (
                                turn.via_node.Nr[(route.no, nodenr)] / originflow)
                        elif turn.from_link in route.path:
                            if turn.from_link.cumflow_at_start[t] == 0.:
                                sumrouteflows_at_start = 0.
                            else:
                                sumrouteflows_at_start += turn.transition_flows[t] * (
                                    turn.from_link.cumrouteflow_at_start[(route.no, t)] /
                                    sumrouteflows)
                    outlink.cumrouteflow_at_start[(route.no, t + time_step)] = outlink.cumrouteflow_at_start[(
                        route.no, t)] + sumrouteflows_at_start

                    # Addition to Yperman. Didn't take into account that at some time cumflows are the same!
                    # ------------------------- START correct delta calc -------------------------
                    tmaxflow = max(route.flow_attime.keys())
                    if abs(outlink.cumrouteflow_at_start[(route.no, t + time_step)] - route.flow_attime[
                        tmaxflow]) < 0.001:
                        if 'link' not in outlink.maxflowdeltas:
                            outlink.maxflowdeltas['link'] = outlink.cumflow_at_start[t + time_step] - \
                                                            outlink.cumflow_at_start[t]
                        if 'route' not in outlink.maxflowdeltas:
                            outlink.maxflowdeltas['route'] = {
                                route.no: outlink.cumrouteflow_at_start[(route.no, t + time_step)] -
                                          outlink.cumrouteflow_at_start[(route.no, t)]}
                        else:
                            if route.no not in outlink.maxflowdeltas['route']:
                                outlink.maxflowdeltas['route'][route.no] = outlink.cumrouteflow_at_start[
                                                                               (route.no, t + time_step)] - \
                                                                           outlink.cumrouteflow_at_start[(route.no, t)]
                                # -------------------------- END correct delta calc --------------------------

        for link in net.links:
            for route in net.routes.get_routes_through_link(link):
                # Addition to Yperman
                # ------------------------- START correct delta calc
                if 'link' in link.maxflowdeltas:
                    linkcumflowdt = link.maxflowdeltas['link']
                else:
                    linkcumflowdt = link.cumflow_at_start[t + time_step] - \
                                    link.cumflow_at_start[t]
                if 'route' in link.maxflowdeltas and route.no in link.maxflowdeltas['route']:
                    routecumflowdt = link.maxflowdeltas['route'][route.no]
                else:
                    routecumflowdt = link.cumrouteflow_at_start[(route.no, t + time_step)] - \
                                     link.cumrouteflow_at_start[(route.no, t)]
                # -------------------------- END correct delta calc --------------------------
                link.cumrouteflow_at_end[(route.no, t + time_step)] = get_inversecumflow(
                    link.cumflow_at_end[t + time_step],
                    link.cumflow_at_start[t],
                    link.cumrouteflow_at_start[(route.no, t)],
                    linkcumflowdt, routecumflowdt)

        # ----------------------------------------------------------------------------------------------
        # ------------------------------------ Macro Characteristics -----------------------------------
        # ----------------------------------------------------------------------------------------------
        # TODO: Make this as separate function and call only for aggregated time steps
        for link in net.links:
            link.vehonlink[t] = link.cumflow_at_start[t] - link.cumflow_at_end[t]
            link.density[t] = link.vehonlink[t] / (link.numlanes * link.length)
            link.flow[t] = link.cumflow_at_start[t + time_step] - link.cumflow_at_start[t]
            link.outflow[t] = link.cumflow_at_end[t + time_step] - link.cumflow_at_end[t]
            if link.density[t] > link.critdensity:
                if link.flow[t] > 0.:
                    linkincurspeed = link.outflow[t] / link.density[t] * (3600. / t)
                else:
                    linkincurspeed = link.v0 * 3600.
                if link.outflow[t] > 0.:
                    linkoutcurspeed = link.outflow[t] / link.density[t] * (3600. / t)
                else:
                    linkoutcurspeed = link.v0 * 3600.
                link.curspeed[t] = (linkincurspeed + linkoutcurspeed) / 2.
            else:
                link.curspeed[t] = link.v0 * 3600.
            link.queue[t] = get_timeinterpolatedcumflow(link.cumflow_at_start,
                                                        max(t + time_step - link.length / link.v0, 0.), t) - \
                            link.cumflow_at_end[t + time_step]
