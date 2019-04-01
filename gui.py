# -*- coding: utf-8 -*-


import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.axes import Axes
from matplotlib import transforms
import numpy as np
from zoompanplot import ZoomPan
# import pandas as pnd


class MapDrawer(object):
    def __init__(self, net, t_horizon, t):
        # type: (netmodel.Net, float, float) -> None
        # TODO: Make aggregation functions for time
        # TODO: Change everything to arrays in numpy and pandas
        self.net = net
        self.t_horizon = t_horizon
        self.t = t
        self.tcur = 0.

        self.links_bars = {}

    def plot_net(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        self.show_bars(ax)

        def onpick(event):
            # type: (matplotlib.backend_bases.PickEvent) -> None
            figobject = event.artist
            if event.mouseevent.button == 1:
                self.plot_link_diags(self.net.links[figobject.get_gid()], self.t)

        def onkeypress(event):
            # type: (matplotlib.backend_bases.KeyEvent) -> None
            # TODO: Add animation
            if event.key == 'w':
                self.tcur += self.t
                self.tcur = min(self.tcur, self.t_horizon)
            elif event.key == 'q':
                self.tcur -= self.t
                self.tcur = max(self.tcur, 0.)
            elif event.key == 'r':
                self.tcur += self.t * 6
                self.tcur = min(self.tcur, self.t_horizon)
            elif event.key == 'e':
                self.tcur -= self.t * 6
                self.tcur = max(self.tcur, 0.)
            self.update_bars(ax)
            fig.canvas.draw()
            fig.canvas.flush_events()

        fig.canvas.mpl_connect('pick_event', onpick)
        fig.canvas.mpl_connect('key_press_event', onkeypress)
        fig.canvas.mpl_connect('motion_notify_event', fig.canvas.onHilite)

        scale = 1.2
        zp = ZoomPan()
        zp.zoom_factory(ax, base_scale=scale)
        zp.pan_factory(ax)

        plt.show()

    def plot_link_diags(self, link, t):
        # TODO: Deactivate elements by clicking on legend elements
        receiving_flows = np.array(sorted([(i[0], i[1] * (3600.0 / t)) for i in link.receiving_flow.items()]))
        sending_flows = np.array(sorted([(i[0], i[1] * (3600.0 / t)) for i in link.sending_flow.items()]))
        inflows = np.array(sorted([(i[0], i[1] * (3600.0 / t)) for i in link.flow.items()]))
        outflows = np.array(sorted([(i[0], i[1] * (3600.0 / t)) for i in link.outflow.items()]))
        cumflows_in = np.array(sorted(link.cumflow_at_start.items()))
        cumflows_out = np.array(sorted(link.cumflow_at_end.items()))
        vehonlink = np.array(sorted(link.vehonlink.items()))
        density = np.array(sorted(link.density.items()))
        flow = np.array(sorted(link.flow.items()))
        queue = np.array(sorted(link.queue.items()))
        curspeed = np.array(sorted(link.curspeed.items()))
        maxjamdensity = np.array(sorted([(time, link.jamdensity) for time, value in link.density.items()]))
        incapa = np.array(sorted([(time, link.capacity) for time, value in link.density.items()]))
        outcapa = np.array(sorted([(time, link.outcapacity) for time, value in link.density.items()]))
        jamstorage = np.array(sorted([(time, link.jamstorage) for time, value in link.density.items()]))
        # plot with various axes scales
        fig = plt.figure()
        fig.suptitle(
            'Link # {} from node # {} to node # {}. L={}. Lanes={}'.format(link.no, link.from_node.no, link.to_node.no,
                                                                           link.length, link.numlanes), fontsize=16)
        gs = gridspec.GridSpec(2, 2)

        ax1 = fig.add_subplot(gs[0, :])
        self.plot_funddiag(*link.get_params_for_hour(t), ax=ax1)

        ax2 = fig.add_subplot(gs[1, 0])
        ax2.plot(*receiving_flows.transpose(), linestyle='--', label='receiving flows', linewidth=3)
        ax2.plot(*sending_flows.transpose(), linestyle='--', label='sending flows', linewidth=3)
        ax2.plot(*cumflows_in.transpose(), label='cumflows IN')
        ax2.plot(*cumflows_out.transpose(), label='cumflows OUT')
        ax2.plot(*inflows.transpose(), label='in flows')
        ax2.plot(*outflows.transpose(), label='out flows')
        ax2.legend(fontsize='small')
        ax2.grid(True)

        ax3 = fig.add_subplot(gs[1, 1])
        ax3.plot(*vehonlink.transpose(), label='vehonlink')
        ax3.plot(*density.transpose(), label='density')
        ax3.plot(*flow.transpose(), label='flow')
        ax3.plot(*queue.transpose(), label='queue')
        ax3.plot(*maxjamdensity.transpose(), label='max jam density')
        ax3.plot(*incapa.transpose(), linestyle='--', label='in capacity', linewidth=3)
        ax3.plot(*outcapa.transpose(), linestyle='--', label='out capacity', linewidth=3)
        ax3.plot(*jamstorage.transpose(), linestyle='--', label='jamstorage', linewidth=3)
        ax3.legend(fontsize='small', loc='upper left')
        ax3.grid(True)
        ax31 = ax3.twinx()
        ax31.plot(*curspeed.transpose(), linestyle='-', color='y', label='curspeed')
        ax31.legend(fontsize='small', loc='upper right')
        ax31.grid(True)

        scale = 1.1
        zp = ZoomPan()
        zp.zoom_factory(ax1, base_scale=scale)
        zp.pan_factory(ax1)
        zp.zoom_factory(ax2, base_scale=scale)
        zp.pan_factory(ax2)
        # TODO: if we have twin axes the zoom and pan don't work!
        zp.zoom_factory(ax3, base_scale=scale)
        zp.pan_factory(ax3)

        fig.show()

    def plot_funddiag(self, v0=60.0, w=-20.0, capa=2000.0, ax=None):
        # type: (float, float, float, matplotlib.axes.Axes) -> None
        """
        Plot fundamental diagram: x=density, y=flow
        :param v0: speed - km per hour
        :param w: backwavespeed - km per hour
        :param capa: capacity - vehicles per hour
        :param ax: axes object to plot on
        :return: None
        """

        def k(qcur, freeflowspeed):
            # k - density, q - flow
            return qcur / freeflowspeed

        def kjam(qcur, capacity, freeflowspeed, backwavespeed):
            b = capacity / freeflowspeed - capacity / backwavespeed  # coefficient for 2nd part of diagram
            return qcur / backwavespeed + b

        q = 0.0

        yvalues = []
        xvalues = []
        while q < capa:
            yvalues.append(q)
            xvalues.append(k(q, v0))
            q += 10.0
            q = min(q, capa)
        while q > 0.0:
            yvalues.append(q)
            xvalues.append(kjam(q, capa, v0, w))
            q -= 10.0

        if ax is None:
            plt.plot(xvalues, yvalues)
            plt.show()
        elif isinstance(ax, Axes):
            ax.plot(xvalues, yvalues, linewidth=3)
            ax.set_title('Fundamental Diagram (v0={:.2f} '
                         'w={:.2f} capa={:.2f} '
                         'kcrit={:.2f} kjam={:.2f})'.format(v0, w, capa, k(capa, v0), kjam(0, capa, v0, w)))
            ax.grid(True)

    def show_bar_label(self, x1, y1, x2, y2):
        # type: (float, float, float, float) -> r1, va, ha, pos
        pos = [(x1 + x2) / 2., (y1 + y2) / 2.]
        v1 = (x2 - x1, y2 - y1)
        v2 = (1., 0.)
        r1 = np.rad2deg(np.arctan2(np.linalg.norm(np.cross(v1, v2)), np.dot(v1, v2)))
        if 0. < r1 < 90.:
            va = 'top'
            ha = 'center'
        elif r1 == 90.:
            va = 'center'
            ha = 'left'
        elif r1 == 0.:
            va = 'top'
            ha = 'center'
        else:
            r1 = 90. - r1
            va = 'bottom'
            ha = 'center'
        return r1, va, ha, pos

    def update_bars(self, ax):
        maxvalue = self.net.links.max_jamstorage()
        # maxwidth and linkwidth etc. should take into account coordinates or use somehow points
        maxwidth = 10.
        maxwidthq = maxwidth * 0.2
        for l in self.net.links:
            curbar = self.links_bars[l.no]
            maxqueue = l.jamstorage
            maxvalue1 = l.capacity
            maxvalue2 = l.outcapacity
            queue = l.queue
            valuetcur = l.vehonlink[self.tcur]
            value1tcur = l.flow[self.tcur]
            value2tcur = l.outflow[self.tcur]
            queuetcur = queue[self.tcur]
            linewidth = valuetcur / maxvalue * maxwidth
            linewidth1 = value1tcur / maxvalue1 * maxwidth
            linewidth2 = value2tcur / maxvalue2 * maxwidth
            rel_queue = queuetcur / maxqueue
            speedtcur = l.curspeed[self.tcur]
            relspeed = speedtcur / (l.v0 * 3600.0)
            if round(relspeed, 2) <= 0.2:
                vcolor = 'r'
            elif round(relspeed, 2) <= 0.4:
                vcolor = 'y'
            elif round(relspeed, 2) <= 0.8:
                vcolor = 'g'
            elif round(relspeed, 2) <= 1.:
                vcolor = 'b'
            else:
                vcolor = 'black'
                print('Speed > 1.: current speed = {}, v0 = {}, link {}'.format(speedtcur, l.v0 * 3600, l.no))
            # TODO: Should do that for each line segment
            x1, y1, x2, y2 = l.geom[0][0], l.geom[0][1], l.geom[1][0], l.geom[1][1]
            doffset = 0.
            dx, dy = self.get_linedxdy(x1, y1, x2, y2)
            verts = self.parallel_bar(x1, x2, y1, y2, dx, dy, doffset, linewidth, linewidth1, linewidth2)
            curbar['parallel'].set_xy(verts)
            curbar['parallel'].set_facecolor(vcolor)
            llx, lly = (verts[2][0] + verts[3][0]) / 2., (verts[2][1] + verts[3][1]) / 2.
            curbar['linelabel'].set_position((llx, lly))
            curbar['linelabel'].set_text('{:.0f}'.format(valuetcur))
            # verts = self.perpendicular_bar(x1, x2, y1, y2, dx, dy, doffset, linewidth, maxwidthq, rel_queue)
            x1, x2, y1, y2 = verts[2][0], verts[3][0], verts[2][1], verts[3][1]
            verts = self.parallel_bar(x1, x2, y1, y2, dx, dy, doffset, rel_queue * maxwidthq, linewidth1, linewidth2)
            curbar['perpendicular'].set_xy(verts)
        ax.get_figure().suptitle(
            'Current time: {0} s; {1} min. Time step: {2}. Time horizon: {3}'.format(self.tcur, self.tcur / 60., self.t, self.t_horizon),
            fontsize=16)

    def show_bars(self, ax):
        """
        :type ax: Axes
        :param ax:
        :return:
        """
        for l in self.net.links:
            self.links_bars[l.no] = {}
            x1, y1, x2, y2 = l.geom[0][0], l.geom[0][1], l.geom[1][0], l.geom[1][1]
            verts = [(0., 0.), (0., 0.), (0., 0.), (0., 0.)]
            line, = ax.plot((x1, x2), (y1, y2), picker=5, gid=l.no, color='black', linewidth=2., alpha=0.5)
            patch = plt.Polygon(verts, lw=1, alpha=0.7)
            ax.add_patch(patch)
            self.links_bars[l.no]['parallel'] = patch
            r1, va, ha, pos = self.show_bar_label(x1, y1, x2, y2)
            dx, dy = self.get_linedxdy(x1, y1, x2, y2)
            offset = 0.
            self.links_bars[l.no]['linelabel'] = ax.text(x=pos[0] + dx * offset, y=pos[1] + dy * offset,
                                                         s='{:.0f}'.format(0.), va=va, ha=ha,
                                                         size=9., rotation=r1)  #, rotation_mode='anchor'
            patch = plt.Polygon(verts, facecolor='m', lw=1, alpha=0.7)
            ax.add_patch(patch)
            self.links_bars[l.no]['perpendicular'] = patch

        ax.axis('square')
        xmin, xmax, ymin, ymax = ax.axis()
        ax.axis(xmin=xmin - 5., xmax=xmax + 5., ymin=ymin - 5., ymax=ymax + 5.)
        ax.get_figure().suptitle(
            'Current time: {0} s; {1} min. Time step: {2}. Time horizon: {3}'.format(self.tcur, self.tcur / 60., self.t, self.t_horizon),
            fontsize=16)

    def get_linedxdy(self, x1, y1, x2, y2):
        a = x1 - x2
        b = y1 - y2
        return a, -b

    def parallel_bar(self, x1, x2, y1, y2, dx, dy, doffset, linewidth, flowinw, flowoutw):
        # xi = (x1 + x2)
        # yi = (y1 + y2)
        verts = [
            (x1 + dy * doffset, y1 + dx * doffset),  # base line P1
            (x2 + dy * doffset, y2 + dx * doffset),  # base line P2
            # (x2 + dy * (doffset + flowoutw), y2 + dx * (doffset + flowoutw)),  # flow out P3
            # (x2 + dy * doffset, y2 + dx * doffset),  # inter flow out vehonlink P4
            # (xi / 1.1 + dy * (doffset + linewidth), yi / 1.1 + dx * (doffset + linewidth)),  # vehonlink P5
            # (xi / 2. + dy * (doffset + linewidth), yi / 2. + dx * (doffset + linewidth)),  # vehonlink P6
            (x2 + dy * (doffset + linewidth), y2 + dx * (doffset + linewidth)),
            (x1 + dy * (doffset + linewidth), y1 + dx * (doffset + linewidth))
            # (x2 + dy * doffset, y2 + dx * doffset),  # inter flow in vehonlink P7
            # (x1 + dy * (doffset + flowinw), y1 + dx * (doffset + flowinw))  # flow in P8
        ]
        return verts

    def perpendicular_bar(self, x1, x2, y1, y2, dx, dy, doffset, linewidth, maxwidthq, rel_queue):
        # if k is None:
        #     dx2 = 0.
        #     dy2 = dy * (y2 - y1) * (1. - rel_queue)
        # else:
        #     dx2 = dy * (x2 - x1) * (1. - rel_queue)
        #     dy2 = k * dy * (x2 - x1) * (1. - rel_queue)
        dx2 = dx * (x1 - x2) * (1. - rel_queue)
        dy2 = -dy * (y1 - y2) * (1. - rel_queue)
        verts = [
            (x2 + dy * (doffset + linewidth), y2 + dx * (doffset + linewidth)),
            (x2 + dy * (doffset + linewidth + maxwidthq), y2 + dx * (doffset + linewidth + maxwidthq)),
            (x1 + dy * (doffset + linewidth + maxwidthq) + dx2,
             y1 + dx * (doffset + linewidth + maxwidthq) + dy2),
            (x1 + dy * (doffset + linewidth) + dx2, y1 + dx * (doffset + linewidth) + dy2)
        ]
        return verts
