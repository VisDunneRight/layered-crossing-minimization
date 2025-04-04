from src import graph
import math
# import altair as alt
import os
try:
    import cairo
    import cairosvg
except OSError:
    print("Need Cairo installed on machine in order to draw graphs")
except ModuleNotFoundError:
    print("Need Cairo installed on machine in order to draw graphs")
import shutil
# from altair_saver import save


def bezier_control_points(g, eid, y1, y2, l1, l2, minl, nxdist, nydist, offset, left_edges, right_edges, straight_edges, long_edges, less_curvy_edges):
    if eid in long_edges:
        les = g.get_long_edges()
        the_le = [v for v in les if v[0] == eid[0] and v[1] == eid[1]][0]
        y2, l2 = g[the_le[-1]].y, g[the_le[-1]].layer
    cp1x = (l1 - 1 - minl) * nxdist + offset + nxdist // 1.5
    cp1y = y1 * nydist + offset
    if eid in left_edges or eid in straight_edges:
        cp1x = (l1 - 1 - minl) * nxdist + offset
    if eid in less_curvy_edges:
        cp1x -= nxdist // 3.5
    cp2x = (l2 - 1 - minl) * nxdist + offset - nxdist // 1.5
    cp2y = y2 * nydist + offset
    if eid in right_edges or eid in straight_edges:
        cp2x = (l2 - 1 - minl) * nxdist + offset
    if eid in less_curvy_edges:
        cp2x += nxdist // 3.5
    endx = (l2 - 1 - minl) * nxdist + offset
    endy = y2 * nydist + offset
    return cp1x, cp1y, cp2x, cp2y, endx, endy


def draw_graph(g: graph.LayeredGraph, svg_name, node_x_distance=150, node_y_distance=100, nested=False, motif=False, groups=None, emphasize_nodes=None, emphasize_edges=None, gravity=False, edge_thickness=False, label_nodes=True, as_png=False, color_scale=None, copies=1, fix_height=-1, remove_witespace=True, straighten_edges=False, node_weight_size=False, dot_group_anchors=False, label_anchors=False, left_straight_edges=None, right_straight_edges=None, full_straight_edges=None, straighten_only_true_edges=False, full_straight_long_edges=None, dont_draw_edges=None, less_curvy_edges=None):
    if left_straight_edges is None:
        left_straight_edges = []
    if right_straight_edges is None:
        right_straight_edges = []
    if full_straight_edges is None:
        full_straight_edges = []
    if full_straight_long_edges is None:
        full_straight_long_edges = []
    if dont_draw_edges is None:
        dont_draw_edges = []
    if less_curvy_edges is None:
        less_curvy_edges = []
    if nested:
        if "Images" not in os.listdir(".."):
            os.mkdir("../Images")
    elif "Images" not in os.listdir():
        os.mkdir("Images")
    offset = 40
    node_radius = 15
    line_width = 4
    font_size = 12
    palette = [(189, 189, 189), (56, 146, 201), (17, 138, 89), (254, 224, 139), (158, 1, 66), (253, 174, 97), (102, 194, 165), (213, 62, 79), (230, 145, 152), (171, 221, 164), (94, 79, 162), (244, 109, 67)]
    palette = [(v[0]/256, v[1]/256, v[2]/256) for v in palette]
    width = (g.n_layers - 1) * node_x_distance + offset * 2
    min_l = min((n.layer for n in g.nodes)) - 1
    min_y = min((n.y for n in g.nodes)) if remove_witespace else 0
    for n in g.nodes:
        n.y -= min_y
    if gravity:
        ndy_save = [nd.y for nd in g.nodes]
        for node_list in g.layers.values():
            min_l_y = min((n.y for n in node_list))
            if min_l_y > min_y:
                for n in node_list:
                    n.y -= min_l_y + min_y
        max_n_nodes = max((len(lay) for lay in g.layers.values()))
        for node_list in g.layers.values():
            for n in node_list:
                n.y += (max_n_nodes - len(node_list)) / 2
    height = max((n.y for n in g.nodes)) * node_y_distance + offset * 2 if fix_height == -1 else fix_height
    if nested:
        # if as_png:
        #     surface = cairo.ImageSurface(f"../Images/{svg_name}.svg", width, height)
        # else:
        surface = cairo.SVGSurface(f"../Images/{svg_name}.svg", width, height)
    elif motif:
        surface = cairo.SVGSurface(f"Images/Crossing-Motifs/{svg_name}.svg", width, height)
    else:
        surface = cairo.SVGSurface(f"Images/{svg_name}.svg", width, height)
    ctx = cairo.Context(surface)
    ctx.set_source_rgb(1, 1, 1)
    ctx.rectangle(0, 0, width, height)
    ctx.fill()
    ctx.set_line_width(line_width)
    for edge in g.edges:  # curve_to(c1x, c1y, c2x, c2y, ex, ey), control points c1, c2, end point e
        if (edge.n1.id, edge.n2.id) in dont_draw_edges:
            continue
        ctx.set_source_rgb(0.8, 0.8, 0.8)
        ctx.move_to((edge.n1.layer - 1 - min_l) * node_x_distance + offset, edge.n1.y * node_y_distance + offset)
        if edge_thickness:
            ctx.set_line_width(edge.weight)
        if emphasize_edges and ((edge.n1.id, edge.n2.id) in emphasize_edges or (edge.n2.id, edge.n1.id) in emphasize_edges):
            ctx.set_source_rgb(17/256, 138/256, 89/256)
        if edge.same_layer_edge:
            ctx.curve_to((edge.n1.layer - 1 - min_l) * node_x_distance + offset + node_x_distance//1.5 - (node_x_distance//2)//(abs(edge.n1.y-edge.n2.y)), edge.n1.y * node_y_distance + offset, (edge.n1.layer - 1 - min_l) * node_x_distance + offset + node_x_distance//1.5 - (node_x_distance//2)//(abs(edge.n1.y-edge.n2.y)), edge.n2.y * node_y_distance + offset, (edge.n1.layer - 1 - min_l) * node_x_distance + offset, edge.n2.y * node_y_distance + offset)
        elif edge.n1.y == edge.n2.y or straighten_edges or (straighten_only_true_edges and not edge.n1.is_anchor_node and not edge.n2.is_anchor_node):
            ctx.line_to((edge.n2.layer - 1 - min_l)*node_x_distance + offset, edge.n2.y*node_y_distance + offset)
        else:
            # ctx.curve_to((edge.n1.layer - 1) * node_x_distance + offset + node_x_distance, edge.n1.y * node_y_distance + offset, (edge.n2.layer - 1) * node_x_distance + offset - node_x_distance, edge.n2.y * node_y_distance + offset, (edge.n2.layer - 1) * node_x_distance + offset, edge.n2.y * node_y_distance + offset)
            p1, p2, p3, p4, p5, p6 = bezier_control_points(g, (edge.n1.id, edge.n2.id), edge.n1.y, edge.n2.y, edge.n1.layer, edge.n2.layer, min_l, node_x_distance, node_y_distance, offset, left_straight_edges, right_straight_edges, full_straight_edges, full_straight_long_edges, less_curvy_edges)
            ctx.curve_to(p1, p2, p3, p4, p5, p6)
        ctx.stroke()
    ctx.set_line_width(line_width)

    ctx.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
    ctx.set_font_size(font_size)
    for i, node in enumerate(g.nodes):  # ctx.arc(2, 1, 0.5, 0, 2 * math.pi), pos (2,1) radius 0.5
        if not node.is_anchor_node or ((groups is not None or "groups" in g.node_data) and not dot_group_anchors):
            if node.stacked or node.fix != 0:
                ctx.set_source_rgb(222/256, 23/256, 56/256)
            elif color_scale is not None:
                max_moves = max(color_scale)
                ctx.set_source_rgb(color_scale[i] / max_moves, 0, 0)
                # ctx.set_source_rgb(color_scale[i] / max_moves, 20 / 256, 120 / 256)
            elif groups is not None:
                ctx.set_source_rgb(palette[groups[i]][0], palette[groups[i]][1], palette[groups[i]][2])
            elif "groups" in g.node_data:
                gp_v = g.node_data["groups"][node.id] + 1 if node.id in g.node_data["groups"] else 0
                ctx.set_source_rgb(palette[gp_v][0], palette[gp_v][1], palette[gp_v][2])
            else:
                ctx.set_source_rgb(163/256, 185/256, 182/256)  # light gray-cyan
            if node.is_anchor_node:
                node_radius_mult = 1 / 3
            elif node_weight_size:
                node_radius_mult = node.weight
            elif emphasize_nodes is not None and emphasize_nodes[node.id]:
                node_radius_mult = 1.7
            else:
                node_radius_mult = 1
            ctx.arc((node.layer - 1 - min_l)*node_x_distance + offset, node.y*node_y_distance + offset, node_radius * node_radius_mult, 0, 2 * math.pi)
            ctx.fill()
            # ctx.set_source_rgb(53 / 256, 83 / 256, 232 / 256)  # blueeeee
            ctx.set_source_rgb(0.1, 0.1, 0.1)
            ctx.arc((node.layer - 1 - min_l) * node_x_distance + offset, node.y * node_y_distance + offset, node_radius * node_radius_mult, 0, 2 * math.pi)
            ctx.stroke()
            ctx.set_source_rgb(0.1, 0.1, 0.1)
            if len(str(node.name)) == 1:
                ctx.move_to((node.layer - 1 - min_l)*node_x_distance + offset - 3, node.y*node_y_distance + offset + 4)
            elif len(str(node.name)) == 2:
                ctx.move_to((node.layer - 1 - min_l) * node_x_distance + offset - 7, node.y * node_y_distance + offset + 4)
            elif len(str(node.name)) == 3:
                ctx.move_to((node.layer - 1 - min_l) * node_x_distance + offset - 11, node.y * node_y_distance + offset + 4)
            else:
                ctx.move_to((node.layer - 1 - min_l) * node_x_distance + offset - 11, node.y * node_y_distance + offset + 4 + (2 * (node.layer % 2) - 1) * 25)
            if (not node.is_anchor_node or label_anchors) and label_nodes:
                ctx.show_text(str(node.name))
        else:
            ctx.set_source_rgb(0.2, 0.2, 0.2)
            ctx.arc((node.layer - 1 - min_l)*node_x_distance + offset, node.y*node_y_distance + offset, line_width//2, 0, 2 * math.pi)
            ctx.fill()
            if label_anchors and label_nodes:
                ctx.move_to((node.layer - 1 - min_l) * node_x_distance + offset + 7, node.y * node_y_distance + offset + 4)
                ctx.show_text(str(node.name))

            # ctx.arc((node.layer - 1 - min_l) * node_x_distance + offset, node.y * node_y_distance + offset,
            #         node_radius // 3, 0, 2 * math.pi)
            # ctx.fill()
            #
            # ctx.arc((node.layer - 1 - min_l) * node_x_distance + offset, node.y * node_y_distance + offset,
            #         node_radius // 3, 0, 2 * math.pi)
            # ctx.stroke()
    surface.finish()
    if gravity:
        for i, nd in enumerate(g.nodes):
            nd.y = ndy_save[i]
    if as_png:
        for i in range(0, copies):
            newname = svg_name.split("_")[0] + "_" + str(int(svg_name.split("_")[1]) + i) if "_" in svg_name else svg_name
            cairosvg.svg2png(url=f"Images/{svg_name}.svg", write_to=f"Images/{newname}.png", output_width=width / 2, output_height=height / 2)
        os.remove(f"Images/{svg_name}.svg")
        # surface.write_to_png(svg_name + ".png")


def draw_altair_scatter(data_points, x_axis, y_axis, color_field, x_title, y_title, chart_name, log_y_scale, plot_loess=False, loess_features=None, opacity=1, xdom=None):
    data = alt.Data(values=data_points)
    chart = alt.Chart(data).mark_circle(size=60).encode(
        x=alt.X(f'{x_axis}:Q', axis=alt.Axis(title=x_title), scale=alt.Scale(domain=xdom) if xdom else alt.Scale()),
        y=alt.Y(f'{y_axis}:Q', scale=alt.Scale(type="log") if log_y_scale else alt.Scale(), axis=alt.Axis(title=y_title)) #,
        # color=alt.Color(f'{color_field}:N', scale=alt.Scale(scheme='dark2'))
    ).configure_mark(opacity=opacity)
    # .facet(column=f'{color_field}:N')
    # if plot_loess:
    #     for feature in loess_features:
    #         chart += chart.transform_filter(alt.FieldEqualPredicate(field=color_field, equal=feature)).transform_loess(
    #             x_axis, y_axis
    #         )
    chart.save(f"charts/{chart_name}.html", embed_options={'renderer': 'svg'})


def draw_altair_line_chart(data_points, x_axis, y_axis, color_field, x_title, y_title, chart_name, log_y_scale):
    data = alt.Data(values=data_points)
    dom = list(set(dp[f"{color_field}"] for dp in data_points))
    # rng = ["#26547C", "#F0567A", "#E09D00", "#4ACB2A"]
    rng = ["#e15759", "#b07aa1", "#9c755f", "#f28e2b", "#ff9da7", "#4e79a7", "#59a14f", "#edc948", "#76b7b2"]
    chart = alt.Chart(data).mark_line(point={"filled": False, "fill": "white"}).encode(
        x=alt.X(f'{x_axis}:Q', axis=alt.Axis(title=x_title)),
        y=alt.Y(f'{y_axis}:Q', scale=alt.Scale(type="log") if log_y_scale else None, axis=alt.Axis(title=y_title)),
        color=alt.Color(f'{color_field}:N', scale=alt.Scale(domain=dom, range=rng))
    )
    chart.save(f"charts/{chart_name}.html", embed_options={'renderer': 'svg'})


def draw_altair_simple_line_chart(data_points, x_axis, y_axis, color_field, x_title, y_title, chart_name, xdom=None, ydom=None, coldom=None, rng2=False):
    data = alt.Data(values=data_points)
    dom = coldom if coldom is not None else sorted(list(set(dp[f"{color_field}"] for dp in data_points)))
    rng = ["#000000", "#e15759", "#4e79a7", "#b07aa1", "#f28e2b", "#9c755f", "#ff9da7", "#59a14f", "#4e79a7", "#ff9da7", "#59a14f"]  # "#edc948", "#76b7b2"]
    # color_palette = ["#ff8a80", "#ff1744", "#d50000", "#ea80fc", "#d500f9", "#aa00ff", "#82b1ff", "#2979ff", "#2962ff", "#ccff90", "#76ff03", "#64dd17", "#ffe57f", "#ffc400", "#ffab00"]
    color_palette = ["#e1bee7", "#ce93d8", "#ab47bc", "#7b1fa2", "#4a148c", "#b3e5fc", "#81d4fa", "#29b6f6", "#0288d1", "#01579b", "#c8e6c9", "#a5d6a7", "#66bb6a", "#388e3c", "#1b5e20", "#ffe0b2", "#ffcc80", "#ffa726", "#f57c00", "#e65100"]
    if rng2:
        rng = ["#000000", "#e15759", "#4e79a7", "#9c755f", "#ff9da7", "#59a14f"]
    chart = alt.Chart(data).mark_line(clip=True).encode(
        x=alt.X(f'{x_axis}:Q', axis=alt.Axis(title=x_title), scale=alt.Scale(domain=xdom) if xdom else alt.Scale()),
        y=alt.Y(f'{y_axis}:Q', axis=alt.Axis(title=y_title), scale=alt.Scale(domain=ydom) if ydom else alt.Scale()),
        # y=alt.Y(f'{y_axis}:Q', axis=alt.Axis(title=y_title, format='%'), scale=alt.Scale(domain=ydom) if ydom else alt.Scale()),
        color=alt.Color(f'{color_field}:N', scale=alt.Scale(domain=dom, range=color_palette)),
        # strokeDash=alt.StrokeDash("Dash:N", sort=["normal", "combined", "optimal"])
    )
    # .properties(
    #     width=800,
    #     height=600,
    #     autosize=alt.AutoSizeParams(
    #         type='fit',
    #         contains='padding'
    #     )
    # )
    chart.save(f"charts/{chart_name}.html", embed_options={'renderer': 'svg'})


def draw_altair_scatter_with_regression_line(data_points, x_axis, y_axis, color_field, x_title, y_title, chart_name):
    data = alt.Data(values=data_points)
    chart = alt.Chart(data).mark_circle(size=60).encode(
        x=alt.X(f'{x_axis}:Q', axis=alt.Axis(title=x_title)),
        y=alt.Y(f'{y_axis}:Q', axis=alt.Axis(title=y_title)),
        color=alt.Color(f'{color_field}:N', scale=alt.Scale(scheme='dark2'))
    )
    chart2 = chart + chart.transform_regression(x_axis, y_axis).mark_line()
    chart2.save(f"charts/{chart_name}.html", embed_options={'renderer': 'svg'})


def draw_altair_scatter_with_custom_line(scatter_data, line_data, x_axis, y_axis, color_field, x_title, y_title, chart_name):
    scatterdata = alt.Data(values=scatter_data)
    scatterchart = alt.Chart(scatterdata).mark_circle(size=60, opacity=0.9).encode(
        x=alt.X(f'{x_axis}:Q', axis=alt.Axis(title=x_title)),
        y=alt.Y(f'{y_axis}:Q', axis=alt.Axis(title=y_title)),
        color=alt.Color(f'{color_field}:N', scale=alt.Scale(scheme='dark2'))
    )
    linedata = alt.Data(values=line_data)
    linechart = alt.Chart(linedata).mark_line(color="black").encode(
        x=alt.X(f'{x_axis}:Q', axis=alt.Axis(title=x_title)),
        y=alt.Y(f'{y_axis}:Q', axis=alt.Axis(title=y_title)),
        # color=alt.Color(f'{color_field}:N', scale=alt.Scale())
    )
    chart = scatterchart + linechart
    chart.save(f"charts/{chart_name}.html", embed_options={'renderer': 'svg'})


def draw_altair_line_compare(data_points, x_axis, y_axis, facet_field, x_title, y_title, chart_name, log_y_scale, experiment_name):
    data = alt.Data(values=data_points)
    # dom = [f"{experiment_name}", "baseline"]
    rng = ["#26547C", "#F0567A"]
    chart = alt.Chart(data).mark_circle(size=60).encode(
        x=alt.X(f'{x_axis}:Q', axis=alt.Axis(title=x_title)),
        y=alt.Y(f'{y_axis}:Q', axis=alt.Axis(title=y_title)),
        # color=alt.Color(f'Technique:N', scale=alt.Scale(domain=dom, range=rng))
    )

    horizline = alt.Chart().mark_rule().encode(
        y='a:Q'
    )

    alt.layer(
        chart, horizline,
        data=data
    ).transform_calculate(
        a="100"
    ).facet(
        column=alt.Column(f'{facet_field}:N', sort=["junger_basic", "vertical_transitivity", "redundancy"])
    ).save(f"charts/{chart_name}.html", embed_options={'renderer': 'svg'})


def draw_altair_colored_line_compare(data_points, x_axis, y_axis, facet_field, color_field, x_title, y_title, chart_name):
    data = alt.Data(values=data_points)
    # dom = [cat1name, cat2name]
    rng = ["#26547C", "#F0567A", "#E09D00"]
    chart = alt.Chart(data).mark_circle(size=60).encode(
        x=alt.X(f'{x_axis}:Q', axis=alt.Axis(title=x_title)),
        y=alt.Y(f'{y_axis}:Q', axis=alt.Axis(title=y_title)),
        color=alt.Color(f'{color_field}:N', scale=alt.Scale(range=rng))
    )

    horizline = alt.Chart().mark_rule().encode(
        y='a:Q'
    )

    alt.layer(
        chart, horizline,
        data=data
    ).transform_calculate(
        a="100"
    ).facet(
        column=alt.Column(f'{facet_field}:N', sort=["junger_basic", "vertical_transitivity"])
    ).save(f"charts/{chart_name}.html", embed_options={'renderer': 'svg'})
