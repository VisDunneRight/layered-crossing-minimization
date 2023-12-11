import cairo
from src import graph
import math
import altair as alt
import os
# from altair_saver import save


def draw_graph(g: graph.LayeredGraph, svg_name, node_x_distance=150, node_y_distance=100, nested=False, motif=False, groups=None, gravity=False, edge_thickness=False, label_nodes=True):
    if nested:
        if "Images" not in os.listdir(".."):
            os.mkdir("../Images")
    elif "Images" not in os.listdir():
        os.mkdir("Images")
    offset = 40
    node_radius = 15
    line_width = 4
    font_size = 12
    palette = [(171/256, 221/256, 164/256), (94/256, 79/256, 162/256), (244/256, 109/256, 67/256), (254/256, 224/256, 139/256), (50/256, 136/256, 189/256), (158/256, 1/256, 66/256), (253/256, 174/256, 97/256), (102/256, 194/256, 165/256), (213/256, 62/256, 79/256), (230/256, 145/256, 152/256)]
    width = (g.n_layers - 1) * node_x_distance + offset * 10
    min_l = min((n.layer for n in g.nodes)) - 1
    min_y = min((n.y for n in g.nodes))
    for n in g.nodes:
        n.y -= min_y
    if gravity:
        for node_list in g.layers.values():
            min_l_y = min((n.y for n in node_list))
            if min_l_y > min_y:
                for n in node_list:
                    n.y -= min_l_y + min_y
        max_n_nodes = max((len(lay) for lay in g.layers.values()))
        for node_list in g.layers.values():
            for n in node_list:
                n.y += (max_n_nodes - len(node_list)) // 2
    height = max((n.y for n in g.nodes)) * node_y_distance + offset * 2
    if nested:
        surface = cairo.SVGSurface(f"../Images/{svg_name}.svg", width, height)
    elif motif:
        surface = cairo.SVGSurface(f"Images/Crossing-Motifs/{svg_name}.svg", width, height)
    else:
        surface = cairo.SVGSurface(f"Images/{svg_name}.svg", width, height)
    ctx = cairo.Context(surface)
    ctx.set_source_rgb(1, 1, 1)
    ctx.rectangle(0, 0, width, height)
    ctx.fill()
    ctx.set_source_rgb(0.2, 0.2, 0.2)
    ctx.set_line_width(line_width)
    for edge in g.edges:  # curve_to(c1x, c1y, c2x, c2y, ex, ey), control points c1, c2, end point e
        ctx.move_to((edge.n1.layer - 1 - min_l) * node_x_distance + offset, edge.n1.y * node_y_distance + offset)
        if edge_thickness:
            ctx.set_line_width(edge.weight)
        if edge.same_layer_edge:
            ctx.curve_to((edge.n1.layer - 1 - min_l) * node_x_distance + offset + node_x_distance//1.5 - (node_x_distance//2)//(abs(edge.n1.y-edge.n2.y)), edge.n1.y * node_y_distance + offset, (edge.n1.layer - 1 - min_l) * node_x_distance + offset + node_x_distance//1.5 - (node_x_distance//2)//(abs(edge.n1.y-edge.n2.y)), edge.n2.y * node_y_distance + offset, (edge.n1.layer - 1 - min_l) * node_x_distance + offset, edge.n2.y * node_y_distance + offset)
        elif edge.n1.y == edge.n2.y:
            ctx.line_to((edge.n2.layer - 1 - min_l)*node_x_distance + offset, edge.n2.y*node_y_distance + offset)
        else:
            # ctx.curve_to((edge.n1.layer - 1) * node_x_distance + offset + node_x_distance, edge.n1.y * node_y_distance + offset, (edge.n2.layer - 1) * node_x_distance + offset - node_x_distance, edge.n2.y * node_y_distance + offset, (edge.n2.layer - 1) * node_x_distance + offset, edge.n2.y * node_y_distance + offset)
            ctx.curve_to((edge.n1.layer - 1 - min_l) * node_x_distance + offset + node_x_distance//1.5, edge.n1.y * node_y_distance + offset, (edge.n2.layer - 1 - min_l) * node_x_distance + offset - node_x_distance//1.5, edge.n2.y * node_y_distance + offset, (edge.n2.layer - 1 - min_l) * node_x_distance + offset, edge.n2.y * node_y_distance + offset)
        ctx.stroke()
    ctx.set_line_width(line_width)

    ctx.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
    ctx.set_font_size(font_size)
    for node in g.nodes:  # ctx.arc(2, 1, 0.5, 0, 2 * math.pi), pos (2,1) radius 0.5
        if not node.is_anchor_node or groups is not None:
            if node.stacked:
                ctx.set_source_rgb(222/256, 23/256, 56/256)
            elif groups is not None:
                ctx.set_source_rgb(palette[groups[node.id]][0], palette[groups[node.id]][1], palette[groups[node.id]][2])
            else:
                ctx.set_source_rgb(163/256, 185/256, 182/256)  # light gray-cyan
            if node.is_anchor_node:
                ctx.arc((node.layer - 1 - min_l)*node_x_distance + offset, node.y*node_y_distance + offset, node_radius//3, 0, 2 * math.pi)
            else:
                ctx.arc((node.layer - 1 - min_l)*node_x_distance + offset, node.y*node_y_distance + offset, node_radius, 0, 2 * math.pi)
            ctx.fill()
            # ctx.set_source_rgb(53 / 256, 83 / 256, 232 / 256)  # blueeeee
            ctx.set_source_rgb(0.1, 0.1, 0.1)
            if node.is_anchor_node:
                ctx.arc((node.layer - 1 - min_l) * node_x_distance + offset, node.y * node_y_distance + offset, node_radius//3, 0, 2 * math.pi)
            else:
                ctx.arc((node.layer - 1 - min_l) * node_x_distance + offset, node.y * node_y_distance + offset, node_radius, 0, 2 * math.pi)
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
            if not node.is_anchor_node and label_nodes:
                ctx.show_text(str(node.name))
        else:
            ctx.set_source_rgb(0.2, 0.2, 0.2)
            ctx.arc((node.layer - 1 - min_l)*node_x_distance + offset, node.y*node_y_distance + offset, line_width//2, 0, 2 * math.pi)
            ctx.fill()

            # ctx.arc((node.layer - 1 - min_l) * node_x_distance + offset, node.y * node_y_distance + offset,
            #         node_radius // 3, 0, 2 * math.pi)
            # ctx.fill()
            #
            # ctx.arc((node.layer - 1 - min_l) * node_x_distance + offset, node.y * node_y_distance + offset,
            #         node_radius // 3, 0, 2 * math.pi)
            # ctx.stroke()


def draw_altair_scatter(data_points, x_axis, y_axis, color_field, x_title, y_title, chart_name, log_y_scale, plot_loess=False, loess_features=None):
    data = alt.Data(values=data_points)
    chart = alt.Chart(data).mark_circle(size=60).encode(
        x=alt.X(f'{x_axis}:Q', axis=alt.Axis(title=x_title)),
        y=alt.Y(f'{y_axis}:Q', scale=alt.Scale(type="log") if log_y_scale else None, axis=alt.Axis(title=y_title)) #,
        # color=alt.Color(f'{color_field}:N', scale=alt.Scale(scheme='dark2'))
    )
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
    dom = coldom if coldom is not None else list(set(dp[f"{color_field}"] for dp in data_points))
    rng = ["#000000", "#e15759", "#4e79a7", "#b07aa1", "#f28e2b", "#9c755f", "#ff9da7", "#59a14f", "#4e79a7", "#ff9da7", "#59a14f"]  # "#edc948", "#76b7b2"]
    if rng2:
        rng = ["#000000", "#e15759", "#4e79a7", "#9c755f", "#ff9da7", "#59a14f"]
    chart = alt.Chart(data).mark_line(clip=True).encode(
        x=alt.X(f'{x_axis}:Q', axis=alt.Axis(title=x_title), scale=alt.Scale(domain=xdom) if xdom else alt.Scale()),
        y=alt.Y(f'{y_axis}:Q', axis=alt.Axis(title=y_title, format='%'), scale=alt.Scale(domain=ydom) if ydom else alt.Scale()),
        color=alt.Color(f'{color_field}:N', scale=alt.Scale(domain=dom, range=rng)),
        strokeDash=alt.StrokeDash("Dash:N", sort=["normal", "combined", "optimal"])
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
