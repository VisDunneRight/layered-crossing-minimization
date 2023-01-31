import cairo
from src import graph
import math
import altair as alt
from altair_saver import save


def draw_graph(g: graph.LayeredGraph, svg_name, node_x_distance=150, node_y_distance=100, nested=False, motif=False, groups=None, gravity=False):
    offset = 40
    node_radius = 15
    line_width = 4
    font_size = 12
    palette = [(171/256, 221/256, 164/256), (94/256, 79/256, 162/256), (244/256, 109/256, 67/256), (254/256, 224/256, 139/256), (50/256, 136/256, 189/256), (158/256, 1/256, 66/256), (253/256, 174/256, 97/256), (102/256, 194/256, 165/256), (213/256, 62/256, 79/256), (230/256, 145/256, 152/256)]
    width = (g.n_layers - 1) * node_x_distance + offset * 2
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
        if edge.same_layer_edge:
            ctx.curve_to((edge.n1.layer - 1 - min_l) * node_x_distance + offset + node_x_distance//1.5 - (node_x_distance//2)//(abs(edge.n1.y-edge.n2.y)), edge.n1.y * node_y_distance + offset, (edge.n1.layer - 1 - min_l) * node_x_distance + offset + node_x_distance//1.5 - (node_x_distance//2)//(abs(edge.n1.y-edge.n2.y)), edge.n2.y * node_y_distance + offset, (edge.n1.layer - 1 - min_l) * node_x_distance + offset, edge.n2.y * node_y_distance + offset)
        elif edge.n1.y == edge.n2.y:
            ctx.line_to((edge.n2.layer - 1 - min_l)*node_x_distance + offset, edge.n2.y*node_y_distance + offset)
        else:
            # ctx.curve_to((edge.n1.layer - 1) * node_x_distance + offset + node_x_distance, edge.n1.y * node_y_distance + offset, (edge.n2.layer - 1) * node_x_distance + offset - node_x_distance, edge.n2.y * node_y_distance + offset, (edge.n2.layer - 1) * node_x_distance + offset, edge.n2.y * node_y_distance + offset)
            ctx.curve_to((edge.n1.layer - 1 - min_l) * node_x_distance + offset + node_x_distance//1.5, edge.n1.y * node_y_distance + offset, (edge.n2.layer - 1 - min_l) * node_x_distance + offset - node_x_distance//1.5, edge.n2.y * node_y_distance + offset, (edge.n2.layer - 1 - min_l) * node_x_distance + offset, edge.n2.y * node_y_distance + offset)
        ctx.stroke()

    ctx.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
    ctx.set_font_size(font_size)
    for node in g.nodes:  # ctx.arc(2, 1, 0.5, 0, 2 * math.pi), pos (2,1) radius 0.5
        if not node.is_anchor_node or groups is not None:
            if node.stacked:
                ctx.set_source_rgb(222/256, 23/256, 56/256)
            elif groups is not None:
                ctx.set_source_rgb(palette[groups[node.name]][0], palette[groups[node.name]][1], palette[groups[node.name]][2])
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
            else:
                ctx.move_to((node.layer - 1 - min_l) * node_x_distance + offset - 7, node.y * node_y_distance + offset + 4)
            if not node.is_anchor_node:
                ctx.show_text(str(node.name))
        else:
            ctx.set_source_rgb(0.2, 0.2, 0.2)
            # ctx.arc((node.layer - 1 - min_l)*node_x_distance + offset, node.y*node_y_distance + offset, line_width//2, 0, 2 * math.pi)
            # ctx.fill()

            ctx.arc((node.layer - 1 - min_l) * node_x_distance + offset, node.y * node_y_distance + offset,
                    node_radius // 3, 0, 2 * math.pi)
            ctx.fill()

            ctx.arc((node.layer - 1 - min_l) * node_x_distance + offset, node.y * node_y_distance + offset,
                    node_radius // 3, 0, 2 * math.pi)
            ctx.stroke()


# data = alt.Data(values=[{'x': 'A', 'y': 5, 'col': 'red'},
#                         {'x': 'B', 'y': 3, 'col': 'red'},
#                         {'x': 'C', 'y': 6, 'col': 'blue'},
#                         {'x': 'D', 'y': 7, 'col': 'blue'},
#                         {'x': 'E', 'y': 2, 'col': 'red'}])

def draw_altair_scatter(data_points, x_axis, y_axis, color_field, x_title, y_title, chart_name, log_y_scale):
    data = alt.Data(values=data_points)
    chart = alt.Chart(data).mark_circle(size=60).encode(
        x=alt.X(f'{x_axis}:Q', axis=alt.Axis(title=x_title)),
        y=alt.Y(f'{y_axis}:Q', scale=alt.Scale(type="log") if log_y_scale else None, axis=alt.Axis(title=y_title)),
        color=alt.Color(f'{color_field}:N', scale=alt.Scale(scheme='dark2'))
    )
    save(chart, f"charts/{chart_name}.svg")
