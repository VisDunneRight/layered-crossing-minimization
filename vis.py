import cairo
import graph
import math


def draw(g: graph.LayeredGraph, svg_name, node_x_distance=100, node_y_distance=100):
    width, height = 500, 420
    offset = 40
    node_radius = 15
    line_width = 4
    font_size = 12
    surface = cairo.SVGSurface(f"Images/{svg_name}.svg", width, height)
    ctx = cairo.Context(surface)
    ctx.set_source_rgb(1, 1, 1)
    ctx.rectangle(0, 0, width, height)
    ctx.fill()

    ctx.set_source_rgb(0.2, 0.2, 0.2)
    ctx.set_line_width(line_width)
    for edge in g.edges:  # curve_to(c1x, c1y, c2x, c2y, ex, ey), control points c1, c2, end point e
        ctx.move_to((edge.n1.layer - 1) * node_x_distance + offset, edge.n1.y * node_y_distance + offset)
        if edge.same_layer_edge:
            ctx.curve_to((edge.n1.layer - 1) * node_x_distance + offset + node_x_distance//2, edge.n1.y * node_y_distance + offset, (edge.n1.layer - 1) * node_x_distance + offset + node_x_distance//2, edge.n2.y * node_y_distance + offset, (edge.n1.layer - 1) * node_x_distance + offset, edge.n2.y * node_y_distance + offset)
        elif edge.n1.y == edge.n2.y:
            ctx.line_to((edge.n2.layer - 1)*node_x_distance + offset, edge.n2.y*node_y_distance + offset)
        else:
            ctx.curve_to((edge.n1.layer - 1) * node_x_distance + offset + node_x_distance, edge.n1.y * node_y_distance + offset, (edge.n2.layer - 1) * node_x_distance + offset - node_x_distance, edge.n2.y * node_y_distance + offset, (edge.n2.layer - 1) * node_x_distance + offset, edge.n2.y * node_y_distance + offset)
        ctx.stroke()

    ctx.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
    ctx.set_font_size(font_size)
    for node in g.nodes:  # ctx.arc(2, 1, 0.5, 0, 2 * math.pi), pos (2,1) radius 0.5
        if not node.is_anchor_node:
            ctx.set_source_rgb(163/256, 185/256, 182/256)  # light gray-cyan
            ctx.arc((node.layer - 1)*node_x_distance + offset, node.y*node_y_distance + offset, node_radius, 0, 2 * math.pi)
            ctx.fill()
            ctx.set_source_rgb(48 / 256, 62 / 256, 63 / 256)  # dark gray-cyan
            ctx.arc((node.layer - 1) * node_x_distance + offset, node.y * node_y_distance + offset, node_radius, 0, 2 * math.pi)
            ctx.stroke()
            ctx.set_source_rgb(0.1, 0.1, 0.1)
            ctx.move_to((node.layer - 1)*node_x_distance + offset - 5, node.y*node_y_distance + offset + 5)
            ctx.show_text(node.name)
        else:
            ctx.set_source_rgb(0.2, 0.2, 0.2)
            ctx.arc((node.layer - 1)*node_x_distance + offset, node.y*node_y_distance + offset, line_width * 2, 0, 2 * math.pi)
            ctx.fill()

    # surface.write_to_png("ex.png")
    # surface.write_to_png(f"./images/{svg_name}.png")


# import math
# import cairo
#
# WIDTH, HEIGHT = 256, 256
# surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, WIDTH, HEIGHT)
# ctx = cairo.Context(surface)
#
# ctx.scale(WIDTH, HEIGHT)  # Normalizing the canvas
#
# ctx.rectangle(0, 0, 1, 1)  # Rectangle(x0, y0, x1, y1)
# ctx.set_source_rgb(1, 1, 1)
# ctx.fill()
#
# ctx.translate(0.1, 0.1)  # Changing the current transformation matrix
#
# ctx.move_to(0, 0)
# # Arc(cx, cy, radius, start_angle, stop_angle)
# ctx.arc(0.2, 0.1, 0.1, -math.pi / 2, 0)
# ctx.line_to(0.5, 0.1)  # Line to (x,y)
# # Curve(x1, y1, x2, y2, x3, y3)
# ctx.curve_to(0.5, 0.2, 0.5, 0.4, 0.2, 0.8)
# ctx.close_path()
#
# ctx.set_source_rgb(0.3, 0.2, 0.5)  # Solid color
# ctx.set_line_width(0.02)
# ctx.stroke()
#
# surface.write_to_png("example.png")
