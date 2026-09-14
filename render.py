import tkinter as tk
from tkinter import ttk
from layered_optimization.graph import LayeredGraph
from layered_optimization.vis import bezier_control_points, edges_skipping_anchors


def cubic_bezier(p0, p1, p2, p3, steps=100):
    """Return list of (x,y) points along cubic Bézier from p0..p3."""
    pts = []
    for i in range(steps + 1):
        t = i / steps
        u = 1 - t
        # Bernstein basis for cubic
        b0 = u*u*u
        b1 = 3 * u*u * t
        b2 = 3 * u * t*t
        b3 = t*t*t
        x = b0 * p0[0] + b1 * p1[0] + b2 * p2[0] + b3 * p3[0]
        y = b0 * p0[1] + b1 * p1[1] + b2 * p2[1] + b3 * p3[1]
        pts.append((x, y))
    return pts


def quad_bezier(p0, p1, p2, steps=100):
    pts = []
    for i in range(steps + 1):
        t = i / steps
        u = 1 - t
        b0 = u*u
        b1 = 2*u*t
        b2 = t*t
        x = b0 * p0[0] + b1 * p1[0] + b2 * p2[0]
        y = b0 * p0[1] + b1 * p1[1] + b2 * p2[1]
        pts.append((x, y))
    return pts


def draw_bezier(canvas, pts, color="black", width=2):
    flat = [coord for p in pts for coord in p]
    canvas.create_line(*flat, fill=color, width=width, smooth=False)


def tkinter_draw(canvas: tk.Canvas, g: LayeredGraph, svg_name, node_x_distance=150, node_y_distance=100, nested=False, motif=False, groups=None, node_outline=False, emphasize_nodes=None, emphasize_edges=None, gravity=False, edge_thickness=False, label_nodes=True, as_png=False, color_scale=None, fix_height=-1, remove_whitespace=True, straighten_edges=False, node_weight_size=False, dot_group_anchors=False, label_anchors=False, ignore_colinear_anchors=True):
    offset = 40
    node_radius = 15
    line_width = 4
    font_size = 12
    min_l = min((n.layer for n in g.nodes)) - 1
    min_y = min((n.y for n in g.nodes)) if remove_whitespace else 0
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

    g_edges = g.edge_ids.keys()
    if ignore_colinear_anchors:
        g_edges, skip_nodes = edges_skipping_anchors(g)

    for edge in g_edges:  # curve_to(c1x, c1y, c2x, c2y, ex, ey), control points c1, c2, end point e
        ctx.set_source_rgb(0.8, 0.8, 0.8)
        ctx.move_to((g[edge[0]].layer - 1 - min_l) * node_x_distance + offset, g[edge[0]].y * node_y_distance + offset)
        if edge_thickness:
            ctx.set_line_width(g.edge_ids[edge].weight)
        if emphasize_edges and (edge in emphasize_edges or (edge[1], edge[0]) in emphasize_edges):
            ctx.set_source_rgb(17/256, 138/256, 89/256)  #008000  #DDDDDD
        if g[edge[0]].layer == g[edge[1]].layer:
            canvas.create_line((g[edge[0]].layer - 1 - min_l) * node_x_distance + offset + node_x_distance//1.5 - (node_x_distance//2)//(abs(g[edge[0]].y-g[edge[1]].y)), g[edge[0]].y * node_y_distance + offset, (g[edge[0]].layer - 1 - min_l) * node_x_distance + offset + node_x_distance//1.5 - (node_x_distance//2)//(abs(g[edge[0]].y-g[edge[1]].y)), g[edge[1]].y * node_y_distance + offset, (g[edge[0]].layer - 1 - min_l) * node_x_distance + offset, g[edge[1]].y * node_y_distance + offset)
        elif g[edge[0]].y == g[edge[1]].y or straighten_edges or (straighten_only_true_edges and not g[edge[0]].is_anchor_node and not g[edge[1]].is_anchor_node):
            ctx.line_to((g[edge[1]].layer - 1 - min_l)*node_x_distance + offset, g[edge[1]].y*node_y_distance + offset)
        else:
            p1, p2, p3, p4, p5, p6 = bezier_control_points(g, edge, g[edge[0]].y, g[edge[1]].y, g[edge[0]].layer, g[edge[1]].layer, min_l, node_x_distance, node_y_distance, offset, None, None, None, None, None)
            ctx.curve_to(p1, p2, p3, p4, p5, p6)
        ctx.stroke()

    ctx.set_line_width(line_width)
    ctx.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)

    if add_text:
        ctx.set_font_size(48)
        ctx.set_source_rgb(0, 0, 0)
        ctx.move_to((g.n_layers - 1) * node_x_distance + offset * 2, offset * 1.5)
        if type(add_text) == list:
            for i, text in enumerate(add_text):
                ctx.move_to((g.n_layers - 1) * node_x_distance + offset * 2, offset * 1.5 + i * 60)
                ctx.show_text(text)
        else:
            ctx.show_text(add_text)

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
            if node_outline:
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
        elif not (ignore_colinear_anchors and node.id in skip_nodes):
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


class ZoomPanCanvas(tk.Canvas):
    def __init__(self, master, **kwargs):
        super().__init__(master, **kwargs)

        # Bindings for panning
        self.bind("<ButtonPress-1>", self.start_pan)
        self.bind("<B1-Motion>", self.do_pan)

        # Bindings for mouse-wheel zoom (cross-platform)
        self.bind("<MouseWheel>", self.zoom)        # Windows, macOS
        self.bind("<Button-4>", self.zoom)         # Linux scroll up
        self.bind("<Button-5>", self.zoom)         # Linux scroll down

        # For keeping track of scaling
        self.scale_factor = 1.0

    # -----------------------------
    # Panning
    # -----------------------------
    def start_pan(self, event):
        self.scan_mark(event.x, event.y)

    def do_pan(self, event):
        self.scan_dragto(event.x, event.y, gain=1)

    # -----------------------------
    # Zooming
    # -----------------------------
    def zoom(self, event):
        # Zoom direction
        if event.num == 5 or event.delta < 0:
            factor = 0.9
        else:
            factor = 1.1

        self.scale_factor *= factor

        # Zoom relative to mouse position
        x = self.canvasx(event.x)
        y = self.canvasy(event.y)

        self.scale("all", x, y, factor, factor)
        self.configure(scrollregion=self.bbox("all"))


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Tkinter Canvas with Side Panel")
        self.geometry("900x600")

        # Main container frame
        container = ttk.Frame(self)
        container.pack(fill="both", expand=True)

        container.columnconfigure(0, weight=3)
        container.columnconfigure(1, weight=1)
        container.rowconfigure(0, weight=1)

        # Canvas frame
        canvas_frame = ttk.Frame(container)
        canvas_frame.grid(row=0, column=0, sticky="nsew")

        # Side panel frame
        panel_frame = ttk.Frame(container, padding=10)
        panel_frame.grid(row=0, column=1, sticky="nsew")

        # Canvas widget
        self.canvas = ZoomPanCanvas(canvas_frame, bg="white")
        self.canvas.pack(fill="both", expand=True)

        # Options panel
        ttk.Label(panel_frame, text="Options", font=("Arial", 14, "bold")).pack(anchor="w", pady=5)

        self.show_grid = tk.BooleanVar()
        ttk.Checkbutton(panel_frame, text="Show Grid", variable=self.show_grid, command=self.toggle_grid).pack(anchor="w")

        self.enable_draw = tk.BooleanVar(value=True)
        ttk.Checkbutton(panel_frame, text="Enable Drawing", variable=self.enable_draw).pack(anchor="w")

        ttk.Label(panel_frame, text="Mode", font=("Arial", 12)).pack(anchor="w", pady=(10, 5))
        self.draw_mode = tk.StringVar(value="line")
        for label, value in [("Line", "line"), ("Rectangle", "rect"), ("Circle", "circle")]:
            ttk.Radiobutton(panel_frame, text=label, variable=self.draw_mode, value=value).pack(anchor="w")

        ttk.Button(panel_frame, text="Clear Canvas", command=self.clear_canvas).pack(pady=20)

        # Draw initial content AFTER the window is actually sized
        self.after(100, self.draw_initial)

    def draw_initial(self):
        w = self.canvas.winfo_width()
        h = self.canvas.winfo_height()

        # Sample graphic
        self.canvas.create_rectangle(w/4, h/4, w/2, h/2, fill="lightblue")
        self.canvas.create_text(w/2.8, h/3, text="Canvas Area", font=("Arial", 14))

    def toggle_grid(self):
        self.canvas.delete("grid")
        if not self.show_grid.get():
            return

        w = self.canvas.winfo_width()
        h = self.canvas.winfo_height()

        for i in range(0, w, 20):
            self.canvas.create_line(i, 0, i, h, fill="#e0e0e0", tags="grid")
        for j in range(0, h, 20):
            self.canvas.create_line(0, j, w, j, fill="#e0e0e0", tags="grid")

    def clear_canvas(self):
        self.canvas.delete("all")


if __name__ == "__main__":
    App().mainloop()