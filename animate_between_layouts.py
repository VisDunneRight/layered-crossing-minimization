from src.optimization import LayeredOptimizer
from src.tabu import tabu
from src.read_data import read
from src.vis import draw_graph
import imageio.v3 as iio
from numpy import stack
import numpy as np
from pygifsicle import optimize


def sigmoid(x):
    # sigmoid function
    # use k to adjust the slope
    s = 1 / (1 + np.exp(-4 * x))
    return s


def draw_chmod():
    gr1 = read("control-flow-graphs/chmod/dbg.main.dot")
    tabu(gr1)
    optimizer = LayeredOptimizer(gr1, symmetry_breaking=True)
    optimizer.just_bendiness_reduction(streamline=False)
    initial_positions = [nd.y for nd in gr1]

    gr2 = read("control-flow-graphs/chmod/dbg.main.dot")
    optimizer = LayeredOptimizer(gr2, symmetry_breaking=True)
    optimizer.optimize_layout()
    optimizer.just_bendiness_reduction(streamline=False)

    frame_count = 50
    pause_frames = 50
    z = np.linspace(-1, 1, frame_count)

    zarr = sigmoid(z)
    zarr[-1] = 1
    height = max(max((n.y for n in gr1.nodes)), max((n.y for n in gr2.nodes))) * 100 + 40 * 2

    draw_graph(gr1, f"chmod/frame_0", as_png=True, label_nodes=False, fix_height=height, copies=pause_frames)
    for i in range(frame_count):
        for nd in gr1.nodes:
            nd.y = zarr[i] * (gr2[nd.id].y - initial_positions[nd.id]) + initial_positions[nd.id]
        draw_graph(gr1, f"chmod/frame_{i + pause_frames}", as_png=True, label_nodes=False, fix_height=height,
                   remove_witespace=False)

    frames = stack([iio.imread(f"Images/chmod/frame_{i}.png") for i in range(frame_count + pause_frames)], axis=0)
    iio.imwrite(f'Images/chmod.gif', frames, format="gif", fps=50)
    optimize(f'Images/chmod.gif')


def animate_two_layouts(g, y_1, y_2, fname, skip_ids=None):
    frame_count = 50
    pause_frames = 50
    z = np.linspace(-1, 1, frame_count)

    zarr = sigmoid(z)
    zarr[-1] = 1
    height = max(y_1 + y_2) * 100 + 40 * 2

    if skip_ids is not None:
        for v in sorted(skip_ids):
            y_1.insert(v, -1)
            y_2.insert(v, -1)
    for i, v in enumerate(y_1):
        if v != -1:
            g[i].y = v
    draw_graph(g, f"{fname}/frame_0", as_png=True, label_nodes=False, fix_height=height, copies=pause_frames, remove_witespace=False)
    for i in range(frame_count):
        for nd in g.nodes:
            nd.y = zarr[i] * (y_2[nd.id] - y_1[nd.id]) + y_1[nd.id]
        draw_graph(g, f"{fname}/frame_{i + pause_frames}", as_png=True, label_nodes=False, fix_height=height, remove_witespace=False)
    for i, v in enumerate(y_2):
        if v != -1:
            g[i].y = v
    draw_graph(g, f"{fname}/frame_100", as_png=True, label_nodes=False, fix_height=height, copies=pause_frames, remove_witespace=False)

    frames = stack([iio.imread(f"Images/{fname}/frame_{i}.png") for i in range(frame_count + 2 * pause_frames)], axis=0)
    iio.imwrite(f'Images/{fname}.gif', frames, format="gif", fps=50)
    optimize(f'Images/{fname}.gif')


if __name__ == '__main__':
    gr = read("./random graphs/networkx2/graph_40_6")
    for to_remove in [18, 12, 28, 29, 22]:
        gr.remove_node(to_remove)
    gr.reindex_nodes()
    gr.relayer()
    opt = LayeredOptimizer(gr)
    opt.optimize_layout(crossing_minimization=True)
    yvs2, yvs2_c = [nd.y for nd in gr.nodes], [nd.y for nd in gr.nodes]
    yvs1 = [max(yvs2) / 2 for _ in range(len(gr.nodes))]
    animate_two_layouts(gr, yvs1, yvs2, "HTDPPT/gr6anim")
    opt.optimize_layout(bendiness_reduction=True, streamline=True, fix_x_vars=True)
    yvs3 = [nd.y for nd in gr.nodes]
    yvs3 = [v + 1 for v in yvs3]
    animate_two_layouts(gr, yvs2_c, yvs3, "HTDPPT/gr6anim2")

    frames = stack([iio.imread(f"Images/HTDPPT/gr6anim/frame_{i}.png") for i in range(150)] + [iio.imread(f"Images/HTDPPT/gr6anim2/frame_{i}.png") for i in range(150)], axis=0)
    iio.imwrite(f'Images/HTDPPT/gr6combo.gif', frames, format="gif", fps=50)
    optimize(f'Images/HTDPPT/gr6combo.gif')
