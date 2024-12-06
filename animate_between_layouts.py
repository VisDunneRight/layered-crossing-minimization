from src.optimization import LayeredOptimizer
from src.tabu import tabu
from src.read_data import read
from src.vis import draw_graph
import imageio.v3 as iio
from numpy import stack
import numpy as np
from pygifsicle import optimize


if __name__ == '__main__':
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

    def sigmoid(x):
        # sigmoid function
        # use k to adjust the slope
        s = 1 / (1 + np.exp(-4*x))
        return s

    zarr = sigmoid(z)
    zarr[-1] = 1
    height = max(max((n.y for n in gr1.nodes)), max((n.y for n in gr2.nodes))) * 100 + 40 * 2

    draw_graph(gr1, f"chmod/frame_0", as_png=True, label_nodes=False, fix_height=height, copies=pause_frames)
    for i in range(frame_count):
        for nd in gr1.nodes:
            nd.y = zarr[i] * (gr2[nd.id].y - initial_positions[nd.id]) + initial_positions[nd.id]
        draw_graph(gr1, f"chmod/frame_{i + pause_frames}", as_png=True, label_nodes=False, fix_height=height, remove_witespace=False)

    frames = stack([iio.imread(f"Images/chmod/frame_{i}.png") for i in range(frame_count + pause_frames)], axis=0)
    iio.imwrite(f'Images/chmod.gif', frames, format="gif", fps=50)
    optimize(f'Images/chmod.gif')
