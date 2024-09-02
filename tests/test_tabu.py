import random
import unittest
import numpy as np
from src import read_data, tabu


class TestTabu(unittest.TestCase):
    def setUp(self) -> None:
        self.g1 = read_data.read("../random graphs/ratio_d3/r1.5k12n8/graph5.lgbin")
        self.g1.y_val_setup()
        self.g2 = read_data.read("../random graphs/ratio_d3/r1.5k12n8/graph10.lgbin")
        self.g2.y_val_setup()

    def test_K_matrix_and_E_values_hard(self):
        """ Validates E array for all nodes by changing its y-value and calling num_edge_crossings() """
        test_g = self.g2
        for nd in test_g:
            nd.y -= 1
        d_adj = test_g.get_double_adj_list()
        rank = np.array([test_g[i].y for i in range(test_g.n_nodes)], np.int_)
        phi = {layer_id: np.array([nd.id for nd in test_g.layers[layer_id]], np.int_) for layer_id in test_g.layers}
        n_cr = test_g.num_edge_crossings()
        for lid in test_g.layers:
            cr_matr, cr_map = tabu.get_crossing_matrix_layer_i(d_adj, rank, test_g.layers[lid])
            for nd in test_g.layers[lid]:
                expected = []
                nd_yv = nd.y
                for i in range(len(test_g.layers[lid])):
                    if i < nd_yv:
                        nd.y = i - 0.5
                        expected.append(n_cr - test_g.num_edge_crossings())
                    elif i == nd_yv:
                        expected.append(0)
                    else:
                        nd.y = i + 0.5
                        expected.append(n_cr - test_g.num_edge_crossings())
                nd.y = nd_yv
                idxs, e_list = tabu.efficient_movement_calculation(cr_matr, cr_map, rank, phi, nd.id, lid)
                print(expected, idxs, e_list)
                self.assertEqual(expected, list(e_list), f"falure at node {nd}")

    def test_repeated_swaps(self):
        """ Validates E array after making multiple random swaps """
        test_g = self.g2
        for nd in test_g:
            nd.y -= 1
        n_random_swaps_per_layer = 10
        d_adj = test_g.get_double_adj_list()
        rank = np.array([test_g[i].y for i in range(test_g.n_nodes)], np.int_)
        phi = {layer_id: np.array([nd.id for nd in test_g.layers[layer_id]], np.int_) for layer_id in test_g.layers}
        n_cr = test_g.num_edge_crossings()
        for lid in test_g.layers:
            cr_matr, cr_map = tabu.get_crossing_matrix_layer_i(d_adj, rank, test_g.layers[lid])
            for _ in range(n_random_swaps_per_layer):
                nd = random.choice(test_g.layers[lid])
                expected = []
                nd_yv = nd.y
                for i in range(len(test_g.layers[lid])):
                    if i < nd_yv:
                        nd.y = i - 0.5
                        expected.append(n_cr - test_g.num_edge_crossings())
                    elif i == nd_yv:
                        expected.append(0)
                    else:
                        nd.y = i + 0.5
                        expected.append(n_cr - test_g.num_edge_crossings())
                nd.y = nd_yv
                idxs, e_list = tabu.efficient_movement_calculation(cr_matr, cr_map, rank, phi, nd.id, lid)
                print(expected, idxs, e_list)
                self.assertEqual(expected, list(e_list), f"falure at node {nd}, rank {rank[nd.id]}")
                tabu.insert_at_position(nd.id, random.randint(0, len(e_list) - 1), rank, phi[lid])
                for nd in test_g.layers[lid]:
                    nd.y = rank[nd.id]
                test_g.layers[lid].sort(key=lambda x: x.y)

    def test_crossing_number_after_swaps(self):
        test_g = self.g2
        for nd in test_g:
            nd.y -= 1
        n_random_swaps_per_layer = 10
        d_adj = test_g.get_double_adj_list()
        rank = np.array([test_g[i].y for i in range(test_g.n_nodes)], np.int_)
        phi = {layer_id: np.array([nd.id for nd in test_g.layers[layer_id]], np.int_) for layer_id in test_g.layers}
        n_cr = test_g.num_edge_crossings()
        for lid in test_g.layers:
            cr_matr, cr_map = tabu.get_crossing_matrix_layer_i(d_adj, rank, test_g.layers[lid])
            # calculate cr from cr_matr
            matr_cr_num = 0
            for nd1 in test_g.layers[lid]:
                for nd2 in test_g.layers[lid]:
                    if nd1.y < nd2.y:
                        matr_cr_num += cr_matr[cr_map[nd1.id]][cr_map[nd2.id]]
                    elif nd1.y == nd2.y and nd1.id != nd2.id:
                        print("un oh")
            for _ in range(n_random_swaps_per_layer):
                nd = random.choice(test_g.layers[lid])
                idxs, e_list = tabu.efficient_movement_calculation(cr_matr, cr_map, rank, phi, nd.id, lid)
                new_pos = random.randint(0, len(e_list)-1)
                print(e_list, nd.y, new_pos, n_cr)
                tabu.insert_at_position(nd.id, new_pos, rank, phi[lid])
                for nd in test_g.layers[lid]:
                    nd.y = rank[nd.id]
                test_g.layers[lid].sort(key=lambda x: x.y)
                new_cr = test_g.num_edge_crossings()

                # calculate cr from cr_matr
                old_matr_cr = matr_cr_num
                matr_cr_num = 0
                for nd1 in test_g.layers[lid]:
                    for nd2 in test_g.layers[lid]:
                        if nd1.y < nd2.y:
                            matr_cr_num += cr_matr[cr_map[nd1.id]][cr_map[nd2.id]]
                        elif nd1.y == nd2.y and nd1.id != nd2.id:
                            print("un oh")

                self.assertEqual(new_cr, n_cr + matr_cr_num - old_matr_cr)
                self.assertEqual(new_cr, n_cr - e_list[new_pos])
                n_cr = new_cr

    def test_swap_easy(self):
        ranking = [1, 2, 3, 0, 4]
        phi_l = [3, 0, 1, 2, 4]
        nid, pos = 3, 3
        tabu.insert_at_position(nid, pos, ranking, phi_l)
        self.assertEqual([0, 1, 2, 3, 4], ranking)
        self.assertEqual([0, 1, 2, 3, 4], phi_l)
        ranking = [1, 2, 3, 0, 4]
        phi_l = [3, 0, 1, 2, 4]
        nid, pos = 4, 1
        tabu.insert_at_position(nid, pos, ranking, phi_l)
        self.assertEqual([2, 3, 4, 0, 1], ranking)
        self.assertEqual([3, 4, 0, 1, 2], phi_l)

    def test_errorcheck_tabu_search(self):
        tabu.tabu(self.g1, maxiter=10000, timelimit=60)
