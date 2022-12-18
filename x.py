from typing import List
import torch
from torch_geometric.data import Data
import torch
from torch_geometric.data import Data
from torch_geometric.datasets import TUDataset

edge_index = torch.tensor([[0, 1, 1, 2], [1, 0, 2, 1]], dtype=torch.long)
x = torch.tensor([[-1], [0], [1]], dtype=torch.float)

data = Data(x=x, edge_index=edge_index)
dataset = TUDataset(root="/tmp/PROTEINS_full", name="PROTEINS_full")


"""
WITH ["hsa:6354", "hsa:3880", "hsa:6372", "hsa:51744", "hsa:958", "cpd:C00768", "cpd:C00180"] AS node_id_list
UNWIND node_id_list AS ids
MATCH (n {kid: ids})
WITH collect(n) AS nds

UNWIND nds AS n1
UNWIND nds AS n2
WITH nds, n1, n2 WHERE id(n1) <> id(n2)
MATCH path = ShortestPath((n1)-[rels*]-(n2))
WHERE all(r IN relationships(path) where type(r) in ["activation", "inhibition", "expression", "repression", "product", "substrate", "enzyme"])
RETURN path ORDER BY length(path) ASC;
"""
