import networkx as nx
from database import graph_driver
import constants as cts

# driver = GraphDatabase.driver('bolt://localhost:7687', auth=("neo4j", "libiao"))

# query = """
# WITH ['hsa:6354', 'hsa:3880', 'hsa:6372', 'hsa:51744', 'hsa:958', 'cpd:C00768', 'cpd:C00180'] AS node_id_list
#     UNWIND node_id_list AS ids
#     MATCH (n {kid: ids})
#     WITH collect(n) AS nds

#     UNWIND nds AS n1
#     UNWIND nds AS n2
#     WITH nds, n1, n2 WHERE id(n1) <> id(n2)
#     MATCH path = allShortestPaths((n1)-[rels*]-(n2))
#     WHERE all(r IN relationships(path) where type(r) in ['activation', 'inhibition', 'expression', 'repression', 'product', 'substrate', 'enzyme'])
#     RETURN path ORDER BY length(path) ASC;
# """

# results = driver.session().run(query)

# G = nx.MultiDiGraph()

# nodes = list(results.graph()._nodes.values())
# for node in nodes:
#     G.add_node(node.id, labels=node['kid'], properties=node._properties)

# rels = list(results.graph()._relationships.values())
# for rel in rels:
#     G.add_edge(rel.start_node.id, rel.end_node.id, key=rel.id, type=rel.type, properties=rel._properties)


def get_nx_graph_with_cypher_query(query: str):
    results = graph_driver.session().run(query)
    G = nx.MultiDiGraph()
    nodes = list(results.graph()._nodes.values())
    for node in nodes:
        G.add_node(node.id, labels=node[cts.DB_PK], properties=node._properties)

    rels = list(results.graph()._relationships.values())
    for rel in rels:
        G.add_edge(
            rel.start_node.id,
            rel.end_node.id,
            key=rel.id,
            type=rel.type,
            properties=rel._properties,
        )
    return G
