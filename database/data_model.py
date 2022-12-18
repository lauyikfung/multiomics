# def create_similarity_relation(node_id_list: list):
#     nodes = []
#     for i in node_id_list:
#         i = i.split(":")[-1]
#         node = graph.nodes.match(id=i).first()
#         if node:
#             nodes.append(node)
#     node_comb_list = list(itertools.product(nodes, nodes))
#     for n1, n2 in node_comb_list:
#         if n1["id"] != n2["id"]:
#             graph.create(Relationship(n1, "similarity", n2))


# def create_similarity_ortholog_with_single_pathway(data):
#     for entry in data["entries"]:
#         if entry["type"] == "ortholog":
#             name_list = entry["name"].split(" ")
#             if len(name_list) > 1:
#                 create_similarity_relation(name_list)
