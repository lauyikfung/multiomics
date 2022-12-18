#!/usr/bin/python
# -*- coding: UTF-8 -*-
from neo4j import GraphDatabase
from py2neo import Node, Relationship, Path
import pandas as pd
import networkx as nx
from py2neo import Graph as pGraph
import py2neo
from matplotlib import pyplot as plt
from igraph import Graph as iGraph
import igraph
import networkx as nx
import time
import sys
# neo4j = pGraph("bolt://localhost:7687", auth=("neo4j", "libiao"))
THRESHOLDS = [float(i) for i in sys.argv[1:]]
def graph_from_cypher(data):
    """Constructs a networkx graph from the results of a neo4j cypher query.
    Example of use:
    >>> result = session.run(query)
    >>> G = graph_from_cypher(result.data())

    Nodes have fields 'labels' (frozenset) and 'properties' (dicts). Node IDs correspond to the neo4j graph.
    Edges have fields 'type_' (string) denoting the type of relation, and 'properties' (dict)."""

    G = nx.MultiDiGraph()
    def add_node(node):
        # Adds node id it hasn't already been added
        u = node.identity
        if G.has_node(u):
            return
        G.add_node(u, labels=node._labels, properties=dict(node))

    def add_edge(relation):
        # Adds edge if it hasn't already been added.
        # Make sure the nodes at both ends are created
        for node in (relation.start_node, relation.end_node):
            add_node(node)
        # Check if edge already exists
        rel_type = str(type(relation))
        u = relation.start_node.identity
        v = relation.end_node.identity
        eid = relation.identity
        # u_type = list(relation.start_node._labels)[0]
        # v_type = list(relation.end_node._labels)[0]
        # u_trend = 
        if G.has_edge(u, v, key=eid):
            return
        # If not, create it
        
        G.add_edge(u, v, key=eid, type_=rel_type[rel_type.find("data.")+5:-2], properties=dict(relation))

    def add_path(entry):
        for rel in entry.relationships:
            add_edge(rel)

    for d in data:
        for entry in d.values():
            # Parse node
            if isinstance(entry, Node):
                add_node(entry)

            # Parse link
            elif isinstance(entry, Relationship):
                add_edge(entry)
            elif isinstance(entry, Path):
                add_path(entry)
            else:
                raise TypeError("Unrecognized object")
    return G
for THRESHOLD in THRESHOLDS:
  start_time = time.time()
  neo4j = pGraph("bolt://localhost:7687", auth=("neo4j", "libiao"))
  dats = neo4j.run("match (n)-[p]-(n2) where (n:reaction OR n:ortholog OR n:gene OR n:compound) \
    AND (n2:reaction OR n2:compound OR n2:ortholog OR n2:gene) return *;")
  end_time = time.time()
  print("Data downloaded")
  print(str(round((end_time-start_time)/60, 2))+ "min")
  G = graph_from_cypher(dats.data())
  print("Graph converted")
  end_time = time.time()
  print(str(round((end_time-start_time)/60, 2))+ "min")
  import copy
  gg = iGraph.from_networkx(G)
  EDGE_TYPE = ['expression', 'gene_expression', 'repression', 'activation', 'inhibition', 'substrate', 'product', 'enzyme']
  for edge in gg.es:
      if edge['type_'] not in EDGE_TYPE:
          gg.delete_edges(edge)
          pass
      # print(edge.source_vertex['properties'])
  for vertice in gg.vs:
      if vertice.degree() == 0:
          gg.delete_vertices(vertice)
  for vertice in gg.vs:
      if vertice.degree() == 0:
          gg.delete_vertices(vertice)
  for vertice in gg.vs:
      v_type = list(vertice['labels'])[0]
      if v_type in ['gene', 'ortholog', 'compound']:
        if vertice["properties"]['trend'] != -2:
          # print(str(vertice["properties"]['log2FC'])+"@"+str(vertice["properties"]['trend']))
          try:
              log2fc = vertice["properties"]['log2FC_covid']
          except:
              vertice["properties"]["trend"] = -2
              continue
          try:
              pvalueFDR = vertice["properties"]['pvalueFDR_covid']
          except:
              vertice["properties"]["trend"] = 0
              continue
          if log2fc > THRESHOLD and pvalueFDR < 0.05:
            vertice["properties"]["trend"] = 1
          elif log2fc < -THRESHOLD and pvalueFDR < 0.05:
            vertice["properties"]["trend"] = -1
          else:
            vertice["properties"]["trend"] = 0
  g = copy.deepcopy(gg)
  g.delete_edges(g.es)
  def adding_edge(e, g):
      if e not in g.es:
          g.add_edge(e.source_vertex, e.target_vertex, type_=e['type_'], properties=e['properties'])
  for vertice in gg.vs:
      v_type = list(vertice['labels'])[0]
      if v_type == "ortholog":
          if vertice["properties"]['trend'] != -2:
              for e in vertice.all_edges():
                  if (e['type_'] in ['expression', 'gene_expression']) and \
                  ([e.source_vertex["properties"]['trend'], e.target_vertex["properties"]['trend']] in [[-1,-1],[1,1]]):
                      adding_edge(e, g)
                  elif (e['type_'] in ['repression']) and \
                  ([e.source_vertex["properties"]['trend'], e.target_vertex["properties"]['trend']] in [[-1,1],[1,-1]]):
                      adding_edge(e, g)
                  elif (e['type_'] in ['activation']) and e.source_vertex == vertice and\
                  ([e.source_vertex["properties"]['trend'], e.target_vertex["properties"]['trend']] in [[-1,-1],[1,1]]):
                      adding_edge(e, g)
                  elif (e['type_'] in ['inhibition']) and e.source_vertex == vertice and \
                  ([e.source_vertex["properties"]['trend'], e.target_vertex["properties"]['trend']] in [[-1,1],[1,-1]]):
                      adding_edge(e, g)
      elif v_type == "reaction":
          d = vertice.all_edges()# 0,1,-2,-1
          d_en = [[], [], [], []]
          d_su = [[], [], [], []]
          d_pr = [[], [], [], []]
          for e in d:
              if e['type_'] == 'enzyme':
                  d_en[e.source_vertex["properties"]['trend']].append(e)
              elif e['type_'] == 'substrate':
                  d_su[e.source_vertex["properties"]['trend']].append(e)
              elif e['type_'] == 'product':
                  d_pr[e.source_vertex["properties"]['trend']].append(e)
          if len(d_en[0]) + len(d_en[1]) + len(d_en[-1]) == 0:
              if len(d_su[1]) > len(d_su[-1]) and len(d_pr[1]) > len(d_pr[-1]):
                  for e in d_su[1] + d_en[-2] + d_pr[1]:
                      adding_edge(e, g)
              elif len(d_su[1]) < len(d_su[-1]) and len(d_pr[1]) < len(d_pr[-1]):
                  for e in d_su[-1] + d_en[-2] + d_pr[-1]:
                      adding_edge(e, g)
          elif len(d_en[1]) > len(d_en[-1]) and len(d_pr[1]) > len(d_pr[-1]):
              for e in d_en[1] + d_pr[1] + d_su[1] + d_su[0] + d_su[-2]:
                  adding_edge(e, g)
          elif len(d_en[1]) < len(d_en[-1]) and len(d_pr[1]) < len(d_pr[-1]):
              for e in d_en[-1] + d_pr[-1] + d_su[-1] + d_su[0] + d_su[-2]:
                  adding_edge(e, g)
          elif len(d_su[1]) > len(d_su[-1]) and len(d_pr[1]) > len(d_pr[-1]):
              for e in d_en[1] + d_pr[1] + d_su[1] + d_en[0] + d_en[-2]:
                  adding_edge(e, g)
          elif len(d_su[1]) < len(d_su[-1]) and len(d_pr[1]) < len(d_pr[-1]):
              for e in d_en[-1] + d_pr[-1] + d_su[-1] + d_en[0] + d_en[-2]:
                  adding_edge(e, g)
      elif vertice['properties']['trend'] == -2:
          d_up = [[], [], [], []]
          d_down = [[], [], [], []]
          for e in vertice.all_edges():
              if e['type_'] in ['enzyme', 'substrate', 'product']:
                  continue
              elif e['type_'] in ['expression', 'gene_expression','activation']:
                  if e.source_vertex == vertice:
                      d_up[e.target_vertex["properties"]['trend']].append(e)
                  else:
                      d_up[e.source_vertex["properties"]['trend']].append(e)
              elif e['type_'] == ['repression', 'inhibition']:
                  if e.source_vertex == vertice:
                      d_down[e.target_vertex["properties"]['trend']].append(e)
                  else:
                      d_down[e.source_vertex["properties"]['trend']].append(e)
          
          if len(d_up[1]) + len(d_down[-1]) >= (len(d_up[-1]) + len(d_down[1])) * 2:
              if len(d_up[1]) > 0 and len(d_down[-1]) > 0:
                  for e in d_up[1] + d_down[-1]:
                      adding_edge(e, g)
          elif len(d_up[-1]) + len(d_down[1]) >= (len(d_up[1]) + len(d_down[-1])) * 2:
              if len(d_up[-1]) > 0 and len(d_down[1]) > 0:
                  for e in d_up[-1] + d_down[1]:
                      adding_edge(e, g)
  edge_clear = True
  while True:
      edge_clear = True
      for vertice in g.vs:
          if vertice.degree() == 0:
              edge_clear = False
              g.delete_vertices(vertice)
      if edge_clear:
          break
  print("Final graph constructed")
  end_time = time.time()
  print(str(round((end_time-start_time)/60, 2))+ "min")
  adj = g.get_adjacency()
  SIDE = adj.shape[0]
  avail = [True for i in range(SIDE)]
  cnt = 0
  comm_cnt = 0
  community_list = []
  community_list_sketch = []
  while cnt < SIDE:
    if not avail[cnt]:
      cnt = cnt + 1
      continue
    stack_n = [cnt]
    k = [[],[],[],[],[]]
    cnt_cpd, cnt_rea, cnt_gen, cnt_orth, cnt_total = 0, 0, 0, 0, 0
    cpd_list = [0,0,0,0]
    gen_list = [0,0,0,0]
    orth_list = [0,0,0,0]
    total_list = [0,0,0,0]
    while len(stack_n) > 0:
      cnt_total = cnt_total + 1
      n = stack_n.pop()
      if list(g.vs["labels"][n])[0] == "compound":
        cnt_cpd = cnt_cpd + 1
        k[0].append((g.vs["properties"][n]['kid'],g.vs["properties"][n]['trend']))
        k[4].append((g.vs["properties"][n]['kid'],g.vs["properties"][n]['trend']))
        total_list[g.vs["properties"][n]['trend']] = total_list[g.vs["properties"][n]['trend']] + 1
        cpd_list[g.vs["properties"][n]['trend']] = cpd_list[g.vs["properties"][n]['trend']] + 1
      elif list(g.vs["labels"][n])[0] == "reaction":
        cnt_rea = cnt_rea + 1
        k[1].append((g.vs["properties"][n]['kid'],-3))
        k[4].append((g.vs["properties"][n]['kid'],-3))
      elif list(g.vs["labels"][n])[0] == "gene":
        cnt_gen = cnt_gen + 1
        k[2].append((g.vs["properties"][n]['kid'],g.vs["properties"][n]['trend']))
        k[4].append((g.vs["properties"][n]['kid'],g.vs["properties"][n]['trend']))
        total_list[g.vs["properties"][n]['trend']] = total_list[g.vs["properties"][n]['trend']] + 1
        gen_list[g.vs["properties"][n]['trend']] = gen_list[g.vs["properties"][n]['trend']] + 1
      elif list(g.vs["labels"][n])[0] == "ortholog":
        cnt_orth = cnt_orth + 1
        k[3].append((g.vs["properties"][n]['kid'],g.vs["properties"][n]['trend']))
        k[4].append((g.vs["properties"][n]['kid'],g.vs["properties"][n]['trend']))
        total_list[g.vs["properties"][n]['trend']] = total_list[g.vs["properties"][n]['trend']] + 1
        orth_list[g.vs["properties"][n]['trend']] = orth_list[g.vs["properties"][n]['trend']] + 1
      avail[n] = False
      for r in range(SIDE):
        if (avail[r] == True) and ((adj[n][r] > 0) or (adj[r][n] > 0))and (r not in stack_n):
          stack_n.append(r)
    community_list.append({"compound": cnt_cpd, "reaction":cnt_rea, "gene":cnt_gen,
    "ortholog": cnt_orth, "total":cnt_total, "kid_list":k})
    community_list_sketch.append({"compound": cnt_cpd, "reaction":cnt_rea, "gene":cnt_gen,
    "ortholog": cnt_orth, "compound_detail":cpd_list,"gen_list":gen_list,"ortholog_list":orth_list,"total":cnt_total})
    comm_cnt = comm_cnt + 1
  def take_total(elem):
    return elem["total"]
  community_list.sort(key=take_total, reverse=True)
  community_list_sketch.sort(key=take_total, reverse=True)
  sorted_cnt = 1
  for item in community_list:
    item["order"] = sorted_cnt
    sorted_cnt = sorted_cnt + 1
  sorted_cnt = 1
  for item in community_list_sketch:
    item["order"] = sorted_cnt
    sorted_cnt = sorted_cnt + 1
  import json
  with open("./free/result_free_"+str(THRESHOLD)+"_.json", "w") as f:
    json.dump(community_list, f)
  with open("./free/result_sketch_free_"+str(THRESHOLD)+"_.json", "w") as f:
    json.dump(community_list_sketch, f)
  print("All done")
  end_time = time.time()#'subGraph_light_'+str(THRESHOLD)+'_.graphml'
  print(str(round((end_time-start_time)/60, 2))+ "min")
  visual_style = {
    # "vertex_size" : 5,
    "bbox": (1000,1000),
    "margin": 60,
    "layout": "kk",
    # "arrow_size": 1,
    "edge_width": 0.3
    }
  node_labels_color = []
  for i in g.vs["properties"]:
      try:
          y = i['trend']
          if y == -1:
              x = "red"
          elif y == 1:
              x = "green"
          elif y == 0:
              x = "black"
          else:
              x = "grey"
      except:
          x = "black"
      node_labels_color.append(x)
  g.vs["label_color"] = node_labels_color
  g.vs["label"] = [i['kid'] for i in g.vs["properties"]]#['_nx_name':<id>, 'labels':类别, 'properties', 'label']
  col = {"compound":"orange", "reaction":"pink", "gene":"brown", "ortholog":"green"}
  siz = {"compound":15, "reaction":5, "gene":15, "ortholog":15}
  g.vs["color"] = [col[list(i)[0]] for i in g.vs["labels"]]
  g.vs["size"] = [siz[list(i)[0]] for i in g.vs["labels"]]
  g.vs["label_dist"] = [1 for i in g.vs["labels"]]
  g.vs["label_size"] = [6 for i in g.vs["labels"]]
  # g.vs["label_color"] = ["blue" for i in g.vs["labels"]]
  #[i['kid'] for i in g.vs["properties"]], [i['name'] for i in g.vs["properties"]]
  #[list(i)[0] for i in g.vs["labels"]]
  g.es['arrow_size'] = [0.5 for i in g.es['type_']]
  e_col = {'expression':"green", 'gene_expression':"green", 'repression':"red", 
  'activation':"green", 'inhibition':"red", 'substrate':"blue", 'product':"magenta", 'enzyme':"pink"}
  g.es['color'] = [e_col[i] for i in g.es['type_']]
  g.es['label'] = [i for i in g.es['type_']]
  g.es['label_color'] = [e_col[i] for i in g.es['type_']]
  g.es["label_size"] = [6 for i in g.es["type_"]]

  igraph.save(g, './free/subGraph_free_'+str(THRESHOLD)+'_.graphml')