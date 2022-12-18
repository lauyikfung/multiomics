from pymongo import MongoClient
from py2neo import Graph, ConnectionUnavailable
from neo4j import GraphDatabase

data_base_dir = "/home/libiao/multiomics_datasets"


# 因为本地和服务器使用可能端口不一样
try:
    graph = Graph("http://localhost:7474", auth=("neo4j", "libiao"))
    graph_driver = GraphDatabase.driver(
        "bolt://localhost:7687", auth=("neo4j", "libiao")
    )
    print("connect 7474")
except ConnectionUnavailable:
    graph = Graph("http://localhost:7475", auth=("neo4j", "libiao"))
    graph_driver = GraphDatabase.driver(
        "bolt://localhost:7687", auth=("neo4j", "libiao")
    )
    print("connect 7475")

# 目前使用默认
mongo_client = MongoClient()
mongo_entity_db = mongo_client["entity"]
mongo_pathway_db = mongo_client["pathway"]
