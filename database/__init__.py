from pymongo import MongoClient
from py2neo import Graph, ConnectionUnavailable
from neo4j import GraphDatabase
from constants.db import YOUR_NEO4J_USER, YOUR_NEO4J_PORT
data_base_dir = "/home/multiomics_datasets"


# 因为本地和服务器使用可能端口不一样
try:
    neo4j = Graph(
        f"bolt://localhost:{YOUR_NEO4J_PORT}", auth=("neo4j", YOUR_NEO4J_USER))
    graph_driver = GraphDatabase.driver(
        f"bolt://localhost:{YOUR_NEO4J_PORT}", auth=("neo4j", YOUR_NEO4J_USER)
    )
    print(f"connect {YOUR_NEO4J_PORT}")
except ConnectionUnavailable:
    graph = Graph(f"http://localhost{YOUR_NEO4J_PORT}",
                  auth=("neo4j", YOUR_NEO4J_USER))
    graph_driver = GraphDatabase.driver(
        f"bolt://localhost:{YOUR_NEO4J_PORT}", auth=("neo4j", YOUR_NEO4J_USER)
    )
    print(f"connect {YOUR_NEO4J_PORT+1}")

# 目前使用默认
mongo_client = MongoClient()
mongo_entity_db = mongo_client["entity"]
mongo_pathway_db = mongo_client["pathway"]
