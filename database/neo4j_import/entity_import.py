import itertools
import re
import tqdm
from typing import List, Tuple
from database import graph, mongo_entity_db, mongo_pathway_db
from py2neo.bulk import create_nodes, create_relationships
from py2neo import *
import constants as cts


MISSING_ENTITY_AUTO_DOWNLOAD = True


def create_node(data: dict):
    node = Node(data["type"], **data)
    graph.create(node)
    return node


def preprocess_entity_data_from_mongodb(origin_data: dict) -> dict:
    target_fields = cts.ENTITY_BASE_PROPS + [
        "SYMBOL",
        "NAME",
        "FORMULA",
        "EQUATION",
        "REL_PATHWAY",
        "DESCRIPTION",
        "organism",
        "REACTION",
        "HISTORY",
        "ENZYME",
        "DEFINITION",
        "MOL_WEIGHT",
        "RCLASS",
        "CLASS",
        "EXACT_MASS",
        "SYSNAME",
    ]

    data = dict()
    for i in target_fields:
        if i in origin_data.keys():
            data[i] = origin_data[i]
    data["id"] = data[cts.DB_PK]

    return data


def download_entity_and_mongo_import(entity_type: str, kid: str):
    from database.kegg_download import download_entity_detail
    from database.mongodb_import import import_single_entity_to_mongodb

    data = download_entity_detail(entity_type, kid)
    data = import_single_entity_to_mongodb(entity_type, kid, data)
    return data


def get_or_create_entity(entity_type: str, kid: str, data: dict = None) -> Node:
    if node := graph.nodes.match(kid=kid).first():
        return node
    print(f"cant find node [{entity_type}]: {kid}")

    if not data:
        data = mongo_entity_db[entity_type].find_one({cts.DB_PK: kid})
        if not data and MISSING_ENTITY_AUTO_DOWNLOAD:
            data = download_entity_and_mongo_import(entity_type, kid)
        else:
            print(f"no data: [{entity_type}]", f"[{kid}]")
            return None
    data = preprocess_entity_data_from_mongodb(data)
    node = create_node(data)
    return node


def get_single_type_entity_list(entity_type: str, kid_list: List[str]) -> List[Node]:
    entity_list = []
    for kid in kid_list:
        if entity := get_or_create_entity(entity_type, kid):
            entity_list.append(entity)
    return entity_list


def get_orthology_list_by_gene_list(gene_id_list: List[str]) -> List[Node]:
    """
    通过 gene 列表，返回其对应的蛋白质列表

    Parameters
    ----------
    gene_id_list : List[str]
        gene_id_list

    Returns
    -------
    List[Node]
        ortholog node list
    """
    data = mongo_entity_db["gene"].find(
        {"kid": {"$in": gene_id_list}}, {"ORTHOLOGY": 1}
    )
    orthology_set = set()
    for i in data:
        if "ORTHOLOGY" in i.keys():
            for o in i["ORTHOLOGY"].keys():
                orthology_set.add(f"ko:{o}")
    return get_single_type_entity_list("ortholog", orthology_set)


def _get_node_data_list_from_mongodb(
    entity_type: str, target_fields: List[str]
) -> List[dict]:
    node_data_list = []
    data = mongo_entity_db[entity_type].find({})
    for i in tqdm.tqdm(data, desc=f"load data [{entity_type}]"):
        node = {field: i.get(field) for field in target_fields}
        node_data_list.append(node)
    return node_data_list


def _batch_create_entity_op(
    entity_type: str, node_data_list: List[dict], batch_size: int = 100
):
    stream = iter(node_data_list)
    epoch = len(node_data_list) // batch_size + 1
    for _ in tqdm.tqdm(range(epoch), desc=f"batch create [{entity_type}]"):
        batch = itertools.islice(stream, batch_size)
        if batch:
            create_nodes(graph.auto(), batch, labels={entity_type})
        else:
            break


def batch_create_entity_with_type(
    entity_type: str, target_fields: List[str], batch_size: int = 100
):
    # 根据传入的entity_type 和fields ，批量在neo4j中创建node
    # 从mongodb获取数据
    node_data_list = _get_node_data_list_from_mongodb(entity_type, target_fields)
    # 批量插入
    _batch_create_entity_op(entity_type, node_data_list, batch_size)


def batch_create_all_entity():
    entity_props_comb = [
        (cts.ENTITY_GENE, cts.ENTITY_BASE_PROPS + ["organism", "SYMBOL"]),
        (cts.ENTITY_ORTHOLOG, cts.ENTITY_BASE_PROPS + ["SYMBOL"]),
        (
            cts.ENTITY_COMPOUND,
            cts.ENTITY_BASE_PROPS + ["EXACT_MASS", "FORMULA", "MOL_WEIGHT"],
        ),
        (cts.ENTITY_ENZYME, cts.ENTITY_BASE_PROPS + ["CLASS", "HISTORY", "SYSNAME"]),
        (
            cts.ENTITY_PATHWAY,
            cts.ENTITY_BASE_PROPS + ["organism", "CLASS", "REL_PATHWAY", "DESCRIPTION"],
        ),
        (
            cts.ENTITY_REACTION,
            cts.ENTITY_BASE_PROPS
            + [
                "RCLASS",
                "NAME",
                "EQUATION",
                "ENZYME",
                "DEFINITION",
            ],
        ),
        (
            cts.ENTITY_REACTION_CLASS,
            cts.ENTITY_BASE_PROPS + ["REACTION", "ENZYME", "DEFINITION"],
        ),
    ]
    for entity_type, props in entity_props_comb:
        batch_create_entity_with_type(entity_type, props)


def update_compound_alias():
    data = mongo_entity_db[cts.ENTITY_COMPOUND].find({}, {"kid": 1, "alias": 1})
    graph.update()
