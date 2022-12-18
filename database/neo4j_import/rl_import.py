import itertools
import re
from typing import List, Set, Tuple

import constants as cts
import tqdm
from database import graph, mongo_entity_db, mongo_pathway_db
from database.neo4j_import.entity_import import (
    get_or_create_entity,
    get_orthology_list_by_gene_list,
    get_single_type_entity_list,
)
from py2neo import *
from py2neo.bulk import create_nodes, create_relationships


def create_relationship(n1: Node, rel_name: str, n2: Node, **rel_props):
    relation_ship = Relationship(n1, rel_name, n2, **rel_props)
    graph.create(relation_ship)
    return relation_ship


def create_relationship_in_gene_ortholog_single(gene_data: dict):
    """
    创建基因对蛋白质的表达关系

    Parameters
    ----------
    gene_data : dict
        gene 的数据，大部分情况下 gene 对 蛋白质 的表达都是一对一的关系，
        有极少数情况存在 gene对蛋白为一对多
    """
    if "ORTHOLOGY" not in gene_data.keys():
        return
    gene = get_or_create_entity("gene", gene_data[cts.DB_PK])

    for k, v in gene_data["ORTHOLOGY"].items():
        ortholog = get_or_create_entity(cts.ENTITY_ORTHOLOG, f"ko:{k}")
        if gene and ortholog:
            create_relationship(
                gene, cts.RL_GENE_EXPRESSION, ortholog, trend=cts.TREND_POSITIVE
            )
        else:
            print(gene_data[cts.DB_PK], k)


def create_relationship_in_gene_ortholog():
    # 创建gene 对 蛋白质的表达关系
    data = mongo_entity_db[cts.ENTITY_GENE].find({}, {cts.DB_PK: 1, "ORTHOLOGY": 1})
    for i in tqdm.tqdm(data):
        create_relationship_in_gene_ortholog_single(i)


def create_relationship_in_ortholog_enzyme_single(enzyme_data: dict):
    """
    创建 蛋白-酶家族 关系

    Parameters
    ----------
    enzyme_data : dict
        酶的数据，8005个酶家族中，有一半 蛋白质和酶家族为一一对应，其他为多对一
    """
    if "ORTHOLOGY" not in enzyme_data.keys():
        return
    enzyme = get_or_create_entity(cts.ENTITY_ENZYME, enzyme_data[cts.DB_PK])

    for k, v in enzyme_data["ORTHOLOGY"].items():
        ortholog = get_or_create_entity(cts.ENTITY_ORTHOLOG, f"ko:{k}")

        if enzyme and ortholog:
            create_relationship(ortholog, cts.RL_ENZYME_FAMILY, enzyme)
        else:
            print(enzyme_data[cts.DB_PK], k)


def create_relationship_in_ortholog_enzyme():
    # 创建 蛋白-酶家族 关系
    data = mongo_entity_db[cts.ENTITY_ENZYME].find({}, {cts.DB_PK: 1, "ORTHOLOGY": 1})
    for i in tqdm.tqdm(data):
        create_relationship_in_ortholog_enzyme_single(i)


def create_relationship_in_rclass_reaction_single(rclass_data: dict):
    """
    创建 rclass 与 reaction的对应

    Parameters
    ----------
    rclass_data : dict
        rclass 的数据, 共有3170个rclass，即reaction的归类。reaction与 rclass为多对一
    """
    if "REACTION" not in rclass_data.keys():
        return
    rclass = get_or_create_entity(cts.ENTITY_REACTION_CLASS, rclass_data[cts.DB_PK])

    for k in rclass_data["REACTION"]:
        reaction = get_or_create_entity(cts.ENTITY_REACTION, f"rn:{k}")
        if rclass and reaction:
            create_relationship(reaction, cts.RL_IN_REACTION_CLASS, rclass)
        else:
            print(rclass_data[cts.DB_PK], k)


def create_relationship_in_rclass_reaction():
    # 创建 rclass 与 reaction的对应
    data = mongo_entity_db[cts.ENTITY_REACTION_CLASS].find(
        {}, {cts.DB_PK: 1, "REACTION": 1}
    )
    for i in tqdm.tqdm(data):
        create_relationship_in_rclass_reaction_single(i)


def process_equation(equation: str) -> Tuple[List[str], List[str]]:
    """
    处理化学方程式，给出化学方程式，返回对应的底物和产物列表

    Parameters
    ----------
    equation : str
        化学方程式的字符串，形如： C19080 + C00002 <=> C19085 + C00013

    Returns
    -------
    Tuple[List[str], List[str]]
        对应的底物和产物列表
    """
    substrate_list, product_list = equation.split("<=>")

    # 去除括号的内容
    substrate_list = re.subn(r"\(.*\)", "", substrate_list)[0]
    product_list = re.subn(r"\(.*\)", "", product_list)[0]

    substrate_list = list(map(lambda x: x.strip(), substrate_list.split("+")))
    product_list = list(map(lambda x: x.strip(), product_list.split("+")))

    substrate_list = list(map(lambda x: x.split(" ")[-1], substrate_list))
    product_list = list(map(lambda x: x.split(" ")[-1], product_list))

    return substrate_list, product_list


def create_relationship_in_reaction():
    # 创建reaction中的关系，包括底物、产物、酶
    col = mongo_entity_db["reaction"]
    for data in tqdm.tqdm(col.find({}, {cts.DB_PK: 1})):
        try:
            create_relationship_in_single_reaction(data[cts.DB_PK])
        except Exception as e:
            continue
            # print(data)
            # print(e)
            # raise e


def create_single_node_in_reaction(entity_id: str):
    if entity_id.startswith("C"):
        node = get_or_create_entity("compound", f"cpd:{entity_id}")
    elif entity_id.startswith("K"):
        node = get_or_create_entity("ortholog", f"ko:{entity_id}")
    elif entity_id.startswith("G"):
        return None
    else:
        node = get_or_create_entity("enzyme", f"ec:{entity_id}")
    return node


def _create_relationship_in_reaction_substrate(reaction, substrate_list):
    substrate_list = [create_single_node_in_reaction(i) for i in substrate_list]
    for substrate in substrate_list:
        if substrate:
            graph.create(Relationship(substrate, cts.RL_SUBSTRATE, reaction))


def _create_relationship_in_reaction_product(reaction, product_list):
    product_list = [create_single_node_in_reaction(i) for i in product_list]
    for product in product_list:
        if product:
            graph.create(Relationship(product, cts.RL_PRODUCT, reaction))


def _create_relationship_in_reaction_ortholog(reaction, ortholog_list):
    ortholog_list = [create_single_node_in_reaction(i) for i in ortholog_list]
    for ortholog in ortholog_list:
        if ortholog:
            graph.create(Relationship(ortholog, cts.RL_ENZYME, reaction))


def _create_relationship_in_reaction_enzyme(reaction: dict, enzyme_list: List):
    """
    在reaction 未指定 ORTHOLOGY 时，会参考 ENZYME 信息。
    通过 ENZYME 寻找 ORTHOLOGY，与 reaction 构建关系，如果 ENZYME没有对应的 ORTHOLOGY，则不构建关系

    Parameters
    ----------
    reaction : dict
        reaction data
    enzyme_list : _type_
        _description_
    """

    def _get_ortholog_list_by_enzyme(_enzyme_id: str) -> Set[str]:
        if not _enzyme_id.startswith("ec"):
            _enzyme_id = f"ec:{_enzyme_id}"
        enzyme = mongo_entity_db[cts.ENTITY_ENZYME].find_one({cts.DB_PK: _enzyme_id})
        return set(enzyme.get("ORTHOLOGY", {}).keys())

    def _create_relationship_in_reaction_ortholog_list(_ortholog_list: List[Node]):
        for ortholog in ortholog_list:
            if ortholog:
                graph.create(Relationship(ortholog, cts.RL_ENZYME, reaction))

    ortholog_set = set()
    for enzyme_id in enzyme_list:
        ortholog_set |= _get_ortholog_list_by_enzyme(enzyme_id)

    ortholog_list = [create_single_node_in_reaction(i) for i in ortholog_set]
    _create_relationship_in_reaction_ortholog_list(ortholog_list)
    # if len(ortholog_list) == 0:
    #     print(reaction["kid"])


def create_relationship_in_single_reaction(reaction_id: str):
    """
    创建在一个reaction中的关系，主要是通过reaction中的EQUATION，构建 底物、产物、酶 的关系。
    需要注意的是，在酶与reaction的关系构建中，我们是使用其对应的蛋白质（即 ORTHOLOGY）与 reaction构建的关系。
    如果 reaction 中有提及 ORTHOLOGY，则直接使用 ORTHOLOGY 构建。
    如果 reaction 中未提及 ORTHOLOGY，则参考 ENZYME，使用 ENZYME 信息中关联的 ORTHOLOGY 与reaction构建关系。
    如果 reaction 中未指定 ENZYME， 或 ENZYME 未指定其对应的 ORTHOLOGY，则不构建此关系。
    Parameters
    ----------
    reaction_id : str
        reaction_id，通过此id在 mongo 中获取详细信息
    """
    data = mongo_entity_db["reaction"].find_one({cts.DB_PK: reaction_id})

    reaction = get_or_create_entity("reaction", data[cts.DB_PK])
    substrate_list, product_list = process_equation(data["EQUATION"])

    _create_relationship_in_reaction_substrate(reaction, substrate_list)
    _create_relationship_in_reaction_product(reaction, product_list)

    # 如果 reaction中指定了ORTHOLOGY
    if "ORTHOLOGY" in data.keys():
        orthology_list = list(data["ORTHOLOGY"].keys())
        _create_relationship_in_reaction_ortholog(reaction, orthology_list)
    # 如果ORTHOLOGY 不在，则探索ENZYME
    elif "ENZYME" in data.keys():
        _create_relationship_in_reaction_enzyme(reaction, data["ENZYME"])


def _create_rel_in_entity_list(
    rel_name: str, entity_list_1: List[Node], entity_list_2: List[Node]
):
    for n1 in entity_list_1:
        for n2 in entity_list_2:
            if n1 and n2:
                create_relationship(n1, rel_name, n2)


def create_ecrel_link(rel: dict, entry1: dict, entry2: dict):
    """
    ecrel-> enzyme-enzyme relation, indicating two enzymes catalyzing successive reaction steps.
    entry1 和 entry2 都为 gene对应的蛋白质。如果是group则跳过

    Parameters
    ----------
    rel : dict
        rel data
    entry1 : dict
        entry1 data
    entry2 : dict
        entry2 data
    """
    if entry1["type"] == "gene" and entry2["type"] == "gene":
        node_list_1 = get_orthology_list_by_gene_list(entry1["name"])
        node_list_2 = get_orthology_list_by_gene_list(entry2["name"])
        _create_rel_in_entity_list(rel["name"], node_list_1, node_list_2)
    else:
        print(rel["link"], rel["name"], entry1["type"], entry2["type"])


def create_pprel_link(rel: dict, entry1: dict, entry2: dict):
    # todo: 确认 pprel 和 ecrel 有啥区别？
    if entry1["type"] == "gene" and entry2["type"] == "gene":
        node_list_1 = get_orthology_list_by_gene_list(entry1["name"])
        node_list_2 = get_orthology_list_by_gene_list(entry2["name"])

        _create_rel_in_entity_list(rel["name"], node_list_1, node_list_2)
    else:
        print(rel["link"], rel["name"], entry1["type"], entry2["type"])


def create_gerel_link(rel: dict, entry1: dict, entry2: dict):
    """
    gene expression interaction, indicating relation of transcription factor and target gene product
    在 gerel 关系中， 起始点为gene对应的蛋白质，终点为gene。
    Parameters
    ----------
    rel : dict
        rel
    entry1 : dict
        entry1
    entry2 : dict
        entry2
    """
    if entry1["type"] == "gene" and entry2["type"] == "gene":
        node_list_1 = get_orthology_list_by_gene_list(entry1["name"])
        node_list_2 = get_single_type_entity_list("gene", entry2["name"])

        _create_rel_in_entity_list(rel["name"], node_list_1, node_list_2)

    else:
        print(rel["link"], rel["name"], entry1["type"], entry2["type"])


def create_pcrel_link(rel: dict, entry1: dict, entry2: dict):
    """
    protein-compound interaction, 即蛋白质与化合物之间的关系。gene会转为蛋白质关联。
    Parameters
    ----------
    rel : dict
        rel
    entry1 : dict
        entry1
    entry2 : dict
        entry2
    """

    if entry1["type"] == "gene":
        node_list_1 = get_orthology_list_by_gene_list(entry1["name"])
    else:
        node_list_1 = get_single_type_entity_list(entry1["type"], entry1["name"])

    if entry2["type"] == "gene":
        node_list_2 = get_orthology_list_by_gene_list(entry2["name"])
    else:
        node_list_2 = get_single_type_entity_list(entry2["type"], entry2["name"])

    _create_rel_in_entity_list(rel["name"], node_list_1, node_list_2)


def create_single_relationship_in_pathway(rel: dict, entry_map: dict):
    rel_name = rel["name"]

    if rel_name in ("compound", "hidden compound"):
        return

    entry1 = entry_map[rel["entry1"]]
    entry2 = entry_map[rel["entry2"]]

    if entry1["type"] == "group" or entry2["type"] == "group":
        return

    if rel["link"] == "ECrel":
        create_ecrel_link(rel, entry1, entry2)
        pass
    elif rel["link"] == "PPrel":
        create_pprel_link(rel, entry1, entry2)
        pass
    elif rel["link"] == "GErel":
        create_gerel_link(rel, entry1, entry2)
        pass
    elif rel["link"] == "PCrel":
        create_pcrel_link(rel, entry1, entry2)
        pass

    # node_list_1 = [get_or_create_entity(entry1["type"], i) for i in entry1["name"]]
    # node_list_2 = [get_or_create_entity(entry2["type"], i) for i in entry2["name"]]


def create_single_pathway_relation(relations: List[dict], entry_map: dict) -> dict:
    target_entity_type = ("compound", "gene", "ortholog")
    for rel in tqdm.tqdm(relations):
        create_single_relationship_in_pathway(rel, entry_map)
        continue
        if rel["name"] in ("compound", "hidden compound"):
            continue
        entry1 = entry_map[rel["entry1"]]
        entry2 = entry_map[rel["entry2"]]

        if (
            entry1["type"] not in target_entity_type
            or entry2["type"] not in target_entity_type
        ):
            continue

        node_list_1 = [get_or_create_entity(entry1["type"], i) for i in entry1["name"]]
        node_list_2 = [get_or_create_entity(entry2["type"], i) for i in entry2["name"]]

        for n1 in node_list_1:
            for n2 in node_list_2:
                if n1 and n2:
                    create_relationship(n1, rel["name"], n2)


def create_single_pathway_entity_relation(pathway_id: str, entry_map: dict):
    target_entity_type = ("compound", "gene", "ortholog")

    pathway = get_or_create_entity("pathway", pathway_id)

    for entry in entry_map.values():
        if entry["type"] not in target_entity_type:
            continue
        for name in entry["name"]:
            entity = get_or_create_entity(entry["type"], name)
            if entity and entity != pathway:
                create_relationship(entity, cts.RL_IN_PATHWAY, pathway)


def create_single_pathway(pathway_id: dict):
    pathway = mongo_pathway_db["pathway"].find_one({cts.DB_PK: pathway_id})

    entry_map = dict()
    for entry in pathway["entries"]:
        entry["name"] = entry["name"].split()
        entry_map[entry["id"]] = entry

    create_single_pathway_relation(pathway["relations"], entry_map)
    create_single_pathway_entity_relation(pathway[cts.DB_PK], entry_map)


def create_all_pathway_relation(organism="hsa"):
    col = mongo_pathway_db["pathway"]
    data = col.find({"organism": organism}, {cts.DB_PK: 1})
    kid_list = [i[cts.DB_PK] for i in data]

    for kid in kid_list:
        print(kid)
        try:
            create_single_pathway(kid)
        except Exception as e:
            print(kid)
            raise e
