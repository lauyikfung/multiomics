from typing import List


def entity_name_to_id(entity_type: str, entity_name: str):
    pass


def generate_subgraph_query(
    target_node_id: List[str],
    target_relationship_list: List[str] = [
        "activation",
        "inhibition",
        "expression",
        "repression",
        "product",
        "substrate",
        "enzyme",
    ],
    full_mode: bool = False,
    max_length: int = None,
):
    """
    根据提供的目标 实体及关注的关系，生产子图的查询语句，需在noe4j中运行的到结果。
    Parameters
    ----------
    target_node_id : List[str]
        关注的实体 对应 kegg的id。目前不支持直接输入名字，因为可能无法对应上，人工评测阶段这样更容易发现问题。
    target_relationship_list : List[str], optional
        关注的关系，默认使用reaction的关系及4个常见关系，可指定更多
        by default [ "activation", "inhibition", "expression", "repression", "product", "substrate", "enzyme", ]
    full_mode : bool, optional
        有两种查询方式，一种为查询一条最短连通路径，一种为查询所有最短连通路径, by default False，即返回一条。
    max_length : int, optional
        可以接受的节点间最长连接长度,默认不指定。 by default None，

    Returns
    -------
    _type_
        _description_
    """
    if max_length:
        max_length = f"..{max_length}"
    else:
        max_length = ""
    if full_mode:
        algo = "allShortestPaths"
    else:
        algo = "ShortestPath"
    query = (
        """WITH $target_node_id$ AS node_id_list
    UNWIND node_id_list AS ids
    MATCH (n {kid: ids})
    WITH collect(n) AS nds

    UNWIND nds AS n1
    UNWIND nds AS n2
    WITH nds, n1, n2 WHERE id(n1) <> id(n2)
    MATCH path = $algo$((n1)-[rels*$max_length$]-(n2))
    WHERE all(r IN relationships(path) where type(r) in $target_relationship_list$)
    RETURN path ORDER BY length(path) ASC;""".replace(
            "$target_node_id$", str(target_node_id)
        )
        .replace("$target_relationship_list$", str(target_relationship_list))
        .replace("$max_length$", max_length)
        .replace("$algo$", algo)
    )
    print(query)
    return query
