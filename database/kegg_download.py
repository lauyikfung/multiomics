import itertools
import json
import os
from typing import List

import constants as cts
import tqdm
from bioservices import KEGG

from database import data_base_dir


def batch_download_entity_detail(_type: str, id_list: List[str]):
    k = KEGG(verbose=True, cache=True)
    res = k.get("+".join(id_list))
    data_list = [k.parse(i.strip()) for i in res.split("\n///") if i.strip()]
    # print(len(id_list), len(data_list))
    for _id, data in zip(id_list, data_list):
        path = os.path.join(data_base_dir, f"./data/kegg/entity/{_type}/{_id}.json")
        with open(path, "w") as fin:
            json.dump(data, fin)


def download_entity_detail(_type: str, _id: str, force: bool = False) -> dict:
    path = os.path.join(data_base_dir, f"./data/kegg/entity/{_type}/{_id}.json")

    # 如果声明强制更新 或者 文件不存在 或者文件为空
    if force or not os.path.exists(path) or not open(path).read():
        k = KEGG(verbose=True, cache=True)
        with open(path, "w") as fin:
            data = k.parse(k.get(_id))
            json.dump(data, fin)
            return data


def batch_download_all_entity_by_type(
    entity_type: str, id_list: List[str], batch_size: int = 10
):
    id_len = len(id_list)

    stream = iter(id_list)
    for _ in tqdm.tqdm(
        range(id_len // batch_size + 1), desc=f"download [{entity_type}]"
    ):
        batch = itertools.islice(stream, batch_size)
        try:
            batch_download_entity_detail(entity_type, list(batch))
        except Exception as e:
            print(batch)
            print(e)
            continue


def download_all_entity_by_type(entity_type: str, ids: List[str], force: bool = False):
    for i in tqdm.tqdm(ids, desc=entity_type):
        download_entity_detail(entity_type, i, force=force)


def download_all_type_entity(force: bool = False):
    k = KEGG(verbose=True, cache=True)

    entity_group = [
        (cts.ENTITY_COMPOUND, k._get_database(cts.ENTITY_COMPOUND)),
        (cts.ENTITY_ENZYME, k._get_database(cts.ENTITY_ENZYME)),
        (cts.ENTITY_REACTION, k._get_database(cts.ENTITY_REACTION)),
        (cts.ENTITY_REACTION_CLASS, k._get_database(cts.ENTITY_REACTION_CLASS)),
        (cts.ENTITY_DRUG, k._get_database(cts.ENTITY_DRUG)),
        (cts.ENTITY_ORTHOLOG, k.koIds),
        (cts.ENTITY_GLYCAN, k._get_database(cts.ENTITY_GLYCAN)),
        (cts.ENTITY_DISEASE, k._get_database(cts.ENTITY_DISEASE)),
        (cts.ENTITY_GENOME, k._get_database(cts.ENTITY_GENOME)),
    ]
    for eg in entity_group:
        download_all_entity_by_type(*eg, force=force)


def download_single_pathway_map(organism: str, _id: str, force: bool = False) -> dict:
    k = KEGG(verbose=True, cache=True)
    path = os.path.join(data_base_dir, f"./data/kegg/pathway/{organism}/{_id}.json")

    if force or not os.path.exists(path):
        with open(path, "w") as fin:
            data = k.parse_kgml_pathway(_id)
            json.dump(data, fin)
            return data


# 下载一个物种所有的pathway
def download_pathway_map_with_organism(organism: str):
    k = KEGG(verbose=True, cache=True)
    k.organism = organism
    for i in tqdm.tqdm(k.pathwayIds, desc="download [pathway]"):

        download_single_pathway_map(organism, i)
        download_entity_detail(cts.ENTITY_PATHWAY, i)


# 批量下载一个物种所有pathway中包含的gene
def batch_download_gene_by_organism(organism: str = "hsa", batch_size: int = 10):
    k = KEGG(verbose=True, cache=True)
    genes = [x.split("\t")[0] for x in k.list(organism).strip().split("\n")]

    id_len = len(genes)
    stream = iter(genes)

    for _ in tqdm.tqdm(range(id_len // batch_size + 1), desc=f"download [gene]"):

        # for i in tqdm.tqdm(range(len(genes) // batch_size + 1)):
        batch = itertools.islice(stream, batch_size)
        # gene_batch = genes[i * batch_size : (i + 1) * batch_size]
        try:
            batch_download_entity_detail(cts.ENTITY_GENE, batch)
        except Exception as e:
            print(batch)
            print(e)
            continue


# 下载一个物种所有pathway中包含的gene
def download_gene_by_organism(organism: str = "hsa"):
    k = KEGG(verbose=True, cache=True)
    genes = [x.split("\t")[0] for x in k.list(organism).strip().split("\n")]
    for gene in tqdm.tqdm(genes):
        download_entity_detail("gene", gene)
