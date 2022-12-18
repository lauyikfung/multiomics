import os
import re
from typing import List
import tqdm
import json
import pymongo
import constants as cts
from database.kegg_download import download_entity_detail, download_single_pathway_map
from database import data_base_dir, mongo_pathway_db, mongo_entity_db


def import_single_pathway_to_mongodb(col, _id, organism, data: dict):
    data[cts.DB_PK] = _id
    data["organism"] = organism
    col.update_one({cts.DB_PK: _id}, {"$set": data}, upsert=True)


def import_pathway_map_to_mongodb(organism):
    base_dir = os.path.join(data_base_dir, f"./data/kegg/pathway/{organism}/")
    col = mongo_pathway_db["pathway"]

    for i in tqdm.tqdm(os.listdir(base_dir), desc=organism):
        try:
            data = json.load(open(os.path.join(base_dir, i)))
            if not data:
                data = download_single_pathway_map(organism, i[:-5], force=True)
            _id = i[:-5]
            import_single_pathway_to_mongodb(col, _id, organism, data)
        except Exception as e:
            print(i)
            raise e


def import_entity_to_mongodb(_type: str, force: bool = True):
    col = mongo_entity_db[_type]
    if force:
        print(mongo_entity_db.drop_collection(_type))
        col.create_index([(cts.DB_PK, pymongo.HASHED)])

    base_dir = os.path.join(data_base_dir, f"data/kegg/entity/{_type}/")
    for i in tqdm.tqdm(os.listdir(base_dir), desc=_type):
        try:
            data = json.load(open(os.path.join(base_dir, i)))
            if not data:
                data = download_entity_detail(_type, i[:-5], force=True)
            _id = i[:-5]
            import_single_entity_to_mongodb(_type, _id, data)
        except Exception as e:
            print(i)
            raise e


def import_all_entity_to_mongodb():
    for _type in [
        cts.ENTITY_PATHWAY,
        cts.ENTITY_GENE,
        cts.ENTITY_ORTHOLOG,
        cts.ENTITY_ENZYME,
        cts.ENTITY_REACTION,
        cts.ENTITY_REACTION_CLASS,
        cts.ENTITY_COMPOUND,
        cts.ENTITY_DRUG,
        cts.ENTITY_GLYCAN,
        cts.ENTITY_GENOME,
        cts.ENTITY_DISEASE,
    ]:
        import_entity_to_mongodb(_type, force=True)


def process_name_list(name_list: List[str]):
    return list(map(lambda x: x.replace(";", ""), name_list))


def import_single_entity_to_mongodb(entity_type: str, kid: str, data: dict):
    col = mongo_entity_db[entity_type]

    data[cts.DB_PK] = kid
    data["type"] = entity_type

    if "NAME" in data:
        data["NAME"] = process_name_list(data["NAME"])
        data["name"] = data["NAME"][0]
    else:
        data["name"] = data[cts.DB_PK]

    if entity_type == cts.ENTITY_GENE:
        data = preprocess_gene_data(data)
    elif entity_type == cts.ENTITY_PATHWAY:
        data = preprocess_pathway_data(data)
    elif entity_type == cts.ENTITY_COMPOUND:
        data = preprocess_compound_data(data)
    elif entity_type == cts.ENTITY_ORTHOLOG:
        data = preprocess_ortholog_data(data)
    elif entity_type == cts.ENTITY_ENZYME:
        data = preprocess_enzyme_data(data)
    col.update_one({cts.DB_PK: kid}, {"$set": data}, upsert=True)
    return data


def preprocess_pathway_data(data: dict):
    data["organism"] = re.findall(r"path\:([a-zA-Z]+)\d+", data[cts.DB_PK])[0]
    data["alias"] = [data["name"]]
    return data


def preprocess_gene_data(data: dict):
    # 去除 gene 名称的无效字段
    data["name"] = data["name"].replace("(RefSeq)", "").strip()
    data["alias"] = [data["name"]]

    if "SYMBOL" in data:
        data["SYMBOL"] = list(map(lambda x: x.strip(), data["SYMBOL"].split(",")))
        data["alias"].extend(data["SYMBOL"])
    data["organism"] = data[cts.DB_PK].split(":")[0]
    return data


def preprocess_compound_data(data: dict):
    if "NAME" in data:
        data["alias"] = [i.strip() for i in data["NAME"]]
    else:
        data["alias"] = [data["name"]]
    return data


def preprocess_ortholog_data(data: dict):
    data["name"] = re.sub(r"\[.*\]", "", data["name"]).strip()
    data["alias"] = [data["name"]]

    if "SYMBOL" in data:
        data["SYMBOL"] = list(map(lambda x: x.strip(), data["SYMBOL"].split(",")))
        data["alias"].extend(data["SYMBOL"])
    return data


def preprocess_enzyme_data(data: dict):
    if "NAME" in data:
        data["alias"] = [i.strip() for i in data["NAME"]]
    else:
        data["alias"] = [data["name"]]
    return data


def update_compound_alias_with_hmdb():
    col = mongo_entity_db["hmdb_metabolites"]
    c_col = mongo_entity_db[cts.ENTITY_COMPOUND]

    _filter = {"kegg_id": {"$ne": None}, "synonyms": {"$ne": None}}

    data = col.find(_filter, {"kegg_id": 1, "synonyms": 1, "name": 1})
    err_count = 0
    for i in tqdm.tqdm(data):
        kegg_id = i["kegg_id"]
        kid = f"cpd:{kegg_id}"
        print(kid)
        synonyms = i["synonyms"]["synonym"]
        if isinstance(synonyms, str):
            synonyms = [synonyms]
        synonyms.append(i["name"])
        c_data = c_col.find_one({"kid": kid}, {"alias": 1})
        if c_data:
            alias = list(set(c_data["alias"] + synonyms))
            c_col.update_one({"kid": kid}, {"$set": {"alias": alias}})
        else:
            err_count += 1
    print(err_count)


def update_compound_alias_with_hmdb_manual():
    import pandas as pd

    df = pd.read_csv("./tmp_data/kegg_compound_map.csv")
    h_col = mongo_entity_db["hmdb_metabolites"]
    c_col = mongo_entity_db[cts.ENTITY_COMPOUND]

    for i, row in df.iterrows():
        hmdb_id = row["HMDB"]
        kegg_id = row["KEGG"]
        if kegg_id == "-":
            continue
        kid = f"cpd:{kegg_id}"
        print(kid, hmdb_id)
        alias = c_col.find_one({"kid": kid}, {"alias": 1})["alias"]
        synonyms = h_col.find_one({"accession": hmdb_id}, {"synonyms": 1})["synonyms"][
            "synonym"
        ]
        alias = list(set(alias + synonyms))
        c_col.update_one({"kid": kid}, {"$set": {"alias": alias}})
