from typing import List
from constants.entities import ENTITY_COMPOUND, ENTITY_GENE, ENTITY_ORTHOLOG
from database import mongo_entity_db
import constants as cts


def entity_name_to_id(entity_type: str, entity_name: str, organism: str = None) -> str:
    entity_name = entity_name.strip()

    def _get_filter_by_name(_name):
        return {"$regex": f"^{_name}$", "$options": "i"}

    _filter = {"alias": _get_filter_by_name(entity_name)}

    if entity_type == cts.ENTITY_GENE:
        assert organism, "未指定 organism"
        _filter["organism"] = organism

    # if entity_type == ENTITY_COMPOUND:
    #     _filter = {"NAME": _get_filter_by_name(entity_name)}
    # elif entity_type == ENTITY_GENE:
    #     assert organism, "未指定 organism"
    #     _filter = {"SYMBOL": _get_filter_by_name(entity_name), "organism": organism}
    # elif entity_type == ENTITY_ORTHOLOG:
    #     _filter = {"SYMBOL": _get_filter_by_name(entity_name)}
    # else:
    #     _filter = {"NAME": _get_filter_by_name(entity_name)}

    if data := mongo_entity_db[entity_type].find_one(_filter):
        return data[cts.DB_PK]
    else:
        if entity_type == ENTITY_COMPOUND:
            if data := mongo_entity_db[entity_type].find_one(
                {"NAME": _get_filter_by_name(f"L-{entity_name}")}
            ):
                return data[cts.DB_PK]
    print(f"can't find {entity_type} entity with {entity_name}")


def batch_entity_name_to_id(entity_type: str, entity_name_list: List[str]):
    results = []
    for entity_name in entity_name_list:
        if entity_id := entity_name_to_id(entity_type, entity_name):
            results.append(entity_id)
    return results


def add_entity_alias(entity_type: str, kid: str, name: str):
    mongo_entity_db[entity_type].update_one(
        {"kid": kid}, {"$addToSet": {"alias": name}}
    )
