import pandas as pd
from typing import List

# metabolomics.xlsx 的第一列是代谢物名称，可转换成KEGG代谢物编号
df_metabolomics = pd.read_excel("./tmp_data/ZIKV课题/metabolomics.xlsx")
# mRNA-file.xlsx 的第一列是基因名，像“X9230116L04Rik”这种的本身不能转化成基因ID，可以不用管，我们只处理能够转化成基因ID的数据。
df_mrna = pd.read_excel("./tmp_data/ZIKV课题/mRNA-file.xlsx")
# protein-file.xlsx的第一列是蛋白名(在KEGG中也是以基因的ID显示出来的)，例如在这个文件中也存在Angel2，但是表示的是基因翻译的蛋白，基因与蛋白质之间大部分是一一对应关系可以由基因ID转化成蛋白ID
df_protein = pd.read_excel("./tmp_data/ZIKV课题/protein-file.xlsx")
# phosphoprotein-file.xlsx的第一列是Protein的NCBI ID，第二列是对应的基因名称Gene Symbol，第三列是磷酸化位点
df_phosphoprotein = pd.read_excel("./tmp_data/ZIKV课题/phosphoprotein-file.xlsx")


def get_trend(
    name_list: List[str],
    p_values: List[float],
    fdr_values: List[float],
    log2fc_values: List[float],
    positive_threshold: float = 1,
    negative_threshold: float = 1,
):
    """
    一般认为校正后的BH P-value或FDR<0.05 认为是在两者之间有显著变化的，
    log2Foldchange>1认为是上调，log2Foldchange<-1认为是下调，其他认为是没有显著变化，
    这个1的阈值也可以适当调整。

    Parameters
    ----------
    name_list : List[str]
        _description_
    p_values : List[float]
        _description_
    fdr_values : List[float]
        _description_
    log2fc_values : List[float]
        _description_
    positive_threshold : float, optional
        _description_, by default 1
    negative_threshold : float, optional
        _description_, by default 1
    """
    df = pd.DataFrame()
    df["name"] = name_list
    df["p_values"] = p_values
    df["fdr_values"] = fdr_values
    df["log2fc_values"] = log2fc_values

    df_target = df[(df["p_values"] < 0.05) | (df["fdr_values"] < 0.05)]
    df_positive = df_target[df_target["log2fc_values"] > positive_threshold]
    df_negative = df_target[df_target["log2fc_values"] < negative_threshold]

    return list(df_positive["name"]), list(df_negative["name"])


positive, negative = get_trend(
    df_protein["Unnamed: 0"],
    df_protein["p_value"],
    df_protein["p_value_fdr"],
    df_protein["log2_FC"],
)


positive, negative = get_trend(
    df_mrna["Unnamed: 0"],
    df_mrna["pvalue"],
    df_mrna["pvalue_fdr"],
    df_mrna["log2FC"],
)
