# 该节点 在某个pathway中
RL_IN_PATHWAY = "in_pathway"

# gene 对 蛋白质的表达
RL_GENE_EXPRESSION = "gene_expression"
# 蛋白质 隶属于 酶家族
RL_ENZYME_FAMILY = "enzyme_family"
# reaction 归属于某个类
RL_IN_REACTION_CLASS = "in_rclass"

# 在reaction中的关系
# reaction的底物
RL_SUBSTRATE = "substrate"
# reaction的产物
RL_PRODUCT = "product"
# 在 reaction 中作为酶的蛋白质
RL_ENZYME = "enzyme"

# 在pathway 中的关系
# 激活作用
RL_ACTIVATION = "activation"
# 抑制作用
RL_INHIBITION = "inhibition"
# 促进表达
RL_EXPRESSION = "expression"
# 抑制表达
RL_REPRESSION = "repression"

# 间接影响，不知细节
RL_INDIRECT_EFFECT = "indirect effect"
# 状态改变
RL_STATE_CHANGE = "state change"
# 关联和分离
RL_BINDING_ASSOCIATION = "binding/association"
RL_DISSOCIATION = "dissociation"
# 由于突变等而缺少相互作用。
RL_MISSING_INTERACTION = "missing interaction"

# 一些分子反应
# 磷酸化作用
RL_PHOSPHORYLATION = "phosphorylation"
# 脱磷酸作用
RL_DEPHOSPHORYLATION = "dephosphorylation"
# 糖基化
RL_GLYCOSYLATION = "glycosylation"
# 泛素化
RL_UBIQUITINATION = "ubiquitination"
# 甲基化作用
RL_METHYLATION = "methylation"

# 对于关系的属性，即对于一条边的start 和 end, start节点丰度变化会如何影响end节点
# positive 即 如果start节点丰度增加，end节点也会呈现增加的趋势，negative相反。
TREND_POSITIVE = "trend_positive"
TREND_NEGATIVE = "trend_negative"

# 具有积极影响的关系列表，即 两点变化关系同向
POSITIVE_RL_LIST = [RL_GENE_EXPRESSION, RL_ACTIVATION, RL_EXPRESSION]

# 具有消积影响的关系列表，即 两点变化关系反向
NEGATIVE_RL_LIST = [RL_INHIBITION, RL_REPRESSION]
