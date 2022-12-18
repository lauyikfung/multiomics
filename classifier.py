import torch
import torch.nn as nn
import torch.nn.functional as F


class ProteinClassifier(nn.Module):
  """
  Only for the proteomics part of Covid-19 data.
  """
  def __init__(self, args, data_size):
    super().__init__()
    self.fc1 = nn.Linear(data_size['protein'], 256)
    self.fc2 = nn.Linear(256, 64)
    self.fc3 = nn.Linear(64, data_size['class_num'])

  def forward(self, inputs):
    outputs = self.fc1(inputs)
    outputs = F.leaky_relu(outputs, 0.2)
    outputs = self.fc2(outputs)
    outputs = torch.sigmoid(outputs)
    outputs = self.fc3(outputs)
    return outputs

  def get_loss(self, source, target_id, reduce=True, **unused):
    outputs = self.forward(source)
    lprobs = F.log_softmax(outputs, dim=-1).view(-1, outputs.size(-1))
    loss = F.nll_loss(
        lprobs,
        target_id.view(-1).long(),
        reduction="sum" if reduce else "none",
    )
    return loss

class MetabolomicsClassifier(nn.Module):
  """
  Only for the metabolomics part of Covid-19 data.
  """
  def __init__(self, args, data_size):
    super().__init__()
    self.fc1 = nn.Linear(data_size['metabolomics'], 256)
    self.fc2 = nn.Linear(256, 64)
    self.fc3 = nn.Linear(64, data_size['class_num'])

  def forward(self, inputs):
    outputs = self.fc1(inputs)
    outputs = F.leaky_relu(outputs, 0.2)
    outputs = self.fc2(outputs)
    outputs = torch.sigmoid(outputs)
    outputs = self.fc3(outputs)
    return outputs

  def get_loss(self, source, target_id, reduce=True, **unused):
    outputs = self.forward(source)
    lprobs = F.log_softmax(outputs, dim=-1).view(-1, outputs.size(-1))
    loss = F.nll_loss(
        lprobs,
        target_id.view(-1).long(),
        reduction="sum" if reduce else "none",
    )
    return loss


class BothClassifier(nn.Module):
  """
  For both proteomics and metabolomics parts of Covid-19 data.
  Not only include two classifiers, but also consider the inaccessible of one part of data.
  The output is the linear combination of two classifiers whose weights are normalized and computed by a third network.
  """
  def __init__(self, args, data_size):
    super().__init__()
    self.data_size = data_size
    self.pro_classifier = ProteinClassifier(args, data_size)
    self.meta_classifier = MetabolomicsClassifier(args, data_size)
    self.fc1 = nn.Linear(data_size['total'], data_size['total'])
    self.fc2 = nn.Linear(data_size['total'], data_size['total'] // 4)
    self.fc3 = nn.Linear(data_size['total'] // 4, 64)
    self.fc4 = nn.Linear(64, data_size['class_num'] * 2)

  def forward(self, inputs, valid):
    inputs_2 = torch.sigmoid(self.fc1(inputs)) + inputs
    w = F.leaky_relu(inputs_2, 0.2)
    w = self.fc2(w)
    w = torch.sigmoid(w)
    w = self.fc3(w)
    w = torch.sigmoid(w)
    w = self.fc4(w)
    p = self.pro_classifier(inputs[:, :self.data_size['protein']])
    m = self.meta_classifier(inputs[:, self.data_size['protein']:])
    w = torch.reshape(w, (-1, self.data_size['class_num'], 2))
    w = F.softmax(w, dim=-1)
    w_p, w_m = w[:, :, 0], w[:, :, 1]
    w_p = torch.mul(w_p.T, valid[:, 0]).T
    w_m = torch.mul(w_m.T, valid[:, 1]).T
    outputs = w_p * p + w_m * m
    return outputs


  def get_loss(self, source, valid, target_id, reduce=True, **unused):
    outputs = self.forward(source, valid)
    lprobs = F.log_softmax(outputs, dim=-1).view(-1, outputs.size(-1))
    loss = F.nll_loss(
        lprobs,
        target_id.view(-1).long(),
        reduction="sum" if reduce else "none",
    )
    return loss
  
  
class BothClassifierWithAdjacencyMatrix(nn.Module):
  """
  For both proteomics and metabolomics parts of Covid-19 data, but with interface of adjacency matrix.
  Not only include two classifiers, but also consider the inaccessible of one part of data.
  The output is the linear combination of two classifiers whose weights are normalized and computed by a third network.
  """
  def __init__(self, args, data_size, device, adj_matrix=torch.ones(1345, 1345)):
    super().__init__()
    self.data_size = data_size
    self.pro_classifier = ProteinClassifier(args, data_size)
    self.meta_classifier = MetabolomicsClassifier(args, data_size)
    self.fc1 = nn.Linear(data_size['total'], data_size['total'])
    self.fc2 = nn.Linear(data_size['total'], data_size['total'] // 4)
    self.fc3 = nn.Linear(data_size['total'] // 4, 64)
    self.fc4 = nn.Linear(64, data_size['class_num'] * 2)
    self.adj_matrix = adj_matrix.to(device)

  def forward(self, inputs, valid):
    inputs_2 = torch.sigmoid(self.fc1(inputs).matmul(self.adj_matrix)) + inputs
    w = F.leaky_relu(inputs_2, 0.2)
    w = self.fc2(w)
    w = torch.sigmoid(w)
    w = self.fc3(w)
    w = torch.sigmoid(w)
    w = self.fc4(w)
    p = self.pro_classifier(inputs[:, :self.data_size['protein']])
    m = self.meta_classifier(inputs[:, self.data_size['protein']:])
    w = torch.reshape(w, (-1, self.data_size['class_num'], 2))
    w = F.softmax(w, dim=-1)
    w_p, w_m = w[:, :, 0], w[:, :, 1]
    w_p = torch.mul(w_p.T, valid[:, 0]).T
    w_m = torch.mul(w_m.T, valid[:, 1]).T
    outputs = w_p * p + w_m * m
    return outputs


  def get_loss(self, source, valid, target_id, reduce=True, **unused):
    outputs = self.forward(source, valid)
    lprobs = F.log_softmax(outputs, dim=-1).view(-1, outputs.size(-1))
    loss = F.nll_loss(
        lprobs,
        target_id.view(-1).long(),
        reduction="sum" if reduce else "none",
    )
    return loss

 
class PnetImprove(nn.Module):
  """
  Improvement of PNET, and waiting for better dataset.
  """
  def __init__(self, args, data_size, device, 
    node_node_matrix = torch.ones(1345, 1345),
    node_module_matrix = torch.ones(1345, 400),
    module_module_matrix = torch.ones(400, 400),
    module_pathway_matrix = torch.ones(400, 300),
    pathway_pathway_matrix = torch.ones(300, 300)):
    super().__init__()
    self.data_size = data_size

    self.node_node = nn.Linear(data_size['node_num'], data_size['node_num'])
    self.node_module = nn.Linear(data_size['node_num'], data_size['module_num'])
    self.module_module = nn.Linear(data_size['module_num'], data_size['module_num'])
    self.module_pathway = nn.Linear(data_size['module_num'], data_size['pathway_num'])
    self.pathway_pathway =nn.Linear(data_size['pathway_num'], data_size['pathway_num'])
    self.projection = nn.Linear(data_size['pathway_num'], data_size['class_num'])

    self.op1 = nn.Linear(data_size['node_num'], data_size['class_num'])
    self.op2 = nn.Linear(data_size['module_num'], data_size['class_num'])
    self.op3 = nn.Linear(data_size['module_num'], data_size['class_num'])
    self.op4 = nn.Linear(data_size['pathway_num'], data_size['class_num'])
    self.op5 = nn.Linear(data_size['pathway_num'], data_size['class_num'])

    self.node_node_matrix = node_node_matrix.to(device)
    self.node_module_matrix = node_module_matrix.to(device)
    self.module_module_matrix = module_module_matrix.to(device)
    self.module_pathway_matrix = module_pathway_matrix.to(device)
    self.pathway_pathway_matrix = pathway_pathway_matrix.to(device)

  def forward(self, inputs, valid):
    x = self.node_node(inputs).matmul(self.node_node_matrix)
    output_1 = torch.sigmoid(self.op1(x))
    x = self.node_module(torch.tanh(x)).matmul(self.node_module_matrix)
    output_2 = torch.sigmoid(self.op2(x))
    x = self.module_module(torch.tanh(x)).matmul(self.module_module_matrix)
    output_3 = torch.sigmoid(self.op3(x))
    x = self.module_pathway(torch.tanh(x)).matmul(self.module_pathway_matrix)
    output_4 = torch.sigmoid(self.op4(x))
    x = self.pathway_pathway(torch.tanh(x)).matmul(self.pathway_pathway_matrix)
    output_5 = torch.sigmoid(self.op5(x))
    x = self.projection(torch.tanh(x))
    output_6 = torch.sigmoid(x)
    outputs = (output_1 + output_2 + output_3 + output_4 + output_5 + output_6) / 6
    return outputs


  def get_loss(self, source, valid, target_id, reduce=True, **unused):
    outputs = self.forward(source, valid)
    lprobs = F.log_softmax(outputs, dim=-1).view(-1, outputs.size(-1))
    loss = F.nll_loss(
        lprobs,
        target_id.view(-1).long(),
        reduction="sum" if reduce else "none",
    )
    return loss
