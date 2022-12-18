import pandas as pd
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch import optim
from torch.utils.data import DataLoader
from tqdm import tqdm
import argparse
import numpy as np
import os
import random
from classifier import ProteinClassifier, MetabolomicsClassifier, BothClassifier, BothClassifierWithAdjacencyMatrix, PnetImprove
from mydataset import ProteinDataset, MetabolomicsDataset, BothDataset

def get_original_data():
  #########################################
  # You should edit the part below to define the method to get the data. 
  # The data in ./covid/covid.xlsx  is an example
  data = pd.read_excel('./covid/covid.xlsx')
  pro_list = data.keys()[7:461]
  data1 = pd.read_excel('./covid/covid.xlsx', sheet_name="Metabolomics")
  meta_list = data1.keys()[7:898]
  #########################################
  pro = [[], [], [], []]
  pro_data = [[], [], [], []]
  for i in range(817):
      if data.iloc[i, 3] == "Healthy":
          pro[0].append(data.iloc[i, 0])
          pro_data[0].append(data.iloc[i, 7:])
      elif data.iloc[i, 3] == "INCOV-T1":
          pro[1].append(data.iloc[i, 0])
          pro_data[1].append(data.iloc[i, 7:])
      elif data.iloc[i, 3] == "INCOV-T2":
          pro[2].append(data.iloc[i, 0])
          pro_data[2].append(data.iloc[i, 7:])
      elif data.iloc[i, 3] == "INCOV-T3":
          pro[3].append(data.iloc[i, 0])
          pro_data[3].append(data.iloc[i, 7:])
  meta = [[], [], [], []]
  meta_data = [[], [], [], []]
  for i in range(683):
      if data1.iloc[i, 3] == "Healthy":
          meta[0].append(data1.iloc[i, 0])
          meta_data[0].append(data1.iloc[i, 7:])
      elif data1.iloc[i, 3] == "INCOV-T1":
          meta[1].append(data1.iloc[i, 0])
          meta_data[1].append(data1.iloc[i, 7:])
      elif data1.iloc[i, 3] == "INCOV-T2":
          meta[2].append(data1.iloc[i, 0])
          meta_data[2].append(data1.iloc[i, 7:])
      elif data1.iloc[i, 3] == "INCOV-T3":
          meta[3].append(data1.iloc[i, 0])
          meta_data[3].append(data1.iloc[i, 7:])
  all_data = []
  for j in range(4):
      for k in range(len(meta[j]) - 1, -1, -1):
          if meta[j][k] in pro[j]:
              idx = pro[j].index(meta[j][k])
              all_data.append((j, (1, 1),list(pro_data[j][idx]) + list(meta_data[j][k])))
              meta[j].pop(k)
              meta_data[j].pop(k)
              pro[j].pop(idx)
              pro_data[j].pop(idx)
  for j in range(4):
      for item in meta_data[j]:
          all_data.append([j, (0, 1),[0 for i in range(len(pro_list))] + list(item)])
  for j in range(4):
      for item in pro_data[j]:
          all_data.append([j, (1, 0),list(item) + [0 for i in range(len(meta_list))]])
  all_data_pro = [i for i in all_data if i[1][0] == 1]
  all_data_meta = [i for i in all_data if i[1][1] == 1]
  return all_data, all_data_pro, all_data_meta    
        

@torch.no_grad()
def evaluate(model, dataset, split="valid", both = False):
    model.eval()
    batch_size = 20
    dataloader = DataLoader(dataset,
                            batch_size=batch_size,
                            collate_fn=dataset.collate_fn)
    accs = []
    losses = []
    for samples in dataloader:
        bsz = len(samples['ids'])
        if both:
            logits = model.forward(samples["source"], samples["valid"])
        else:
            logits = model.forward(samples["source"])
        lprobs = F.log_softmax(logits, dim=-1).view(-1, logits.size(-1))
        entropy = F.nll_loss(lprobs.float(),
                             samples["target_id"].view(-1).long(),
                             reduction="none").view(bsz, -1)
        acc = (samples["target_id"].view(-1) == torch.max(lprobs.float(), 1)[1])
        accs.extend(acc.tolist())
        losses.append(entropy.mean().item())
    print("%s: loss: %.3f, acc: %.3f" %
          (split, np.mean(losses), np.mean(accs)))
    return np.mean(losses), np.mean(accs)


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--batch-size", default=4, type=int)
    parser.add_argument("--lr", default=1e-4, type=float)
    parser.add_argument("--weight-decay", default=1e-3, type=float)
    parser.add_argument("--num-epoch", default=30, type=int)
    parser.add_argument("--save-interval", default=1, type=int)
    parser.add_argument("--gpunum", default=0, type=int)
    parser.add_argument("--seed", default=2022, type=int)
    parser.add_argument("--save-dir", default="./models")
    parser.add_argument("--type",
                        help='Choose Data Type',
                        default="both",
                        choices=["protein", "metabolomics", "both", "with_adj_matrix", "pnet_improve"])
    parser.add_argument("--scheduler",
                        default="None",
                        choices=["OneCycleLR", "CosineAnnealingLR", "MultiStepLR", "None"])
    args = parser.parse_args()
    return args


def train(args):
    args.save_dir += "_" + args.type
    all_data, all_data_pro, all_data_meta = get_original_data()
    os.makedirs(args.save_dir, exist_ok=True)
    device = "cuda:" + str(args.gpunum) if torch.cuda.is_available() else "cpu"
    print(device)
    random.seed(args.seed)
    random.shuffle(all_data_pro)
    random.shuffle(all_data_meta)
    random.shuffle(all_data)

    #####################################################################
    # You should edit the following arguments if use other datasets
    data_size = {
        'total' : 1345,
        'protein' : 454,
        'metabolomics' : 891,
        'class_num' : 4,
        'node_num' : 1345,
        'module_num': 400,
        'pathway_num' : 300,
    }
    adj_matrix=torch.ones(1345, 1345)# default, you should edit it for adaption of your own dataset, same as below
    node_node_matrix = torch.ones(1345, 1345)# default
    node_module_matrix = torch.ones(1345, 400)# default
    module_module_matrix = torch.ones(400, 400)# default
    module_pathway_matrix = torch.ones(400, 300)# default
    pathway_pathway_matrix = torch.ones(300, 300)# default
    if args.type == 'protein':
        train_set = ProteinDataset(all_data_pro[:500], device=device)#806
        valid_set = ProteinDataset(all_data_pro[500:630], device=device)
        test_set = ProteinDataset(all_data_pro[630:], device=device)
        model = ProteinClassifier(args, data_size).to(device)
    elif args.type == 'metabolomics':
        train_set = MetabolomicsDataset(all_data_meta[:420], device=device)
        valid_set = MetabolomicsDataset(all_data_meta[420:530], device=device)
        test_set = MetabolomicsDataset(all_data_meta[530:], device=device)#672
        model = MetabolomicsClassifier(args, data_size).to(device)
    elif args.type == "both":
        train_set = BothDataset(all_data[:600], device=device)
        valid_set = BothDataset(all_data[600:750], device=device)
        test_set = BothDataset(all_data[750:], device=device)#959
        model = BothClassifier(args, data_size).to(device)
    elif args.type == "with_adj_matrix":
        train_set = BothDataset(all_data[:600], device=device)
        valid_set = BothDataset(all_data[600:750], device=device)
        test_set = BothDataset(all_data[750:], device=device)#959
        model = BothClassifierWithAdjacencyMatrix(args, data_size, device, adj_matrix).to(device)
    elif args.type == "pnet_improve":
        train_set = BothDataset(all_data[:600], device=device)
        valid_set = BothDataset(all_data[600:750], device=device)
        test_set = BothDataset(all_data[750:], device=device)#959
        model = PnetImprove(args, data_size, device, 
        node_node_matrix = node_node_matrix,
        node_module_matrix = node_module_matrix,
        module_module_matrix = module_module_matrix,
        module_pathway_matrix = module_pathway_matrix,
        pathway_pathway_matrix = pathway_pathway_matrix).to(device)
    #####################################################################

    optimizer = optim.Adam(model.parameters(),
                           lr=args.lr,
                           weight_decay=args.weight_decay)
    train_loader = DataLoader(train_set,
                                  batch_size=args.batch_size,
                                  collate_fn=train_set.collate_fn,
                                  shuffle=False)
    if args.scheduler == "CosineAnnealingLR":
        train_scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
            optimizer, args.num_epoch, eta_min=0, last_epoch=- 1, verbose=False)
    if args.scheduler == "MultiStepLR":
        train_scheduler = torch.optim.lr_scheduler.MultiStepLR(optimizer, [
                                                               10, 20], 0.5)
    if args.scheduler == "OneCycleLR":
        train_scheduler = torch.optim.lr_scheduler.OneCycleLR(
            optimizer, max_lr=2 * args.lr, epochs=args.num_epoch, steps_per_epoch=len(train_loader), three_phase=False)    
    evaluate(model, valid_set, both = (args.type == 'both') or (args.type == "with_adj_matrix"))
    valid_losses, valid_accs, test_losses, test_accs = [], [], [], []
    for epoch in range(args.num_epoch):
        model.train()
        with tqdm(train_loader, desc="training") as pbar:
            losses = []
            for samples in pbar:
                optimizer.zero_grad()
                loss = model.get_loss(**samples)
                loss.backward()
                optimizer.step()
                losses.append(loss.item())
                pbar.set_description("Epoch: %d, Loss: %0.8f, lr: %0.6f" %
                                     (epoch + 1, np.mean(losses) / args.batch_size ,
                                      optimizer.param_groups[0]['lr']))
                if args.scheduler == "OneCycleLR":
                    train_scheduler.step()
            if args.scheduler in ["CosineAnnealingLR", "MultiStepLR"]:
                train_scheduler.step()
        if epoch % args.save_interval == 0:
            torch.save(
                model,
                args.save_dir + "/{}.pt".format(epoch + 1))
        v_loss, v_acc = evaluate(model, valid_set, both = (args.type == 'both') or (args.type == "with_adj_matrix"))
        valid_losses.append(v_loss)
        valid_accs.append(v_acc)
        t_loss, t_acc = evaluate(model, test_set, split="test", both = (args.type == 'both') or (args.type == "with_adj_matrix"))
        test_losses.append(t_loss)
        test_accs.append(t_acc)
    print(min(valid_losses), max(valid_accs), min(test_losses), max(test_accs))


if __name__ == "__main__":
    args = get_args()
    train(args)
