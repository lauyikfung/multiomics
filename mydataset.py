from torch.utils.data import Dataset
import torch

class ProteinDataset(Dataset):
    """
    Process of the proteomics data from the csv files.
    """
    def __init__(self,
                 data,
                 device="cpu"):

        self.data = data
        self.device = device

    def __len__(self):
        return len(self.data)

    # @profile
    def __getitem__(self, index):

        content = self.data[index]
        return {"id": index, "target":content[0], "valid":content[1], "source": content[2][:454]}

    # @profile
    def collate_fn(self, samples):
        bsz = len(samples)
        ids = []
        source = torch.Tensor(bsz, 454)
        target = torch.Tensor(bsz, 4) * 0.0
        target_id = torch.Tensor(bsz, 1) * 0.0
        for idx, sample in enumerate(samples):
            ids.append(sample['id'])
            source[idx, :] = torch.Tensor(sample['source'])
            target[idx, sample['target']] = target[idx, sample['target']] + 1.0
            target_id[idx, 0] = target_id[idx, 0] + sample['target']
        
        return {
            "ids": torch.Tensor(ids).to(self.device),
            "source": source.to(self.device),
            "target": target.to(self.device),
            "target_id": target_id.to(self.device),
        }


class MetabolomicsDataset(Dataset):
    """
    Process of the metabolomics data from the csv files.
    """
    def __init__(self,
                 data,
                 device="cpu"):

        self.data = data
        self.device = device

    def __len__(self):
        return len(self.data)

    # @profile
    def __getitem__(self, index):

        content = self.data[index]
        return {"id": index, "target":content[0], "valid":content[1], "source": content[2][454:]}

    # @profile
    def collate_fn(self, samples):
        bsz = len(samples)
        ids = []
        source = torch.Tensor(bsz, 891)
        target = torch.Tensor(bsz, 4) * 0.0
        target_id = torch.Tensor(bsz, 1) * 0.0
        for idx, sample in enumerate(samples):
            ids.append(sample['id'])
            source[idx, :] = torch.Tensor(sample['source'])
            target[idx, sample['target']] = target[idx, sample['target']] + 1.0
            target_id[idx, 0] = target_id[idx, 0] + sample['target']
        
        return {
            "ids": torch.Tensor(ids).to(self.device),
            "source": source.to(self.device),
            "target": target.to(self.device),
            "target_id": target_id.to(self.device),
        }


class BothDataset(Dataset):
    """
    Process of data of two omics from the csv files.
    """
    def __init__(self,
                 data,
                 device="cpu"):

        self.data = data
        self.device = device

    def __len__(self):
        return len(self.data)

    # @profile
    def __getitem__(self, index):

        content = self.data[index]
        return {"id": index, "target":content[0], "valid":content[1], "source": content[2]}

    # @profile
    def collate_fn(self, samples):
        bsz = len(samples)
        ids = []
        source = torch.Tensor(bsz, 1345)
        target = torch.Tensor(bsz, 4) * 0.0
        target_id = torch.Tensor(bsz, 1) * 0.0
        valid = torch.Tensor(bsz, 2) * 0.0
        for idx, sample in enumerate(samples):
            ids.append(sample['id'])
            source[idx, :] = torch.Tensor(sample['source'])
            valid[idx, 0] = valid[idx, 0] + sample['valid'][0]
            valid[idx, 1] = valid[idx, 1] + sample['valid'][1]
            target[idx, sample['target']] = target[idx, sample['target']] + 1.0
            target_id[idx, 0] = target_id[idx, 0] + sample['target']
        
        return {
            "ids": torch.Tensor(ids).to(self.device),
            "source": source.to(self.device),
            "valid": valid.to(self.device),
            "target": target.to(self.device),
            "target_id": target_id.to(self.device),
        }
