
from transformers import AutoTokenizer, AutoModelForMaskedLM, TrainingArguments, AutoModelForSequenceClassification, AutoConfig ,AutoModelForCausalLM
import torch

from sklearn.metrics import matthews_corrcoef, f1_score
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
import numpy as np

model_path = "/scratch/lb4489/LLM"

tokenizer = AutoTokenizer.from_pretrained("PoetschLab/GROVER", cache_dir=model_path)
model = AutoModelForMaskedLM.from_pretrained("PoetschLab/GROVER", cache_dir=model_path)

num_labels_promoter = 6

model = AutoModelForSequenceClassification.from_pretrained("PoetschLab/GROVER", num_labels=num_labels_promoter, cache_dir=model_path,trust_remote_code=True)


import pandas as pd

# 读取 TSV 文件
df = pd.read_csv("/scratch/lb4489/project/mttRNA/Aund_trainset/trainset.csv", sep=",")

# 显示前几行
print(df.head())

df_selected = df[['seq', 'group']].copy()

print(df_selected["group"].value_counts())

from sklearn.model_selection import train_test_split
from datasets import load_dataset, Dataset

df_selected.columns = ["data", "label"]  
df_selected["label"] = df_selected["label"].astype(int)

train_data, val_data, train_labels, val_labels = train_test_split(
    df_selected["data"], df_selected["label"], test_size=0.2, random_state=42
)

train_df = pd.DataFrame({
    "data": train_data,
    "label": train_labels
})

max_count = train_df['label'].value_counts().max()
max_count = max_count// 3

df_upsampled = pd.concat([
    train_df,
    *[
        train_df[train_df['label'] == label].sample(max_count - len(train_df[train_df['label'] == label]), replace=True)
        for label in train_df['label'].unique() if len(train_df[train_df['label'] == label]) < max_count
    ]
])

print(df_upsampled['label'].value_counts())

train_data = df_upsampled['data']

train_labels = df_upsampled['label']

ds_train = Dataset.from_dict({"data": train_data.tolist(), "label": train_labels.tolist()})
ds_validation = Dataset.from_dict({"data": val_data.tolist(), "label": val_labels.tolist()})

print(ds_train)
print(ds_validation)

def tokenize_function(examples):
    outputs = tokenizer(examples["data"])
    return outputs

tokenized_datasets_train = ds_train.map(
    tokenize_function,
    batched=True,
    remove_columns=["data"],
)
tokenized_datasets_validation = ds_validation.map(
    tokenize_function,
    batched=True,
    remove_columns=["data"],
)

import numpy as np
import torch
import torch.nn.functional as F
from sklearn.metrics import accuracy_score, f1_score, recall_score
from transformers import TrainingArguments, Trainer

def compute_metrics(eval_pred):
    logits, labels = eval_pred
    predictions = np.argmax(logits, axis=-1)

    # 计算 Accuracy, F1-score, Recall
    acc = accuracy_score(labels, predictions)
    f1 = f1_score(labels, predictions, average="macro")

    return {
        "accuracy": acc,
        "f1_score": f1
    }

device = torch.device("cuda")
model = model.to(device)

batch_size = 128
model_name='GROVER'

args_promoter = TrainingArguments(
    f"{model_name}-finetuned",
    remove_unused_columns=False,
    evaluation_strategy="steps",
    save_strategy="steps",
    learning_rate=8e-5,
    per_device_train_batch_size=batch_size,
    gradient_accumulation_steps=1,
    per_device_eval_batch_size=256,
    logging_steps=500,
    load_best_model_at_end=True,
    metric_for_best_model="f1_score",
    label_names=["labels"],
    dataloader_drop_last=True,
    max_steps=30000,   
    weight_decay=0.01,
    lr_scheduler_type="cosine_with_restarts",
    bf16=True,
    warmup_ratio=0.05, 
)

trainer = Trainer(
    model.to(device),
    args=args_promoter,
    train_dataset=tokenized_datasets_train,
    eval_dataset=tokenized_datasets_validation,
    tokenizer=tokenizer,
    compute_metrics=compute_metrics,  # 使用新的 metrics 计算函数
)

trainer.train() 

curve_evaluation_f1_score = [[a['step'], a['eval_f1_score']] for a in trainer.state.log_history if 'eval_f1_score' in a]

curve_evaluation_accuracy = [[a['step'], a['eval_accuracy']] for a in trainer.state.log_history if 'eval_accuracy' in a]


eval_f1_score = [c[1] for c in curve_evaluation_f1_score]
eval_accuracy = [c[1] for c in curve_evaluation_accuracy]

steps = [c[0] for c in curve_evaluation_f1_score] 

plt.figure(figsize=(10, 5))

plt.plot(steps, eval_f1_score, 'b', marker="o", linestyle="-", label='Validation F1 Score')

plt.plot(steps, eval_accuracy, 'g', marker="o", linestyle="--", label='Validation Accuracy')

plt.title('Evaluation Metrics Over Training Steps')
plt.xlabel('Number of Training Steps')
plt.ylabel('Score')
plt.legend()
plt.grid(True)

#plt.show()

plt.savefig('/scratch/lb4489/project/mttRNA/Aund_trainset/GROVER_evaluation_metrics.png', dpi=300, bbox_inches='tight')


from datasets import Dataset
import pandas as pd

dfv = pd.read_csv("/scratch/lb4489/project/mttRNA/Aund_trainset/valinset.csv", sep=",")
dfv = dfv[["seq", "group"]]

#label_mapping = {0:0,1:1,2:2,3: 3, 4: 3, 5: 3}  # 重新映射小类
#dfv["group"] = dfv["group"].map(label_mapping)
true_labels = dfv["group"].astype(int).tolist()


ds_validation = Dataset.from_dict({
    "data": dfv["seq"].tolist(),
})

def tokenize_function(examples):
    outputs = tokenizer(examples["data"])
    return outputs

tokenized_datasets_vali = ds_validation.map(
    tokenize_function,
    batched=True,
    remove_columns=["data"],
)

from transformers import Trainer

trainer = Trainer(model=model, tokenizer=tokenizer)

predictions = trainer.predict(tokenized_datasets_vali)

import numpy as np
import torch
import torch.nn.functional as F
from sklearn.metrics import accuracy_score, f1_score, recall_score
from transformers import TrainingArguments, Trainer

def compute_metrics(eval_pred):
    logits, labels = eval_pred
    predictions = np.argmax(logits, axis=-1)

    # 计算 Accuracy, F1-score, Recall
    acc = accuracy_score(labels, predictions)
    f1 = f1_score(labels, predictions, average="macro")

    return {
        "accuracy": acc,
        "f1_score": f1
    }

logits = predictions.predictions

results = compute_metrics((logits, true_labels))
print(results)

from sklearn.metrics import confusion_matrix
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.special import softmax

probs = softmax(logits, axis=-1)  # 进行 softmax 转换
preds = np.argmax(probs, axis=-1)

cm = confusion_matrix(true_labels, preds)

plt.figure(figsize=(8, 6))
sns.heatmap(cm, annot=True, fmt="d", cmap="Blues")
plt.xlabel("Predicted Label")
plt.ylabel("True Label")
plt.title("Confusion Matrix")
plt.show()

plt.savefig('/scratch/lb4489/project/mttRNA/Aund_trainset/GROVER_Confusion.png', dpi=300, bbox_inches='tight')

import pandas as pd

label_mapping = {0:0,1:1,2:2,3:3, 4:4, 5:5}

# 先转换为 Pandas Series，再进行映射
true_labels = pd.Series(true_labels).map(label_mapping).tolist()

print(true_labels[:10])

import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
from scipy.special import softmax

probs = softmax(logits, axis=-1)  

from sklearn.preprocessing import label_binarize
num_classes = probs.shape[1] 
true_labels_onehot = label_binarize(true_labels, classes=np.arange(num_classes))

plt.figure(figsize=(8, 6))
for i in range(num_classes):
    fpr, tpr, _ = roc_curve(true_labels_onehot[:, i], probs[:, i]) 
    roc_auc = auc(fpr, tpr) 
    
    plt.plot(fpr, tpr, lw=2, label=f'Class {i} (AUC = {roc_auc:.2f})')

plt.plot([0, 1], [0, 1], 'k--', lw=1)  
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("ROC Curve for Each Class")
plt.legend(loc="lower right")
plt.show()

plt.savefig('/scratch/lb4489/project/mttRNA/Aund_trainset/GROVER_AUC.png', dpi=300, bbox_inches='tight')

from sklearn.metrics import precision_recall_curve, auc

num_classes = probs.shape[1] 
true_labels_onehot = label_binarize(true_labels, classes=np.arange(num_classes))

plt.figure(figsize=(8, 6))
for i in range(num_classes):
    precision, recall, _ = precision_recall_curve(true_labels_onehot[:, i], probs[:, i])  # 计算 Precision-Recall
    pr_auc = auc(recall, precision)  

    plt.plot(recall, precision, lw=2, label=f'Class {i} (PR-AUC = {pr_auc:.2f})')

plt.plot([0, 1], [0.5, 0.5], 'k--', lw=1)  
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel("Recall")
plt.ylabel("Precision")
plt.title("Precision-Recall Curve for Each Class")
plt.legend(loc="lower left")
plt.show()

plt.savefig('/scratch/lb4489/project/mttRNA/Aund_trainset/GROVER_PRAUC.png', dpi=300, bbox_inches='tight')

