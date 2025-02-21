from collections import Counter
import pandas as pd
from sklearn import metrics
import matplotlib.pyplot as plt

#Input file paths for your VSEARCH .cd-hit file, training data fasta, and testing data fasta
VSEARCHresults = 'temp/members.cd-hit'
train_fasta= 'RDP/RDPtraining.fasta'
test_fasta= 'RDP/RDPtest.fasta'

#Function to assign labels based on majority vote
def assign_labels(cluster, train_dict, test_dict, level):
    cluster_labels = []
    sub_cluster_tr = []
    sub_cluster_te = []
    most_common_label = ''

    for member in cluster:
        if member in train_dict:
            sub_cluster_tr.append(member)
        elif member in test_dict:
            sub_cluster_te.append(member)

    for i in sub_cluster_tr:
        label = train_dict[i]
        cluster_labels.append(label.split(';')[level])
    
    if cluster_labels:
        counter = Counter(cluster_labels)
        most_common_label = counter.most_common(1)[0][0]

    for j in sub_cluster_te:
        labeled_cluster_dict[j] = most_common_label

def parse_cdhit_clstr(file_path):
    cluster_data = {'representative': [], 'member': []}
    
    current_cluster = None
    representative_sequence = None

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()

            if line.startswith('>Cluster'):
                current_cluster = int(line.split()[1])
                representative_sequence = None
            elif line.endswith('*'):
                representative_sequence = line[:-1].split('>')[1].split('...')[0].strip()
                cluster_data['representative'].append(representative_sequence)
                cluster_data['member'].append(representative_sequence)  
            elif current_cluster is not None:
                member_sequence = line.split('>')[1].split('...')[0].strip()
                cluster_data['representative'].append(representative_sequence)
                cluster_data['member'].append(member_sequence)

    return pd.DataFrame(cluster_data)

df = parse_cdhit_clstr(VSEARCHresults)

#Turn the training and testing FASTA files into dictionaries
train_dict = {}
test_dict = {}

def clean_label(label):
    return label.strip()

def process_fasta(data, data_dict):
    sequences = data.split(">")
    sequences.pop(0)
    for sequence in sequences:
        lines = sequence.strip().split("\n")
        if lines:
            seq_info = lines[0].split(None, 1)
            seq_id = seq_info[0].strip()
            label = clean_label(seq_info[1]) if len(seq_info) > 1 else ""
            data_dict[seq_id] = label

with open(train_fasta) as f:
    process_fasta(f.read(), train_dict)

with open(test_fasta) as f:
    process_fasta(f.read(), test_dict)

results = {}
num_clusters = 0
predacclist = []
unlabacclist = []

#Loop to run calculations for every taxonomic depth
for level in range(2, 7):
    labeled_cluster_dict = {}
    unlabeled_cluster_dict = {}
    
    labeled_clusters = []
    unlabeled_clusters = []
    labeled_singletons = []
    unlabeled_singletons = []  
    cluster = []
    singletons= []
    singnum = 0
    
    for index, row in df.iterrows():
        rep_str = row["representative"]
        mem_str = row["member"]
        
        if mem_str == rep_str:
           if cluster:
               if len(cluster) == 1: 
                   singletons.append(cluster[0])  
                   singnum += 1
               elif any(member in train_dict for member in cluster):
                   labeled_clusters.append(cluster)
               else:
                   unlabeled_clusters.append(cluster)
           cluster = [mem_str]
        elif cluster:
           cluster.append(mem_str)

    num_clusters = len(labeled_clusters) + len(unlabeled_clusters)

    for cluster in labeled_clusters:
        assign_labels(cluster, train_dict, test_dict, level)

    if unlabeled_clusters:
        for cluster in unlabeled_clusters:
            testlabels = [test_dict[member].split(';')[level] for member in cluster if member in test_dict]

            if testlabels:
                counter = Counter(testlabels)
                most_common_labeltest = counter.most_common(1)[0][0]
                
                for member in cluster:
                    unlabeled_cluster_dict[member] = most_common_labeltest
    
    predicted_labels = {**labeled_cluster_dict, **unlabeled_cluster_dict}

    labeled_correct_predictions = sum(1 for key in labeled_cluster_dict 
                                      if test_dict.get(key, "").split(';')[level] == labeled_cluster_dict[key])
    labeled_total_predictions = len(labeled_cluster_dict)

    unlabeled_correct_predictions = sum(1 for key in unlabeled_cluster_dict 
                                        if test_dict.get(key, "").split(';')[level] == unlabeled_cluster_dict[key])
    unlabeled_total_predictions = len(unlabeled_cluster_dict)

    labeled_accuracy = labeled_correct_predictions / labeled_total_predictions if labeled_total_predictions else 0
    results[level] = {
        "Predicted Accuracy": labeled_accuracy,
        "Clusters": num_clusters,
        "Singletons": singnum
    }

    if unlabeled_total_predictions:
        unlabeled_accuracy = unlabeled_correct_predictions / unlabeled_total_predictions
        results[level]["Unlabeled Accuracy"] = unlabeled_accuracy
        unlabacclist.append(unlabeled_accuracy)

    predacclist.append(labeled_accuracy)
  
    labels_true = []
    labels_pred = []
    testing_dict = {key: value.split(';') for key, value in test_dict.items()}

    for key, value in predicted_labels.items():
        if key in testing_dict:
            labels_pred.append(value)
            labels_true.append(testing_dict[key][level])

    homogeneity, completeness, v_measure = metrics.homogeneity_completeness_v_measure(labels_true, labels_pred)
    results[level]["Completeness"] = completeness
    results[level]["Homogeneity"] = homogeneity
    results[level]["V-measure"] = v_measure

for level, result in results.items():
    print(f"Results for Level {level}:")
    print("Predicted Accuracy:", result["Predicted Accuracy"])
    print("Clusters:", result["Clusters"])
    print("Singletons:", result["Singletons"])
    if "Unlabeled Accuracy" in result:
        print("Unlabeled Accuracy:", result["Unlabeled Accuracy"])
    print("Completeness:", result["Completeness"])
    print("Homogeneity:", result["Homogeneity"])
    print("V-measure:", result["V-measure"])
    print()

#Optional: Plot cluster distribution
#clustersizelist = labeled_clusters  
#cluster_sizes = [len(cluster) for cluster in clustersizelist]
#plt.hist(cluster_sizes, bins=range(1, max(cluster_sizes) + 2), align="left", edgecolor="black")
#plt.xlabel("Cluster Size")
#plt.ylabel("Frequency")
#plt.title("Distribution of Labeled Cluster Sizes")
#plt.show()
