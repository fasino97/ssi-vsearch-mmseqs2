from collections import Counter
import pandas as pd
from sklearn import metrics
import matplotlib.pyplot as plt

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

#Read Output tsv file and assign column names
df = pd.read_csv("MMSeqs2Silva/silva2024/inc/75.tsv", sep='\t', header=None)
df.columns = ['representative', 'member']

#Turn the training and testing FASTA files into dictionaries
train_dict = {}
test_dict = {}

with open('MMSeqs2Silva/silva2024/timingtesting/silvatraining2024batch12345.fasta') as f:
    trainingdata = f.read()
    trsequences = trainingdata.split(">")
    trsequences.pop(0)
    for sequence in trsequences:
        lines = sequence.split("\n")
        seq_id = lines[0].split(" ", 1)[0]
        label = " ".join(lines[0].split(" ")[1:])
        seq = "".join(lines[1:])
        
        train_dict[seq_id] = label
with open('MMSeqs2Silva/silva2024/timingtesting/silvatesting2024.fasta') as f:
    testingdata = f.read()
    tesequences = testingdata.split(">")
    tesequences.pop(0)
    for sequence in tesequences:
        lines = sequence.split("\n")
        seq_id = lines[0].split(" ", 1)[0]
        label = " ".join(lines[0].split(" ")[1:])
        seq = "".join(lines[1:])
        
        test_dict[seq_id] = label


results = {}

#Initialize counts for clusters and singletons
num_clusters = 0
predacclist = []
unlabacclist = []

#Loop through taxonomic depth levels 2 to 6
for level in range(2, 7):
    #Initialize separate dictionaries for each loop
    labeled_cluster_dict = {}
    unlabeled_cluster_dict = {}
    
    #Initialize separate lists for labeled and unlabeled clusters for each taxonomic depth level
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
        
        #Section results into clusters; adds single data point clusters to the singletons list
        if mem_str == rep_str:
           if cluster:
               if len(cluster) == 1: 
                   singletons.append(cluster[0])
                   singnum = singnum + 1
               elif any(member in train_dict for member in cluster):
                   labeled_clusters.append(cluster)
               else:
                   unlabeled_clusters.append(cluster)
           cluster = [mem_str]
        elif cluster:
           cluster.append(mem_str)

    #Update counts for clusters and singletons
    num_clusters = len(labeled_clusters) + len(unlabeled_clusters)
    
    #Give the labeled clusters labels based on majority vote
    member_clusters = []
    
    for cluster in labeled_clusters:
        assign_labels(cluster, train_dict, test_dict, level)

    #Assign the most common label from testing data to each unlabeled cluster
    if len(unlabeled_clusters) != 0:
        for cluster in unlabeled_clusters:
            testlabels = []

            for member in cluster:
                labeltest = test_dict.get(member) 
                if labeltest:
                    testlabels.append(labeltest.split(';')[level])

            if testlabels:
                counter = Counter(testlabels)
                most_common_labeltest = counter.most_common(1)[0][0]
                
                for member in cluster:
                    unlabeled_cluster_dict[member] = most_common_labeltest
    
    
    #Combine the labeled and unlabeled cluster labels
    predicted_labels = {}
    predicted_labels.update(labeled_cluster_dict)
    predicted_labels.update(unlabeled_cluster_dict)

    #Calculate labeled and unlabeled accuracy
    labeled_correct_predictions = 0
    labeled_total_predictions = 0
    unlabeled_correct_predictions = 0
    unlabeled_total_predictions = 0
    
    for key in predicted_labels:
        if key in labeled_cluster_dict:
            labeled_total_predictions += 1
            if predicted_labels[key] == test_dict.get(key, "").split(';')[level]:  # Use .get() to handle missing keys
                labeled_correct_predictions += 1
        else:
            if len(unlabeled_clusters) != 0 and len(unlabeled_clusters) != singnum:
                if key in test_dict:
                    unlabeled_total_predictions += 1
                    if predicted_labels[key] == test_dict[key].split(';')[level]:
                        unlabeled_correct_predictions += 1

    labeled_accuracy = labeled_correct_predictions / labeled_total_predictions
    results[level] = {
        "Predicted Accuracy": labeled_accuracy,
        "Clusters": num_clusters,
        "Singletons": singnum
    }

    if len(unlabeled_clusters) != 0:
        if unlabeled_total_predictions != 0:
            unlabeled_accuracy = unlabeled_correct_predictions / unlabeled_total_predictions
            results[level]["Unlabeled Accuracy"] = unlabeled_accuracy
            unlabacclist.append(unlabeled_accuracy)

    predacclist.append(labeled_accuracy)
  
    labels_true = []
    labels_pred = []
    testing_dict = {}

    for key, value in test_dict.items():
        label_list = value.split(';')
        testing_dict[key] = label_list

    for key, value in predicted_labels.items():
        if key in testing_dict:
            labels_pred.append(value)
            labels_true.append(testing_dict[key][level])

    #Compute comleteness, homogeneity, and V-measure
    homogeneity, completeness, v_measure = metrics.homogeneity_completeness_v_measure(labels_true, labels_pred)
    results[level]["Completeness"] = completeness
    results[level]["Homogeneity"] = homogeneity
    results[level]["V-measure"] = v_measure

#Print the results for each taxonomic depth level
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

#Code below is optional; provides cluster distribution plot

#clustersizelist = labeled_clusters  
#cluster_sizes = [len(cluster) for cluster in clustersizelist]
#plt.hist(cluster_sizes, bins=range(1, max(cluster_sizes) + 2), align="left", edgecolor="black")
#plt.xlabel("Cluster Size")
#plt.ylabel("Frequency")
#plt.title("Distribution of Labeled Cluster Sizes")
#plt.show()