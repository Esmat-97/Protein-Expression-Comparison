import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# تحميل البيانات
protein_data = pd.read_csv(
    "https://raw.githubusercontent.com/ProteintechLab/Python-for-Biologists/8d513c5ecd83b76a261b8e1bdbfcaa160b1aacf8/Plotting%20and%20Visualisation/mice_protein_expression.csv"
)

# --- الجزء الأول: مقارنة بروتينين ---
def compare_expression(a, b):
    if a not in protein_data.columns or b not in protein_data.columns:
        raise ValueError("One or both proteins not found in the dataset.")
    subset = protein_data[[a, b]].dropna()
    diff = subset[a] - subset[b]
    print(f"Average difference between {a} and {b}: {np.mean(diff):.2f}")
    return diff

diff = compare_expression("DYRK1A_N", "CDK5_N")

plt.figure(figsize=(6, 4))
plt.bar(range(len(diff)), diff)
plt.xlabel("Mouse")
plt.ylabel("Difference")
plt.title("Expression Difference: DYRK1A vs CDK5")
plt.savefig("difference_bar.png")
plt.show()





# --- الجزء الثاني: مقارنة مجموعة بروتينات ---
proteins = ["DYRK1A_N", "CDK5_N", "NR1_N", "NR2A_N", "BDNF_N"]
subset = protein_data[proteins].dropna()

diff_matrix = pd.DataFrame(index=proteins, columns=proteins)
for p1 in proteins:
    for p2 in proteins:
        diff_matrix.loc[p1, p2] = np.mean(subset[p1] - subset[p2])

diff_matrix = diff_matrix.astype(float)

plt.figure(figsize=(8, 6))
sns.heatmap(diff_matrix, annot=True, cmap="coolwarm", fmt=".2f")
plt.title("Average Expression Differences Between Proteins")
plt.savefig("heatmap.png")
plt.show()




# اختيار مجموعة بروتينات للتحليل
proteins = ["DYRK1A_N", "CDK5_N", "NR1_N", "NR2A_N", "BDNF_N"]

# فلترة البيانات
subset = protein_data[proteins].dropna()

# حساب مصفوفة الارتباط
corr_matrix = subset.corr()

# رسم Heatmap للارتباطات
plt.figure(figsize=(8, 6))
sns.heatmap(corr_matrix, annot=True, cmap="coolwarm", fmt=".2f")
plt.title("Correlation Heatmap of Protein Expressions")
plt.savefig("correlation_heatmap.png")  # حفظ الصورة في فولدر النتائج
plt.show()





# اختيار بروتين معين للمقارنة
protein = "BDNF_N"

# رسم Boxplot لمقارنة التعبير حسب الـ Genotype
plt.figure(figsize=(8, 6))
sns.boxplot(x="Genotype", y=protein, data=protein_data)
plt.title(f"Expression of {protein} by Genotype")
plt.savefig("group_comparison_genotype.png")
plt.show()

# رسم Boxplot لمقارنة التعبير حسب الـ Treatment
plt.figure(figsize=(8, 6))
sns.boxplot(x="Treatment", y=protein, data=protein_data)
plt.title(f"Expression of {protein} by Treatment")
plt.savefig("group_comparison_treatment.png")
plt.show()




from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler


# اختيار مجموعة بروتينات للتحليل
proteins = ["DYRK1A_N", "CDK5_N", "NR1_N", "NR2A_N", "BDNF_N"]

# فلترة البيانات
subset = protein_data[proteins].dropna()

# توحيد القيم (Scaling)
scaler = StandardScaler()
scaled_data = scaler.fit_transform(subset)

# --- PCA ---
pca = PCA(n_components=2)
pca_result = pca.fit_transform(scaled_data)

plt.figure(figsize=(8, 6))
plt.scatter(pca_result[:, 0], pca_result[:, 1], alpha=0.7)
plt.xlabel("PCA Component 1")
plt.ylabel("PCA Component 2")
plt.title("PCA Projection of Protein Expression")
plt.savefig("pca_projection.png")
plt.show()

# --- Clustering باستخدام K-means ---
kmeans = KMeans(n_clusters=3, random_state=42)
clusters = kmeans.fit_predict(scaled_data)

plt.figure(figsize=(8, 6))
plt.scatter(pca_result[:, 0], pca_result[:, 1], c=clusters, cmap="viridis", alpha=0.7)
plt.xlabel("PCA Component 1")
plt.ylabel("PCA Component 2")
plt.title("K-means Clustering on PCA Projection")
plt.savefig("pca_kmeans.png")
plt.show()
