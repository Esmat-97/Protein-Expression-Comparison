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

