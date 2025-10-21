import matplotlib.pyplot as plt
import numpy as np

# Dataset names
datasets = ['PRJNA544129', 'PRJNA62742', 'PRJNA43443']

# QC and Mapping percentages
qc_pass = [81.74, 99.77, 99.44]
mapping_rate = [82.82, 98.04, 97.70]

# Bar chart setup
x = np.arange(len(datasets))
width = 0.35

fig, ax = plt.subplots(figsize=(10, 6))
bars1 = ax.bar(x - width/2, qc_pass, width, label='QC Passed (%)', color='#4CAF50')
bars2 = ax.bar(x + width/2, mapping_rate, width, label='Mapped (%)', color='#2196F3')

# Labels and title
ax.set_ylabel('Percentage (%)')
ax.set_title('Quality Control and Mapping Rates by Dataset')
ax.set_xticks(x)
ax.set_xticklabels(datasets)
ax.set_ylim(0, 110)
ax.legend()

# Add text labels on bars
for bar in bars1 + bars2:
    height = bar.get_height()
    ax.annotate(f'{height:.2f}%',
                xy=(bar.get_x() + bar.get_width() / 2, height),
                xytext=(0, 3),
                textcoords="offset points",
                ha='center', va='bottom', fontsize=9)

plt.tight_layout()
plt.show()
