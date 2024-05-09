
def alpha_plot():
    import os
    import matplotlib.pyplot as plt
    import re
    import numpy as np

    alpha_diversity_output = '/scratch365/zhuang8/R16S_data/workspace/alpha_diversity_output'

    diversity_data = {'BP': [], 'F': [], 'ISi': [], 'Sh': [], 'Si': []}
    sample_names = []

    custom_legends = {
        'F': "Fisher's index",
        'BP': "Berger-parker's diversity",
        'ISi': "Simpson's Reciprocal Index",
        'Sh': "Shannon's diversity",
        'Si': "Simpson's index of diversity"
    }
    custom_colors = {
        'BP': '#5996d2',
        'F': '#d3d764',
        'ISi': '#92c6e1',
        'Sh': '#f0ccd5',
        'Si': '#dfa166'
    }

    metric_order = ['Si', 'BP', 'F', 'Sh', 'ISi']

    for file in os.listdir(alpha_diversity_output):
        file_path = os.path.join(alpha_diversity_output, file)
        if file.endswith('.txt'):
            metric = file.split('_')[-1].split('.')[0]
            base_name = file.split('_')[0]
            if base_name not in sample_names:
                sample_names.append(base_name)
            with open(file_path, 'r') as f:
                content = f.read()
                match = re.search(r": (\d+\.\d+)", content)
                if match:
                    value = float(match.group(1))
                    diversity_data[metric].append(value)

    num_samples = len(sample_names)
    num_metrics = len(diversity_data)
    indices = np.arange(num_samples)
    width = 0.15
    space = 0.03
    fig, ax = plt.subplots(figsize=(14, 8))


    for i, metric in enumerate(metric_order):
        values = diversity_data[metric]
        ax.bar(indices + i * (width + space), values, width, label=custom_legends[metric], color=custom_colors[metric])

    ax.set_ylabel('Alpha Diversity Value')
    ax.set_title('Alpha Diversity Across Samples for Different Metrics')
    ax.set_xticks(indices + (width + space) * (num_metrics - 1) / 2)
    ax.set_xticklabels(sample_names)
    ax.legend()

    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig('/scratch365/zhuang8/R16S_data/workspace/plot_data/alpha_diversity_plot.png')


def beta_plot():
    import pandas as pd
    import numpy as np
    import seaborn as sns
    import io
    import matplotlib.pyplot as plt

    file_path = '/scratch365/zhuang8/R16S_data/workspace/beta_diversity_output/output.txt'
    sample_names = []


    with open(file_path, 'r') as file:
        lines = file.readlines()
        data_lines = []
        for line in lines:
            if line.startswith('#'):
                file_name = line.split('/')[-1].split('.')[0]
                sample_names.append(file_name)
            else:
                data_lines.append(line)

    data = pd.DataFrame(data_lines)
    data = pd.read_csv(io.StringIO(''.join(data_lines)), sep='\s+', header=0, index_col=0)

    data.replace('x.xxx', np.nan, inplace=True)
    data = data.infer_objects() 
    data = data.astype(float)

    plt.figure(figsize=(12, 12))
    heatmap = sns.heatmap(data, annot=True, cmap='viridis', fmt=".3f", xticklabels=sample_names, yticklabels=sample_names)
    heatmap.set_title('Beta Diversity Heatmap', fontdict={'fontsize':12}, pad=12)
    plt.savefig("/scratch365/zhuang8/R16S_data/workspace/plot_data/beta_diversity_plot.png")

def comparison():
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    from skbio.stats.distance import permanova
    from skbio import DistanceMatrix

    file_path = '/scratch365/zhuang8/R16S_data/workspace/beta_diversity_output/output.txt'

    with open(file_path, 'r') as f:
        lines = f.readlines()

    data = [line.strip().split('\t')[1:] for line in lines[11:]]
    size = len(data)
    distance_matrix = np.zeros((size, size))
    for i in range(size):
        for j in range(size):
            if data[i][j] != 'x.xxx':
                distance_matrix[i, j] = float(data[i][j])
            if i != j and data[j][i] != 'x.xxx':
                distance_matrix[i, j] = float(data[j][i])

    print(distance_matrix)
    groups = ['Group1' if 'B1' in line else 'Group2' if 'B3' in line else 'Group3' for line in lines[1:11]]
    group_ids = [str(i) for i in range(len(groups))] 
    groups = pd.Series(groups, index=group_ids, name='Group') 
    dm = DistanceMatrix(distance_matrix, ids=group_ids)

    result = permanova(dm, groups)
    df = pd.DataFrame(distance_matrix, columns=groups, index=groups)
    print(result)

    plt.figure(figsize=(10, 8))
    sns.boxplot(data=df.melt(var_name='Group', value_name='Distance'), x='Group', y='Distance', palette="Set3")
    plt.title('PERMANOVA Results')
    plt.savefig("/scratch365/zhuang8/R16S_data/workspace/plot_data/permanova_plot.png")
    plt.show()

if __name__ == "__main__":
    comparison()
