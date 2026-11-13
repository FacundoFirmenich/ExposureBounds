from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root = Path(__file__).resolve().parent
for lang in ('EN', 'ES'):
    es = lang == 'ES'
    fig, ax = plt.subplots(figsize=(7.2, 3.4), layout='constrained')
    ax.axhspan(7/8, 9/8, color='#83a9e3', alpha=.36, label='Brazo A' if es else 'Arm A')
    ax.axhspan(15/16, 17/16, color='#a97bd0', alpha=.52, label='Brazo B' if es else 'Arm B')
    ax.plot([0, 1], [1, 1], color='#333064', linewidth=1.5, linestyle='--',
            label='Entrada nominal' if es else 'Nominal input')
    for y, c in [(7/8, '#375fa1'), (9/8, '#375fa1'), (15/16, '#7947a2'), (17/16, '#7947a2')]:
        ax.hlines(y, 0, 1, color=c, linewidth=1.1)
        ax.scatter([0], [y], color=c, s=24, zorder=3)
    ax.set(xlim=(-.045, 1.025), ylim=(.84, 1.16), xticks=[0, .25, .5, .75, 1],
           xlabel='Tiempo normalizado' if es else 'Normalised time',
           ylabel='Exposición compatible' if es else 'Compatible exposure')
    ax.set_yticks([7/8, 15/16, 1, 17/16, 9/8], ['7/8', '15/16', '1', '17/16', '9/8'])
    ax.legend(loc='upper center', bbox_to_anchor=(.5, 1.18), ncols=3, frameon=False)
    ax.spines[['top', 'right']].set_visible(False)
    fig.savefig(root / f'ARM_ENVELOPES_{lang}.pdf')
    plt.close(fig)
