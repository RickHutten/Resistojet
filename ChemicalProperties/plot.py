import matplotlib.pyplot as plt
from reader import FileReader, Files


def plot_specific_heat_ratio_isobaric(filename, plot_cp=True):
    data = FileReader(filename)

    # Plot 1
    fig, ax1 = plt.subplots()
    if plot_cp:
        color = 'b'
    else:
        color = 'k'
    ax1.plot(data.temperature, data.gamma, color, label='Specific Heat Ratio')
    ax1.set_xlabel(data.titles.temperature)
    ax1.set_ylabel(data.titles.gamma)

    if plot_cp:
        # Plot 2
        ax2 = ax1.twinx()
        ax2.plot(data.temperature, data.cp, 'r', label='Cp')
        ax2.set_ylabel(data.titles.cp)
        fig.legend(loc="upper center", bbox_to_anchor=(0.5, 0.9))

    plt.title('Specific heat ratio of water (g) at 6 bar')
    fig.tight_layout()
    plt.show()


plot_specific_heat_ratio_isobaric(Files.isobaric_water_6_bar, True)
