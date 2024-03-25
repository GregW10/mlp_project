import matplotlib.pyplot as plt
import math

# plt.rcParams["text.usetex"] = True
# plt.rcParams["font.family"] = "Helvetica"

path = "aggregated_NAppf1bs100L13.csv"


def main():
    f = open(path, "r")
    losses = dict()
    for line in f.readlines()[1:]:
        l = line.split(',')
        losses[float(l[0])] = float(l[1])
    f.close()
    fig, ax = plt.subplots()
    title = f"Mean Prediction Error vs Time into Future"
    ax.set_title(title)
    ax.plot(losses.keys(), losses.values(), label="Test set")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Mean Error")
    ax.set_xlim([0, None])
    ax.set_ylim([0, None])
    ax.grid()
    ax.legend()
    plt.show()
    fig.savefig(f"{path.split(sep='.')[0]}.pdf")


if __name__ == "__main__":
    main()
