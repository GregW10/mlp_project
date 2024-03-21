import matplotlib.pyplot as plt
import math

plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "Helvetica"

tpath = "train_losses.csv"
vpath = "val_losses.csv"


def main():
    tf = open(tpath, "r")
    vf = open(vpath, "r")
    tlosses = dict()
    for line in tf.readlines()[1:]:
        l = line.split(',')
        tlosses[int(l[0])] = float(l[1])
    vlosses = dict()
    for line in vf.readlines()[1:]:
        l = line.split(',')
        vlosses[int(l[0])] = float(l[1])
    tf.close()
    vf.close()
    fig, ax = plt.subplots()
    title = f"Training and Validation Losses vs Epoch Number"
    ax.set_title(title)
    ax.plot(tlosses.keys(), tlosses.values(), label="Training set")
    ax.plot(vlosses.keys(), vlosses.values(), label="Val. set")
    ax.set_xlabel("Epoch Number")
    ax.set_ylabel("Loss")
    # ax.set_xlim([0, None])
    ax.set_ylim([0, None])
    ax.grid()
    ax.legend()
    plt.show()
    fig.savefig(f"figure.png")


if __name__ == "__main__":
    main()
