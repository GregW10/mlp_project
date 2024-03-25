import matplotlib.pyplot as plt
import math
import sys

# plt.rcParams["text.usetex"] = True
# plt.rcParams["font.family"] = "Helvetica"

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
    title = f"Training and Validation Errors vs Epoch Number"
    ax.set_title(title)
    last_val = 50
    ax.plot(list(tlosses.keys())[:last_val], list(tlosses.values())[:last_val], label="Training set")
    ax.plot(list(vlosses.keys())[:last_val], list(vlosses.values())[:last_val], label="Val. set")
    ax.set_xlabel("Epoch Number")
    ax.set_ylabel("Mean Square Error")
    # ax.set_xlim([0, None])
    ax.set_ylim([0, None])
    ax.grid()
    ax.legend()
    plt.show()
    if len(sys.argv) == 2:
        fig.savefig(sys.argv[1])
    else:
        fig.savefig(f"figure.pdf")


if __name__ == "__main__":
    main()
