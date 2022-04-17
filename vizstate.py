import sys
import imageio
import os
import matplotlib.pyplot as plt

def get_fname(sdir, itr):
    return sdir + "state" + str(itr) + ".txt"
def get_pname(pdir, itr):
    return pdir + "state" + str(itr) + ".png"

def make_gif(numitrs, pdir):

    with imageio.get_writer('output.gif', mode='I') as writer:
        for i in range(numitrs):
            image = imageio.imread(get_pname(pdir, i))
            writer.append_data(image)
            os.remove(get_pname(pdir, i))

statedir = "./states/"
plotdir = "./plots/"
numitrs = int(sys.argv[1])

if not os.path.exists(plotdir):
    os.makedirs(plotdir)

for i in range(numitrs):

    file = get_fname(statedir, i)
    x, y, f2v = [], [], []
    ball_pos, ball_radius = [], 0.0
    status = 0

    with open(file, 'r') as f:

        for line in f:
            line = line.strip().rstrip(')').split('(')[-1]

            if line=="pos": continue
            elif line=="f2v": status = 1
            elif line=="ball_pos": status = 2
            elif line=="ball_radius": status = 3

            else:

                arr = line.split()
                arr = [float(th) for th in arr]
                if status==0:
                    x.append(arr[0])
                    y.append(arr[1])
                elif status==1: f2v.append(arr)
                elif status==2: ball_pos.append(arr[0])
                else: ball_radius = arr[0]

    plt.clf()
    plt.triplot(x, y, f2v)

    ax = plt.gca()
    c1 = plt.Circle(ball_pos, ball_radius, color = 'k')
    ax.add_artist(c1)

    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.savefig(get_pname(plotdir,i))

make_gif(numitrs, plotdir)
os.rmdir(plotdir)
