# このプログラムは、天頂角分布に従うμ粒子がある距離だけ離れた板の間を移動しどのくらい衝突するかを
# MonteCarlo法によって計算するためのプログラムです。
# 乱数はn個だけ用意して使用する(イベントはn回発生させる)一方、それぞれのイベントの発生確率は
# μ粒子の運動ベクトルの射影面積にも依存するので、後で重み付けをしています。

import random
import numpy as np
import matplotlib.pyplot as plt

n = 10000
# the information of the situation such as the distance between the two boards
length = 30
upperx = 25  # default:25
uppery = 18  # default:18

# create the zenith angle distribution by the accept-reject method
theta = []
dist = []
i = 0
while i < n:
    thetatmp = np.pi * random.random() / 2
    disttmp = random.random()
    if disttmp < np.cos(thetatmp) ** 2:
        theta = np.append(theta, thetatmp)
        dist = np.append(dist, disttmp)
        i += 1

randx = []
randy = []
phi = []
for i in range(n):
    # generate random numbers
    randx.append(random.random())
    randy.append(random.random())
    phi.append(random.random() * 2 * np.pi)

j = 0
angles = []
signal = []
effeventsarr = []
effnarr = []
actharr = []
probarr = []
while j < 91:
    degangle = j
    angle = degangle * np.pi / 180
    boardvec = [0, np.sin(angle), -np.cos(angle)]
    acctheta = []
    inproarr = []
    effevents = 0
    effn = 0
    for i in range(n):
        muvec = [np.sin(theta[i]) * np.cos(phi[i]), np.sin(theta[i]) * np.sin(phi[i]), -np.cos(theta[i])]
        upperco = [upperx * randx[i], uppery * randy[i] * np.cos(angle), uppery * randy[i] * np.sin(angle)]
        upperproj = [upperx * randx[i], uppery * randy[i], 0]
        inpro = 0
        for k in range(3):
            inpro += boardvec[k] * muvec[k]
        # the inner product also means the probability of the event
        inpro = np.absolute(inpro)
        reachco = []
        effn = effn + inpro
        for k in range(3):
            reachco.append(upperco[k] + length * (-inpro * boardvec[k] + muvec[k]))
            # this is the condition that each cosmic ray penetrate both boards.
        if 0 < reachco[0] < upperco[0] and 0 < reachco[1] < upperco[1]:
            acctheta.append(theta[i])
            inproarr.append(inpro)
            effevents += inpro
    angles = np.append(angles, degangle)
    signal = np.append(signal, len(acctheta))
    effeventsarr = np.append(effeventsarr, effevents)
    effnarr = np.append(effnarr, effn)
    actharr.append(acctheta)
    probarr.append(inproarr)
    j = j + 1

divarr = [effeventsarr[i] / effnarr[i] for i in range(91)]

# convert the angle to degree
degacctheta = [acctheta[i] * 180 / np.pi for i in range(len(acctheta))]


def cos2(x):
    return np.cos(x * np.pi / 180) ** 2


x = np.linspace(0, 90)
# plot a histogram of the angles of the penetrated particles
# the probability of the events should be taken into account as the weights.
fig1 = plt.figure(figsize=(8, 6))
ax = fig1.add_subplot(111)
plt.hist(actharr[0], bins=20, weights=probarr[0])
plt.savefig("hist.png", bbox_inches="tight")

# this graph means the dependence of the effective signals on the slope angle.
fig2 = plt.figure(figsize=(8, 6))
ax2 = fig2.add_subplot(111)
plt.plot(angles, effeventsarr, "o")
plt.savefig("effevents.png", bbox_inches="tight")

# this graph means the zenith angle distribution
fig3 = plt.figure(figsize=(8, 6))
ax3 = fig3.add_subplot(111)
plt.plot(theta, dist, "o")
plt.savefig("dist.png", bbox_inches="tight")

fig4 = plt.figure(figsize=(8, 6))
ax4 = fig4.add_subplot(111)
plt.plot(angles, effnarr, "o")
'''plt.plot(x,effnarr[0]*np.cos(x*np.pi/180))'''
plt.savefig("effn.png", bbox_inches="tight")

fig5 = plt.figure(figsize=(8, 6))
ax5 = fig5.add_subplot(111)
plt.plot(angles, divarr, "o")
plt.savefig("proportion.png", bbox_inches="tight")
