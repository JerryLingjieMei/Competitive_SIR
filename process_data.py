import numpy as np
import json
import matplotlib.pyplot as plt

# result_prob = json.load(open("immunization_prob.json", "r"))
# result_prob = json.load(open("immunization_power.json", "r"))
# result_prob = json.load(open("immunization_small.json", "r"))
# result_prob = json.load(open("immunization_erdos.json", "r"))
# result_prob = json.load(open("immunization_fb.json", "r"))
# result_prob = json.load(open("immunization_fb1.json", "r"))
# result_prob = json.load(open("immunization_prob_1.json", "r"))
# result_prob = json.load(open("immunization_power_1.json", "r"))
# result_prob = json.load(open("immunization_small_1.json", "r"))
# result_prob = json.load(open("immunization_erdos_1.json", "r"))
# result_prob = json.load(open("immunization_power_2.json", "r"))
result_prob = json.load(open("immunization_small_2.json", "r"))

x = np.linspace(0, 1, 100)
u = np.zeros((50, 100))
v = np.zeros((50, 100))

for i, (result1, result2) in enumerate(result_prob):
    u[:, 99 - i] = np.log(np.array(result1)) / np.log(7057)
    v[:, 99 - i] = np.log(np.array(result2)) / np.log(7057)

y = np.zeros(100)
z = np.zeros(100)

for i in range(100):
    u0 = u[:, i]
    # u0 = u[:, max(0, i - 1):min(99, i + 2)]
    y[i] = np.mean(np.extract(u0 > .2, u0))
    v0 = v[:, i]
    # v0 = v[:, max(0, i - 1):min(99, i + 2)]
    z[i] = np.mean(np.extract(v0 > .2, v0))

plt.plot(x, y, "vr", label="Red Epidemic")
plt.plot(x, z, "^b", label="Blue Epidemic")
plt.xlabel(r"$\varphi$")
plt.ylabel(r'$\log n_\mathrm{infected}/\log_n$')

plt.xlim(.15, 1.0)

# plt.axvline(c="g", x=0.2529182916786208)
# plt.axvline(c="g", x=0.34297965187745055)
# plt.axvline(c="g", x=0.7810228345567294)
# plt.axvline(c="g", x=0.7776641044074543)
# plt.axvline(c="g", x=0.1894480522443251)
# plt.axvline(c="g", x=0.757792208977301)
# plt.axvline(c="g", x=0.4835549032247652)
# plt.axvline(c="g", x=0.6835549032247652)
# plt.axvline(c="g", x=0.4955694436766997)
# plt.axvline(c="g", x=0.48989768744529805)
# plt.axvline(c="g", x=0.5443274028526786)
plt.axvline(c="g", x=0.7861974824401493)
plt.legend()
# plt.savefig("prob.png")
# plt.savefig("power.png")
# plt.savefig("small.png")
# plt.savefig("erdos.png")
# plt.savefig("fb.png")
# plt.savefig("fb1.png")
# plt.savefig("prob1.png")
# plt.savefig("power1.png")
# plt.savefig("small1.png")
# plt.savefig("erdos1.png")
# plt.savefig("power2.png")
plt.savefig("small2.png")
