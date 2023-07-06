import matplotlib.pyplot as plt
import numpy as np



from numpy import sqrt, sin, cos, tan, pi

def ellipse(x0, y0, a, b):
    X0, Y0 = x0, y0 - a
    X1, Y0 = x0 + a, y0 + a
    X = np.linspace(X0, X1, 1000)
    Y = np.concatenate([[x + sqrt(sqrt(x - x ** 6)) for x in X], [x - sqrt(sqrt(x - x ** 6)) for x in X[::-1]]])
    X = np.concatenate([np.linspace(X0, X1, 1000), np.linspace(X1, X0, 1000)])
    return X, Y

def periodic(x0, y0, T):
    X =  np.linspace(x0 - T, x0 +  T, 1000)
    Y =  [y0 + sin((x - x0)/T * 2 * pi) ** 10 for x in X]
    return X, Y

X1, Y1 = ellipse(0, 0, 1, 1)
plt.figure(figsize=(10, 6))

ax = plt.gca()
#ax.spines['left'].set_position('center')
#ax.spines['bottom'].set_position('center')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.title("Либрации", fontsize = 16)
plt.plot(X1, Y1 , linewidth = 2, color = 'red')
plt.xlabel("q", fontsize = 16)
plt.ylabel("p", fontsize = 16)
plt.show()

plt.figure(figsize=(10, 6))

ax = plt.gca()
#ax.spines['left'].set_position('center')
#ax.spines['bottom'].set_position('center')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

X, Y = periodic(0, 2, 1)

plt.title("Периодичность импульса", fontsize = 16)
plt.plot(X, Y, linewidth = 2, color = 'red')
plt.ylim((1, 3.5))
plt.xlabel("q", fontsize = 16)
plt.ylabel("p", fontsize = 16)
plt.show()
