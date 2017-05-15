import pylab as plt

fig = plt.figure()
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

def get_axis_limits(ax, x_scale=0, y_scale=1.02):
    return ax.get_xlim()[1]*x_scale, ax.get_ylim()[1]*y_scale

ax1.annotate('(a)', xy=get_axis_limits(ax1), annotation_clip=False)
ax2.annotate('(b)', xy=get_axis_limits(ax2), annotation_clip=False)
plt.show()
