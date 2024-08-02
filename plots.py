import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from ptolemy import Sun, Moon, Mercury, Venus, Mars, Jupiter, Saturn

# optional progress bar
try:
    from tqdm import tqdm
except ImportError:
    def tqdm(x):
        return x

def altaz_plot(bodies, alts, azs):
    alts = np.radians(alts)
    azs = np.radians(azs)
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    ax.set_rmax(1.0)
    for j, body in enumerate(bodies):
        name = body.__class__.__name__
        r = np.cos(alts[j])
        theta = -azs[j]
        r[alts[j] < 0] = np.nan
        ax.plot(theta, r, label=name)
    return fig, ax


def animation_alt_az(bodies, alts, azs, step):
    alts = np.radians(alts)
    azs = np.radians(azs)
    n = len(alts[0])
    print(f"Animating {n} frames")
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(10, 6))
    ax.set_rmax(1.0)
    lines = []
    colors = {
        "Sun": "yellow",
        "Moon": "grey",
        "Mercury": "orange",
        "Venus": "pink",
        "Mars": "red",
        "Jupiter": "brown",
        "Saturn": "blue"
    }
    
    for j, body in enumerate(bodies):
        name = body.name
        r = np.cos(alts[j, 0])
        theta = -azs[j, 0]
        line, = ax.plot(theta, r, '.', color=colors[name], label=name, markersize=10)
        lines.append(line)
    ax.legend(loc=(1.0, 0), ncol=2)
    pbar = tqdm(total=n)

    def update(i):
        # brightness following sun's height above the horizon
        sun_alt = alts[0][i] # this is in radians
        night = sun_alt < -20
        day = sun_alt >= -20
        if night:
            brightness = 0
        else:
            brightness = (np.sin(sun_alt) + 1) / 2
        color = plt.cm.Greys_r(brightness)
        ax.set_facecolor(color)
        for j, line in enumerate(lines):
            if day and (j > 1) and (bodies[j].name not in ["Sun", "Moon"]):
                r = np.nan
                theta = np.nan
            else:
                r = np.cos(alts[j, i])
                theta = -azs[j, i]
                if alts[j, i] < 0:
                    r = np.nan
                    theta = np.nan
            line.set_xdata(theta)
            line.set_ydata(r)
        pbar.update(1)
        return lines

    print(n)
    ani = FuncAnimation(fig, update, frames=n)
    ani.save("animation.gif", writer='ffmpeg', fps=10)
    pbar.close()


def animation_3d(bodies, x, y, z):
    trails = {
        "Moon": 10,
        "Mercury": 20,
        "Mars": 200,
        "Jupiter": 300,
        "Saturn": 500
    }
    colors = {
        "Earth": "green",
        "Sun": "yellow",
        "Moon": "grey",
        "Mercury": "orange",
        "Venus": "pink",
        "Mars": "red",
        "Jupiter": "brown",
        "Saturn": "blue"
    }
    n = len(x[0]) // 2
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, figsize=(12,12))
    points = []
    lines = []
    ax.plot([0], [0], [0], 'o', color=colors["Earth"], label="Earth", markersize=10)
    for j, body in enumerate(bodies):
        point, = ax.plot(x[j, 0], y[j, 0], '.', markersize=10, label=body.name, color=colors[body.name])
        line, = ax.plot(x[j, 0], y[j, 0], '-', color=colors[body.name])
        points.append(point)
        lines.append(line)
    plt.xlim(x.min(), x.max())
    plt.ylim(y.min(), y.max())
    # plt.axis('off')
    ax.set_aspect('equal',adjustable='box')
    plt.legend()
    pbar = tqdm.tqdm(total=n)
    def update(i):
        for j, (point, line) in enumerate(zip(points, lines)):
            point.set_xdata(x[j, i])
            point.set_ydata(y[j, i])
            point.set_3d_properties(z[j, i])
            trail = trails.get(bodies[j].name, 40)
            line.set_xdata(x[j, max(0, i-trail):i])
            line.set_ydata(y[j, max(0, i-trail):i])
            line.set_3d_properties(z[j, max(0, i-trail):i])
        pbar.update(1)
        return lines + points

    print(n)
    ani = FuncAnimation(fig, update, frames=n, cache_frame_data=False)
    ani.save("animation_3d.mp4", fps=40)
    pbar.close()


def main():
    sun = Sun()
    bodies = [
        sun,
        Moon(),
        Mercury(sun),
        Venus(sun),
        Mars(sun),
        Jupiter(),
        Saturn()
    ]
    nbody = len(bodies)
    # start of egyptian calendar epoch.
    # julian days start at noon, change to midnight here
    jd = 1448637.9044675925 + 0.25
    day = 1.
    year = day * 365
    hour = day / 24
    minute = hour / 60

    step = 2 * day
    n = int(50 * year / step)

    jds = jd + np.arange(n) * step
    sun.altaz_at_date(jds)

    alts = np.zeros((nbody, n))
    azs = np.zeros((nbody, n))
    x = np.zeros((nbody, n))
    y = np.zeros((nbody, n))
    z = np.zeros((nbody, n))
    for j, body in enumerate(bodies):
        x[j], y[j], z[j] = body.xyz_at_date(jds)
        alts[j], azs[j] = body.xyz_to_altaz(jds, x[j], y[j], z[j])
    # animation_3d(bodies, x, y, z)
    animation_alt_az(bodies, alts, azs, step)



if __name__ == "__main__":
    main()

