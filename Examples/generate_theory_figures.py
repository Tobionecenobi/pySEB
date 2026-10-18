"""Reproduce the Theory page's SVGs offline; site builds use the saved images.

Run with the documentation Python environment. Use --preview-dir for PNG previews.
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, Polygon, Patch
import numpy as np

OUTPUT = Path(__file__).resolve().parents[1] / 'docs/figures'
INK = '#172f46'
BLUE = '#176b96'
ORANGE = '#b95017'
RED = '#a4333f'
GREY = '#8d9caa'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 12,
                     'svg.fonttype': 'none', 'svg.hashsalt': 'pyseb-theory'})


def save(fig, name, previews):
    path = OUTPUT / f'{name}.svg'
    fig.savefig(path, facecolor='white', metadata={'Date': None})
    path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines()) + '\n')
    if previews:
        previews.mkdir(parents=True, exist_ok=True)
        fig.savefig(previews / f'{name}.png', dpi=150, facecolor='white')
    plt.close(fig)


SHAPE_COLORS = ('#c58b19', '#318eaf', '#9b62b5', '#4c9660')
POSES = ((0., 0., 0., 0.), (0., 1.1, -.9, .7), (0., np.pi, 0., np.pi))


def shape_definitions():
    triangle = np.array([[0., 0.], [1.2, 0.], [.6, .75]])
    triangle_points = np.array([(i*triangle[0]+j*triangle[1]+(5-i-j)*triangle[2])/5
                                for i in range(6) for j in range(6-i)])
    rng = np.random.default_rng(17)
    sphere = rng.normal(size=(36, 3))
    sphere /= np.linalg.norm(sphere, axis=1)[:, None]
    sphere *= (.85*rng.random(36)**(1/3))[:, None]
    sphere += np.array([.85, 0., 0.])
    square = np.array([[0., -.55], [1.1, -.55], [1.1, .55], [0., .55]])
    rectangle = np.array([[0., -.35], [1.7, -.35], [1.7, .35], [0., .35]])
    square_points = np.array([(x, y) for x in np.linspace(0, 1.1, 5)
                              for y in np.linspace(-.55, .55, 5)])
    rectangle_points = np.array([(x, y) for x in np.linspace(0, 1.7, 6)
                                 for y in np.linspace(-.35, .35, 4)])
    flat = lambda points: np.column_stack([points, np.zeros(len(points))])
    return [('triangle', 1.2, triangle, flat(triangle_points)),
            ('sphere', 1.7, None, sphere),
            ('square', 1.1, square, flat(square_points)),
            ('rectangle', 1.7, rectangle, flat(rectangle_points))]


def chain_pose(angles):
    anchor = np.zeros(2)
    pieces = []
    for (name, length, outline, points), angle, color in zip(shape_definitions(), angles, SHAPE_COLORS):
        rotation = np.array([[np.cos(angle), -np.sin(angle)],
                             [np.sin(angle), np.cos(angle)]])
        cloud = points.copy()
        cloud[:, :2] = points[:, :2]@rotation.T+anchor
        center = np.array([.85, 0.])@rotation.T+anchor if name == 'sphere' else None
        boundary = outline@rotation.T+anchor if outline is not None else None
        end = np.array([length, 0.])@rotation.T+anchor
        pieces.append(dict(name=name, color=color, points=cloud, outline=boundary,
                           center=center, start=anchor.copy(), end=end.copy()))
        anchor = end
    bounds = []
    for piece in pieces:
        if piece['name'] == 'sphere':
            bounds.extend([piece['center']-.85, piece['center']+.85])
        else:
            bounds.extend(piece['outline'])
    bounds = np.asarray(bounds)
    shift = (bounds.min(axis=0)+bounds.max(axis=0))/2
    for piece in pieces:
        piece['points'][:, :2] -= shift
        for field in ('start', 'end', 'center', 'outline'):
            if piece[field] is not None:
                piece[field] -= shift
    return pieces


def draw_chain(ax, pieces):
    # Transparent fills let each shape remain visible inside the other shapes.
    for piece in pieces:
        if piece['name'] == 'sphere':
            ax.add_patch(Circle(piece['center'], .85, facecolor=piece['color'],
                                edgecolor='none', alpha=.24, zorder=1))
            ax.add_patch(Circle(piece['center']+[-.23, .28], .13,
                                facecolor='white', edgecolor='none', alpha=.65, zorder=2))
        else:
            ax.add_patch(Polygon(piece['outline'], closed=True, facecolor=piece['color'],
                                 edgecolor='none', alpha=.28, zorder=1))
    for piece in pieces:
        if piece['name'] == 'sphere':
            ax.add_patch(Circle(piece['center'], .85, fill=False, edgecolor=piece['color'], lw=2, zorder=3))
        else:
            ax.add_patch(Polygon(piece['outline'], closed=True, fill=False,
                                 edgecolor=piece['color'], lw=2, zorder=3))
    joints = np.array([piece['start'] for piece in pieces]+[pieces[-1]['end']])
    ax.scatter(*joints.T, s=12, color=INK, zorder=4)
    ax.set(xlim=(-3.15, 3.15), ylim=(-1.5, 1.5), aspect='equal')
    ax.axis('off')


def assumptions(previews):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), facecolor='white')
    for ax, angles, title in zip(axes, (POSES[0], POSES[2]), ('Spread out', 'Inside one another')):
        draw_chain(ax, chain_pose(angles))
        ax.set_title(title, color=INK, weight='bold', fontsize=14, pad=10)
    handles = [Patch(facecolor=color, edgecolor=color, alpha=.55) for color in SHAPE_COLORS]
    fig.legend(handles, ['Triangle', 'Sphere', 'Square', 'Rectangle'],
               loc='lower center', bbox_to_anchor=(.5, .12), ncol=4, frameon=False, fontsize=11)
    fig.text(.5, .055, 'Same subunits. Same connections. They can pass through one another like ghosts.',
             ha='center', color=INK, fontsize=12)
    fig.subplots_adjust(left=.025, right=.975, top=.84, bottom=.26, wspace=.16)
    save(fig, 'theory-ghost-subunits', previews)


def no_loops(previews):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), facecolor='white')
    pieces = chain_pose((0., np.pi/2, np.pi, -3*np.pi/4))
    start, end = pieces[0]['start'], pieces[-1]['end']
    for ax, title in zip(axes, ('Allowed: open chain', 'Not allowed: closing the loop')):
        draw_chain(ax, pieces)
        ax.set(xlim=(-2.2, 2.2), ylim=(-1.7, 1.7))
        ax.scatter(*np.array([start, end]).T, s=65, facecolors='white',
                   edgecolors=INK, linewidths=1.8, zorder=5)
        ax.set_title(title, color=INK, weight='bold', fontsize=14, pad=12)
    axes[0].text(-1.8, -1.5, 'Ends stay separate', color=INK, fontsize=12)
    axes[1].plot(*np.array([start, end]).T, '--', color=RED, lw=2.5, zorder=6)
    midpoint = (start+end)/2
    axes[1].text(*midpoint, '×', color=RED, fontsize=32, weight='bold',
                 ha='center', va='center', zorder=7,
                 bbox=dict(facecolor='white', edgecolor='none', pad=0))
    axes[1].text(-1.8, -1.5, 'Do not add this connection', color=RED, fontsize=12)
    fig.text(.5, .065, 'Subunits may overlap. Their connections must not close a loop.',
             ha='center', color=INK, fontsize=12)
    fig.subplots_adjust(left=.035, right=.965, top=.84, bottom=.19, wspace=.18)
    save(fig, 'theory-no-loops', previews)


def ensemble_curves():
    chains = [chain_pose(angles) for angles in POSES]
    q = np.linspace(0., 6., 241)
    curves = []
    for pieces in chains:
        points = np.concatenate([piece['points'] for piece in pieces])
        distances = np.linalg.norm(points[:, None, :]-points[None, :, :], axis=-1)
        curves.append(np.sinc(q[:, None, None]*distances/np.pi).mean(axis=(1, 2)))
    return chains, q, np.asarray(curves)


def ensemble(previews):
    chains, q, curves = ensemble_curves()
    assert np.allclose(curves[:, 0], 1.)
    assert np.all(curves >= -1e-12)
    fig = plt.figure(figsize=(11, 6.2), facecolor='white')
    grid = fig.add_gridspec(2, 3, height_ratios=(1, 1.15))
    for index, (pieces, title) in enumerate(zip(chains, ('Stretched', 'Bent', 'Folded'))):
        ax = fig.add_subplot(grid[0, index])
        draw_chain(ax, pieces)
        ax.set_title(title, fontsize=13, color=INK)
    ax = fig.add_subplot(grid[1, :])
    for index, curve in enumerate(curves):
        ax.plot(q, curve, color=GREY, lw=1.3, alpha=.8,
                label='Example arrangements' if index == 0 else None)
    ax.plot(q, curves.mean(axis=0), color=INK, lw=2.7, label='Their average')
    ax.set(xlabel=r'$q$ (example units)', ylabel=r'$F(q)$', xlim=(0, 6), ylim=(0, 1.05))
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(labelsize=11)
    ax.legend(fontsize=11, frameon=False)
    fig.text(.5, .53, '↓ Average their scattering', ha='center', color=INK, fontsize=13)
    fig.text(.5, .03, 'The same four shapes in three arrangements. Curves use point samples inside the shapes.',
             ha='center', color=INK, fontsize=11)
    fig.subplots_adjust(left=.07, right=.97, top=.91, bottom=.14, hspace=.6, wspace=.16)
    save(fig, 'theory-conformation-average', previews)


def intensity_curves():
    x = np.linspace(0., 8., 321)  # x = q R for a uniform sphere.
    amplitude = np.ones_like(x)
    nonzero = x != 0
    amplitude[nonzero] = 3*(np.sin(x[nonzero])-x[nonzero]*np.cos(x[nonzero]))/x[nonzero]**3
    form = amplitude**2
    beta_squared, density, background = 4., .5, .15
    return x, form, beta_squared*form, density*beta_squared*form+background, background


def intensity(previews):
    x, form, single, measured, background = intensity_curves()
    assert np.isclose(form[0], 1.)
    assert np.allclose(measured-background, .5*single)
    fig = plt.figure(figsize=(11, 3.8), facecolor='white')
    grid = fig.add_gridspec(1, 5, width_ratios=(1, .22, 1, .22, 1))
    for position, y, title, color in ((0, form, r'$F(q)$', BLUE),
                                      (2, single, r'$\beta^2 F(q)$', ORANGE),
                                      (4, measured, r'$n\beta^2 F(q)+B$', INK)):
        ax = fig.add_subplot(grid[0, position])
        ax.plot(x, y, color=color, lw=2.3)
        ax.set_title(title, color=INK, fontsize=14, pad=12)
        ax.set(xlabel=r'$qR$', xlim=(0, 8), ylim=(0, y.max()*1.1))
        ax.spines[['top', 'right']].set_visible(False)
        ax.tick_params(labelsize=10)
        if position == 4:
            ax.axhline(background, color=GREY, ls='--', lw=1.2)
            ax.text(6.5, background+.13, r'$B$', color=GREY, ha='center')
    for position, label in ((1, r'$\times\,\beta^2$'), (3, r'$\times\,n$'+'\n'+r'$+\,B$')):
        ax = fig.add_subplot(grid[0, position]); ax.axis('off')
        ax.add_patch(FancyArrowPatch((0, .5), (1, .5), transform=ax.transAxes,
                                    arrowstyle='-|>', mutation_scale=16, lw=1.8, color=INK))
        ax.text(.5, .68, label, ha='center', va='center', transform=ax.transAxes, fontsize=12, color=INK)
    fig.subplots_adjust(left=.04, right=.98, top=.77, bottom=.24, wspace=.45)
    fig.text(.5, .05, 'Illustrative sphere curves; arbitrary consistent units. Scaling and background are added separately.',
             ha='center', color=INK, fontsize=11)
    save(fig, 'theory-form-factor-to-intensity', previews)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--preview-dir', type=Path)
    previews = parser.parse_args().preview_dir
    assumptions(previews)
    no_loops(previews)
    ensemble(previews)
    intensity(previews)
