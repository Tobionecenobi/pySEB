"""Compact F127 galleries without discarding the model/unfolded controls."""
from collections import deque

import plotly.graph_objects as go
import pyseb
from pyseb._viewer import figure_class


def compact_scene(traces):
    geometry = [trace for trace in traces
                if (trace.meta or {}).get('pyseb', {}).get('role') == 'geometry']
    adjacency = {trace.meta['pyseb']['instance']: trace.meta['pyseb'].get('neighbors', [])
                 for trace in geometry}
    partitions = {}
    for start in adjacency:
        if start in partitions:
            continue
        partitions[start] = 0
        pending = deque([start])
        while pending:
            node = pending.popleft()
            for neighbor in adjacency.get(node, []):
                if neighbor in adjacency and neighbor not in partitions:
                    partitions[neighbor] = 1-partitions[node]
                    pending.append(neighbor)
    batches, membership = {}, {}
    for trace in geometry:
        instance = trace.meta['pyseb']['instance']
        group = instance if instance in ('core', 'anchor') else f'Arm segments:Group {partitions[instance]+1}'
        color = trace.color if trace.type == 'mesh3d' else trace.line.color
        key = (group, trace.type, color, trace.opacity)
        batch = batches.setdefault(key, dict(x=[], y=[], z=[], i=[], j=[], k=[], members=set()))
        batch['members'].add(instance)
        membership.setdefault(instance, set()).add(group)
        offset = len(batch['x'])
        for axis in ('x', 'y', 'z'):
            batch[axis].extend(getattr(trace, axis))
            if trace.type != 'mesh3d':
                batch[axis].append(None)
        if trace.type == 'mesh3d':
            for axis in ('i', 'j', 'k'):
                batch[axis].extend(int(index)+offset for index in getattr(trace, axis))
    result = []
    for (group, kind, color, opacity), batch in batches.items():
        members = batch.pop('members')
        neighbors = {target for member in members for neighbor in adjacency[member]
                     for target in membership.get(neighbor, {neighbor}) if target != group}
        common = dict(name=group, opacity=opacity, showlegend=False, hoverinfo='skip',
                      meta={'pyseb': {'instance': group, 'role': 'geometry',
                                      'neighbors': sorted(neighbors)}})
        if kind == 'mesh3d':
            trace = go.Mesh3d(**batch, color=color, **common)
        else:
            trace = go.Scatter3d(x=batch['x'], y=batch['y'], z=batch['z'],
                                mode='lines', line={'color': color, 'width': 5}, **common)
        result.append(trace)
    return result


def compact_structure_view(world, parameters, options, title):
    source = pyseb.structure_figure(world, 'micelle', parameters, options, height=420,
                                    layouts=['random', 'readable'])
    figure = figure_class()(compact_scene(source.data))
    figure.update_layout(source.layout)
    figure.update_layout(title_text=title)
    figure._pyseb_variants = [
        {**variant, 'data': [trace.to_plotly_json() for trace in compact_scene(
            go.Figure(variant['data']).data)]}
        for variant in source._pyseb_variants]
    return figure
