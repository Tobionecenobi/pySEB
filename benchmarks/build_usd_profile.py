"""Build a native USD stage benchmark without production profiling hooks.

An instrumented source copy lives in the chosen /tmp build directory. Exact
anchors fail loudly if the exporter changes, rather than silently mistiming a
stage. The regular library supplies all other objects; its USD object is not
extracted because the benchmark copy already defines the export symbols.
"""
import argparse
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path('seb/USDExport.cpp'))
    parser.add_argument('--build-dir', type=Path, required=True)
    parser.add_argument('--library-dir', type=Path, default=Path('build-python'))
    args = parser.parse_args()
    args.build_dir.mkdir(parents=True, exist_ok=True)
    source = args.source.read_text()

    def replace(old, new):
        nonlocal source
        if source.count(old) != 1:
            raise ValueError(f'Expected one profiling anchor: {old!r}')
        source = source.replace(old, new)

    # An envelope timer is nested inside serialization, so its duration is
    # subtracted when reporting the remaining geometry/serialization cost.
    prelude = '''#include <chrono>
extern void usdProfileStage(const char*);
extern void usdProfileEnvelope(double);
namespace { struct USDEnvelopeTimer {
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    ~USDEnvelopeTimer() { usdProfileEnvelope(std::chrono::duration<double>(
        std::chrono::steady_clock::now()-start).count()); }
}; }
'''
    replace('    const auto& scatterers=cloud.getScatterers();\n    TriangleMesh mesh;',
            '    USDEnvelopeTimer usdEnvelopeTimer;\n    const auto& scatterers=cloud.getScatterers();\n    TriangleMesh mesh;')
    replace('    const bool validExtension =', '    usdProfileStage("prepare");\n    const bool validExtension =')
    replace('    std::sort(sceneLinks.begin(),sceneLinks.end(),',
            '    usdProfileStage("layout");\n    std::sort(sceneLinks.begin(),sceneLinks.end(),')
    # Match both the original exporter and an implementation using stringstream.
    marker = '    std::ostringstream output;'
    if marker not in source:
        marker = '    std::stringstream output;'
    replace(marker, '    usdProfileStage("serialize");\n'+marker)
    replace('    const std::string temporary=path+".tmp";',
            '    usdProfileStage("write");\n    const std::string temporary=path+".tmp";')
    marker = 'throw SEBException("unable to commit USD output","World::ExportUSD");}\n}'
    replace(marker, marker[:-1]+'    usdProfileStage("");\n}')
    generated = args.build_dir/'USDExport.profile.cpp'
    generated.write_text(prelude+source)
    # Use the same optimization, PIC and language flags as the Release library.
    executable = args.build_dir/'profile-usd'
    subprocess.run(['g++', '-O3', '-DNDEBUG', '-fPIC', '-std=c++14',
                    '-DSPECIALFUNCTIONS_IMPL', '-I.', '-Iseb', '-Iseb-symbolic',
                    'benchmarks/profile_usd.cpp', str(generated),
                    str(args.library_dir/'libseb.a'), str(args.library_dir/'libseb-symbolic.a'),
                    '-lgsl', '-lgslcblas', '-lpthread', '-o', str(executable)], check=True)
    print(executable)


if __name__ == '__main__':
    main()
