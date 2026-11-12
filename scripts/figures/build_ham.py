"""Usage: python scripts/figures/build_ham.py CONFIG.json NEW_OUTPUT"""
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'src'))
import ham
import ham_maps


def main():
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    config_path = Path(sys.argv[1]).resolve()
    config = json.loads(config_path.read_text())

    def resolve(path):
        return (config_path.parent / path).resolve()

    inputs = {}
    for name, entry in config['inputs'].items():
        path = resolve(entry['path'])
        if ham.sha256_file(path) != entry['sha256']:
            raise ValueError(f'Input hash mismatch: {path}')
        inputs[name] = path
    cotton = resolve(config['cotton_images'])
    maize = resolve(config['maize_images'])
    for name, expected in config['image_sha256'].items():
        path = (cotton if '_337_03_' in name else maize) / name
        if ham.sha256_file(path) != expected:
            raise ValueError(f'Image hash mismatch: {path}')
    sys.argv = [sys.argv[0], '--output-dir', sys.argv[2],
                '--cotton-images', str(cotton), '--maize-images', str(maize)]
    for name, path in inputs.items():
        sys.argv += ['--' + name.replace('_', '-'), str(path)]
    ham_maps.main()


if __name__ == '__main__':
    main()
