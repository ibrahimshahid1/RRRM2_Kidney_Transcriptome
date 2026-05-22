#!/usr/bin/env python3
import re
from pathlib import Path

root = Path('.')
main = root / 'manuscript_v8_results_supplement.tex'
text = main.read_text()

pattern = re.compile(r"\\input\{([^}]+)\}")

# Replace repeatedly until no changes
while True:
    def repl(m):
        path = root / m.group(1)
        if path.exists() and path.is_file():
            return path.read_text()
        return m.group(0)
    new_text = pattern.sub(repl, text)
    if new_text == text:
        break
    text = new_text

out = root / 'manuscript_v8_results_supplement_fully_expanded.tex'
out.write_text(text)
print('Wrote', out)
