"""Pre-strip cards that Node Runner doesn't render, edit, or otherwise use.

Big aerospace decks are often dominated (by raw card count) by analysis-
only cards - nodal temperatures, dynamic load functions, response
spectra, optimization variables, output sets - that Node Runner has no
UI for. v3.1.x sent those cards through pyNastran's parser anyway,
which on a 4.6M-TEMP-card deck cost many minutes of wall clock for
zero user-visible benefit.

v3.2.2 strips these cards from the inlined buffer BEFORE pyNastran
sees them. The raw lines are stashed verbatim and written back on
export so round-trip is preserved.

The skip list is a hardcoded frozenset - no user toggle. The rationale
is that a toggle would be a UX footgun: a user who flipped it off
would just get a slow import with no UI clue why. When (if) we later
add a viewer or editor for any of these card types, that type gets
removed from the frozenset and starts going through the full parser.

If you're adding a new card type that we DON'T render: add it here.
If you're adding UI for a card type that's currently here: remove it.
"""

from __future__ import annotations

from collections import Counter
from typing import Tuple


# Card types Node Runner does not currently render, edit, or display.
# These are parsed into pyNastran's model dicts but no UI looks at them.
# Stripping them from the inline buffer is information-preserving (raw
# lines are stashed and rewritten on export) and gives an order-of-
# magnitude import speedup on decks dominated by these cards.
ANALYSIS_ONLY_CARDS = frozenset({
    # Thermal loads (nodal temperatures, etc.)
    'TEMP', 'TEMPD', 'TEMPF', 'TEMPB', 'TEMPRB', 'TEMPAX', 'TEMPP1',
    'TEMPP2', 'TEMPP3', 'TEMPP4',
    # Transient / frequency-response loads
    'TLOAD1', 'TLOAD2', 'RLOAD1', 'RLOAD2',
    'DAREA', 'DPHASE', 'DELAY',
    # Function tables (analysis input, not geometry)
    'TABLED1', 'TABLED2', 'TABLED3', 'TABLED4',
    'TABLEM1', 'TABLEM2', 'TABLEM3', 'TABLEM4',
    'TABLES1', 'TABRND1', 'TABRNDG', 'TABDMP1',
    # Optimization
    'DCONSTR', 'DESVAR', 'DLINK',
    'DRESP1', 'DRESP2', 'DRESP3',
    'DVPREL1', 'DVPREL2', 'DVCREL1', 'DVCREL2', 'DVMREL1', 'DVMREL2',
    'DOPTPRM',
    # Output-control sets (for analysis case control only)
    'SET1', 'SET2', 'SET3', 'BOUTPUT',
    # Eigenvalue control (analysis-only; we don't run modal analyses)
    'EIGRL', 'EIGR', 'EIGB', 'EIGC', 'EIGP',
})


def _card_name(line: str) -> str:
    """Return the upper-case card name on this line, or '' if the line
    isn't a card start (blank, comment, or continuation).

    Handles fixed-format (`TEMP    ...` in cols 1-8), free-field
    comma-delimited (`TEMP,...`), and the 16-char `TEMP*   ...` form.
    """
    if not line:
        return ''
    stripped = line.lstrip()
    if not stripped or stripped.startswith('$'):
        return ''
    # Continuation lines.
    if line[:1] in ('+', '*'):
        return ''
    if line[:1] == ' ' and stripped and stripped[0].isdigit():
        return ''
    # Free-field with comma.
    if ',' in stripped[:16]:
        head = stripped.split(',', 1)[0]
    else:
        # Fixed format: card name lives in cols 1-8.
        head = line[:8]
    return head.strip().upper().rstrip('*')


def strip_analysis_only_cards(flat_text: str) -> Tuple[str, str, Counter]:
    """Walk the inlined buffer line by line. Move any card that starts
    a name in ANALYSIS_ONLY_CARDS (plus its continuation lines) into a
    stash buffer. Everything else stays in the kept buffer.

    Returns ``(kept_text, skipped_raw_text, counts_by_card_name)``.

    Algorithm: single pass over splitlines(keepends=True). Track a
    state variable ``in_skip_block`` that's True while we're inside a
    skipped card (including its continuation lines). New card-start
    lines flip the state based on whether the card name is in the
    skip set. Comment and blank lines are appended to whichever
    buffer matches the current state - they're treated as part of the
    surrounding card. This keeps the visual layout intact on export.
    """
    kept_parts = []
    skipped_parts = []
    counts = Counter()
    in_skip_block = False

    # Use splitlines(keepends=True) so we preserve line endings exactly
    # - critical for round-trip on Windows decks that may use \r\n.
    for raw in flat_text.splitlines(keepends=True):
        name = _card_name(raw)
        if name:
            # New card start: decide which bucket.
            if name in ANALYSIS_ONLY_CARDS:
                in_skip_block = True
                counts[name] += 1
                skipped_parts.append(raw)
            else:
                in_skip_block = False
                kept_parts.append(raw)
        else:
            # Continuation / comment / blank: stays with current state.
            if in_skip_block:
                skipped_parts.append(raw)
            else:
                kept_parts.append(raw)

    kept_text = ''.join(kept_parts)
    skipped_text = ''.join(skipped_parts)
    return kept_text, skipped_text, counts
