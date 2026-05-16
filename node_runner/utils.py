import os

# Project root is one level up from this file's directory (node_runner/)
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def get_entity_title_from_comment(comment, entity_type, entity_id):
    """
    Intelligently parses a Nastran comment string to find a title.
    It looks for a colon as a delimiter, typical of professional-tool-generated titles.
    If no valid title format is found, it returns a default name.

    Handles multi-line comments by extracting the last non-empty comment
    line (which is typically the meaningful title line closest to the card).
    """
    default_title = f"{entity_type} {entity_id}"
    if not comment:
        return default_title

    # Split into lines and find the last non-empty comment line
    lines = comment.strip().splitlines()
    meaningful = ''
    for line in reversed(lines):
        stripped = line.strip().lstrip('$').strip()
        if stripped:
            meaningful = stripped
            break

    if not meaningful:
        return default_title

    # Look for a colon, which is the key delimiter for professional titles.
    if ':' in meaningful:
        # Split at the first colon and take the second part as the title.
        parts = meaningful.split(':', 1)
        if len(parts) > 1:
            title = parts[1].strip()
            # Return the title only if it's not empty after stripping.
            return title if title else default_title

    # If no colon but the comment text is non-empty, use it as the title.
    # This handles plain user-entered titles like "$ Aluminum".
    return meaningful
