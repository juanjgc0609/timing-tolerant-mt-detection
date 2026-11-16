"""Append amendment 6 to the ATLAS-0 card, from amendment-6-draft.md, with the date of anchoring.

    python append_amendment.py "16 November 2026"            # prints the text; add --write to append it to the card

The heading, the owner's approval and everything from "Task ATLAS-0, amendment 6" on are copied unchanged; only the
date placeholder is filled and the note on the draft is dropped. Run it once, after amendment 5 is in the card.
"""
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
CARD = HERE.parent / "ATLAS-0" / "card.md"
PLACEHOLDER = "⟨date of anchoring⟩"


def text(date):
    s = (HERE / "amendment-6-draft.md").read_text(encoding="utf-8")
    head, rest = s.split("\n", 1)
    approval = rest[rest.index("Owner approval."):rest.index("Task ATLAS-0, amendment 6")]
    body = rest[rest.index("Task ATLAS-0, amendment 6"):]
    if body.count(PLACEHOLDER) != 1:
        raise SystemExit("the draft must hold the date placeholder exactly once")
    return head.replace("draft, to be anchored", "to be anchored") + "\n" + body.replace(PLACEHOLDER, date).rstrip() \
        + "\n\n" + approval.rstrip() + "\n"


def main(argv):
    if not argv or argv[0].startswith("-"):
        raise SystemExit(__doc__)
    out, card = text(argv[0]), CARD.read_text(encoding="utf-8")
    if "--write" not in argv:
        print(out)
        return 0
    if "## Amendment 6" in card:
        raise SystemExit("the card already holds amendment 6")
    if "## Amendment 5" not in card:
        raise SystemExit("amendment 5 is not in the card yet: anchor it first")
    CARD.write_text(card.rstrip("\n") + "\n\n" + out, encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
