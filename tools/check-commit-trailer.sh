#!/bin/sh
# T085 + T086: enforce the per-commit invariants documented in CLAUDE.md:
#
#   - No `Co-Authored-By:` trailers (rule rolled back; AI-Assisted is now
#     the canonical attribution).
#   - No "Generated with [Claude Code]" marketing tag.
#   - Every commit on this branch (that isn't on main yet) ends with the
#     `AI-Assisted: Claude (Spec-Driven Development, TDD methodology)`
#     trailer.
#
# Run in CI's `lint` job. The trailer rule applies only to commits that
# are on the working branch but not on the base (main). Merge commits
# from the base branch are exempt.

set -eu

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

base_ref="${1:-main}"
if ! git rev-parse --verify "$base_ref" >/dev/null 2>&1; then
    echo "SKIP: base ref '$base_ref' not available — running outside a checkout that has it"
    exit 0
fi

merge_base=$(git merge-base HEAD "$base_ref")
range="${merge_base}..HEAD"
commits=$(git rev-list "$range")
if [ -z "$commits" ]; then
    echo "PASS: no commits on the branch ahead of $base_ref"
    exit 0
fi

violations=0
for sha in $commits; do
    msg=$(git log -1 --pretty=%B "$sha")
    short=$(git log -1 --pretty='%h %s' "$sha")

    # `Co-Authored-By:` must not appear as a line-anchored trailer. We
    # check for the pattern at the start of any line so it isn't fooled
    # by quoted prose discussing the rule itself.
    if echo "$msg" | grep -qiE '^Co-Authored-By:'; then
        echo "FAIL: $short has a forbidden Co-Authored-By: trailer" >&2
        violations=$((violations + 1))
    fi

    # The marketing tag emitted by older Claude Code versions is
    # `🤖 Generated with [Claude Code](https://claude.com/claude-code)`,
    # appearing as its own line. We anchor the regex to line-start so it
    # doesn't trip on prose that mentions the forbidden form (e.g. the
    # commit message that introduces this very check).
    if echo "$msg" | grep -qE '^(.*🤖 )?Generated with \[Claude Code\]'; then
        echo "FAIL: $short includes the 'Generated with [Claude Code]' marketing tag" >&2
        violations=$((violations + 1))
    fi

    # AI-Assisted trailer must appear (anywhere in the message — but
    # conventionally at the end).
    if ! echo "$msg" | grep -qE '^AI-Assisted: Claude '; then
        echo "FAIL: $short is missing the required AI-Assisted trailer" >&2
        violations=$((violations + 1))
    fi
done

if [ "$violations" -gt 0 ]; then
    echo "" >&2
    echo "Required trailer (CLAUDE.md): AI-Assisted: Claude (Spec-Driven Development, TDD methodology)" >&2
    exit 1
fi

echo "PASS: T085+T086 commit trailers OK for $(echo "$commits" | wc -l | tr -d ' ') commit(s)"
