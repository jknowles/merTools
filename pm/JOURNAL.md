# Project journal

Append-only, newest first. **Entries are never edited** — the value of this file is
that it records what was believed at the time, including the parts that turned out
wrong. Where things stand *today* is in `STATUS.md`, which is generated.

Four lines per entry. The analysis belongs in the issue or the decision record; this
file carries the reasoning and the pointers.

- **Why** — the driver. The one line git cannot reconstruct later.
- **Obligates** — issues this change created elsewhere. Numbers, not prose.
- **Refs** — commits, issues, decision records.

---

## 2026-11-06 · quality · compass adopted

**Why:** merTools is GitHub-only, and it is the first repository compass manages there. It is in maintenance after 1.0.0, so reviews stay off to keep it out of the local review queue.
**Obligates:** nothing filed. `TODO.md` is left as it was, from the 1.0.0 submission.
**Refs:** #138, #139 (setup test issues, closed), #140 (board)
