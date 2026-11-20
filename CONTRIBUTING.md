# Contributing

PEBBLE is a student entry for CERN Beamline for Schools 2027, released under the GNU General Public License v3.0.

## Licensing of contributions

By submitting a pull request, patch or any other change, you agree that your contribution is licensed under GPL-3.0, the same license as the project. This is the usual inbound = outbound arrangement: you keep copyright in what you write, and you grant the project permission to distribute it under GPL-3.0.

Every commit needs a Developer Certificate of Origin sign-off. Add it automatically:

```bash
git commit -s
```

which appends a line to the commit message:

```
Signed-off-by: Your Name <you@example.com>
```

The sign-off states that you wrote the change, or have the right to submit it, and that it may be distributed under the project license. Do not sign off on work that came from someone else without their permission.

## Practical notes

- Keep changes small and focused, one topic per pull request.
- Run `python run_pipeline.py --smoke-test` before opening a pull request, and mention in the description whether it passed.
- Do not commit `build/`, `results/raw_csv/` or anything containing credentials. They are covered by `.gitignore`.
- Open an issue before large changes so the approach can be agreed first.
