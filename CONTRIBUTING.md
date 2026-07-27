# Contributing

Open an issue before broad API or wire-format changes. Keep changes narrowly
owned, add JVM and Scala.js tests in shared sources, and run:

```sh
./tools/publish-gale-local.sh   # once per machine / Gale pin change
sbt compileAll testAll
sbt smokeCheck                  # publishedLocal consumer graph
```

For numerical behavior, include an analytic law, adversarial case, or
independent differential oracle. Preserve existing schema identifiers unless a
new version and migration path are part of the change.

Gale is consumed by Maven coordinate. Do not reintroduce a Git `ProjectRef`
dependency: published multivar POMs must declare a resolvable `gale-core`
artifact. Until Gale is on Maven Central, `tools/publish-gale-local.sh`
installs the pinned revision under a revision-qualified version.
