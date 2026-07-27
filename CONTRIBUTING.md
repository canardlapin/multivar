# Contributing

Open an issue before broad API or wire-format changes. Keep changes narrowly
owned, add JVM and Scala.js tests in shared sources, and run:

```sh
sbt compileAll testAll
```

For numerical behavior, include an analytic law, adversarial case, or
independent differential oracle. Preserve existing schema identifiers unless a
new version and migration path are part of the change.
