# API reference

`docsCheck` generates Scaladoc for both public modules and adds it to the
rendered guide. These APIs come from the shared sources that compile and run on
the JVM and Scala.js.

- <a href="../api/core/index.html">multivar-core API</a>
- <a href="../api/ir/index.html">multivar-ir API</a>

`multivar-core` contains analysis methods, semantic matrix geometry, solver
adapters, fitted capabilities, and workflow declarations. Most applications
need this module.

`multivar-ir` contains schemas, JSON codecs, hashes, validation, and
conformance documents for portable programs and evidence.

Generate the site and API reference from the repository root:

```sh
sbt docsCheck
```

The HTML site is written under `site/target/docs/site`. Generated files are
build output and are not committed.
